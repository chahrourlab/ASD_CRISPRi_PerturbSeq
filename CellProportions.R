#Load libraries
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(speckle)
library(limma)
library(reticulate)
library(monocle3)
library(ggridges)
library(ggrepel)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)

#Set working directory
setwd("/working_directory")

#Read in files
scObj <- readRDS("After_guide_and_region_filtering.rds")
cells <- rownames(scObj@meta.data)
pd <- import("pandas")
df <- pd$read_pickle("After_guide_and_region_filtering.pkl")
df <- t(df)
df <- as.data.frame(df)
annotation <- read.csv(file="annotation.tsv",sep="\t")
metadata <- scObj@meta.data

#Set up library id in metadata
lib_map <- setNames(
  paste0("LW", 558:581),
  1:24                   
)
meta_base <- scObj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(
    lib_num = str_extract(cell_id, "(?<=-)\\d+"),
    library = lib_map[lib_num]
  ) %>%
  select(-lib_num)

#NTC cells
ntc_guides <- annotation[annotation$intended_target == "non", 1]
ntc_matrix <- df[, colnames(df) %in% ntc_guides, drop = FALSE]
ntc_cells <- rownames(ntc_matrix)[rowSums(ntc_matrix) > 0]

#Unique target genes
genes <- unique(annotation$intended_target)
genes <- setdiff(genes, c("non")) 

#Clusters
clusters <- unique(meta_base$merged_clusters)

#CMH test loop (Iterating over target genes and clusters)
cmh_results <- map_dfr(genes, function(roi) {
  #Identify guides matching the current target gene
  guides <- annotation[annotation$intended_target == roi, 1]
  #Subset guide matrix for this ROI and identify positive cells
  grna_matrix <- df[, colnames(df) %in% guides, drop = FALSE]
  if (is.null(dim(grna_matrix))) return(NULL)
  pert_cells <- rownames(grna_matrix)[rowSums(grna_matrix) > 0]
  #Filter cells to only compare ROI-perturbed cells vs. NTC cells
  sub_meta <- meta_base %>%
    filter(cell_id %in% c(pert_cells, ntc_cells)) %>%
    mutate(
      is_pert = factor(
        ifelse(cell_id %in% pert_cells, "Perturbed", "NTC"),
        levels = c("Perturbed", "NTC")
      )
    )
  #Skip if we don't have at least 50 cells for either group
  if (sum(sub_meta$is_pert == "Perturbed") < 50) return(NULL)
  #Test across each cluster
  map_dfr(clusters, function(cls) {
    sub_cluster_meta <- sub_meta %>%
      mutate(
        is_cluster = factor(
          ifelse(merged_clusters == cls, "In_Cluster", "Out_Cluster"),
          levels = c("In_Cluster", "Out_Cluster")
        )
      )
    #Construct 3D Contingency Table: [Perturbation Status x Cluster Status x Stratification Factor]
    tab <- table(
      sub_cluster_meta$is_pert,
      sub_cluster_meta$is_cluster,
      sub_cluster_meta$library 
    )
    #Run CMH test
    test_res <- tryCatch({
      res <- mantelhaen.test(tab)
      data.frame(
        target = roi,
        cluster = cls,
        n_pert = length(pert_cells),
        pval = res$p.value,
        odds_ratio = as.numeric(res$estimate),
        log_odds_ratio = log(as.numeric(res$estimate))
      )
    }, error = function(e) {
      data.frame(
        target = roi,
        cluster = cls,
        n_pert = length(pert_cells),
        pval = NA_real_,
        odds_ratio = NA_real_,
        log_odds_ratio = NA_real_
      )
    })
    
    return(test_res)
  })
})

#FDR correction
cmh_results_clean <- cmh_results %>%
  filter(!is.na(pval)) %>%
  mutate(
    padj = p.adjust(pval, method = "BH"),
    log_padj = -log10(padj),
    log_padj_capped = pmin(log_padj, 5), 
    signed_log_padj = sign(log_odds_ratio) * log_padj_capped
  )
write.table(cmh_results_clean, file = "cmh_stats_percluster.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Empirical permutation test against background
set.seed(42)
n_perms <- 100 

#Build background log odds ratios across shuffled gRNA labels
null_dist <- map_dfr(1:n_perms, function(i) {
  #Permute row names (cells) relative to gRNA matrix
  shuffled_df <- df
  rownames(shuffled_df) <- sample(rownames(df))
  #Pick a random set of 6 guides to act as a "pseudo target"
  random_guides <- sample(colnames(shuffled_df), size = 6) 
  pseudo_matrix <- shuffled_df[, random_guides, drop = FALSE]
  pseudo_pert <- rownames(pseudo_matrix)[rowSums(pseudo_matrix) > 0]
  #Create null metadata subset
  sub_null <- meta_base %>%
    filter(cell_id %in% c(pseudo_pert, ntc_cells)) %>%
    mutate(is_pert = factor(
      ifelse(cell_id %in% pseudo_pert, "Perturbed", "NTC"),
      levels = c("Perturbed", "NTC")
    ))
  #Skip if insufficient cells in pseudo-target
  if (sum(sub_null$is_pert == "Perturbed") < 10) return(NULL)
  #Test across each cluster
  map_dfr(clusters, function(cls) {
    sub_cluster_null <- sub_null %>%
      mutate(is_cluster = factor(
        ifelse(merged_clusters == cls, "In_Cluster", "Out_Cluster"),
        levels = c("In_Cluster", "Out_Cluster") # Col 1 = In_Cluster, Col 2 = Out_Cluster
      ))
    tab_null <- table(
      sub_cluster_null$is_pert,
      sub_cluster_null$is_cluster,
      sub_cluster_null$library
    )
    tryCatch({
      res <- mantelhaen.test(tab_null)
      data.frame(
        perm_id = i,
        cluster = cls,
        null_log_or = log(as.numeric(res$estimate))
      )
    }, error = function(e) {
      NULL
    })
  })
}, .progress = TRUE) %>%
  filter(!is.na(null_log_or) & !is.infinite(null_log_or))
write.table(null_dist, file="nulldist.txt",sep="\t")

#Summarize background distribution (Mean and SD) per cluster
null_summary <- null_dist %>%
  group_by(cluster) %>%
  summarise(
    null_mean = mean(null_log_or),
    null_sd = sd(null_log_or),
    n_valid_perms = n(),
    .groups = "drop"
  )

#Compare real CMH results against empirical background distribution
cmh_empirical <- cmh_results_clean %>%
  inner_join(null_summary, by = "cluster") %>%
  mutate(
    #Score observed log OR relative to null noise
    t_stat = (log_odds_ratio - null_mean) / null_sd,
    #Two-tailed p-value using t-distribution 
    emp_pval = 2 * pt(-abs(t_stat), df = n_valid_perms - 1),
    #Adjust p-values for FDR
    emp_padj = p.adjust(emp_pval, method = "BH"),
    emp_log_padj = -log10(pmax(emp_padj, 1e-300)),
    signed_emp_log_padj = sign(log_odds_ratio) * pmin(emp_log_padj, 5) # Capped at 5
  )
write.table(cmh_empirical, file = "cmh_stats_empirical.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Significant hits
clusters_order <- c("NPC", "mixture", "neuron")
plot_base_df <- cmh_empirical %>%
  mutate(
    cluster = factor(cluster, levels = clusters_order),
    padj_safe = pmax(padj, 1e-300),
    emp_padj_safe = pmax(emp_padj, 1e-300),
    minus_log10_padj = -log10(padj_safe),
    plot_size_val = pmin(minus_log10_padj, 50),
    is_dual_sig = (padj < 0.05 & emp_padj < 0.05)
  )
dual_hits <- plot_base_df %>%
  filter(is_dual_sig)
write.table(
  dual_hits, 
  file = "cmh_stats_dual_fdr_filtered.txt", 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)

#Volcano plot per cluster
for (cls in clusters_order) {
  cls_data <- plot_base_df %>% 
    filter(cluster == cls) %>%
    mutate(
      significance_cat = case_when(
        is_dual_sig & log_odds_ratio > 0 ~ "Enriched",
        is_dual_sig & log_odds_ratio < 0 ~ "Depleted",
        TRUE ~ "Not Significant"
      ),
      significance_cat = factor(
        significance_cat, 
        levels = c("Enriched", "Depleted", "Not Significant")
      )
    )
  if (nrow(cls_data) == 0) next
  top_enriched_labels <- cls_data %>%
    filter(significance_cat == "Enriched") %>%
    slice_min(order_by = padj_safe, n = 5, with_ties = FALSE)
  top_depleted_labels <- cls_data %>%
    filter(significance_cat == "Depleted") %>%
    slice_min(order_by = padj_safe, n = 5, with_ties = FALSE)
  top_labels <- bind_rows(top_enriched_labels, top_depleted_labels)
  
  p_volcano <- ggplot(cls_data, aes(x = log_odds_ratio, y = minus_log10_padj)) +
    geom_point(
      aes(color = significance_cat, alpha = significance_cat), 
      size = 2
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.6) +
    scale_color_manual(
      values = c(
        "Enriched" = "#a50026",      
        "Depleted" = "#313695",      
        "Not Significant" = "grey70" 
      ),
      drop = FALSE
    ) +
    scale_alpha_manual(
      values = c(
        "Enriched" = 0.9,
        "Depleted" = 0.9,
        "Not Significant" = 0.35
      ),
      drop = FALSE
    ) +
    geom_text_repel(
      data = top_labels,
      aes(label = target),
      size = 3.5,
      fontface = "italic",
      max.overlaps = 15,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "grey50",
      segment.size = 0.2
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 10, face = "bold"),
      panel.grid.major = element_line(color = "grey92")
    ) +
    labs(
      title = paste0("Volcano Plot: ", cls),
      subtitle = "Red = Enriched, Blue = Depleted | Labeled: Top 5 Enriched & Top 5 Depleted",
      x = "Log Odds Ratio (Depleted < 0 < Enriched)",
      y = expression(-log[10]~"(CMH "*padj*")")
    )
  clean_cls_name <- gsub("[^A-Za-z0-9_]", "_", as.character(cls))
  ggsave(paste0("volcano_", clean_cls_name, ".pdf"), plot = p_volcano, width = 6.5, height = 5.5)
}