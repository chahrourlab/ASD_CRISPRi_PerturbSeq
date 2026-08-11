setwd("/working_directory")
library(tidyverse)
library(pheatmap)
library(reshape2)
library(tidygraph)
library(ggraph)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggnewscale)

#Load the significant annotated DEGs
df <- read.delim("Significant_annotated_DEGs.txt", sep="\t", stringsAsFactors = FALSE)
#Read in perturbed targets that are not SFARI genes
nonsfari <- read.delim("nonsfari.txt", sep="\t", stringsAsFactors = FALSE)
nonsfari <- c(nonsfari$Gene)

#Create the Signed Gene identifier
sig_degs_signed <- df %>%
  mutate(direction = ifelse(log2fc > 0, "UP", "DOWN"),
         signed_gene = paste0(DEG, "_", direction)) %>%
  select(target, signed_gene)
signed_gene_sets <- split(sig_degs_signed$signed_gene, sig_degs_signed$target)

#Initialize the Jaccard matrix
pert_names <- names(signed_gene_sets)
n <- length(pert_names)
signed_jaccard_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(pert_names, pert_names))

#Calculate Jaccard for the signed sets
for (i in 1:n) {
  for (j in i:n) {
    set1 <- signed_gene_sets[[i]]
    set2 <- signed_gene_sets[[j]]
    
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    
    j_index <- if (union > 0) intersection / union else 0
    
    signed_jaccard_matrix[i, j] <- j_index
    signed_jaccard_matrix[j, i] <- j_index
  }
}

#Prepare links with filtered weight
links <- melt(signed_jaccard_matrix) %>%
  rename(from = Var1, to = Var2, weight = value) %>%
  filter(from != to, weight > 0.1) 
write.csv(links, file="jaccard_links_afterFDR.csv")

#Prepare Nodes
nodes <- data.frame(name = unique(c(as.character(links$from), as.character(links$to))))

#Update the graph data
graph <- tbl_graph(nodes = nodes, edges = links, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(
    community = as.factor(group_louvain()),
    degree = centrality_degree(),
    short_label = sub("_.*", "", name),
    is_nonsfari = name %in% nonsfari
  )

#Plot
p <- ggraph(graph, layout = "nicely") + 
  geom_edge_link(aes(edge_alpha = weight, edge_width = weight), 
                 color = "royalblue") + 
  scale_edge_width(name = "Jaccard Index", range = c(0.2, 2)) +
  scale_edge_alpha(name = "Jaccard Index", range = c(0.1, 1)) +
  geom_node_point(aes(color = community, size = degree)) +
  scale_color_brewer(palette = "Set1") +
  new_scale_color() + 
geom_node_text(aes(label = short_label, color = is_nonsfari), 
               repel = TRUE, 
               size = 3, 
               fontface = "bold", 
               family = "Arial",
               max.overlaps = 15) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), 
                     guide = "none") + 
  
  theme_void() +
  theme(
    text = element_text(family = "Arial", face = "bold"),
    plot.title = element_text(size = 14, margin = margin(b = 5)),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Jaccard Similarity Network",
    color = "Cluster",
    size = "Centrality"
  )

ggsave("Jaccard_Network_Plot.pdf", plot = p, device = cairo_pdf, width = 11, height = 9)

#Read in per target summary file listing the number of significant DEGs in each category
df <- read.csv(file="Summary_per_target.txt", sep="\t")
df$Target <- as.character(df$Target)

#Barplot function
save_formatted_plot <- function(plot_obj, filename) {
  ggsave(filename, plot = plot_obj, width = 5, height = 6, dpi = 300)
}
plain_theme <- theme_classic() + 
  theme(
    text = element_text(family = "Arial", face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

#Total DEGs
plot1_data <- df %>% arrange(desc(Count.of.DEGs)) %>% head(10)
p1 <- ggplot(plot1_data, aes(x = reorder(Target, -Count.of.DEGs), y = Count.of.DEGs)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  plain_theme +
  labs(title = "Top 10 Perturbations: Total DEG Count", x = "Target Gene", y = "Total DEGs")
save_formatted_plot(p1, "Top10_Total_DEGs.png")

#SFARI DEGs
plot2_data <- df %>% arrange(desc(Count.of.SFARI.genes)) %>% head(10)
p2 <- ggplot(plot2_data, aes(x = reorder(Target, -Count.of.SFARI.genes), y = Count.of.SFARI.genes)) +
  geom_bar(stat = "identity", fill = "coral") +
  plain_theme +
  labs(title = "Top 10 Significant SFARI Enrichments", x = "Target Gene", y = "Count of SFARI DEGs")
save_formatted_plot(p2, "Top10_SFARI_DEGs.png")

#TF DEGs
plot3_data <- df %>% arrange(desc(Count.of.TFs)) %>% head(10)
p3 <- ggplot(plot3_data, aes(x = reorder(Target, -Count.of.TFs), y = Count.of.TFs)) +
  geom_bar(stat = "identity", fill = "mediumseagreen") +
  plain_theme +
  labs(title = "Top 10 Significant TF Enrichments", x = "Target Gene", y = "Count of TF DEGs")
save_formatted_plot(p3, "Top10_TF_DEGs.png")

#NDD DEGs
plot4_data <- df %>% arrange(desc(Count.of.BrainspanNDD)) %>% head(10)
p4 <- ggplot(plot4_data, aes(x = reorder(Target, -Count.of.BrainspanNDD), y = Count.of.BrainspanNDD)) +
  geom_bar(stat = "identity", fill = "orchid") +
  plain_theme +
  labs(title = "Top 10 Significant Brainspan NDD Enrichments", x = "Target Gene", y = "Count of NDD DEGs")
save_formatted_plot(p4, "Top10_NDD_DEGs.png")

#lncRNA DEGs
plot5_data <- df %>% arrange(desc(Count.of.lncRNAs)) %>% head(10)
p5 <- ggplot(plot5_data, aes(x = reorder(Target, -Count.of.lncRNAs), y = Count.of.lncRNAs)) +
  geom_bar(stat = "identity", fill = "gold") +
  plain_theme +
  labs(title = "Top 10 Significant lncRNA Enrichments", x = "Target Gene", y = "Count of lncRNA DEGs")
save_formatted_plot(p5, "Top10_lncRNA_DEGs.png")

#Read in per DEG summary file listing the number of times a DEG occurs across pertubations 
sens_df <- read.delim("sensitive_genes.txt", sep="\t", stringsAsFactors = FALSE)

#Barplot function
plain_theme_sens <- theme_classic() + 
  theme(
    text = element_text(family = "Arial", face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

#Global sensitive DEGs 
p1_data <- sens_df %>%
  select(DEG, count) %>%
  filter(!is.na(DEG) & DEG != "") %>%
  arrange(desc(count)) %>%
  head(10)
p1_sens <- ggplot(p1_data, aes(x = reorder(DEG, -count), y = count)) +
  geom_bar(stat = "identity", fill = "#a6cbe3") + 
  plain_theme_sens +
  labs(title = "Top 10 Sensitive DEGs", x = "Gene Symbol", y = "Perturbation Count")
ggsave("Top10_Sensitive_Total.png", p1_sens, width = 5, height = 6)

#Sensitive SFARI genes (Cols 3 & 4)
p2_data <- sens_df %>%
  select(sfari, sfari_count) %>%
  filter(!is.na(sfari) & sfari != "" & sfari != "None") %>%
  arrange(desc(sfari_count)) %>%
  head(10)
p2_sens <- ggplot(p2_data, aes(x = reorder(sfari, -sfari_count), y = sfari_count)) +
  geom_bar(stat = "identity", fill = "#f8b195") +
  plain_theme_sens +
  labs(title = "Top 10 Sensitive SFARI Genes", x = "SFARI Gene", y = "Perturbation Count")
ggsave("Top10_Sensitive_SFARI.png", p2_sens, width = 5, height = 6)

#Sensitive Brainspan NDD genes (Cols 5 & 6)
p3_data <- sens_df %>%
  select(brainspan, brainspan_count) %>%
  filter(!is.na(brainspan) & brainspan != "" & brainspan != "None") %>%
  arrange(desc(brainspan_count)) %>%
  head(10)
p3_sens <- ggplot(p3_data, aes(x = reorder(brainspan, -brainspan_count), y = brainspan_count)) +
  geom_bar(stat = "identity", fill = "#e1bee7") +
  plain_theme_sens +
  labs(title = "Top 10 Sensitive NDD Genes", x = "NDD Gene", y = "Perturbation Count")
ggsave("Top10_Sensitive_NDD.png", p3_sens, width = 5, height = 6)

#Sensitive lncRNAs (Cols 7 & 8)
p4_data <- sens_df %>%
  select(lncrna, lncrna_count) %>%
  filter(!is.na(lncrna) & lncrna != "" & lncrna != "None") %>%
  arrange(desc(lncrna_count)) %>%
  head(10)
p4_sens <- ggplot(p4_data, aes(x = reorder(lncrna, -lncrna_count), y = lncrna_count)) +
  geom_bar(stat = "identity", fill = "#fff176") +
  plain_theme_sens +
  labs(title = "Top 10 Sensitive lncRNAs", x = "lncRNA", y = "Perturbation Count")
ggsave("Top10_Sensitive_lncRNA.png", p4_sens, width = 5, height = 6)

#Sensitive transcription factors (Cols 9 & 10) ---
p5_data <- sens_df %>%
  select(TF, TF_count) %>%
  filter(!is.na(TF) & TF != "" & TF != "None") %>%
  arrange(desc(TF_count)) %>%
  head(10)
p5_sens <- ggplot(p5_data, aes(x = reorder(TF, -TF_count), y = TF_count)) +
  geom_bar(stat = "identity", fill = "#b7e4c7") +
  plain_theme_sens +
  labs(title = "Top 10 Sensitive Transcription Factors", x = "TF Gene", y = "Perturbation Count")
ggsave("Top10_Sensitive_TF.png", p5_sens, width = 5, height = 6)

#Create the data frame for most frequently occuring lncRNA-PCG pairs
df <- data.frame(
  Pairs = c("LINC00261;FOXA2", "MIR137HG;DPYD", "XACT;AMOT", "LINC00491;SLCO6A1", 
            "LINC01090;TFPI", "DLX6-AS1;DLX5", "H19;IGF2", "DANT1;PLS3", 
            "LINC00645;FOXG1", "MEG3;DLK1", "TEX41;ARHGAP15"),
  Occurrences = c(17, 14, 10, 9, 9, 7, 7, 6, 6, 5, 5)
)
#Generate the barplot
p<- ggplot(df, aes(x = reorder(Pairs, -Occurrences), y = Occurrences)) +
  geom_bar(stat = "identity", fill = "#fff176", color = "black") +
  theme_classic() +
  labs(
    title = "Number of Occurrences per Pair",
    x = "Pairs",
    y = "Number of Occurrences"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10) # Rotates long labels
  )
p
ggsave("Top_freq_dysregulated_pairs.png", p, width = 8, height = 6)

#Create the data frame for most frequently dysregulated TE families
df <- data.frame(
  Element = c("LTR", "DNA", "LINE", "SINE", "Satellite", "Retroposon"),
  Occurrences = c(2662, 1035, 429, 194, 60, 25)
)
#Generate the barplot
p <- ggplot(df, aes(x = reorder(Element, -Occurrences), y = Occurrences)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  theme_classic() +
  labs(
    title = "Number of Occurrences by Element Type",
    x = "Element Type",
    y = "Number of Occurrences"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10) # Keeps labels readable
  )
p
ggsave("Top_TE_fams.png", p, width = 5, height = 6)


#Create the data frame for the top 10 targets that dysregulate TEs
df <- data.frame(
  Gene = c("CHAMP1", "CSDE1", "SIN3B", "MYH9", "BTAF1", "POGZ", "CNOT3", "ZMYM2", "C11orf30", "HIRA"),
  Occurrences = c(132, 104, 92, 70, 46, 40, 38, 36, 35, 32)
)
#Generate the barplot
p <- ggplot(df, aes(x = reorder(Gene, -Occurrences), y = Occurrences)) +
  geom_bar(stat = "identity", fill = "darkgrey", color = "black") +
  theme_classic() +
  labs(
    title = "Number of Occurrences by Gene",
    x = "Gene",
    y = "Number of Occurrences"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10) # Rotates labels for clarity
  )
p
ggsave("Top_TE_perts.png", p, width = 5, height = 6)
