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
library(pheatmap)
library(gridExtra)
library(grid)

#Set working directory
setwd("/working_directory")

#Read in files
scObj <- readRDS("After_guide_and_region_filtering.rds")
metadata <- scObj@meta.data
pd <- reticulate::import("pandas")
df <- pd$read_pickle("After_guide_and_region_filtering.pkl")
df <- t(df)
source <- read.csv(file="annotation.tsv",sep="\t")

#Create cds object
Idents(scObj)<-"merged_clusters"
expr_data <- GetAssayData(scObj, slot="counts")
cell_data <- as.data.frame(scObj@meta.data)
gene_data <- as.data.frame(x=row.names(expr_data), row.names=row.names(expr_data))
colnames(gene_data) <- "gene_short_name"
cds <- new_cell_data_set(expr_data, cell_metadata=cell_data, gene_metadata=gene_data)
list_cluster<-scObj@active.ident
cds@clusters$UMAP$clusters<-list_cluster
recreate.partition<-c(rep(1,length(cds@colData@rownames)))
names(recreate.partition)<-cds@colData@rownames
cds@clusters$UMAP$partitions <- recreate.partition
recreate.partition<-as.factor(recreate.partition)
cds@int_colData@listData$reducedDims$UMAP<-scObj@reductions$UMAP@cell.embeddings
q<-plot_cells(cds,color_cells_by="cluster",label_groups_by_cluster = FALSE,group_label_size = 5)+theme(legend.position = "right")
q
q<-plot_cells(cds,color_cells_by="partition",label_groups_by_cluster = FALSE,group_label_size = 5)+theme(legend.position = "right")
q

#Learn trajectory and order cells
cds<-learn_graph(cds,use_partition=FALSE, verbose=TRUE)
q<- plot_cells(cds,color_cells_by="cluster",label_groups_by_cluster = FALSE,label_branch_points=FALSE,label_roots=FALSE,label_leaves=FALSE,group_label_size = 5)
q
cds <- order_cells(cds)
max_wt <- max(cds@principal_graph_aux@listData$UMAP$pseudotime)
cds@colData$Scaled_Pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime/max_wt
q <- plot_cells(cds,color_cells_by="Scaled_Pseudotime",label_branch_points=FALSE,label_roots=FALSE,label_leaves=FALSE,show_trajectory_graph = FALSE)
q <- q + scale_color_viridis_c("Scaled pseudotime",option="plasma") +theme(aspect.ratio=1)+xlim(-10,20)+ylim(-10,20)+xlab("UMAP_1")+ylab("UMAP_2")
q
pdf(file="allcells_pseudotime.pdf", height=5, width=5)
print(q)
dev.off()

#Save data
scObj <- AddMetaData(scObj,metadata = cds@colData$Scaled_Pseudotime,col.name = "Monocle_pseudo_wt")
meta <- scObj@meta.data
write.table(meta, file="monocle_allcells_meta.txt", sep="\t")
saveRDS(cds, file = "monocle3_allcells.rds")

#Make boxplot for NTC median pseudotime and PRGRIP1L and FEZF2 examples
metadata <- read.csv(file="monocle_allcells_meta.txt", sep="\t")
pseudotime_ntc<-data.frame(cell=rownames(metadata),pseudotime=metadata$Monocle_pseudo_wt)
pseudotime_ntc$treatment<-"NTC"
roi <- "non"
guides <-  c(source[source$intended_target==roi,1])
target <- roi
grna_matrix <- df[,colnames(df) %in% guides]
grna_matrix <- as.data.frame(grna_matrix)
grna_matrix$all <- rowSums(grna_matrix)
region_positive_df <- grna_matrix[grna_matrix$all > 0,]
region_positive_cells <- rownames(region_positive_df)
pseudotime_ntc <- pseudotime_ntc[pseudotime_ntc$cell %in% region_positive_cells,c(2,3)]
pseudotime_rpgrip1l<-data.frame(cell=rownames(metadata),pseudotime=metadata$Monocle_pseudo_wt)
pseudotime_rpgrip1l$treatment<-"RPGRIP1L"
roi <- "RPGRIP1L"
guides <-  c(source[source$intended_target==roi,1])
target <- roi
grna_matrix <- df[,colnames(df) %in% guides]
grna_matrix <- as.data.frame(grna_matrix)
grna_matrix$all <- rowSums(grna_matrix)
region_positive_df <- grna_matrix[grna_matrix$all > 0,]
region_positive_cells <- rownames(region_positive_df)
pseudotime_rpgrip1l <- pseudotime_rpgrip1l[pseudotime_rpgrip1l$cell %in% region_positive_cells,c(2,3)]
pseudotime_fezf2<-data.frame(cell=rownames(metadata),pseudotime=metadata$Monocle_pseudo_wt)
pseudotime_fezf2$treatment<-"FEZF2"
roi <- "FEZF2"
guides <-  c(source[source$intended_target==roi,1])
target <- roi
grna_matrix <- df[,colnames(df) %in% guides]
grna_matrix <- as.data.frame(grna_matrix)
grna_matrix$all <- rowSums(grna_matrix)
region_positive_df <- grna_matrix[grna_matrix$all > 0,]
region_positive_cells <- rownames(region_positive_df)
pseudotime_fezf2 <- pseudotime_fezf2[pseudotime_fezf2$cell %in% region_positive_cells,c(2,3)]
pseudotime_df<-rbind(pseudotime_ntc,pseudotime_rpgrip1l,pseudotime_fezf2)
desired_order <- c("NTC", "RPGRIP1L", "FEZF2")
pseudotime_df$treatment <- factor(pseudotime_df$treatment, levels = desired_order)
p <- ggplot(pseudotime_df, aes(x = treatment, y = pseudotime, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1.5) +
  coord_flip() +
  scale_fill_manual(values = c("lightgreen","#F8766D","#00BFC4")) +
  labs(x = "", y = "Scaled Pseudotime") +
  theme_classic() +
  theme(
    axis.text    = element_text(face = "bold", size = 12, colour = "black"),
    axis.title   = element_text(face = "bold", size = 12, colour = "black"),
    legend.position = "none", # Hide legend since treatments are labeled on x-axis
    aspect.ratio = 1
  )
p
pdf(file="boxplot.pdf", height=5, width=5)
print(p)
dev.off()

#Generate UMAPS per cell proportion significant target colored by pseudotime 
output_dir <- "pseudotime_loop_allsigtargets"
#Read in file with target names
targets_file <- "target_genes.txt"  

#Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#Set up metadata having pseudotime values
metadata <- read.csv(file = "monocle_allcells_meta.txt", sep = "\t")
scObj@meta.data <- metadata
scObj$Monocle_pseudo_wt <- as.numeric(scObj$Monocle_pseudo_wt)

#Read targets
additional_targets <- readLines(targets_file)
additional_targets <- trimws(additional_targets)
additional_targets <- additional_targets[additional_targets != ""]
all_targets <- unique(c("non", additional_targets))

#Color palette generation
if (length(all_targets) > 1) {
  target_colors <- c("lightgreen", scales::hue_pal()(length(all_targets) - 1))
} else {
  target_colors <- "lightgreen"
}
names(target_colors) <- all_targets

#Loop over targets
pseudotime_list <- list()
for (roi in all_targets) {
  message(paste0("Processing target: ", roi))
  target_lookup <- if (roi == "non") "non" else roi
  guides <- source[source$intended_target == target_lookup, 1]
  if (length(guides) == 0) {
    warning(paste("No guides found for target:", roi, "- skipping."))
    next
  }
  grna_matrix <- df[, colnames(df) %in% guides, drop = FALSE]
  grna_matrix <- as.data.frame(grna_matrix)
  if (ncol(grna_matrix) > 0) {
    grna_matrix$all <- rowSums(grna_matrix)
    region_positive_cells <- rownames(grna_matrix[grna_matrix$all > 0, ])
  } else {
    region_positive_cells <- character(0)
  }
  #Extract pseudotime data for heatmap
  pt_sub <- data.frame(
    cell = rownames(metadata),
    pseudotime = metadata$Monocle_pseudo_wt,
    stringsAsFactors = FALSE
  )
  pt_sub <- pt_sub[pt_sub$cell %in% region_positive_cells, c("pseudotime"), drop = FALSE]
  if (nrow(pt_sub) > 0) {
    pt_sub$treatment <- roi
    pseudotime_list[[roi]] <- pt_sub
  }
  scObj$perturbed <- "No"
  scObj$perturbed[rownames(scObj@meta.data) %in% region_positive_cells] <- roi
  Idents(scObj) <- "perturbed"
  if (roi %in% levels(Idents(scObj))) {
    scObj_sub <- subset(scObj, idents = roi)
    p_feat <- FeaturePlot(scObj_sub, features = "Monocle_pseudo_wt") + 
      scale_color_viridis_c("Scaled pseudotime", option = "plasma") + 
      theme(aspect.ratio = 1) + 
      xlim(-10, 20) + 
      ylim(-10, 20)
    pdf_filename <- file.path(output_dir, paste0("pseudotime_", tolower(roi), "_allcells.pdf"))
    pdf(file = pdf_filename, height = 4, width = 4)
    print(p_feat)
    dev.off()
  }
}

#Save pseudotime data per target
pseudotime_df <- do.call(rbind, pseudotime_list)
pseudotime_df$treatment <- factor(pseudotime_df$treatment, levels = all_targets)
write.table(pseudotime_df , file="pseudotime_values_for_summaryplots.txt",sep="\t")

#Heatmap
pseudotime_df <- read.csv(file="pseudotime_values_for_summaryplots.txt",sep="\t")
num_targets <- length(unique(pseudotime_df$treatment))
#Bin pseudotime into 3 intervals
pseudotime_df$pt_bin <- cut(
  pseudotime_df$pseudotime, 
  breaks = seq(0, 1, length.out = 4), 
  include.lowest = TRUE,
  labels = c("Early","Mid","Late")
)
mat <- table(pseudotime_df$treatment, pseudotime_df$pt_bin)
mat_prop <- prop.table(mat, margin = 1)
ntc_row <- mat_prop["non", , drop = FALSE]
targets_mat <- mat_prop[rownames(mat_prop) != "non", ]
#Order target rows descending based on the "Early" bin column
target_order <- order(targets_mat[, "Early"], decreasing = TRUE)
sorted_targets_mat <- targets_mat[target_order, , drop = FALSE]
ordered_matrix <- rbind(ntc_row, sorted_targets_mat)
annotation_row <- data.frame(
  Control = ifelse(rownames(ordered_matrix) == "non", "NTC (Baseline)", "Target")
)
rownames(annotation_row) <- rownames(ordered_matrix)
anno_colors <- list(
  Control = c("NTC (Baseline)" = "lightgreen", "Target" = "grey85")
)
color_palette <- colorRampPalette(c("#0571b0", "lightyellow", "#ca0020"))(100) # Blue -> White -> Red
genes_to_show <- c("FEZF2", "RPGRIP1L", "non")
custom_row_labels <- ifelse(
  rownames(ordered_matrix) == "non", "NTC",
  ifelse(rownames(ordered_matrix) %in% c("FEZF2", "RPGRIP1L"), rownames(ordered_matrix), "")
)
p <- pheatmap(
  ordered_matrix,
  cluster_cols = FALSE,       
  cluster_rows = FALSE,       
  show_rownames = TRUE,       
  labels_row = custom_row_labels, 
  annotation_row = annotation_row,
  annotation_colors = anno_colors,
  color = color_palette,
  legend_breaks = seq(0, 1, by = 0.2),
  main = "Cell Distribution Across Pseudotime Bins (NTC at Top)"
)
p
pdf(file = "heatmap_pseudotime_6x4.pdf", height=6,width=4)
print(p)
dev.off()

#NTC density ridge
ntc_df <- pseudotime_df %>% filter(treatment == "non")
p_ntc <- ggplot(ntc_df, aes(x = pseudotime)) +
  geom_density(fill = "lightgreen", color = "black", alpha = 0.8) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "NTC Baseline Density", y = "Density") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title   = element_text(face = "bold", size = 11, hjust = 0.5)
  )
pdf(file="ntc_ridge_heatmaptop.pdf",height=2,width=4)
print(p_ntc)
dev.off()