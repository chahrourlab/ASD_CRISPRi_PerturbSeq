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

#Set working directory
setwd("/working_directory")

#Read in an rds file containing processed single-cell data, including a PCA embedding, converted from h5ad format
scObj_full <- readRDS(file="data.rds")

#Read in the gRNA dataframe in pkl format
pd <- import("pandas")
df <- pd$read_pickle("sgRNA_dataframe.pkl")
df <- t(df)

#Read in the gRNA annotations
annotation <- read.csv(file="annotations.tsv",sep="\t")

#Read in good gRNAs - ie non-outliers from energy distance output
non_outlier_targeting <- read.csv(file="targeting_outlier_table.csv")
non_outlier_targeting <- non_outlier_targeting[non_outlier_targeting$pval_outlier > 0.05,]
non_outlier_nontargeting <- read.csv(file="non_targeting_outlier_table.csv")
non_outlier_nontargeting <- non_outlier_nontargeting[non_outlier_nontargeting$pval_outlier > 0.05,]
non_outlier <- rbind(non_outlier_targeting,non_outlier_nontargeting)
df <- df[, colnames(df) %in% non_outlier$X]
df <- as.data.frame(df)

#Get rid of cells that have 0 guides assigned now
df$row_sum <- rowSums(df)
subset_df <- subset(df, row_sum > 0)
subset_df <- subset_df[, !(names(subset_df) == "row_sum")]

#get rid of cells that have too many unique regions targeted 
barplot <- read.csv(file="unique_regions_per_cell.txt",sep="\t",header=FALSE)
barplot_summary <- barplot %>%
  group_by(V2) %>%
  summarize(Cell_Count = n())
p <- ggplot(barplot_summary, aes(x = V2, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +  
  labs(title = "Number of Cells vs. Number of uniquely targeted regions",
       x = "Number of uniquely targeted regions",
       y = "Number of Cells") +
  theme_bw()  
p
pdf(file="unique_regions_per_cell.pdf")
print(p)
dev.off()
quantiles <- quantile(barplot$V2, probs = c(0.25, 0.75))
#To detect too many targeted cells: ((1.5*(Q3-Q1)) + Q3)
#Q1 is 2 and Q3 is 5. 1.5x(5-2) + 5  = 9.5
barplot_subset <- barplot[barplot$V2 <10,]
subset_df <- subset_df[rownames(subset_df) %in% barplot_subset$V1,]
cells <- as.data.frame(rownames(subset_df))
guides <- as.data.frame(colnames(subset_df))
write.table(cells,file="cells.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(guides,file="guides.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

#Create a filtered Seurat object
scObj <- subset(scObj_full,cells=c(cells$`rownames(subset_df)`))
Idents(scObj)<-"louvain"
p<- DimPlot(scObj, reduction = "UMAP", label=FALSE, raster=FALSE)+theme(aspect.ratio=1)+xlim(-10,20)+ylim(-10,20)+xlab("UMAP_1")+ylab("UMAP_2")
p
pdf(file="UMAP_after_guide_and_region_filtering.pdf")
print(p)
dev.off()
saveRDS(scObj, file = "After_guide_and_region_filtering.rds")
#Subset the sgRNA dataframe pkl file to include the same cells and guides in cells.txt and guides.txt and save as After_guide_and_region_filtering.pkl


#Merge clusters and annotate
scObj <- readRDS("After_guide_and_region_filtering.rds")
metadata<-scObj@meta.data
f_cortex<-c("PAX6","OTX2","DCX","TUBB3")
p<-DotPlot(scObj, features=f_cortex, dot.min=0.1,cols=c("gray95", "purple")) + theme(axis.text.x = element_text(angle = 90, hjust=1))
p
pdf(file="Unmerged_dotplot.pdf")
print(p)
dev.off()
metadata$merged_clusters <- metadata$louvain
metadata["merged_clusters"][metadata["merged_clusters"] == 2] <- "NPC"
metadata["merged_clusters"][metadata["merged_clusters"] == 4] <- "NPC"
metadata["merged_clusters"][metadata["merged_clusters"] == 6] <- "NPC"
metadata["merged_clusters"][metadata["merged_clusters"] == 0] <- "mixture"
metadata["merged_clusters"][metadata["merged_clusters"] == 3] <- "mixture"
metadata["merged_clusters"][metadata["merged_clusters"] == 5] <- "mixture"
metadata["merged_clusters"][metadata["merged_clusters"] == 7] <- "mixture"
metadata["merged_clusters"][metadata["merged_clusters"] == 8] <- "mixture"
metadata["merged_clusters"][metadata["merged_clusters"] == 1] <- "neuron"
scObj@meta.data<-metadata
clusters <- metadata %>% group_by(merged_clusters) %>% tally()
#Visualize merged clustering
Idents(scObj)<-"merged_clusters"
our_colors=c("#118941","steelblue","#D55E00")
p<- DimPlot(scObj, reduction = "UMAP", label=FALSE, raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-10,20)+ylim(-10,20)+xlab("UMAP_1")+ylab("UMAP_2")
p
pdf(file="merged_umap.pdf")
print(p)
dev.off()
f_cortex<-c("PAX6","OTX2","DCX","TUBB3")
levels(scObj) <- c("NPC","mixture","neuron")
p<-DotPlot(scObj, features=f_cortex, dot.min=0.1,cols=c("gray95", "purple")) + theme(axis.text.x = element_text(angle = 90, hjust=1))
p
pdf(file="merged_dotplot.pdf")
print(p)
dev.off()
p<-FeaturePlot(scObj, features=f_cortex)& WhiteBackground() & theme(aspect.ratio=1) &xlim(-10,20)&ylim(-10,20) &xlab("UMAP_1")&ylab("UMAP_2")
p
pdf(file="featureplot.pdf")
print(p)
dev.off()
saveRDS(scObj, file = "After_guide_and_region_filtering.rds")

#NTC plot
pd <- import("pandas")
df <- pd$read_pickle("After_guide_and_region_filtering.pkl")
df <- t(df)
df <- as.data.frame(df)
annotation <- read.csv(file="annotation.tsv",sep="\t")
goi <- colnames(df)[startsWith(colnames(df), "non-targeting")]
#get the number of cells that have these guides and the barcodes
grna_matrix <- df[,colnames(df) %in% goi]
grna_matrix <- as.data.frame(grna_matrix)
grna_matrix$all <- rowSums(grna_matrix)
region_positive_df <- grna_matrix[grna_matrix$all > 0,]
region_positive_cells <- rownames(region_positive_df)
#Visualize the specific perturbation on the UMAP
metadata <- scObj@meta.data
metadata$perturbed <- "Yes"
metadata$perturbed[rownames(metadata) %in% region_positive_cells] <- "NTC"
scObj@meta.data<-metadata
Idents(scObj)<-"perturbed"
perturbed <- metadata %>% group_by(perturbed) %>% tally()
p<- DimPlot(scObj, reduction = "UMAP", label=FALSE, raster=FALSE,cols=c("grey","black"))+theme(aspect.ratio=1)+xlim(-10,20)+ylim(-10,20)+xlab("UMAP_1")+ylab("UMAP_2")
p
pdf(file="NTC_UMAP.pdf")
print(p)
dev.off()