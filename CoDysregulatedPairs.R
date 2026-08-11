#Load libraries
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(dplyr)

#Set working directory
setwd("/working_directory")

#Load the lncRNA coordinates file downloaded from BioMart
lncRNA_df <- read.csv(file="lncrna_coordinates.txt",sep="\t")
lncRNA_gr <- GRanges(
  seqnames = lncRNA_df$Chromosome.scaffold.name,
  ranges = IRanges(start = lncRNA_df$Gene.start..bp., end = lncRNA_df$Gene.end..bp.),
  strand = lncRNA_df$Strand,
  lncRNA_name = lncRNA_df$lncRNA_gene_name
)
edb <- EnsDb.Hsapiens.v86
all_ensemble_genes <- genes(edb)
all_ensemble_genes <- genes(edb, filter = GeneBiotypeFilter("protein_coding"))
all_ensemble_genes <- genes(edb, filter = GeneBiotypeFilter(value = "lincRNA", condition = "!="))
seqlevelsStyle(lncRNA_gr) <- "Ensembl"
nearest_genes_hits <- distanceToNearest(lncRNA_gr, all_ensemble_genes)
lncRNA_indices <- queryHits(nearest_genes_hits)
gene_indices <- subjectHits(nearest_genes_hits)
distances <- mcols(nearest_genes_hits)$distance
closest_gene_gr <- all_ensemble_genes[gene_indices]
results_df <- data.frame(
  lncRNA_name = lncRNA_df$lncRNA_gene_name[lncRNA_indices],
  lncRNA_chr = lncRNA_df$Chromosome.scaffold.name[lncRNA_indices],
  lncRNA_start = lncRNA_df$Gene.start..bp.[lncRNA_indices],
  lncRNA_stop = lncRNA_df$Gene.end..bp.[lncRNA_indices],
  closest_gene_ensembl_id = names(closest_gene_gr),
  closest_gene_symbol = closest_gene_gr$gene_name, # EnsDb usually provides this
  closest_gene_biotype = closest_gene_gr$gene_biotype, # Here's the biotype!
  closest_gene_chr = as.character(seqnames(closest_gene_gr)),
  closest_gene_start = start(closest_gene_gr),
  closest_gene_stop = end(closest_gene_gr),
  distance_to_gene = distances
)
write.table(results_df,file="lncRNA_closestPCG_annotated.txt",sep="\t")

#Load the significant annotated DEGs and lncRNA-PCG file generated above 
deg_df <- read.csv(file="annotated_DEGs.txt",sep="\t")
lncrna_annotation_df <- read.csv(file="lncRNA_closestPCG_annotated.txt",sep="\t")

#Find co-dysregulated pairs
lncrna_degs <- deg_df %>%
  dplyr::filter(lncRNAs == TRUE) %>%
  dplyr::select(gene_names, region, gene, log2fc, Significance_score, abs.sig,key) %>%
  dplyr::rename(
    lncRNA_DEG_name = gene_names,
    lncRNA_DEG_log2fc = log2fc,
    lncRNA_DEG_Significance_score = Significance_score,
    lncRNA_DEG_abs_sig = abs.sig
  )
merged_lncrna_pcg <- lncrna_degs %>%
  inner_join(lncrna_annotation_df %>% dplyr::select(lncRNA_name, closest_gene_symbol,distance_to_gene),
             by = c("lncRNA_DEG_name" = "lncRNA_name"))
pcg_degs <- deg_df %>%
  dplyr::filter(lncRNAs == FALSE) %>%
  dplyr::select(gene_names, region, gene, log2fc, Significance_score, abs.sig,key) %>%
  dplyr::rename(
    PCG_DEG_name = gene_names,
    PCG_DEG_log2fc = log2fc,
    PCG_DEG_Significance_score = Significance_score,
    PCG_DEG_abs_sig = abs.sig
  )

final_output <- merged_lncrna_pcg %>%
  inner_join(pcg_degs,
             by = c("key", "closest_gene_symbol" = "PCG_DEG_name"))

final_output_df <- final_output %>%
  dplyr::select(
    key,
    region = region.x,
    target = gene.x,
    lncRNA_DEG_name,
    lncRNA_DEG_log2fc,
    lncRNA_DEG_Significance_score,
    closest_gene_symbol,
    PCG_DEG_log2fc,
    PCG_DEG_Significance_score,
    distance_to_gene
  )
write.table(final_output_df,file="co-dysregulated_pairs.txt",sep="\t")