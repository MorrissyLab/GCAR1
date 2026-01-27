library(Seurat)
library(ggplot2)
library(patchwork)
library(STdeconvolve)
library(WGCNA)
library(data.table)
library(ComplexHeatmap)

setwd("/gpnmb_car/")

primary_no_addon_path <- paste("/gpnmb_car/datasets/SPSQ_primary_HD/", "outs/", sep = '/')
biopsy_s1_path <- paste("/gpnmb_car/datasets/SPSQ_biopsy_S1_DS_RA_M/", "outs/", sep = '/')
biopsy_s2_path <- paste("/gpnmb_car/datasets/SPSQ_biopsy_S2_DS_RA_M/", "outs/", sep = '/')

# Load Samples ----
primary_no_addon <- Load10X_Spatial(data.dir = primary_no_addon_path, bin.size = c(24))
biopsy_s1 <- Load10X_Spatial(data.dir = biopsy_s1_path, bin.size = c(24))
biopsy_s2 <- Load10X_Spatial(data.dir = biopsy_s2_path, bin.size = c(24))

primary_no_addon <- subset(x = primary_no_addon, subset = nCount_Spatial.024um > 0, slot = "counts")
biopsy_s1 <- subset(x = biopsy_s1, subset = nCount_Spatial.024um > 0, slot = "counts")
biopsy_s2 <- subset(x = biopsy_s2, subset = nCount_Spatial.024um > 0, slot = "counts")

combined <- merge(primary_no_addon, y=c(biopsy_s1, biopsy_s2), 
                  add.cell.ids = c("pNA", "bS1", "bS2"), 
                  project = "combined_samples")

# Create the modify the metadata file ----
combined_meta <- combined@meta.data
combined_meta$orig.ident <- rownames(combined_meta)

combined_meta$sample <- paste(sapply(strsplit(row.names(combined_meta), "_"), "[", 1))
colnames(combined_meta)[1] <- "Spots"

# Write a tsv for the resulting sample_metadata
fwrite(combined_meta, file = "primary_biopsy_DS_RA_M_metadata_24um.tsv", row.names=FALSE, sep="\t", quote=TRUE)

# Create the cell x gene count table ----
combined_layers <- JoinLayers(combined)
combined_count_table <- combined_layers[["Spatial.024um"]]$counts

# Transpose the matrix
combined_count_table <- as.data.frame(t(as.matrix(combined_count_table)))

# Remove rows containing all zeros and then remove columns containing all zeros.
combined_count_table_filtered <- combined_count_table[rowSums(combined_count_table)>0,]
combined_count_table_filtered <- combined_count_table_filtered[,colSums(combined_count_table_filtered[])>0]

fwrite(combined_count_table_filtered, file = "primary_biopsy_DS_RA_M_count_filtered_024um.tsv", row.names=TRUE, sep="\t")

# Get the OD gene list using STdeconvolve ----
pNA_count_table <- primary_no_addon[["Spatial.024um"]]$counts
bS1_count_table <- biopsy_s1[["Spatial.024um"]]$counts
bS2_count_table <- biopsy_s2[["Spatial.024um"]]$counts

pNA_corpus <- restrictCorpus(pNA_count_table, removeAbove = 1.0, removeBelow = 0.01, alpha = 0.05,
                             plot = FALSE, verbose = TRUE, nTopOD = NA)
bS1_corpus <- restrictCorpus(bS1_count_table, removeAbove = 1.0, removeBelow = 0.01, alpha = 0.05,
                             plot = FALSE, verbose = TRUE, nTopOD = NA)
bS2_corpus <- restrictCorpus(bS2_count_table, removeAbove = 1.0, removeBelow = 0.01, alpha = 0.05,
                             plot = FALSE, verbose = TRUE, nTopOD = NA)

pNA_ODgenes_list <- rownames(pNA_corpus)
bS1_ODgenes_list <- rownames(bS1_corpus)
bS2_ODgenes_list <- rownames(bS2_corpus)

combined_ODgenes_list <- Reduce(union, c(pNA_ODgenes_list, bS1_ODgenes_list, bS2_ODgenes_list)) 
combined_ODgenes_list <- as.matrix(data.frame(combined_ODgenes_list))

write.table(combined_ODgenes_list, file = "primary_biopsy_DS_RA_M_ODgenes_list_024um.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
