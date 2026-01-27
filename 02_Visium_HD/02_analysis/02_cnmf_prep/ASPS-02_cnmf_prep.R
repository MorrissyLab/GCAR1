library(Seurat)
library(ggplot2)
library(patchwork)
library(STdeconvolve)
library(WGCNA)
library(data.table)
library(ComplexHeatmap)
setwd("/gpnmb_car/")

# 1. Load the sample metadata ----
sample_metadata <- read.csv("analysis/SPSC/samples.csv")
colnames(sample_metadata)[1] <- "ID"

# 2. Load Samples ----
for (sample in sample_metadata$ID){
  sample_path <- paste("/work/morrissy_lab/gurveer/gpnmb_car/analysis/seurat/",
                       sample,
                       "/024/PCA_50_Res_1.2/Cutoff_1/vHD_024_PCA_50_Res_1.2_post_QC_1_2.rds",
                       sep = "")
  sample_so <- readRDS(sample_path)

  assign(sample, sample_so, envir = .GlobalEnv)
}
rm(sample_so, sample_path, sample)

# 3. Create a combined Seurat object ----
so_list <- sort(ls()[grep("sample_metadata", ls(), invert = TRUE)])

# See: https://github.com/satijalab/seurat/issues/10170
SPSC_biopsy1 <- UpdateSeuratObject(SPSC_biopsy1)
SPSC_biopsy2 <- UpdateSeuratObject(SPSC_biopsy2)

combined <- merge(SPSC_biopsy1, y = c(SPSC_biopsy2),
                  add.cell.ids = c(so_list),
                  project = "combined_samples")

# 4. Create the modify the metadata file ----
combined_meta <- combined@meta.data
combined_meta$orig.ident <- rownames(combined_meta)
pattern_regex <- paste(sample_metadata$ID, collapse = "|")
combined_meta$sample <- gsub(paste0("(", pattern_regex, ").*"), "\\1", rownames(combined_meta))
colnames(combined_meta)[1] <- "Spots"

# 5. Remove low quality clusters per sample ----
cluster_info <- combined_meta %>%
  group_by(sample, seurat_clusters) %>%
  summarize(counts_avg = mean(nCount_Spatial.024um),
            counts_std = sd(nCount_Spatial.024um),
            features_avg = mean(nFeature_Spatial.024um),
            features_std = sd(nFeature_Spatial.024um)) %>%
  ungroup() %>%
  mutate_at("seurat_clusters", as.numeric) %>%
  arrange(seurat_clusters)

SPSC_biopsy1 <- subset(SPSC_biopsy1, idents = c(0,1,3,5,7,8))
SPSC_biopsy2 <- subset(SPSC_biopsy2, idents = c(2:4,6,9,10))

# 6. Recreate a combined Seurat object ----
so_list <- sort(ls()[grep("SPSC", ls())])
combined <- merge(SPSC_biopsy1, y = c(SPSC_biopsy2),
                  add.cell.ids = c(so_list),
                  project = "combined_samples")
combined_layers <- JoinLayers(combined)

# 7. Write a tsv for the resulting sample_metadata ----
combined_meta <- combined_layers@meta.data
combined_meta$orig.ident <- rownames(combined_meta)
pattern_regex <- paste(sample_metadata$ID, collapse = "|")
combined_meta$sample <- gsub(paste0("(", pattern_regex, ").*"), "\\1", rownames(combined_meta))
colnames(combined_meta)[1] <- "Spots"
fwrite(combined_meta, file = "SPSC_metadata.tsv", row.names = FALSE, sep = "\t", quote = TRUE)

# 8. Create the cell x gene count table ----
combined_count_table <- combined_layers[["Spatial.024um"]]$counts

# Transpose the matrix
combined_count_table <- as.data.frame(t(as.matrix(combined_count_table)))

# Remove rows containing all zeros and then remove columns containing all zeros.
combined_count_table_filtered <- combined_count_table[rowSums(combined_count_table)>0,]
combined_count_table_filtered <- combined_count_table_filtered[,colSums(combined_count_table_filtered[])>0]

fwrite(combined_count_table_filtered, file = "SPSC_count_filtered_024um.tsv", row.names=TRUE, sep="\t")

# 9. Get the OD gene list using STdeconvolve ----
OD_gene_list <- list()
for (sample in c(sample_metadata$ID, "combined_layers")){
  sample_so <- eval(as.name(sample))
  sample_count_table <- sample_so[["Spatial.024um"]]$counts
  sample_corpus <- restrictCorpus(sample_count_table, removeAbove = 1.0, removeBelow = 0.01, alpha = 0.05,
                                  plot = FALSE, verbose = TRUE, nTopOD = NA)
  OD_gene_list[[sample]] <- rownames(sample_corpus)
}

combined_OD_genes_list <- unique(unlist(OD_gene_list))
combined_OD_genes_list <- as.matrix(data.frame(combined_OD_genes_list))

write.table(combined_OD_genes_list, file = "SPSC_ODgenes_list_024um.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
