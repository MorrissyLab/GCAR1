#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(glmGamPoi)
  library(future)
  library(reshape2)
})


#  GCAR1 PT02 scRNA-seq preprocessing, QC, and integration
# ============================================================

samples <- c(
  "apheresis", "product", "D8", "D15_SOC",
  "D17_SOC", "D22", "D31", "D46"
)

# QC thresholds
min_features <- 200
max_features <- 6000
max_mito     <- 20


# Parallelization
plan("multicore", workers = 16)


# Directories
# -----------------------------
data_dir <- Sys.getenv("GCAR_DATA_DIR")
stopifnot(
  nzchar(data_dir),
  dir.exists(data_dir)
)

out_base <- "results"
qc_dir   <- file.path(out_base, "qc")
rds_dir  <- file.path(out_base, "seurat_objects")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)


# read + QC + filter
# -----------------------------
process_sample <- function(sample, data_dir) {

  message("Processing sample: ", sample)

  tenx_path <- file.path(
    data_dir,
    paste0("GCAR1_PT02_", sample),
    "outs/per_sample_outs",
    paste0("GCAR1_PT02_", sample),
    "count/sample_filtered_feature_bc_matrix"
  )

  sobj <- CreateSeuratObject(
    counts = Read10X(tenx_path),
    project = sample
  )

  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

  # QC plot (pre-filter)
  pdf(file.path(qc_dir, paste0(sample, "_QC_violin.pdf")),
      width = 7, height = 5)
  print(
    VlnPlot(
      sobj,
      features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
      ncol = 3,
      log = TRUE
    )
  )
  dev.off()

  # Filtering
  mt_thresh <- ifelse(sample == "product", 20, 10)

  sobj_filtered <- subset(
    sobj,
    subset =
      nCount_RNA > 5000 &
      percent.mt < mt_thresh
  )
}


seurat_list <- lapply(samples, process_sample, data_dir = data_dir)
names(seurat_list) <- samples

# -----------------------------
# Combined QC plots
# -----------------------------
meta_all <- bind_rows(
  lapply(names(seurat_list), function(nm) {
    cbind(sample = nm, seurat_list[[nm]]@meta.data)
  })
)

meta_long <- melt(
  meta_all,
  id.vars = "sample",
  measure.vars = c("nCount_RNA", "nFeature_RNA", "percent.mt")
)

pdf(file.path(qc_dir, "AllSamples_QC_violin.pdf"),
    width = 10, height = 5)
ggplot(meta_long, aes(x = sample, y = value, fill = sample)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(width = 0.1, outlier.size = 0.3) +
  facet_wrap(~variable, scales = "free_y", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "QC metrics across GCAR1 PT02 samples")
dev.off()


# SCTransform + PCA (per sample)
# -----------------------------

# Integration parameters
nfeatures_integration <- 3000


seurat_list <- lapply(seurat_list, function(obj) {
  obj <- SCTransform(
    obj,
    vars.to.regress = "percent.mt",
    method = "glmGamPoi",
    verbose = FALSE
  )
  RunPCA(obj, verbose = FALSE)
})


# RPCA integration
# -----------------------------
features <- SelectIntegrationFeatures(
  object.list = seurat_list,
  nfeatures = nfeatures_integration
)

seurat_list <- PrepSCTIntegration(
  object.list = seurat_list,
  anchor.features = features
)

anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca"
)

integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)


# Dimensional reduction + clustering
# -----------------------------
dims_use              <- 1:30
cluster_resolutions   <- c(0.5, 0.8, 1, 1.5, 2, 2.4)

DefaultAssay(integrated) <- "integrated"

integrated <- integrated |>
  RunPCA(verbose = FALSE) |>
  FindNeighbors(dims = dims_use, nn.method = "hnsw")

for (res in cluster_resolutions) {
  message("Clustering at resolution ", res)
  integrated <- FindClusters(integrated, resolution = res)
}

integrated <- RunUMAP(integrated, dims = dims_use)


# Save final object
# -----------------------------
saveRDS(
  integrated,
  file = file.path(
    rds_dir,
    "GCAR1_PT02_integrated_SCT_RPCA_UMAP.rds"
  )
)
