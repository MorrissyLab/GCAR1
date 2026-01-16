t#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggalluvial)
  library(tidyr)
})

# GCAR1 PT02  manuscript plots (Fig 4)
# ============================================================

# Parameters

sample_order <- c(
  "apheresis", "harvest", "D8", "D15_SOC",
  "D17_SOC", "D22", "D31", "D46"
)

col_gcar_pos <- "#be0119"
col_gcar_neg <- "#d8dcd6"

out_dir <- "results/figures/Fig4"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load integrated object

integrated <- readRDS("results/seurat_objects/GCAR1_PT02_integrated_SCT_RPCA_UMAP.rds")
meta <- integrated@meta.data


## Fig 4  GCAR⁺ proportion over time (bar + flow)
# ============================================================

df_gcar <- meta %>%
  group_by(orig.ident) %>%
  summarise(
    total_cells = n(),
    gcar_pos_cells = sum(GCAR_pos, na.rm = TRUE),
    prop = gcar_pos_cells / total_cells,
    .groups = "drop"
  ) %>%
  mutate(
    orig.ident = factor(orig.ident, levels = sample_order),
    x_label = paste0(orig.ident, "\n(", format(total_cells, big.mark = ","), ")"),
    label = ifelse(
      prop < 0.01,
      sprintf("%.2f%%\n(%s)", prop * 100, gcar_pos_cells),
      sprintf("%.1f%%\n(%s)", prop * 100, gcar_pos_cells)
    )
  )

p_fig4d <- ggplot(df_gcar, aes(x = orig.ident, y = prop, group = 1)) +
  geom_flow(
    aes(alluvium = 1, stratum = orig.ident),
    fill = col_gcar_pos, alpha = 0.2, width = 0.5
  ) +
  geom_col(fill = col_gcar_pos, width = 0.8) +
  geom_text(aes(label = label), vjust = -0.3, size = 5) +
  scale_y_continuous(
    limits = c(0, 0.20),
    labels = percent_format(),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_discrete(labels = df_gcar$x_label) +
  labs(x = "Sample / Timepoint", y = "GCAR1⁺ cells") +
  theme_classic(base_size = 13)

ggsave(
  file.path(out_dir, "Fig4d_GCARpos_proportion_bar_alluvial.pdf"),
  p_fig4d, width = 10, height = 6
)


## Fig  UMAP GCAR⁺ vs GCAR⁻
# ============================================================
umap_df <- FetchData(integrated, vars = "GCAR_pos") %>%
  cbind(Embeddings(integrated, "umap"))

p_fig4e <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  geom_point(
    data = subset(umap_df, GCAR_pos == FALSE),
    color = col_gcar_neg, size = 0.6
  ) +
  geom_point(
    data = subset(umap_df, GCAR_pos == TRUE),
    color = col_gcar_pos, size = 0.6
  ) +
  theme_classic(base_size = 14) +
  labs(title = "GCAR1 expression")

ggsave(
  file.path(out_dir, "Fig4e_UMAP_GCARpos_vs_neg.pdf"),
  p_fig4e, width = 8, height = 6
)

## Fig 4  GCAR+ cells composition by cluster level annotation
# ============================================================

# Remove cell-level / cluster-level clashes for T-cell only
# Cell-level (SCimilarity / hint-based)
tcell_labels_hint <- c(
  "CD4", "CD4cm", "CD4em", "CD4helper", "CD4naive",
  "CD8", "CD8cytotoxic", "CD8em", "CD8em_td", "CD8naive",
  "MAIT", "matureNKT", "T_cell", "Treg", "lymphocyte"
)

# Cluster-level (Seurat Final_Annotation)
tcell_labels_final <- c(
  "CD4", "CD4em", "CD4naive",
  "CD8em", "CD8em_td", "CD8naive",
  "MAIT", "Tcell_mix", "Treg"
)

#Assign T-cell status at both levels
df <- df %>%
  mutate(
    Tcell_status_cell =
      ifelse(celltype_hint_Final_Annotation_v2 %in% tcell_labels_hint,
             "Tcell", "nonTcell"),

    Tcell_status_cluster =
      ifelse(Final_Annotation %in% tcell_labels_final,
             "Tcell", "nonTcell")
  )
table(df$Tcell_status_cell, df$Tcell_status_cluster)

# Keep only concordant cells
df_matched <- df %>%
  filter(Tcell_status_cell == Tcell_status_cluster)

#Cluster-level composition of GCAR pos T cells
sample_order <- c(
  "product", "D8", "D15_SOC",
  "D17_SOC", "D22", "D31", "D46"
)

annotation_colors <- c(
  classical_monocyte       = "#0000FF",
  "non-classical_monocyte"   = "#36648B",
  "classical_monocyte_CD14+" = "#27408B",

  CD4naive = "#FF7F00",
  CD4em    = "#6A33C2",
  CD4      = "#FF0000",

  CD8e = "#FFD700",
  CD8em = "#FB6496",
  CD8em_td = "#FDBF6F",
  CD8naive = "#C68642",

  Treg      = "#33A02C",
  Tcell_mix = "#B15928",

  MAIT = "#C814FA",

  NK_cells_CD16pos_CD56neg = "#00CED1",
  NK_cells_CD16neg_CD56pos = "#20B2AA",
  NK_cells                 = "#00E2E5",

  cDC            = "#EEE685",
  pDC            = "#BDB76B",
  DC_progenitors = "#9ACD32",

  B = "#1F78C8"
)

df_prop <- df_matched %>%
  filter(
    GCAR_pos == TRUE,
    orig.ident %in% sample_order
  ) %>%
  count(orig.ident, Final_Annotation) %>%
  group_by(orig.ident) %>%
  mutate(
    prop = n / sum(n),
    total_cells = sum(n)
  ) %>%
  ungroup()

df_prop$orig.ident <- factor(df_prop$orig.ident, levels = sample_order)

p <- ggplot(df_prop, aes(
  x = orig.ident,
  y = prop,
  fill = Final_Annotation
)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = annotation_colors, drop = TRUE) +
  labs(
    y = "Proportion of GCAR pos cells",
    fill = "Cluster annotation"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90))

pdf("GCARpos_FinalAnnotation_clusterlevel_matched.pdf", 6, 6)
print(p)
dev.off()

########################### CD4/CD8 only ######################## 

cd4_types <- c("CD4naive", "CD4em", "CD4", "Treg")

meta <- integrated@meta.data %>%
  filter(GCAR_pos == TRUE, orig.ident %in% sample_order)


  cd4_counts <- meta %>%
  filter(Final_Annotation %in% cd4_types) %>%
  count(orig.ident, Final_Annotation)

totals <- meta %>%
  count(orig.ident, name = "total_cells")

df_cd4 <- cd4_counts %>%
  left_join(totals, by = "orig.ident") %>%
  mutate(
    prop = n / total_cells,
    orig.ident = factor(orig.ident, levels = sample_order)
  )

 ggplot(df_cd4,
       aes(
         x = orig.ident,
         y = prop,
         stratum = Final_Annotation,
         alluvium = Final_Annotation,
         fill = Final_Annotation
       )) +
  geom_flow(alpha = 0.3, width = 0.5) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = annotation_colors[cd4_types]) +
  labs(
    x = "Timepoint",
    y = "Proportion of all GCAR⁺ cells",
    fill = "CD4 T-cell subtype"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

pdf("GCARpos_CD4_subtypes_proportion.pdf", 7, 6)
print(last_plot())
dev.off()

### same for CD8 cd8_types <- c("CD8em_td", "CD8e", "CD8em", "CD8naive", "MAIT")


## Fig 4 GCAR pos T-cell dotplot
###########################
dict_marker_genes <- list(
  "stem-like_memory" = c("HNRNPLL","CCR7","SELL","LEF1","IL7R","TCF7"),
  "activation_cytotoxicity_effector" = c(
    "CD28","CD27","GZMK","GZMA","PRF1","CCL4","GZMM",
    "CCL5","NKG7","GZMB","GNLY","LYAR","CXCR3","CXCR4","TXNIP"
  ),
  "pre-exhaustion_exhaustion" = c(
    "TIGIT","ENTPD1","HAVCR2","CTLA4","LAG3","PDCD1",
    "EOMES","TOX","PMCH","BATF","PRDM1","CXCR1"
  ),
  "resident" = c("MCM5","TNFRSF9","ITGAE","NR4A2")
)
genes_use <- unlist(dict_marker_genes)


## T-cell subtypes to keep
Tcell_types <- c(
  "CD4naive", "CD4em", "CD4", "Treg",
  "CD8em_td", "CD8e", "CD8em", "CD8naive",
  "MAIT", "Tcell_mix"
)


## Sample order (top → bottom in dotplot)
sample_order <- c(
  "apheresis", "product", "D8", "D15_SOC",
  "D17_SOC", "D22", "D31", "D46"
)

## Subset GCAR⁺ T cells
obj <- subset(
  integrated,
  subset = GCAR_pos == TRUE & Final_Annotation %in% Tcell_types
)
DefaultAssay(obj) <- "RNA"
expr <- as.matrix(GetAssayData(obj, layer = "data"))
expr_df <- as.data.frame(expr)
expr_df <- tibble::rownames_to_column(expr_df, var = "gene")

## Pivot to long format
expr_long <- expr_df %>%
  pivot_longer(
    cols = -gene,
    names_to = "cell",
    values_to = "expr"
  )

## Add metadata (timepoints / samples)
meta_df <- obj@meta.data %>%
  tibble::rownames_to_column(var = "cell") %>%
  select(cell, orig.ident)

expr_long <- expr_long %>%
  left_join(meta_df, by = "cell")

## Filter to marker genes
expr_long <- expr_long %>%
  filter(gene %in% genes_use)

## Summarize: % cells expressing + avg expression
dotplot_df <- expr_long %>%
  group_by(orig.ident, gene) %>%
  summarise(
    pct_expr = mean(expr > 0) * 100,
    avg_expr = mean(expr),
    .groups = "drop"
  )

## Factor ordering for dotplot
dotplot_df <- dotplot_df %>%
  mutate(
    gene = factor(gene, levels = genes_use),
    orig.ident = factor(orig.ident, levels = rev(sample_order))
  )

## Dotplot
p <- ggplot(dotplot_df, aes(x = gene, y = orig.ident)) +
  geom_point(aes(size = pct_expr, color = avg_expr)) +
  scale_color_gradientn(
    colors = c("#FED976","#FD8D3C","#f62929ff","#ab090cff"),
    name = "Mean expression\nin group"
  ) +
  scale_size(range = c(0,6), name = "% cells") +
  labs(x = "Gene", y = "Timepoint") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    panel.grid = element_blank()
  )

pdf("GCARpos_Tcells_marker_dotplot.pdf", width = 12, height = 5)
print(p)
dev.off()
