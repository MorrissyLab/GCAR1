library(ActivePathways)
library(GSEABase)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(circlize)
library(tidyr)
library(gprofiler2)
library(purrr)
library(tibble)
library(lsa)
library(RColorBrewer)
library(GSVA)
library(Seurat)
library(msigdbr)
library(stringr)
library(parallel)
library(gridExtra)
library(fastcluster)
library(cowplot)
library(ggrepel)
library(data.table)
library(future)
library(future.apply)
library(doParallel)
library(rstatix)
library(paletteer)
library(igraph)
library(ggraph)
library(ALDEx2)
library(BiocParallel)
library(scico)

args <- commandArgs()

setwd("~/Documents/GitHub/GCAR1/")
SAVE_FOLDER <- "02_Visium_HD/04_results/ASPS-02/"


PROGRAM_ANNOTATION <- NULL
RANK <- 15
GLOBAL_THRESHOLD <- 0.05
sample_metadata <- data.frame(sample = c("SPSC_biopsy1", "SPSC_biopsy2"))
sample_metadata$TREATMENT <- c("Biopsy_1", "Biopsy_2")
sample_metadata$LOCATION <- c(rep("Lung", 2))
sample_metadata$COLS <- paletteer_d("ggsci::dark_uchicago")[6:7]
sample_cols <- setNames(sample_metadata$COLS, sample_metadata$sample)
pattern_regex <- paste(sample_metadata$sample, collapse = "|")

# MGS Parameters and Switches ----
SPECTRA <- "02_Visium_HD/03_datasets/02_cnmf_outputs/ASPS-02/cNMF_5_100_5_F.gene_spectra_score.k_15.dt_0_1.txt"
GMT_PATH <- "/bulk/morrissy_bulk/GMT/"

GENE_RANK_PLOTTING <- T
COMPUTE_MGS <- T
TOP_GENES_SPECTRA <- T
gPROFILER_RESULTS <- T

# Program Analysis Parameters and Switches ----
usage <- "02_Visium_HD/03_datasets/02_cnmf_outputs/ASPS-02/cNMF_5_100_5_F.usages.k_15.dt_0_1.consensus.txt"
usage <- data.frame(fread(file = usage, header = TRUE, sep = "\t"), row.names = 1)
colnames(usage) <- paste("GEP", seq(1, length(usage)), sep = "")
usage_norm <- t(apply(usage, 1, function (x) x/sum(x)))
usage_norm <- usage_norm[complete.cases(usage_norm),]
program_annotation <- sprintf(paste0("P%0", str_length(RANK), "d"), c(1:RANK))
program_annotation_list <- c("pneumocyte", "T_metabolism", "T", "T",
                             "CD8", "T", "T", "macrophage", "VSMC_pericyte",
                             "cell_cycle", "CAF_ECM", "endothelial",
                             "lung_goblet_cell", "lung_ciliated_cell", "MES")
program_annotations <- program_annotation

PROGRAM_ENRICHMENT <- T
PROGRAM_GENE_EXPRESSION <- T

# Neighbourhood Analysis Parameters and Switches ----
NEIGHBORHOOD_ANALYSIS_PROGRAMS <- T
NEIGHBOURHOOD_ANALYSIS_GCAR_CELLS <- T
NEIGHBORHOOD_ANALYSIS_SUMMARY <- T
NEIGHBORHOOD_ANALYSIS_SUMMARY_PROGRAMS <- c("5")

# Gene Expression Analysis Parameters and Switches ----
GENE_EXPRESSION_LIST <- c("GCAR", "GPNMB", "PDCD1",
                          "CXCL9", "CXCL10", "CCL5", "PPBP", "GBP1",
                          "GBP5", "IFI30", "CD74", "TAP1", "TIMP1",
                          "SPP1", "ENPP2", "IDO1", "FN1", "SERPINE1",
                          "CD274", "PVR", "NECTIN2", "TREM1", "TREM2", "SPP1")
GENE_EXPRESSION_NORMALIZE <- list("1" = c("GPNMB"),
                                  "2" = c("CD3D", "CD3E", "CD3G", "CD4", "CD8A"))

INSIDE_OUTSIDE_T_CELL <- T
INSIDE_OUTSIDE_NICHE <- T
INSIDE_OUTSIDE_GENES <- data.frame(category = c(rep("Pair_1", 2),
                                                rep("Pair_2", 3),
                                                rep("Pair_3", 5),
                                                rep("Pair_4", 2),
                                                rep("Pair_5", 2),
                                                rep("Pair_6", 3),
                                                rep("Pair_7", 3)),
                                   gene = c("PDCD1", "CD274",
                                            "CTLA4", "CD80", "CD86",
                                            "TIGIT", "PVR", "NECTIN2", "CD226", "CD96",
                                            "KLRB1", "CLEC2D",
                                            "FAS", "FASLG",
                                            "HAVCR2", "LGALS9", "CEACAM1",
                                            "LAG3", "FGL1", "LGALS3"))

PROGRAM_EXPRESSION_ASSOCIATIONS <- T
BIN_FRACTION <- T

# Save global variable names ----
program_var_names <- c(ls(), "program_var_names")

# General Functions ----
create_folder <- function(path){
  if (file.exists(path)){
  } else {
    dir.create(file.path(path), recursive = TRUE)
  }
}

gmt_to_list <- function(gmt_obj){
  geneset_list <- list()

  for (pathway in gmt_obj){
    geneset_name <- pathway$id
    geneset_genes <- pathway$genes

    if (geneset_name %in% names(geneset_list)){
      geneset_name_repeat_suffix <- length(which(names(geneset_list) == geneset_name)) + 1
      geneset_name <- paste(geneset_name, geneset_name_repeat_suffix, sep = ".")
    }

    geneset_list[[geneset_name]] <- geneset_genes
  }

  return(geneset_list)
}

# MGS Functions ----
load_spectra <- function(spectra_file){
  gene_score <- read.table(spectra_file, sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  gene_score <- t(as.matrix(gene_score))
}

ranked_spectra <- function(spectra_df, num_genes = "all", program_annotation = NULL){
  gene_score <- spectra_df

  # Set the max genes to include for each GEP
  if (num_genes == "all"){
    num_genes <- length(gene_score[,1])
  }

  # Add the top genes from each program from the selected rank to the GEP_df
  GEP_df <- data.frame(V1 = 1:num_genes)
  for (program in 1:ncol(gene_score)){
    sorted_spectra_df <- data.frame(genes = rownames(gene_score), value = gene_score[, program])
    sorted_spectra_df <- sorted_spectra_df[order(match(sorted_spectra_df$value, sort(sorted_spectra_df$value, decreasing = TRUE))),]

    top_genes <- sorted_spectra_df[1:num_genes, 1]
    GEP_df[program] <- top_genes
  }

  # Rename the columns for the GEP_df
  if (is.null(program_annotation)){
    program_annotation <- seq(1, ncol(gene_score))
    program_annotation <- ifelse(program_annotation <= 9, paste("GEP", program_annotation, sep = "_0"), paste("GEP", program_annotation, sep = "_"))
  }
  colnames(GEP_df) <- program_annotation

  return(GEP_df)
}

marker_gene_score <- function(gep_df, gene_set, scale = FALSE, detailed = FALSE, filter = FALSE,
                              species = "Human", rank = 30){
  gene_score <- gep_df
  rownames(gene_score) <- paste("Usage_", seq(1, nrow(gene_score)), sep="")
  gene_score <- t(gene_score)
  gene_score <- data.frame(gene_score)

  final_df <- data.frame()

  if (detailed == FALSE){
    # For each pathway, get usage of the pathway genes in the GEP.
    if (!detailed){
      for (pathway in 1:length(geneset)){
        pathway_name <- geneset[[pathway]]$id
        marker_genes <- geneset[[pathway_name]]$genes

        # The human gene IDs are all upper case, but mouse ones are in title case.
        if (species == "Human"){
          marker_genes <- tolower(marker_genes)
          marker_genes <- toupper(marker_genes)
        }else if (species == "Mouse"){
          marker_genes <- tolower(marker_genes)
          marker_genes <- tools::toTitleCase(marker_genes)
        }
        marker_genes <- gsub(" ", "", unique(marker_genes))

        pathway_gene_df <- data.frame(cell_type = pathway_name, n_markers = length(marker_genes))
        sub1 <- data.frame(Index = 1)

        for (GEP in 1:nrow(gene_score)){
          int_N <- length(intersect(gep_df[, GEP], marker_genes))
          ranks <- 1/(which(gep_df[, GEP] %in% intersect(gep_df[, GEP], marker_genes)))
          ranks_sum <- sum(ranks)

          marker_score <- int_N*ranks_sum
          marker_score <- marker_score/log10(length(marker_genes))

          sub2 <- data.frame(topN_marker_score = marker_score)
          colnames(sub2) <- paste("K", rank, "_GEP", GEP, "_", colnames(sub2), sep="")
          colnames(sub2) <- gsub("topN", paste("top", dim(gep_df)[1], sep=""), colnames(sub2))

          sub1 <- cbind(sub1, sub2)
        }
        sub1 <- cbind(pathway_gene_df, sub1)
        final_df <- rbind(final_df, sub1)
      }

      rownames_names <- names(geneset)
      rownames_names <- make.unique(rownames_names)
      rownames(final_df) <- rownames_names
    }
  } else {
    for (pathway in 1:length(geneset)){
      dataset_name = geneset[[pathway]]$name
      pathway_name <- geneset[[pathway]]$id

      marker_genes <- geneset[[pathway_name]]$genes

      # The human gene IDs are all upper case, but mouse ones are in title case.
      if (species != "Human"){
        marker_genes <- tolower(marker_genes)
        marker_genes <- tools::toTitleCase(marker_genes)
      }
      marker_genes <- gsub(" ", "", unique(marker_genes))

      marker_genes <- setdiff(marker_genes, "")

      pathway_gene_df <- data.frame(dataset = dataset_name, cell_type = pathway_name, n_markers = length(marker_genes))
      sub1 <- data.frame(Index = 1)

      for (GEP in 1:nrow(gene_score)){
        int_N <- length(intersect(gep_df[, GEP], marker_genes))
        ranks <- 1/(which(gep_df[, GEP] %in% intersect(gep_df[, GEP], marker_genes)))
        ranks_sum <- sum(ranks)

        marker_score <- int_N*ranks_sum
        marker_score <- marker_score/log10(length(marker_genes))

        intersect_genes = paste(intersect(gep_df[, GEP], marker_genes ), collapse="|")
        intersect_genes_rank = paste(which(gep_df[, GEP] %in% intersect(gep_df[,GEP], marker_genes )), collapse="|")

        sub2 <- data.frame(n_Int = int_N, genes_Int = intersect_genes, ranks_Int = intersect_genes_rank, topN_marker_score = marker_score)
        colnames(sub2) <- paste("K", rank, "_GEP", GEP, "_", colnames(sub2), sep="")
        colnames(sub2) <- gsub("topN", paste("top", dim(gep_df)[1], sep=""), colnames(sub2))

        pathway_gene_df <- cbind(pathway_gene_df, sub2)
      }
      max_score <- max(pathway_gene_df[,grep("_marker_score", colnames(pathway_gene_df), value = TRUE)])
      pathway_gene_df <- pathway_gene_df %>% add_column(maxScore = max_score, ".after" = "n_markers")
      final_df <- rbind(final_df, pathway_gene_df)
    }
  }

  if (filter) {
    if (detailed == FALSE){
      # If there's a cell_type, that is not present in the GEPs, remove the row
      final_df <- final_df[!apply(final_df[, 4:ncol(final_df)], 1, function(row) {
        all(is.na(row) | is.infinite(row) | row == 0 | is.nan(row))
      }), ]
    }
    else if (detailed == TRUE){
      # If there's a cell_type, that is not present in the GEPs, remove the row
      final_df <- final_df[!is.nan(final_df$maxScore),]
      final_df <- final_df[final_df$maxScore != 0,]
    }
  }

  if (scale){
    final_df <- scale(t(final_df))
  }

  return(final_df)
}

gene_position <- function(gene_by_position_df, gene_list){
  # Gets the rank of select genes across every column (program/cluster)
  gene_position_df <- data.frame(gene_name = NA, row = NA, col = NA)

  for (gene in gene_list){
    temp_locations <- which(gep_df == gene, arr.ind = TRUE)
    temp_locations <- as.data.frame(temp_locations)
    temp_locations$gene_name <- gene
    gene_position_df <- rbind(gene_position_df, temp_locations)
  }

  gene_position_df <- gene_position_df[-1,]

  gene_position_df$row_inverse <- (gene_position_df$row ** -1)
  gene_position_df$row_inverse_logged <- log(gene_position_df$row_inverse)

  return(gene_position_df)
}

find_element_indices <- function(df, element){
  sapply(df, function(column) {
    index <- which(column == element)
    if (length(index) == 0) NA else index
  })
}

# Plotting specific functions
gmt_heatmap <- function(final_df, type = 'cnmf',
                        file_name = "", file_height = 10, file_width = 25){
  # Remove the text based columns, empty rows, and empty columns
  plot_df <- final_df[-c(1,3)]
  plot_df <- plot_df[complete.cases(plot_df),]
  plot_df <- plot_df[apply(plot_df, 1, sum) != 0,]
  plot_df <- plot_df[apply(plot_df, 2, sum) != 0]

  col_markers <- colorRamp2(breaks = c(0, max(plot_df$n_markers)), hcl_palette = "Prgn")
  haR <- HeatmapAnnotation(GMT_nr_marker = plot_df$n_markers, which = 'row',
                           col = list(GMT_nr_marker = col_markers))

  plot_df <- plot_df[-1]
  plot_df <- as.matrix(plot_df)

  col_matrix <- colorRamp2(breaks = c(0, max(plot_df)),
                           hcl_palette = "Blues", reverse = TRUE)

  if (type == "cnmf"){
    ht_clustered <- Heatmap(plot_df, right_annotation = haR, column_title_rot = 45,
                            cluster_columns = TRUE, cluster_rows = TRUE, col = col_matrix)

    ht_col_clustered <- Heatmap(plot_df, right_annotation = haR, column_title_rot = 45,
                                cluster_columns = TRUE, cluster_rows = FALSE, col = col_matrix)

    ht_row_clustered <- Heatmap(plot_df, right_annotation = haR, column_title_rot = 45,
                                cluster_columns = FALSE, cluster_rows = TRUE, col = col_matrix)

    ht_unclustered <- Heatmap(plot_df, right_annotation = haR, column_title_rot = 45,
                              cluster_columns = FALSE, cluster_rows = FALSE, col = col_matrix)
  }

  pdf(file_name, height = file_height, width = file_width)
  draw(ht_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_col_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_row_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_unclustered, padding = unit(c(3, 3, 3, 3), "cm"))
  dev.off()
}

gene_position_linegraph <- function(gep_gene_position_df,
                                    plot_title = "Temp Plot Title (K5)",
                                    plot_x_title = "Gene Expression Program",
                                    file_name = "temp.pdf",
                                    file_height = 10, file_width = 25){

  gep_gene_position <- gep_gene_position_df

  gg_r <- ggplot(gep_gene_position, aes(x=col, y=row)) +
    geom_line() + scale_x_continuous(breaks=seq(0, max(gep_gene_position$col), 4)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme_minimal() +
    labs(title = plot_title, x = plot_x_title, y = "Gene Rank") +
    facet_wrap(. ~ gene_name, nrow = length(unique(gep_gene_position$gene_name)))
  gg_ri <- ggplot(gep_gene_position, aes(x=col, y=row_inverse)) +
    geom_line() + scale_x_continuous(breaks=seq(0, max(gep_gene_position$col), 4)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme_minimal() +
    labs(title = plot_title, x = plot_x_title, y = "Inverse Gene Rank") +
    facet_wrap(. ~ gene_name, nrow = length(unique(gep_gene_position$gene_name)))
  gg_ril <- ggplot(gep_gene_position, aes(x=col, y=row_inverse_logged)) +
    geom_line() + scale_x_continuous(breaks=seq(0, max(gep_gene_position$col), 4)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme_minimal() +
    labs(title = plot_title, x = plot_x_title, y = "Natural Log Inverse of Gene Rank") +
    facet_wrap(. ~ gene_name, nrow = length(unique(gep_gene_position$gene_name)))

  pdf(file_name, height = file_height, width = file_width)
  plot(gg_r)
  plot(gg_ri)
  plot(gg_ril)
  dev.off()
}

gene_position_heatmap <- function(gep_gene_position_df, heatmap_col = "row_inverse",
                                  file_name = "", file_height = 10, file_width = 25,
                                  program_annotations = NULL){
  gep_gene_position_wide <- gep_gene_position_df
  gep_gene_position_wide <- gep_gene_position_wide %>%
    pivot_wider(names_from = col, values_from = heatmap_col, id_cols = gene_name)
  gep_gene_position_wide <- as.data.frame(gep_gene_position_wide)
  row.names(gep_gene_position_wide) <- gep_gene_position_wide$gene_name
  gep_gene_position_wide <- gep_gene_position_wide[,-1]
  gep_gene_position_wide <- as.matrix(gep_gene_position_wide)

  scaled_mat <- t(scale(t(gep_gene_position_wide)))

  if (is.null(program_annotations)){
    ht_clustered <- Heatmap(scaled_mat, cluster_columns = TRUE, cluster_rows = TRUE)
    ht_col_clustered <- Heatmap(scaled_mat, cluster_columns = TRUE, cluster_rows = FALSE)
    ht_row_clustered <- Heatmap(scaled_mat, cluster_columns = FALSE, cluster_rows = TRUE)
    ht_unclustered <- Heatmap(scaled_mat, cluster_columns = FALSE, cluster_rows = FALSE)

    pdf(file_name, height = file_height, width = file_width)
    draw(ht_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_col_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_row_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_unclustered, padding = unit(c(3, 3, 3, 3), "cm"))
    dev.off()
  }

  if (!is.null(program_annotations)){
    colnames(scaled_mat) <- program_annotations

    program_level_one_annotations <- str_to_title(program_annotation_list_label_level_one)
    label_colors <- paletteer_d("colorBlindness::paletteMartin")
    names(label_colors) <- unique(program_level_one_annotations)

    haT <- HeatmapAnnotation(Label = program_level_one_annotations,
                             col = list(Label = label_colors),
                             which = "col")

    ht_clustered <- Heatmap(scaled_mat, top_annotation = haT,
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            column_names_rot = 45,
                            name = "Scaled Inverse Gene Rank",
                            row_title = "Genes", column_title = "cNMF Programs")
    ht_col_clustered <- Heatmap(scaled_mat, top_annotation = haT,
                                cluster_columns = TRUE, cluster_rows = FALSE,
                                column_names_rot = 45,
                                name = "Scaled Inverse Gene Rank",
                                row_title = "Genes", column_title = "cNMF Programs")
    ht_row_clustered <- Heatmap(scaled_mat, top_annotation = haT,
                                cluster_columns = FALSE, cluster_rows = TRUE,
                                column_names_rot = 45,
                                name = "Scaled Inverse Gene Rank",
                                row_title = "Genes", column_title = "cNMF Programs")
    ht_unclustered <- Heatmap(scaled_mat, top_annotation = haT,
                              cluster_columns = FALSE, cluster_rows = FALSE,
                              column_names_rot = 45,
                              name = "Scaled Inverse Gene Rank",
                              row_title = "Genes", column_title = "cNMF Programs")

    pdf(file_name, height = file_height, width = file_width)
    draw(ht_clustered, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "Scaled Inverse Gene Rank Of Genes Involved in Immune Resistance Mechanisms")
    draw(ht_col_clustered, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "Scaled Inverse Gene Rank Of Genes Involved in Immune Resistance Mechanisms")
    draw(ht_row_clustered, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "Scaled Inverse Gene Rank Of Genes Involved in Immune Resistance Mechanisms")
    draw(ht_unclustered, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "Scaled Inverse Gene Rank Of Genes Involved in Immune Resistance Mechanisms")
    dev.off()
  }
}

gProfiler_results <- function(GEP_df, species = "hsapiens", source_filter = TRUE, term_filter = TRUE){
  gep_df_list <- as.list(GEP_df)

  # Run gProfiler
  gprofiler_output <- gost(gep_df_list, organism = species, ordered_query = TRUE, multi_query = TRUE)

  # Convert the pVals from lists to their own dataframe with columns as programs
  gprofiler_pVal <- t(as.data.frame(gprofiler_output$result$p_values))
  gprofiler_pVal <- as.data.frame(gprofiler_pVal)

  # Set the column names to the program name inputs, -log10 of the pVals
  colnames(gprofiler_pVal) <- names(gep_df_list)
  gprofiler_pVal <- -log10(gprofiler_pVal)
  gprofiler_pVal <- as.data.frame(gprofiler_pVal)

  # Merge the pathway information with the pVals
  gprofiler_pVal_info <- data.frame(source = gprofiler_output$result$source,
                                    term_name = gprofiler_output$result$term_name,
                                    term_size = gprofiler_output$result$term_size)
  gprofiler_merged_pVal_info <- cbind(gprofiler_pVal_info, gprofiler_pVal)

  # Filter the output data frame to only GO:BP
  if (source_filter){
    gprofiler_merged_pVal_info <- gprofiler_merged_pVal_info %>% filter(source == "GO:BP")
  }
  # Filter the output data frame to only pathways with > 10 and < 2000 terms
  if (term_filter){
    gprofiler_merged_pVal_info <- gprofiler_merged_pVal_info %>% filter(term_size > 10) %>% filter(term_size < 2000)
  }

  # Set the rownames to be the term names and make sure they're unique
  rownames(gprofiler_merged_pVal_info) <- make.unique(gprofiler_merged_pVal_info$term_name)
  gprofiler_merged_pVal_info$term_name <- NULL

  # For each term, identify which program has the greatest significance for it
  gProf_res_sub2 <- gprofiler_merged_pVal_info
  gProf_res_sub2 <- data.frame(gProf_res_sub2) %>% add_column("group" = "", .before = colnames(gProf_res_sub2)[1])
  for(i in 1:nrow(gProf_res_sub2)){
    gProf_res_sub2$group[i] <- names(which.max(gProf_res_sub2[i,][4:ncol(gProf_res_sub2)]))
  }

  # Create a vector/df to reorder the terms based on the most significant term
  groups <- unique(gProf_res_sub2$group)
  Terms_order <- data.frame()
  for(i in 1:length(groups)){
    gProf_res_sub3 <- gProf_res_sub2 %>% filter(group == groups[i])
    gProf_res_sub3 <- gProf_res_sub3[which(colnames(gProf_res_sub3) %in% c("group", groups[i]))]
    gProf_res_sub3 <- gProf_res_sub3[order(gProf_res_sub3[,2], decreasing = TRUE),]
    Terms_order <- rbind(Terms_order, data.frame(V1 = rownames(gProf_res_sub3), V2 = groups[i]))
  }

  # Reorder the terms in the dataframe using the vector/df
  gprofiler_merged_pVal_info_reordered <- gprofiler_merged_pVal_info[Terms_order$V1,]

  # Remove the term information
  gprofiler_merged_pVal_info_reordered <- as.data.frame(gprofiler_merged_pVal_info_reordered)
  gprofiler_merged_pVal_info_reordered_data <- gprofiler_merged_pVal_info_reordered[,-c(1,2)]
  gprofiler_merged_pVal_info_reordered_data <- as.matrix(gprofiler_merged_pVal_info_reordered_data)

  # Assign certain variables to the GlobalEnv, so they can be accessed for other functions later on
  assign("gProfiler_output", gprofiler_output, envir = .GlobalEnv)
  assign("gProfiler_merged_pVal_info", gprofiler_merged_pVal_info, envir = .GlobalEnv)
  assign("gProfiler_merged_pVal_info_reordered", gprofiler_merged_pVal_info_reordered, envir = .GlobalEnv)

  return(gprofiler_merged_pVal_info_reordered_data)
}

gProfiler_results_plots <- function(gProfiler_merged_pVal_info_reordered_df, pathway = "all", save_path){
  ht_opt$message = FALSE

  if (pathway != "all"){
    pathway_sub <- gsub(":", "", pathway)
    save_folder <- paste(save_path, pathway_sub, "/", sep = "")
    create_folder(save_folder)
  }
  else{
    pathway_sub <- "combinedPathways"
    save_folder <- save_path
  }

  base_file_name <- paste(save_folder, "K", RANK, "_gProfiler_", sep = "")
  output_file_name <- paste(base_file_name, pathway_sub, "_top1000.pdf", sep = "")
  capped_output_file_name <- paste(base_file_name, "capped_", pathway_sub, "_top1000.pdf", sep = "")
  capped_output_diff_file_name <- paste(base_file_name, "capped_diff_", pathway_sub, "_top1000.pdf", sep = "")

  # Retrieve the specific output ----
  if (pathway != "all"){
    output <- gProfiler_merged_pVal_info_reordered %>% filter(source == pathway)
  }
  else{
    output <- gProfiler_merged_pVal_info_reordered
  }
  output <- output[,-c(1,2)]

  capped_output <- output
  capped_output[capped_output > 10] <- 10

  # Remove the pathways in the output that are the same across all the programs
  capped_output_diff <- as.data.frame(capped_output)
  capped_output_diff$total <- rowSums(capped_output_diff)
  capped_output_diff <- capped_output_diff[capped_output_diff$total != max(capped_output_diff),]
  capped_output_diff <- capped_output_diff[,-length(capped_output_diff)]

  # Convert all inputs to a matrix
  final_outputs <- list()
  final_outputs[[1]] <- as.matrix(output)
  final_outputs[[2]] <- as.matrix(capped_output)
  final_outputs[[3]] <- as.matrix(capped_output_diff)

  # Plotting ----
  final_outputs_save_folders <- c(output_file_name, capped_output_file_name, capped_output_diff_file_name)
  col_fun <- brewer.pal(name = "Blues", n = 9)

  for (output_index in seq(1:length(final_outputs))){
    temp <- final_outputs[[output_index]]
    temp_save_name <- final_outputs_save_folders[output_index]

    ht_clustered <- Heatmap(temp, show_row_names = TRUE,
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            col = col_fun, height = unit(0.5, "cm")*nrow(temp))
    ht_col_clustered <- Heatmap(temp, show_row_names = TRUE,
                                cluster_columns = TRUE, cluster_rows = FALSE,
                                col = col_fun, height = unit(0.5, "cm")*nrow(temp))
    ht_row_clustered <- Heatmap(temp, show_row_names = TRUE,
                                cluster_columns = FALSE, cluster_rows = TRUE,
                                col = col_fun, height = unit(0.5, "cm")*nrow(temp))
    ht_unclustered <- Heatmap(temp, show_row_names = TRUE,
                              cluster_columns = FALSE, cluster_rows = FALSE,
                              col = col_fun, height = unit(0.5, "cm")*nrow(temp))
    pdf(temp_save_name, width = 25, height = (0.1969*nrow(capped_output)) + 1.3386)
    draw(ht_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_col_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_row_clustered, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_unclustered, padding = unit(c(3, 3, 3, 3), "cm"))
    dev.off()
  }
}

GEP_correlation_singlerank <- function(GS_dataset1, GS_dataset2){
  GS1 <- GS_dataset1
  GS2 <- GS_dataset2

  correlation_mat <- cor(GS1, GS2, method = "spearman")
  return(correlation_mat)
}

# Program Analysis Functions ----
program_enrichment <- function(usage_norm, threshold = 0.1, plot = "barplot", file_name, file_height = 6, file_width = 8, program_annotation = NULL){
  cutoff <- as.data.frame(usage_norm)
  cutoff$sample <- gsub(paste0("(", pattern_regex, ").*"), "\\1", rownames(cutoff))

  if (plot == "barplot"){
    # Get a df that contains the sample, # of total sample spots, # of the sample spots above the threshold per program (cpc), and a proportion.
    # The proportion is the cpc / # of total sample spots.
    p1_df <- cutoff %>%
      group_by(sample) %>%
      summarise(across(1:(ncol(cutoff) - 1), ~ sum(. > threshold), .names = "GEP{col}"),
                total_spots = n()) %>%
      pivot_longer(cols = starts_with("GEP"), names_to = "program", values_to = "cpc") %>%
      mutate(program = as.numeric(gsub("GEP", "", program)),
             proportion = cpc / total_spots) %>%
      ungroup()

    # Do the same as above but, this time, combine the samples by the TREATMENT column.
    p2_df <- merge(sample_metadata[, c("sample", "TREATMENT")], cutoff,
                   by.x = "sample", by.y = "sample", all.x = TRUE)
    p2_df <- p2_df[,-c(1)]; colnames(p2_df)[1] <- "sample"
    p2_df <- p2_df %>%
      group_by(sample) %>%
      summarise(across(1:(ncol(p2_df) - 1), ~ sum(. > threshold), .names = "GEP{col}"),
                total_spots = n()) %>%
      pivot_longer(cols = starts_with("GEP"), names_to = "program", values_to = "cpc") %>%
      mutate(program = as.numeric(gsub("GEP", "", program)),
             proportion = cpc / total_spots) %>%
      ungroup()

    # Do the same as above but, this time, combine the samples by the summarized response column.
    p3_df <- merge(sample_metadata[, c("sample", "LOCATION")], cutoff,
                   by.x = "sample", by.y = "sample", all.x = TRUE)
    p3_df <- p3_df[,-c(1)]; colnames(p3_df)[1] <- "sample"
    p3_df <- p3_df %>%
      group_by(sample) %>%
      summarise(across(1:(ncol(p3_df) - 1), ~ sum(. > threshold), .names = "GEP{col}"),
                total_spots = n()) %>%
      pivot_longer(cols = starts_with("GEP"), names_to = "program", values_to = "cpc") %>%
      mutate(program = as.numeric(gsub("GEP", "", program)),
             proportion = cpc / total_spots) %>%
      ungroup()

    # Plotting ----
    # If the programs have been annotated, then change the program names from GEP... to the annotated names.
    if (!is.null(program_annotation)){
      p1_df$program <- rep(program_annotations, length(unique(p1_df$sample)))
      p2_df$program <- rep(program_annotations, length(unique(p2_df$sample)))
      p3_df$program <- rep(program_annotations, length(unique(p3_df$sample)))
    } else{
      p1_df$program <- ifelse(p1_df$program <= 9, paste("GEP", p1_df$program, sep = "0"), paste("GEP", p1_df$program, sep = ""))
      p2_df$program <- ifelse(p2_df$program <= 9, paste("GEP", p2_df$program, sep = "0"), paste("GEP", p2_df$program, sep = ""))
      p3_df$program <- ifelse(p3_df$program <= 9, paste("GEP", p3_df$program, sep = "0"), paste("GEP", p3_df$program, sep = ""))
    }

    p1_df <- as.data.frame(p1_df); p2_df <- as.data.frame(p2_df); p3_df <- as.data.frame(p3_df)
    colnames(p1_df)[1] <- "Sample"; colnames(p2_df)[1] <- "Sample"; colnames(p3_df)[1] <- "Sample"

    # Plot options
    plot_title <- paste("K", ncol(usage_norm), "Program Enrichment", sep = " ")

    pdf(file = file_name, height = file_height, width = file_width)
    p1 <- ggplot(p1_df, aes(x = program, y = cpc, fill = Sample)) +
      geom_col(colour = "black", position = "fill") + theme_classic2() +
      labs(title = plot_title, x = "Gene Expression Program", y = "Proportion of Spots", color = "Sample",
           subtitle = paste("Threshold - ", threshold, sep = "")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour="black"),
            legend.position = "top", plot.margin = unit(c(1,10,1,1), "pt"), text = element_text(size = 20)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = sample_cols) + facet_grid(~ program, scales = "free")
    p2 <- ggplot(p1_df, aes(x = program, y = proportion, fill = Sample)) +
      geom_bar(stat="identity", color = "black", position = "dodge", width = .8) + theme_classic2() +
      labs(title = plot_title, x = "Gene Expression Program", y = "Proportion of Spots", color = "Sample",
           subtitle = paste("Threshold - ", threshold, sep = "")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour="black"),
            legend.position = "top", plot.margin = unit(c(1,10,1,1), "pt"), text = element_text(size = 20)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = sample_cols) + facet_wrap(~ program, nrow = 5, scales = "free")
    p3 <- ggplot(p2_df, aes(x = program, y = cpc, fill = Sample)) +
      geom_col(colour = "black", position = "fill") + theme_classic2() +
      labs(title = plot_title, x = "Gene Expression Program", y = "Proportion of Spots", color = "Sample",
           subtitle = paste("Threshold - ", threshold, sep = "")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour="black"),
            legend.position = "top", plot.margin = unit(c(1,10,1,1), "pt"), text = element_text(size = 20)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = unname(sample_cols[1:4]))
    p4 <- ggplot(p2_df, aes(x = program, y = proportion, fill = Sample)) +
      geom_bar(stat="identity", color = "black", position = "dodge", width = .8) + theme_classic2() +
      labs(title = plot_title, x = "Gene Expression Program", y = "Proportion of Spots", color = "Sample",
           subtitle = paste("Threshold - ", threshold, sep = "")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour="black"),
            legend.position = "top", plot.margin = unit(c(1,10,1,1), "pt"), text = element_text(size = 20)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = unname(sample_cols[1:4])) + facet_wrap(~ program, nrow = 5, scales = "free")
    p5 <- ggplot(p3_df, aes(x = program, y = cpc, fill = Sample)) +
      geom_col(colour = "black", position = "fill") + theme_classic2() +
      labs(title = plot_title, x = "Gene Expression Program", y = "Proportion of Spots", color = "Sample",
           subtitle = paste("Threshold - ", threshold, sep = "")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour="black"),
            legend.position = "top", plot.margin = unit(c(1,10,1,1), "pt"), text = element_text(size = 20)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = unname(sample_cols[1:2]))
    p6 <- ggplot(p3_df, aes(x = program, y = proportion, fill = Sample)) +
      geom_bar(stat="identity", color = "black", position = "dodge", width = .8) + theme_classic2() +
      labs(title = plot_title, x = "Gene Expression Program", y = "Proportion of Spots", color = "Sample",
           subtitle = paste("Threshold - ", threshold, sep = "")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour="black"),
            legend.position = "top", plot.margin = unit(c(1,10,1,1), "pt"), text = element_text(size = 20)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = unname(sample_cols[1:2])) + facet_wrap(~ program, nrow = 5, scales = "free")
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    print(p6)
    dev.off()

    # This plot is for a version of p2_df with the same y-axis
    temp_p2_df <- p2_df
    temp_p2_df$Sample <- factor(temp_p2_df$Sample, levels = c("Biopsy_1", "Biopsy_2"))
    pdf(file = sub(".pdf", "_flat.pdf", file_name), height = 6.816, width = 14)
    p7 <- ggplot(temp_p2_df, aes(x = program, y = proportion, fill = Sample)) +
      geom_bar(stat="identity", color = "black", position = "dodge", width = .8) + theme_classic2() +
      labs(title = plot_title, x = "Gene Expression Program", y = "Proportion of Spots", color = "Sample",
           subtitle = paste("Threshold - ", threshold, sep = "")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour="black"),
            legend.position = "top", theme(plot.margin = unit(c(1,1,1,1), "cm")), text = element_text(size = 20)) +
      scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = unname(sample_cols[1:4]))
    print(p7)
    dev.off()
  }
}

# Neighbourhood Analysis Functions ----
seurat_object_coords <- function(){
  union_spots <- data.frame()

  for (sample in so_list){
    sample_so <- eval(as.name(sample))

    # Get spot centroids from each sample
    sample_spots <- as.data.frame(sample_so@images$slice1.024um@boundaries$centroids@coords)

    # Set the rownames to the cell IDs
    rownames(sample_spots) <- sample_so@images$slice1.024um@boundaries$centroids@cells

    # Add a prefix to the rownames, so they match what's in the merged seurat object
    rownames(sample_spots) <- paste(sample, rownames(sample_spots), sep = "_")

    union_spots <- rbind(union_spots, sample_spots)
  }

  # Take the union of the spots, rename the columns, and make a column to specify the sample
  colnames(union_spots) <- c("row", "col")
  union_spots$sample <- gsub(paste0("(", pattern_regex, ").*"), "\\1", rownames(union_spots))

  return(union_spots)
}

program_neighbours <- function(program_spots, cell_centroids, radius = 75, max_radius = 3){
  program_spots <- program_spots
  distances <- seq(0, max_radius) * radius

  # Convert cell_centroids to a matrix for faster indexing
  # Store sample column separately for indexing
  cell_centroids_matrix <- as.matrix(cell_centroids[, c("col", "row")])
  cell_centroid_samples <- cell_centroids$sample

  # Pre-compute bounds for each radius level to avoid repetitive calculations
  bounds <- lapply(distances, function(k) list(lower = -k, upper = k))

  # Parallelize across `program_spots` if possible
  # Use available cores minus one
  cl <- makeCluster(parallelly::availableCores() - 1)
  clusterExport(cl, c("program_spots", "cell_centroids_matrix", "cell_centroid_samples",
                      "bounds", "distances", "radius", "max_radius", "pattern_regex"),
                envir=environment())
  clusterEvalQ(cl, library(data.table))

  final_results <- parLapply(cl, seq_along(program_spots), function(i) {
    spot <- cell_centroids_matrix[program_spots[i], ]
    sample <- gsub(paste0("(", pattern_regex, ").*"), "\\1", rownames(cell_centroids)[i])

    union_spots_subset <- cell_centroids_matrix[cell_centroid_samples == sample, , drop = FALSE]
    subset_indices <- which(cell_centroid_samples == sample)

    # Store intermediate results for each radius
    spot_results <- vector("list", length(bounds))

    for (j in seq_along(bounds)) {
      k_bounds <- bounds[[j]]

      # Apply bounds to subset neighbors
      mask_x <- (union_spots_subset[, "col"] >= (spot["col"] + k_bounds$lower)) &
        (union_spots_subset[, "col"] <= (spot["col"] + k_bounds$upper))
      mask_y <- (union_spots_subset[, "row"] >= (spot["row"] + k_bounds$lower)) &
        (union_spots_subset[, "row"] <= (spot["row"] + k_bounds$upper))

      neighbors <- subset_indices[mask_x & mask_y]

      # Record result
      spot_results[[j]] <- data.frame(
        sample = sample,
        spot_of_interest = program_spots[i],
        r = distances[j],
        type = "circular",
        neighbours = I(list(rownames(cell_centroids)[neighbors])),
        nr_neighbours = length(neighbors)
      )
    }
    rbindlist(spot_results)
  })

  stopCluster(cl)  # Stop the parallel cluster

  # Combine all parallelized results
  final_neigh_df <- rbindlist(final_results)

  return(final_neigh_df)
}

remove_overlapping_neighbours <- function(spot_neighbours_df){
  # Convert to data.table for efficient updates
  spot_neighbours_dt <- as.data.table(spot_neighbours_df)

  # Ensure radii are unique and sorted
  radii <- sort(unique(as.numeric(spot_neighbours_dt$r)))

  # Set up parallel cluster backend
  num_cores <- parallelly::availableCores() - 1  # Leave one core free
  plan(multicore, workers = num_cores)

  # Split unique spots into chunks for parallel processing
  unique_spots <- unique(spot_neighbours_dt$spot_of_interest)
  spot_chunks <- split(unique_spots, cut(seq_along(unique_spots), num_cores, labels = FALSE))

  # Define function to process each chunk
  process_chunk <- function(spots_chunk) {
    # Initialize a list to store neighbors encountered for each spot of interest
    neighbors_encountered <- vector("list", length(spots_chunk))
    names(neighbors_encountered) <- spots_chunk

    # Process each spot of interest in the chunk
    for (soi in spots_chunk) {
      # Initialize encountered neighbors
      neighbors_encountered[[soi]] <- character(0)

      # Loop over each radius
      for (rad in radii) {
        # Logical indexing to select rows by spot of interest and radius
        soi_rows <- spot_neighbours_dt[spot_of_interest == soi & r == rad]

        # Apply unique and setdiff to the neighbors
        current_neighbors <- unique(unlist(soi_rows$neighbours))
        unique_neighbors <- setdiff(current_neighbors, neighbors_encountered[[soi]])

        # Update encountered neighbors list
        neighbors_encountered[[soi]] <- c(neighbors_encountered[[soi]], unique_neighbors)

        # Update the data.table with unique neighbors and count
        spot_neighbours_dt[spot_of_interest == soi & r == rad, `:=`(
          neighbours = list(unique_neighbors),
          nr_neighbours = length(unique_neighbors)
        )]

        # Debug statement: Check if unique_neighbors is being correctly assigned
        if (length(unique_neighbors) == 0) {
          message("Empty unique_neighbors for spot_of_interest = ", soi, ", radius = ", rad)
        }
      }
    }
    return(spot_neighbours_dt[spot_of_interest %in% spots_chunk, ])
  }

  # Apply processing in parallel using future_lapply with cluster strategy
  results <- future_lapply(spot_chunks, process_chunk)

  # Combine results into a single data.table
  final_result <- rbindlist(results)

  # Return final result as a data.frame if needed
  return(as.data.frame(final_result))
}

# 1. MGS ----
# 1.1 Make the gene rank graphs ----
if (GENE_RANK_PLOTTING){
  # 1.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "annotation/K", RANK, "/genes_ranked/", sep = "")
  create_folder(save_path)
  create_folder(paste(save_path, "T_Cell_Subunits/", sep = ""))
  create_folder(paste(save_path, "PD1_PDL1/", sep = ""))
  create_folder(paste(save_path, "resistance_genes/", sep = ""))
  create_folder(paste(save_path, "spectral_flow/", sep = ""))
  create_folder(paste(save_path, "fas/", sep = ""))
  create_folder(paste(save_path, "dnam/", sep = ""))
  create_folder(paste(save_path, "kevin/", sep = ""))
  create_folder(paste(save_path, "myeloid_lineage/", sep = ""))
  create_folder(paste(save_path, "myeloid_lineage_kev/", sep = ""))
  create_folder(paste(save_path, "franz/", sep = ""))
  create_folder(paste(save_path, "top_per_program/", sep = ""))

  # 1.1 Calculate all the gene ranks per program from the spectra ----
  spectra_df <- load_spectra(SPECTRA)
  gep_df <- ranked_spectra(spectra_df)

  indices <- data.frame("PD1" = find_element_indices(gep_df, "PDCD1"),
                        "PDL1" = find_element_indices(gep_df, "CD274"),
                        "CTLA4" = find_element_indices(gep_df, "CTLA4"),
                        "CD80" = find_element_indices(gep_df, "CD80"),
                        "CD86" = find_element_indices(gep_df, "CD86"),
                        "TIGIT" = find_element_indices(gep_df, "TIGIT"),
                        "PVR" = find_element_indices(gep_df, "PVR"),
                        "NECTIN2" = find_element_indices(gep_df, "NECTIN2"),
                        "KLRB1" = find_element_indices(gep_df, "KLRB1"),
                        "CLEC2D" = find_element_indices(gep_df, "CLEC2D"))

  # 1.2 Retrieve specific gene positions ----
  gene_position_t_cell <- gene_position(gep_df, c("CD3D", "CD3E", "CD3G", "CD8A", "CD4"))
  gene_position_resistance <- gene_position(gep_df, c("PDCD1", "CD274"))
  gene_position_resistance_extened <- gene_position(gep_df, c("PDCD1", "CD274",
                                                              "CTLA4", "CD80","CD86",
                                                              "TIGIT", "PVR", "NECTIN2", "CD226", "CD96",
                                                              "KLRB1", "CLEC2D",
                                                              "FAS", "FASLG",
                                                              "HAVCR2", "LGALS9", "CEACAM1",
                                                              "LAG3", "FGL1", "LGALS3"))
  gene_position_spectral_flow <- gene_position(gep_df, c("PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "FOXP3", "CD19", "KLRB1",
                                                         "GZMB", "CD38", "CD69", "ITGAL", "SELL", "CXCR5", "ITGA4", "PDCD1", "MKI67",
                                                         "ITGAM", "CX3CR1", "LY6G6C", "ITGAX", "MPO", "CD40",
                                                         "TEK", "CCR2"))
  gene_position_FAS <- gene_position(gep_df, c("FASLG"))
  gene_position_DNAM1 <- gene_position(gep_df, c("CD226", "TIGIT"))
  gene_position_kevin <- gene_position(gep_df, c("DLL3"))
  gene_position_myeloid <- gene_position(gep_df, c("CD14", "FCGR3A",
                                                   "CD68", "CSF1R", "ITGAM", "FCGR1A", "CX3CR1",
                                                   "CD1A", "THBD", "CD1C",
                                                   "ANPEP", "CCR3", "CD33", "IL5RA",
                                                   "MRC1", "CD163"))
  gene_position_myeloid_kev <- gene_position(gep_df, c("CD68", "CD163", "MRC1", "CD274"))
  gene_position_franz <- gene_position(gep_df, c("TIGIT", "PVR", "NECTIN2", "CD226", "CD96"))
  gene_position_tpp <- gene_position(gep_df, unique(unlist(gep_df[c(1:10),])))

  # 1.3 Save the ranked genes per gep data frame ----
  write.table(gep_df, paste(save_path, "K", RANK, "_ranked_genes.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)

  # 1.4 Plot the gene ranks ----
  gene_position_linegraph(gene_position_t_cell, plot_title = paste("T-Cell Subunits (K", RANK, ")", sep = ""),
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "T_Cell_Subunits/K", RANK, "_CD3_subunits_spike_graph.pdf", sep = ""),
                          file_height = 10, file_width = 25)
  gene_position_heatmap(gene_position_t_cell, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "T_Cell_Subunits/K", RANK, "_CD3_subunits_heatmap.pdf", sep = ""),
                        file_height = 10, file_width = 25)

  # For CD274 and PDCD1
  gene_position_linegraph(gene_position_resistance, plot_title = paste("K", RANK, " - PDCD1 and CD274 Spike Graph", sep = ""),
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "PD1_PDL1/K", RANK, "_PDCD1_CD274_spike_graph.pdf", sep = ""),
                          file_height = 10, file_width = 25)
  gene_position_heatmap(gene_position_resistance, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "PD1_PDL1/K", RANK, "_PDCD1_CD274_heatmap.pdf", sep = ""),
                        file_height = 10, file_width = 25)

  # For all the resistance genes
  write.table(indices, file = paste(save_path, "/resistance_genes/resistance_gene_ranks.tsv", sep = ""), quote = FALSE, sep = "\t")
  gene_position_linegraph(gene_position_resistance_extened, plot_title = paste("K", RANK, " - Resistance Genes Spike Graph", sep = ""),
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "resistance_genes/K", RANK, "_resistance_genes_spike_graph.pdf", sep = ""),
                          file_height = 10, file_width = 25)
  gene_position_heatmap(gene_position_resistance_extened, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "resistance_genes/K", RANK, "_resistance_genes_heatmap.pdf", sep = ""),
                        file_height = 10, file_width = 25)

  # For the spectral flow panel
  gene_position_linegraph(gene_position_spectral_flow, plot_title = paste("K", RANK, " - Spectral Flow Spike Graph", sep = ""),
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "spectral_flow/K", RANK, "_spectral_flow_genes_spike_graph.pdf", sep = ""),
                          file_height = 30, file_width = 25)
  gene_position_heatmap(gene_position_spectral_flow, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "spectral_flow/K", RANK, "_spectral_flow_genes_heatmap.pdf", sep = ""),
                        file_height = 30, file_width = 25)

  # For the FASL position
  gene_position_linegraph(gene_position_FAS, plot_title = paste("FASL (K", RANK, ")", sep = ""),
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "fas/K", RANK, "_fas_spike_graph.pdf", sep = ""),
                          file_height = 10, file_width = 25)
  gene_position_heatmap(gene_position_FAS, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "fas/K", RANK, "_fas_heatmap.pdf", sep = ""),
                        file_height = 10, file_width = 25)

  # For the DNAM position
  gene_position_linegraph(gene_position_DNAM1, plot_title = paste("DNAM (K", RANK, ")", sep = ""),
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "dnam/K", RANK, "_dnam_spike_graph.pdf", sep = ""),
                          file_height = 10, file_width = 25)
  gene_position_heatmap(gene_position_DNAM1, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "dnam/K", RANK, "_dnam_heatmap.pdf", sep = ""),
                        file_height = 10, file_width = 25)

  # For Kevins genes of interest
  gene_position_linegraph(gene_position_kevin, plot_title = paste("Kevin Genes of Interest (K", RANK, ")", sep = ""),
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "kevin/K", RANK, "_kevin_spike_graph.pdf", sep = ""),
                          file_height = 10, file_width = 25)
  gene_position_heatmap(gene_position_kevin, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "kevin/K", RANK, "_kevin_heatmap.pdf", sep = ""),
                        file_height = 10, file_width = 25)

  # For Myeloid genes of interest
  gene_position_linegraph(gene_position_myeloid, plot_title = paste("Myeloid Genes of Interest (K", RANK, ")", sep = ""),
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "myeloid_lineage/K", RANK, "_myeloid_spike_graph.pdf", sep = ""),
                          file_height = 10, file_width = 25)
  gene_position_heatmap(gene_position_myeloid, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "myeloid_lineage/K", RANK, "_myeloid_heatmap.pdf", sep = ""),
                        file_height = 10, file_width = 25)

  # For Myeloid genes of interest
  gene_position_linegraph(gene_position_myeloid_kev, plot_title = paste("Myeloid Genes of Interest (K", RANK, ")", sep = ""),
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "myeloid_lineage_kev/K", RANK, "_myeloid_spike_graph.pdf", sep = ""),
                          file_height = 10, file_width = 25)
  gene_position_heatmap(gene_position_myeloid_kev, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "myeloid_lineage_kev/K", RANK, "_myeloid_heatmap.pdf", sep = ""),
                        file_height = 10, file_width = 25)

  # For Franz genes of interest
  gene_position_linegraph(gene_position_franz, plot_title = paste("Franz Genes of Interest (K", RANK, ")", sep = ""),
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "franz/K", RANK, "_franz_spike_graph.pdf", sep = ""),
                          file_height = 10, file_width = 25)
  gene_position_heatmap(gene_position_franz, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "franz/K", RANK, "_franz_heatmap.pdf", sep = ""),
                        file_height = 10, file_width = 25)

  # For the top 10 genes per program
  gene_position_linegraph(gene_position_tpp, plot_title = paste("Top Genes Per Program at (K", RANK, ")", sep = ""),
                          plot_x_title = "Gene Expression Program",
                          file_name = paste(save_path, "top_per_program/K", RANK, "_t10pp_spike_graph.pdf", sep = ""),
                          file_height = 500, file_width = 50)
  gene_position_heatmap(gene_position_tpp, heatmap_col = "row_inverse",
                        file_name = paste(save_path, "top_per_program/K", RANK, "_t10pp_heatmap.pdf", sep = ""),
                        file_height = 125, file_width = 25)
}

# 1.2 Compute MGS for the top 1000 genes ----
if (COMPUTE_MGS){
  # 2.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "annotation/K", RANK, "/MarkerGeneScore/", sep = "")
  create_folder(save_path)

  # 2.1 Perform GMT-specific MGS ----
  # Set up parallel backend
  cl <- makeCluster(parallelly::availableCores() - 1 )
  registerDoParallel(cl)

  # List of files to process
  gmt_files <- list.files(GMT_PATH)

  # Parallel processing
  foreach(file = gmt_files, .packages = c("stringr", "ActivePathways",
                                          "tibble", "ComplexHeatmap",
                                          "circlize")) %dopar% {
                                            gmt_alias <- noquote(str_split_1(file, ".gmt")[1])
                                            create_folder(paste(save_path, gmt_alias, sep = ""))

                                            # 2.1.1 Calculate all the gene ranks per program from the spectra ----
                                            spectra_df <- load_spectra(SPECTRA)
                                            gep_df <- ranked_spectra(spectra_df, num_genes = 1000)

                                            # 2.1.2 Load the GMT file ----
                                            geneset <- read.GMT(paste0("/bulk/morrissy_bulk/GMT/", file))

                                            # 2.1.3 Compute the regular and detailed marker gene scores ----
                                            final_df_simple <- marker_gene_score(gep_df, geneset, scale = FALSE,
                                                                                 detailed = FALSE, filter = TRUE,
                                                                                 species = "Human", rank = RANK)
                                            final_df_detailed <- marker_gene_score(gep_df, geneset, scale = FALSE,
                                                                                   detailed = TRUE, filter = TRUE,
                                                                                   species = "Human", rank = RANK)

                                            # 2.1.4 Save the output tables ----
                                            write.table(final_df_simple, paste0(save_path, gmt_alias, "/K", RANK, "_MGS_top1000_genes.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
                                            write.table(final_df_detailed, paste0(save_path, gmt_alias, "/K", RANK, "_MGS_top1000_genes_detailed.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

                                            # Plot a heatmap of the GMT gene scores. This should be ran with the top 1000 genes.
                                            gmt_heatmap(final_df_simple,
                                                        file_name = paste0(save_path, gmt_alias, "/K", RANK, "_MGS_top1000_genes.pdf"),
                                                        file_height = round(length(geneset) / 4), file_width = 2*ncol(final_df_simple))
                                          }

  # Stop the parallel cluster
  stopCluster(cl)

  # 2.2 Merge the detailed GMT results ----
  merged_detailed_results <- data.frame()
  for (file in list.files(save_path)) {
    gmt_alias <- noquote(file)
    gmt_alias_folder <- list.files(paste0(save_path, gmt_alias))

    # Make sure the detailed output file exists
    if (length(list.files(paste0(save_path, gmt_alias))[grepl("detailed", list.files(paste0(save_path, gmt_alias)))]) > 0){
      gmt_alias_detailed <- list.files(paste0(save_path, gmt_alias))[grepl("detailed", list.files(paste0(save_path, gmt_alias)))]

      gmt_detailed_file <- fread(paste0(save_path, gmt_alias, "/", gmt_alias_detailed))

      if (ncol(merged_detailed_results) == 0) {
        merged_detailed_results <- gmt_detailed_file
      } else{
        merged_detailed_results <- rbind(merged_detailed_results, gmt_detailed_file)
      }
    }
  }
  write.table(merged_detailed_results, paste0(save_path, "/K", RANK, "_MGS_top1000_genes_detailed_merged.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}

# 1.3 Plot the top 10 genes per program based on their spectra score ----
if (TOP_GENES_SPECTRA){
  # 3.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "annotation/K", RANK, "/genes_ranked/", sep = "")
  save_path <- paste0(save_path, "custom/")
  create_folder(save_path)

  # 3.1 Calculate all the gene ranks per program from the spectra ----
  spectra_df <- load_spectra(SPECTRA)
  gep_df <- ranked_spectra(spectra_df, num_genes = 10, program_annotation = PROGRAM_ANNOTATION)

  # 3.2 Get the top 10 genes per program ranked based on the spectra ----
  top_genes <- gep_df[c(1:10),]
  top_genes <- top_genes %>%
    pivot_longer(everything(), names_to = "gep", values_to = "gene") %>%
    mutate(program = as.numeric(gsub("GEP_", "", gep))) %>%
    arrange(program)
  top_gene_scores <- t(spectra_df[c(top_genes$gene),]); top_gene_scores <- scale(top_gene_scores)

  top_gene_scores_annot <- as.numeric(rownames(top_gene_scores))
  haT <- HeatmapAnnotation(Program = top_genes$program,
                           col = list(Program = colorRamp2(range(1:length(top_gene_scores_annot)),
                                                           hcl_palette = "Inferno", reverse = TRUE)),
                           which = "column")
  haR <- HeatmapAnnotation(Program = top_gene_scores_annot,
                           col = list(Program = colorRamp2(range(top_gene_scores_annot),
                                                           hcl_palette = "Blue-Red", reverse = TRUE)),
                           which = "row")

  col_fun <- colorRamp2(breaks = range(top_gene_scores), hcl_palette = "RdYlBu", reverse = TRUE)
  ht <- Heatmap(top_gene_scores,
                right_annotation = haR, top_annotation = haT,
                cluster_columns = FALSE, cluster_rows = FALSE,
                col = col_fun, column_names_rot = 45,
                name = "Scaled Spectra Scores for \nTop 10 Scoring Genes Per Program",
                row_title = "Program", column_title = "Genes")
  ht_row_clust <- Heatmap(top_gene_scores,
                          right_annotation = haR, top_annotation = haT,
                          cluster_columns = FALSE, cluster_rows = TRUE,
                          col = col_fun, column_names_rot = 45,
                          name = "Scaled Spectra Scores for \nTop 10 Scoring Genes Per Program",
                          row_title = "Program", column_title = "Genes")
  ht_col_clust <- Heatmap(top_gene_scores,
                          right_annotation = haR, top_annotation = haT,
                          cluster_columns = TRUE, cluster_rows = FALSE,
                          col = col_fun, column_names_rot = 45,
                          name = "Scaled Spectra Scores for \nTop 10 Scoring Genes Per Program",
                          row_title = "Program", column_title = "Genes")
  ht_clust <- Heatmap(top_gene_scores,
                      right_annotation = haR, top_annotation = haT,
                      cluster_columns = TRUE, cluster_rows = TRUE,
                      col = col_fun, column_names_rot = 45,
                      name = "Scaled Spectra Scores for \nTop 10 Scoring Genes Per Program",
                      row_title = "Program", column_title = "Genes")
  pdf(paste(save_path, "K", RANK, "_Top_10_Scoring_Genes.pdf", sep = ""), width = (0.5*ncol(top_gene_scores)) + 1.3386, height = 20)
  draw(ht, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "Scaled Spectra Scores for Top 10 Scoring Genes Per Program")
  draw(ht_row_clust, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "Scaled Spectra Scores for Top 10 Scoring Genes Per Program")
  draw(ht_col_clust, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "Scaled Spectra Scores for Top 10 Scoring Genes Per Program")
  draw(ht_clust, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "Scaled Spectra Scores for Top 10 Scoring Genes Per Program")
  dev.off()
}

# 1.4 Generate gProfiler results ----
if (gPROFILER_RESULTS){
  # 4.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "annotation/K", RANK, "/gProfiler/", sep = "")
  create_folder(save_path)

  # 4.1 Calculate all the gene ranks per program from the spectra ----
  spectra_df <- load_spectra(SPECTRA)
  gep_df <- ranked_spectra(spectra_df, num_genes = 1000, program_annotation = PROGRAM_ANNOTATION)

  # 4.2 Run gProfiler on the gep_df ----
  output <- gProfiler_results(gep_df, source_filter = FALSE)

  # 4.3 Create a gost plot of the results ----
  pdf(paste(save_path, "K", RANK, "_gostplot_top1000.pdf", sep = ""), height = 50, width = 7)
  gostplot(gProfiler_output, capped = TRUE, interactive = FALSE)
  dev.off()

  # 4.4 Create the pathways plot of the results ----
  pathways <- unique(gProfiler_merged_pVal_info_reordered$source)
  pathways <- c("all", pathways)
  for (pathway in pathways){
    print(pathway)
    gProfiler_results_plots(gProfiler_merged_pVal_info_reordered,
                            pathway = pathway,
                            save_path)
  }

  # 4.5 Make the selected ranks annotation file ----
  sra_file <- data.frame(Rank = RANK, GEP = seq(1:RANK),
                         Label = "", Annotation_notes_SM = "",
                         Annotation = "",
                         Comment = "", GO_BP = "")

  sra_file$GO_BP <- apply(output, 2, function(col) {
    top_500 <- order(col, decreasing = TRUE)[1:500]
    top_500_pathways <- rownames(output)[top_500]
    paste(top_500_pathways, collapse = " | ")
  })
  write.table(sra_file, file = paste(SAVE_FOLDER, "annotation/K", RANK, "/SelectedRankAnnotation.txt", sep = ""),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # 4.6 Write the whole matrix of the gProfiler results ----
  write.table(gProfiler_merged_pVal_info_reordered, file = paste0(save_path, "gProfiler_results.csv"), sep = ",")

  save(list = ls(environment()), file = paste0(save_path, "gprofiler_object_save.RData"))
}

# Load the Seurat Objects ----
rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))

# 1. Classical loading (This is how it was loaded on an HPC, but for reviewer ease, use the RDS loading)
# name_prefix_select <- 1
# for (sample in sample_metadata$sample){
#   sample_path <- paste0("/work/morrissy_lab/gurveer/gpnmb_car/datasets/", sample, "/outs/")
#   sample_so <- Load10X_Spatial(data.dir = sample_path, bin.size = c(24))

#   sample_so@meta.data$sample <- sample
#   assign(sample, sample_so, envir = .GlobalEnv)

#   name_prefix_select <- name_prefix_select + 1
# }
# rm(sample_so, sample_path)
# so_list <- sort(ls()[grepl("SPSC", ls())])

# 2. Reviewer loading
SPSC_biopsy1 <- readRDS("02_Visium_HD/03_datasets/01_seurat_objects/ASPS-02/biopsy_sample_1_seurat_object.rds")
SPSC_biopsy2 <- readRDS("02_Visium_HD/03_datasets/01_seurat_objects/ASPS-02/biopsy_sample_2_seurat_object.rds")

# Update the objects because of Seurat backwards compatability issues
SPSC_biopsy1 <- UpdateSeuratObject(SPSC_biopsy1)
SPSC_biopsy2 <- UpdateSeuratObject(SPSC_biopsy2)

# 2.2 Add the sample name to the metadata table
SPSC_biopsy1@meta.data["Sample"] <- "SPSC_biopsy1"
SPSC_biopsy2@meta.data["Sample"] <- "SPSC_biopsy2"

so_list <- sort(ls()[grepl("SPSC", ls())])

# 3. Merge the samples
combined <- merge(SPSC_biopsy1, y = c(SPSC_biopsy2),
                  add.cell.ids = c(so_list),
                  project = "combined_samples")

combined_layers <- JoinLayers(combined)
combined_layers <- subset(combined_layers, nCount_Spatial.024um > 0)
combined_layers <- subset(combined_layers, cells = rownames(usage_norm))

program_var_names <- c(ls(), "program_var_names")

# 2. Program Analysis ----
# 2.1 Program Enrichment ----
if (PROGRAM_ENRICHMENT){
  # 1.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "K", RANK, "/program_enrichment/", sep = "")
  create_folder(save_path)

  # 1.1 Basic program enrichment across samples ----
  program_enrichment(usage_norm, plot = "barplot", threshold = GLOBAL_THRESHOLD,
                     file_name = paste(save_path, "/K", RANK, "_Program_Enrichment.pdf", sep = ""),
                     file_height = 15, file_width = 20 + (RANK / 5),
                     program_annotation = program_annotation)

}

# 2.2 Panel Gene Expression  ----
if (PROGRAM_GENE_EXPRESSION){
  # 2.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "K", RANK, "/program_gene_expression/", sep = "")
  franz_panel_path <- paste(save_path, "/franz_panel/", sep = "")
  singh_panel_path <- paste(save_path, "/singh_panel/", sep = "")
  resistance_panel_path <- paste(save_path, "/resistance_panel/", sep = "")
  fas_panel_path <- paste(save_path, "/fas_panel/", sep = "")
  resistance_we_gpnmb_panel_path <- paste(save_path, "/resistance_panel/wo_gpnmb/", sep = "")
  kevin_panel_path <- paste(save_path, "/kevin_panel/", sep = "")
  franz_resistance_panel_path <- paste(save_path, "/franz_resistance_panel/", sep = "")
  gpnmb_only_path <- paste(save_path, "/gpnmb_only/", sep = "")
  create_folder(save_path)
  create_folder(franz_panel_path)
  create_folder(singh_panel_path)
  create_folder(resistance_panel_path)
  create_folder(fas_panel_path)
  create_folder(resistance_we_gpnmb_panel_path)
  create_folder(kevin_panel_path)
  create_folder(franz_resistance_panel_path)
  create_folder(gpnmb_only_path)

  # 2.1 Make a plot showcasing gene expression across programs for Franz's flow panel----
  flow_gene_panel <- data.frame(category = c("some monocyte and some NK", "some monocyte and some NK",
                                             "NK and NKT", "all T cells but Effector Memory",
                                             "Treg - yes but also way better markers than this (Tox, FoxP3, Helois, etc)",
                                             "Treg and recently activated T cells - yes but also way better markers than this",
                                             "Naïve or Central Memory", "All lymphocytes", "B", "B",
                                             "some monocytes", "B and all T cells but less so on Effector Memory",
                                             "T", "T", "T", "T", "T", "T",
                                             "CAR", "CD4/CD8 T", "CD4/CD8 T", "CD4/CD8 T",
                                             "CD4/CD8 T", "All T cells except Naïve", "Naïve or Central Memory"),
                                gene = c("FCGR3A", "FCGR3B", "NCAM1", "IL7R", "IL2RA", "ISG20",
                                         "CCR7", "PTPRC", "CD19", "IGHD", "CD14", "CD27",
                                         "CD8A", "CD8B", "CD4", "CD3D", "CD3E", "CD3G",
                                         "MYC", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "FAS",
                                         "SELL"))
  gene_panel <- flow_gene_panel
  genes_in_data <- rownames(combined_layers@assays$Spatial.024um$counts)[rownames(combined_layers@assays$Spatial.024um$counts) %in% gene_panel$gene]
  combined_count_table <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[genes_in_data ,])))

  combined_count_table_subset <- combined_count_table[, colnames(combined_count_table) %in% gene_panel$gene]
  combined_count_table_subset <- combined_count_table_subset[rownames(combined_count_table_subset) %in% rownames(usage_norm),]

  corr_matrix <- cor(combined_count_table_subset, usage_norm)
  corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset, usage_norm)))
  corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

  gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
  corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
  corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

  split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                  levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

  ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_split = split, row_title_rot = 0)
  ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = TRUE,
                    row_split = split, row_title_rot = 0)
  ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = TRUE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0)
  ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                     cluster_columns = TRUE, cluster_rows = TRUE,
                     row_split = split, row_title_rot = 0)
  ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       row_split = split, row_title_rot = 0)
  ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = TRUE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = TRUE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            row_split = split, row_title_rot = 0)

  p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
  p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
  p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
  p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

  pdf(paste(save_path, "/franz_panel/", "K", RANK, "_gene_expression_program_usage_correlation.pdf", sep = ""), width = 50, height = 10)
  print(plot_grid(p1, p2))
  print(plot_grid(p3, p4))
  print(plot_grid(p5, p6))
  print(plot_grid(p7, p8))
  dev.off()

  # 2.2 Make a plot showcasing gene expression across programs for Singh's panel----
  singh_gene_panel <- data.frame(category = c("Endothelial Cells", rep("Glioma/Tumor Cells", 3),
                                              rep("T Cells and Subsets", 6),
                                              rep("Monocytes/Macrophages", 5), "Microglia",
                                              rep("Pro-Inflammatory and Cytotoxic Factors", 3),
                                              rep("Myeloid Cells/Antigen Presentation", 11),
                                              "Cell Death Markers"),
                                 gene = c("CD31", "GFAP", "SOX2", "GPNMB", "CD4", "CD8A",
                                          "CD8B", "FOXP3", "PTPRC", "GZMB", "CD68", "CD163",
                                          "ITGAX", "MARCO", "MRC1", "P2RY12", "IFNG", "IFNGR1",
                                          "IFNGR2", "HLA-DRA", "HLA-DRB1", "HLA-DRB2", "HLA-DRB3", "HLA-DRB4",
                                          "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8", "HLA-DRB9", "CD74", "CASP3"))
  gene_panel <- singh_gene_panel
  genes_in_data <- rownames(combined_layers@assays$Spatial.024um$counts)[rownames(combined_layers@assays$Spatial.024um$counts) %in% gene_panel$gene]
  combined_count_table <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[genes_in_data ,])))

  combined_count_table_subset <- combined_count_table[, colnames(combined_count_table) %in% gene_panel$gene]
  combined_count_table_subset <- combined_count_table_subset[rownames(combined_count_table_subset) %in% rownames(usage_norm),]
  combined_count_table_subset <- combined_count_table_subset[,colSums(combined_count_table_subset) != 0]

  corr_matrix <- cor(combined_count_table_subset, usage_norm)
  corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset, usage_norm)))
  corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

  gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
  corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
  corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

  split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                  levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

  ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_split = split, row_title_rot = 0)
  ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = TRUE,
                    row_split = split, row_title_rot = 0)
  ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = TRUE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0)
  ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                     cluster_columns = TRUE, cluster_rows = TRUE,
                     row_split = split, row_title_rot = 0)
  ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       row_split = split, row_title_rot = 0)
  ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = TRUE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = TRUE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            row_split = split, row_title_rot = 0)

  p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
  p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
  p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
  p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

  pdf(paste(save_path, "/singh_panel/" , "K", RANK, "_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
  print(plot_grid(p1, p2))
  print(plot_grid(p3, p4))
  print(plot_grid(p5, p6))
  print(plot_grid(p7, p8))
  dev.off()

  # 2.3 Make a plot showcasing gene expression across programs for resistance genes ----
  resistance_gene_panel <- data.frame(category = c("ASPS",
                                                   rep("Pair_1", 2),
                                                   rep("Pair_2", 3),
                                                   rep("Pair_3", 5),
                                                   rep("Pair_4", 2),
                                                   rep("Pair_5", 2),
                                                   rep("Pair_6", 3),
                                                   rep("Pair_7", 3)),
                                      gene = c("GPNMB",
                                               "PDCD1", "CD274",
                                               "CTLA4", "CD80", "CD86",
                                               "TIGIT", "PVR", "NECTIN2", "CD226", "CD96",
                                               "KLRB1", "CLEC2D",
                                               "FAS", "FASLG",
                                               "HAVCR2", "LGALS9", "CEACAM1",
                                               "LAG3", "FGL1", "LGALS3"))
  gene_panel <- resistance_gene_panel
  genes_in_data <- rownames(combined_layers@assays$Spatial.024um$counts)[rownames(combined_layers@assays$Spatial.024um$counts) %in% gene_panel$gene]
  combined_count_table <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[genes_in_data ,])))

  combined_count_table_subset <- combined_count_table[, colnames(combined_count_table) %in% gene_panel$gene]
  combined_count_table_subset <- combined_count_table_subset[rownames(combined_count_table_subset) %in% rownames(usage_norm),]

  # 2.3.1 Make a plot that is tissue agnostic ----
  corr_matrix <- cor(combined_count_table_subset, usage_norm)
  corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset, usage_norm)))
  corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

  gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
  corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
  corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

  split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                  levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

  ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_split = split, row_title_rot = 0)
  ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = TRUE,
                    row_split = split, row_title_rot = 0)
  ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = TRUE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0)
  ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                     cluster_columns = TRUE, cluster_rows = TRUE,
                     row_split = split, row_title_rot = 0)
  ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       row_split = split, row_title_rot = 0)
  ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = TRUE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = TRUE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            row_split = split, row_title_rot = 0)

  p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
  p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
  p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
  p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

  pdf(paste(save_path, "/resistance_panel/", "K", RANK, "_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
  print(plot_grid(p1, p2))
  print(plot_grid(p3, p4))
  print(plot_grid(p5, p6))
  print(plot_grid(p7, p8))
  dev.off()

  # 2.3.2 Make a plot that is tissue specific ----
  tissues_list <- unique(sub(paste0("(", pattern_regex, ").*"), "\\1", rownames(usage_norm)))

  for (tissue in tissues_list){
    combined_count_table_subset_tissue <- combined_count_table_subset[grep(tissue, rownames(combined_count_table_subset)), ]
    combined_count_table_subset_tissue <- combined_count_table_subset_tissue[,  colSums(combined_count_table_subset_tissue) != 0,]
    usage_norm_subset <- usage_norm[grep(tissue, rownames(usage_norm)), ]

    corr_matrix <- cor(combined_count_table_subset_tissue, usage_norm_subset)
    corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset_tissue, usage_norm_subset)))

    corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

    gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
    corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]

    corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

    split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                    levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

    ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                  cluster_columns = FALSE, cluster_rows = FALSE,
                  row_split = split, row_title_rot = 0)
    ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = FALSE, cluster_rows = TRUE,
                      row_split = split, row_title_rot = 0)
    ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = TRUE, cluster_rows = FALSE,
                      row_split = split, row_title_rot = 0)
    ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                       cluster_columns = TRUE, cluster_rows = TRUE,
                       row_split = split, row_title_rot = 0)
    ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                         cluster_columns = FALSE, cluster_rows = FALSE,
                         row_split = split, row_title_rot = 0)
    ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = FALSE, cluster_rows = TRUE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = TRUE, cluster_rows = FALSE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                              cluster_columns = TRUE, cluster_rows = TRUE,
                              row_split = split, row_title_rot = 0)

    p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
    p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
    p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
    p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

    pdf(paste(save_path, "/resistance_panel/", "K", RANK, "_", tissue, "_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
    print(plot_grid(p1, p2))
    print(plot_grid(p3, p4))
    print(plot_grid(p5, p6))
    print(plot_grid(p7, p8))
    dev.off()
  }

  # 2.4 For FASLG ----
  resistance_fas_gene_panel <- data.frame(category = c("FASLG_priority"),
                                          gene = c("FASLG", "FAS"))

  gene_panel <- resistance_fas_gene_panel
  genes_in_data <- rownames(combined_layers@assays$Spatial.024um$counts)[rownames(combined_layers@assays$Spatial.024um$counts) %in% gene_panel$gene]
  combined_count_table <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[genes_in_data ,])))

  combined_count_table_subset <- combined_count_table[, colnames(combined_count_table) %in% gene_panel$gene]
  combined_count_table_subset <- combined_count_table_subset[rownames(combined_count_table_subset) %in% rownames(usage_norm),]

  # 2.4.1 Make a plot that is tissue agnostic ----
  corr_matrix <- cor(combined_count_table_subset, usage_norm)
  corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset, usage_norm)))
  corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

  gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
  corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
  corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

  split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                  levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

  ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_split = split, row_title_rot = 0)
  ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = TRUE,
                    row_split = split, row_title_rot = 0)
  ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = TRUE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0)
  ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                     cluster_columns = TRUE, cluster_rows = TRUE,
                     row_split = split, row_title_rot = 0)
  ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       row_split = split, row_title_rot = 0)
  ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = TRUE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = TRUE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            row_split = split, row_title_rot = 0)

  p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
  p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
  p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
  p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

  pdf(paste(save_path, "/fas_panel/", "K", RANK, "_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
  print(plot_grid(p1, p2))
  print(plot_grid(p3, p4))
  print(plot_grid(p5, p6))
  print(plot_grid(p7, p8))
  dev.off()

  # 2.4.2 Make a plot that is tissue specific ----
  tissues_list <- unique(sub(paste0("(", pattern_regex, ").*"), "\\1", rownames(usage_norm)))

  for (tissue in tissues_list){
    combined_count_table_subset_tissue <- combined_count_table_subset[grep(tissue, rownames(combined_count_table_subset)), ]
    combined_count_table_subset_tissue <- combined_count_table_subset_tissue[,  colSums(combined_count_table_subset_tissue) != 0,]
    if (!is.null(dim(combined_count_table_subset_tissue))){
      usage_norm_subset <- usage_norm[grep(tissue, rownames(usage_norm)), ]

      corr_matrix <- cor(combined_count_table_subset_tissue, usage_norm_subset)
      corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset_tissue, usage_norm_subset)))

      corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

      gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
      corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]

      corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

      split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                      levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

      ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0)
      ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                        cluster_columns = FALSE, cluster_rows = TRUE,
                        row_split = split, row_title_rot = 0)
      ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                        cluster_columns = TRUE, cluster_rows = FALSE,
                        row_split = split, row_title_rot = 0)
      ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                         cluster_columns = TRUE, cluster_rows = TRUE,
                         row_split = split, row_title_rot = 0)
      ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0)
      ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                               cluster_columns = FALSE, cluster_rows = TRUE,
                               row_split = split, row_title_rot = 0)
      ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                               cluster_columns = TRUE, cluster_rows = FALSE,
                               row_split = split, row_title_rot = 0)
      ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                                cluster_columns = TRUE, cluster_rows = TRUE,
                                row_split = split, row_title_rot = 0)

      p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
      p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
      p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
      p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

      pdf(paste(save_path, "/fas_panel/", "K", RANK, "_", tissue, "_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
      print(plot_grid(p1, p2))
      print(plot_grid(p3, p4))
      print(plot_grid(p5, p6))
      print(plot_grid(p7, p8))
      dev.off()
    }
  }
  # 2.5 For resistance w/o GPNMB ----
  resistance_gene_panel <- data.frame(category = c(rep("Pair_1", 2),
                                                   rep("Pair_2", 3),
                                                   rep("Pair_3", 5),
                                                   rep("Pair_4", 2),
                                                   rep("Pair_5", 2),
                                                   rep("Pair_6", 3),
                                                   rep("Pair_7", 3)),
                                      gene = c("PDCD1", "CD274",
                                               "CTLA4", "CD80", "CD86",
                                               "TIGIT", "PVR", "NECTIN2", "CD226", "CD96",
                                               "KLRB1", "CLEC2D",
                                               "FAS", "FASLG",
                                               "HAVCR2", "LGALS9", "CEACAM1",
                                               "LAG3", "FGL1", "LGALS3"))

  gene_panel <- resistance_gene_panel
  genes_in_data <- rownames(combined_layers@assays$Spatial.024um$counts)[rownames(combined_layers@assays$Spatial.024um$counts) %in% gene_panel$gene]
  combined_count_table <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[genes_in_data ,])))

  combined_count_table_subset <- combined_count_table[, colnames(combined_count_table) %in% gene_panel$gene]
  combined_count_table_subset <- combined_count_table_subset[rownames(combined_count_table_subset) %in% rownames(usage_norm),]

  # 2.5.1 Make a plot that is tissue agnostic ----
  corr_matrix <- cor(combined_count_table_subset, usage_norm)
  corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset, usage_norm)))
  corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

  gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
  corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
  corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

  split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                  levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

  ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_split = split, row_title_rot = 0)
  ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = TRUE,
                    row_split = split, row_title_rot = 0)
  ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = TRUE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0)
  ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                     cluster_columns = TRUE, cluster_rows = TRUE,
                     row_split = split, row_title_rot = 0)
  ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       row_split = split, row_title_rot = 0)
  ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = TRUE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = TRUE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            row_split = split, row_title_rot = 0)

  p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
  p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
  p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
  p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

  pdf(paste(save_path, "/resistance_panel/wo_gpnmb/", "K", RANK, "_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
  print(plot_grid(p1, p2))
  print(plot_grid(p3, p4))
  print(plot_grid(p5, p6))
  print(plot_grid(p7, p8))
  dev.off()

  # 2.5.2 Make a plot that is tissue specific ----
  tissues_list <- unique(sub(paste0("(", pattern_regex, ").*"), "\\1", rownames(usage_norm)))

  for (tissue in tissues_list){
    combined_count_table_subset_tissue <- combined_count_table_subset[grep(tissue, rownames(combined_count_table_subset)), ]
    combined_count_table_subset_tissue <- combined_count_table_subset_tissue[,  colSums(combined_count_table_subset_tissue) != 0,]
    usage_norm_subset <- usage_norm[grep(tissue, rownames(usage_norm)), ]

    corr_matrix <- cor(combined_count_table_subset_tissue, usage_norm_subset)
    corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset_tissue, usage_norm_subset)))

    corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

    gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
    corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]

    corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

    split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                    levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

    ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                  cluster_columns = FALSE, cluster_rows = FALSE,
                  row_split = split, row_title_rot = 0)
    ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = FALSE, cluster_rows = TRUE,
                      row_split = split, row_title_rot = 0)
    ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = TRUE, cluster_rows = FALSE,
                      row_split = split, row_title_rot = 0)
    ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                       cluster_columns = TRUE, cluster_rows = TRUE,
                       row_split = split, row_title_rot = 0)
    ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                         cluster_columns = FALSE, cluster_rows = FALSE,
                         row_split = split, row_title_rot = 0)
    ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = FALSE, cluster_rows = TRUE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = TRUE, cluster_rows = FALSE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                              cluster_columns = TRUE, cluster_rows = TRUE,
                              row_split = split, row_title_rot = 0)

    p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
    p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
    p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
    p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

    pdf(paste(save_path, "/resistance_panel/wo_gpnmb/", "K", RANK, "_", tissue, "_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
    print(plot_grid(p1, p2))
    print(plot_grid(p3, p4))
    print(plot_grid(p5, p6))
    print(plot_grid(p7, p8))
    dev.off()
  }
  # 2.6 Make a plot showcasing gene expression across programs for Kevin's genes of interest ----
  kevin_gene_panel <- data.frame(category = c(rep("molecular_subtypes", 5),
                                              rep("proliferation_cell_cycle", 5),
                                              "CAR_Target"),
                                 gene = c("IDH1", "IDH2", "ATRX", "CIC", "FUBP1",
                                          "MKI67", "CDKN2A", "CDKN2B", "CCND1", "CCND2",
                                          "DLL3"))
  gene_panel <- kevin_gene_panel
  genes_in_data <- rownames(combined_layers@assays$Spatial.024um$counts)[rownames(combined_layers@assays$Spatial.024um$counts) %in% gene_panel$gene]
  combined_count_table <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[genes_in_data ,])))

  combined_count_table_subset <- combined_count_table[, colnames(combined_count_table) %in% gene_panel$gene]
  combined_count_table_subset <- combined_count_table_subset[rownames(combined_count_table_subset) %in% rownames(usage_norm),]

  # 2.6.1 Make a plot that is tissue agnostic ----
  corr_matrix <- cor(combined_count_table_subset, usage_norm)
  corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset, usage_norm)))
  corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

  gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
  corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
  corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

  split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                  levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

  ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_split = split, row_title_rot = 0)
  ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = TRUE,
                    row_split = split, row_title_rot = 0)
  ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = TRUE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0)
  ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                     cluster_columns = TRUE, cluster_rows = TRUE,
                     row_split = split, row_title_rot = 0)
  ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       row_split = split, row_title_rot = 0)
  ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = TRUE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = TRUE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            row_split = split, row_title_rot = 0)

  p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
  p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
  p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
  p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

  pdf(paste(save_path, "/kevin_panel/", "K", RANK, "_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
  print(plot_grid(p1, p2))
  print(plot_grid(p3, p4))
  print(plot_grid(p5, p6))
  print(plot_grid(p7, p8))
  dev.off()

  # 2.6.2 Make a plot that is tissue specific ----
  tissues_list <- unique(sub(paste0("(", pattern_regex, ").*"), "\\1", rownames(usage_norm)))

  for (tissue in tissues_list){
    combined_count_table_subset_tissue <- combined_count_table_subset[grep(tissue, rownames(combined_count_table_subset)), ]
    combined_count_table_subset_tissue <- combined_count_table_subset_tissue[,  colSums(combined_count_table_subset_tissue) != 0,]
    usage_norm_subset <- usage_norm[grep(tissue, rownames(usage_norm)), ]

    corr_matrix <- cor(combined_count_table_subset_tissue, usage_norm_subset)
    corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset_tissue, usage_norm_subset)))

    corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

    gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
    corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]

    corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

    split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                    levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

    ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                  cluster_columns = FALSE, cluster_rows = FALSE,
                  row_split = split, row_title_rot = 0)
    ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = FALSE, cluster_rows = TRUE,
                      row_split = split, row_title_rot = 0)
    ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = TRUE, cluster_rows = FALSE,
                      row_split = split, row_title_rot = 0)
    ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                       cluster_columns = TRUE, cluster_rows = TRUE,
                       row_split = split, row_title_rot = 0)
    ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                         cluster_columns = FALSE, cluster_rows = FALSE,
                         row_split = split, row_title_rot = 0)
    ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = FALSE, cluster_rows = TRUE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = TRUE, cluster_rows = FALSE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                              cluster_columns = TRUE, cluster_rows = TRUE,
                              row_split = split, row_title_rot = 0)

    p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
    p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
    p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
    p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

    pdf(paste(save_path, "/kevin_panel/", "K", RANK, "_", tissue, "_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
    print(plot_grid(p1, p2))
    print(plot_grid(p3, p4))
    print(plot_grid(p5, p6))
    print(plot_grid(p7, p8))
    dev.off()
  }

  # 2.7 For Franz Resistance ----
  resistance_franz_gene_panel <- data.frame(category = c("franz_resistance"),
                                            gene = c("TIGIT", "PVR", "NECTIN2", "CD226", "CD96"))

  gene_panel <- resistance_franz_gene_panel
  genes_in_data <- rownames(combined_layers@assays$Spatial.024um$counts)[rownames(combined_layers@assays$Spatial.024um$counts) %in% gene_panel$gene]
  combined_count_table <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[genes_in_data ,])))

  combined_count_table_subset <- combined_count_table[, colnames(combined_count_table) %in% gene_panel$gene]
  combined_count_table_subset <- combined_count_table_subset[rownames(combined_count_table_subset) %in% rownames(usage_norm),]

  # 2.7.1 Make a plot that is tissue agnostic ----
  corr_matrix <- cor(combined_count_table_subset, usage_norm)
  corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset, usage_norm)))
  corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

  gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
  corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
  corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

  split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                  levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

  ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_split = split, row_title_rot = 0)
  ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = TRUE,
                    row_split = split, row_title_rot = 0)
  ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = TRUE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0)
  ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                     cluster_columns = TRUE, cluster_rows = TRUE,
                     row_split = split, row_title_rot = 0)
  ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       row_split = split, row_title_rot = 0)
  ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = TRUE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = TRUE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            row_split = split, row_title_rot = 0)

  p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
  p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
  p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
  p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

  pdf(paste(franz_resistance_panel_path, "K", RANK, "_franz_resistance_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
  print(plot_grid(p1, p2))
  print(plot_grid(p3, p4))
  print(plot_grid(p5, p6))
  print(plot_grid(p7, p8))
  dev.off()

  # 2.7.2 Make a plot that is tissue specific ----
  tissues_list <- unique(sub(paste0("(", pattern_regex, ").*"), "\\1", rownames(usage_norm)))

  for (tissue in tissues_list){
    combined_count_table_subset_tissue <- combined_count_table_subset[grep(tissue, rownames(combined_count_table_subset)), ]
    combined_count_table_subset_tissue <- combined_count_table_subset_tissue[,  colSums(combined_count_table_subset_tissue) != 0,]
    usage_norm_subset <- usage_norm[grep(tissue, rownames(usage_norm)), ]

    corr_matrix <- cor(combined_count_table_subset_tissue, usage_norm_subset)
    corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset_tissue, usage_norm_subset)))

    corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

    gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
    corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]

    corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

    split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                    levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

    ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                  cluster_columns = FALSE, cluster_rows = FALSE,
                  row_split = split, row_title_rot = 0)
    ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = FALSE, cluster_rows = TRUE,
                      row_split = split, row_title_rot = 0)
    ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = TRUE, cluster_rows = FALSE,
                      row_split = split, row_title_rot = 0)
    ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                       cluster_columns = TRUE, cluster_rows = TRUE,
                       row_split = split, row_title_rot = 0)
    ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                         cluster_columns = FALSE, cluster_rows = FALSE,
                         row_split = split, row_title_rot = 0)
    ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = FALSE, cluster_rows = TRUE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = TRUE, cluster_rows = FALSE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                              cluster_columns = TRUE, cluster_rows = TRUE,
                              row_split = split, row_title_rot = 0)

    p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
    p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
    p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
    p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

    pdf(paste(franz_resistance_panel_path, "K", RANK, "_", tissue, "_franz_resistance_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
    print(plot_grid(p1, p2))
    print(plot_grid(p3, p4))
    print(plot_grid(p5, p6))
    print(plot_grid(p7, p8))
    dev.off()
  }
  # 2.8 For GPNMB ----
  gpnmb_gene_panel <- data.frame(category = c("ASPS", "ASPS"),
                                 gene = c("GPNMB", "ANGPTL2"))

  gene_panel <- gpnmb_gene_panel
  genes_in_data <- rownames(combined_layers@assays$Spatial.024um$counts)[rownames(combined_layers@assays$Spatial.024um$counts) %in% gene_panel$gene]
  combined_count_table <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[genes_in_data ,])))

  combined_count_table_subset <- combined_count_table[, colnames(combined_count_table) %in% gene_panel$gene]
  combined_count_table_subset <- combined_count_table_subset[rownames(combined_count_table_subset) %in% rownames(usage_norm), ]
  combined_count_table_subset <- combined_count_table_subset[, "GPNMB", drop = FALSE]

  # 2.8.1 Make a plot that is tissue agnostic ----
  corr_matrix <- cor(combined_count_table_subset, usage_norm)
  corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset, usage_norm)))
  corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

  corr_matrix <- as.matrix(corr_matrix)
  corr_matrix_cosine <- t(as.matrix(corr_matrix_cosine)); rownames(corr_matrix_cosine) <- "GPNMB"

  split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                  levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

  ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_split = split, row_title_rot = 0)
  ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = TRUE,
                    row_split = split, row_title_rot = 0)
  ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = TRUE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0)
  ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                     cluster_columns = TRUE, cluster_rows = TRUE,
                     row_split = split, row_title_rot = 0)
  ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       row_split = split, row_title_rot = 0)
  ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = TRUE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = TRUE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0)
  ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            row_split = split, row_title_rot = 0)

  p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
  p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
  p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
  p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

  pdf(paste(save_path, "/gpnmb_only/", "K", RANK, "_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
  print(plot_grid(p1, p2))
  print(plot_grid(p3, p4))
  print(plot_grid(p5, p6))
  print(plot_grid(p7, p8))
  dev.off()

  # 2.8.2 Make a plot that is tissue specific ----
  tissues_list <- unique(sub(paste0("(", pattern_regex, ").*"), "\\1", rownames(usage_norm)))

  tissue_df_cor <- list()
  tissue_df_cor_cosine <- list()

  for (tissue in tissues_list){
    combined_count_table_subset_tissue <- combined_count_table_subset[grep(tissue, rownames(combined_count_table_subset)), ]
    # combined_count_table_subset_tissue <- combined_count_table_subset_tissue[, colSums(combined_count_table_subset_tissue) != 0,]
    usage_norm_subset <- usage_norm[grep(tissue, rownames(usage_norm)), ]

    corr_matrix <- cor(combined_count_table_subset_tissue, usage_norm_subset)
    corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset_tissue, usage_norm_subset)))
    corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

    corr_matrix <- as.matrix(corr_matrix); rownames(corr_matrix) <- "GPNMB"
    corr_matrix_cosine <- t(as.matrix(corr_matrix_cosine)); rownames(corr_matrix_cosine) <- "GPNMB"

    tissue_df_cor[[tissue]] <- as.data.frame(corr_matrix)
    tissue_df_cor_cosine[[tissue]] <- as.data.frame(corr_matrix_cosine)

    split <- factor(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)],
                    levels = unique(gene_panel$category[gene_panel$gene %in% rownames(corr_matrix)]))

    ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                  cluster_columns = FALSE, cluster_rows = FALSE,
                  row_split = split, row_title_rot = 0)
    ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = FALSE, cluster_rows = TRUE,
                      row_split = split, row_title_rot = 0)
    ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                      cluster_columns = TRUE, cluster_rows = FALSE,
                      row_split = split, row_title_rot = 0)
    ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                       cluster_columns = TRUE, cluster_rows = TRUE,
                       row_split = split, row_title_rot = 0)
    ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                         cluster_columns = FALSE, cluster_rows = FALSE,
                         row_split = split, row_title_rot = 0)
    ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = FALSE, cluster_rows = TRUE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                             cluster_columns = TRUE, cluster_rows = FALSE,
                             row_split = split, row_title_rot = 0)
    ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                              cluster_columns = TRUE, cluster_rows = TRUE,
                              row_split = split, row_title_rot = 0)

    p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
    p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
    p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
    p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

    pdf(paste(save_path, "/gpnmb_only/", "K", RANK, "_", tissue, "_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
    print(plot_grid(p1, p2))
    print(plot_grid(p3, p4))
    print(plot_grid(p5, p6))
    print(plot_grid(p7, p8))
    dev.off()
  }

  tissue_df_cor <- rbindlist(tissue_df_cor, idcol = "ID")
  tissue_df_cor_cosine <- rbindlist(tissue_df_cor_cosine, idcol = "ID")

  rn_cor <- tissue_df_cor$ID
  rn_cos <- tissue_df_cor_cosine$ID

  tissue_df_cor <- tissue_df_cor[, !"ID", with = FALSE]
  tissue_df_cor_cosine <- tissue_df_cor_cosine[, !"ID", with = FALSE]

  tissue_df_cor <- as.matrix(tissue_df_cor)
  tissue_df_cor_cosine <- as.matrix(tissue_df_cor_cosine)

  rownames(tissue_df_cor) <- rn_cor
  rownames(tissue_df_cor_cosine) <- rn_cos

  ht <- Heatmap(tissue_df_cor, name = "GPNMB Gene Expression and Program Usage Pearson Correlation",
                cluster_columns = FALSE, cluster_rows = FALSE, row_title_rot = 0)
  ht_row <- Heatmap(tissue_df_cor, name = "GPNMB Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = TRUE, row_title_rot = 0)
  ht_col <- Heatmap(tissue_df_cor, name = "GPNMB Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = TRUE, cluster_rows = FALSE, row_title_rot = 0)
  ht_both <- Heatmap(tissue_df_cor, name = "GPNMB Gene Expression and Program Usage Pearson Correlation",
                     cluster_columns = TRUE, cluster_rows = TRUE, row_title_rot = 0)
  ht_cosine <- Heatmap(tissue_df_cor_cosine, name = "GPNMB Gene Expression and Program Usage Cosine Correlation",
                       cluster_columns = FALSE, cluster_rows = FALSE, row_title_rot = 0)
  ht_cosine_row <- Heatmap(tissue_df_cor_cosine, name = "GPNMB Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = TRUE, row_title_rot = 0)
  ht_cosine_col <- Heatmap(tissue_df_cor_cosine, name = "GPNMB Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = TRUE, cluster_rows = FALSE, row_title_rot = 0)
  ht_cosine_both <- Heatmap(tissue_df_cor_cosine, name = "GPNMB Gene Expression and Program Usage Cosine Correlation",
                            cluster_columns = TRUE, cluster_rows = TRUE, row_title_rot = 0)

  p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
  p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
  p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
  p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

  pdf(paste(save_path, "/gpnmb_only/", "K", RANK, "_tissue_combined_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 10)
  print(plot_grid(p1, p2))
  print(plot_grid(p3, p4))
  print(plot_grid(p5, p6))
  print(plot_grid(p7, p8))
  dev.off()
}

# Neighbourhood Analysis ----
# Create or just load the spot_neighbours df for all neighborhood analysis
if (NEIGHBORHOOD_ANALYSIS_PROGRAMS || NEIGHBOURHOOD_ANALYSIS_GCAR_CELLS ||
    NEIGHBORHOOD_ANALYSIS_SUMMARY ||
    INSIDE_OUTSIDE_T_CELL || INSIDE_OUTSIDE_NICHE){
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))

  union_cell_coordinate <- seurat_object_coords()
  combined_metadata_table <- as.data.frame(combined_layers@meta.data)
  options(future.globals.maxSize = 1e20)

  # Check if the spot neighbours file exists, if not, then make it and save it
  if (!file.exists(paste("02_Visium_HD/03_datasets/04_bin_neighbours/ASPS-02/", "spot_neighbours_df.rds", sep = ""))){
    spot_neighbours_df <- data.frame()

    for (sample in so_list) {
      sample_so <- eval(as.name(sample))

      # Get spot radii from each sample
      sample_spot_radii <- sample_so@images[["slice1.024um"]]@boundaries[["centroids"]]@radius

      sample_cell_coordinates <- union_cell_coordinate[union_cell_coordinate$sample == sample,]
      sample_spot_neighbours_df <- program_neighbours(rownames(sample_cell_coordinates),
                                                      sample_cell_coordinates,
                                                      radius = (sample_spot_radii+1))

      spot_neighbours_df <- rbind(spot_neighbours_df, sample_spot_neighbours_df)
    }

    saveRDS(spot_neighbours_df, file = paste(SAVE_FOLDER, "spot_neighbours_df.rds", sep = ""))
  } else{spot_neighbours_df <- readRDS(paste("02_Visium_HD/03_datasets/04_bin_neighbours/ASPS-02/", "spot_neighbours_df.rds", sep = ""))}
  print("Spot Neighbours Dataframe Created and Saved:")
  print(Sys.time())
  if (!file.exists(paste("02_Visium_HD/03_datasets/04_bin_neighbours/ASPS-02/", "spot_neighbours_df_wo_overlaps.rds", sep = ""))){
    spot_neighbours_df_wo_overlaps <- data.frame()

    for (i in 1:length(so_list)){
      # Subset the spot_neighbours_df to a specific sample
      spot_neighbours_subset_df <- spot_neighbours_df[spot_neighbours_df$sample == so_list[i],]

      # Remove the overlapping neighbours for a specific sample
      spot_neighbours_df_wo_overlaps_subset <- remove_overlapping_neighbours(spot_neighbours_subset_df)

      # Rbind the results to the final spot_neighbours_df_wo_overlaps
      spot_neighbours_df_wo_overlaps <- rbind(spot_neighbours_df_wo_overlaps, spot_neighbours_df_wo_overlaps_subset)
    }

    # spot_neighbours_df_wo_overlaps <- spot_neighbours_df_wo_overlaps[spot_neighbours_df_wo_overlaps$nr_neighbours != 0,]
    saveRDS(spot_neighbours_df_wo_overlaps, file = paste(SAVE_FOLDER, "spot_neighbours_df_wo_overlaps.rds", sep = ""))
  } else{spot_neighbours_df_wo_overlaps <- readRDS(paste("02_Visium_HD/03_datasets/04_bin_neighbours/ASPS-02/", "spot_neighbours_df_wo_overlaps.rds", sep = ""))}
  print("Spot Neighbours Without Overlaps Dataframe Created and Saved:")
  print(Sys.time())
  
  combined_metadata_table <- combined_metadata_table[rownames(combined_metadata_table) %in% rownames(usage_norm),]
  spot_neighbours_df <- spot_neighbours_df[spot_neighbours_df$spot_of_interest %in% rownames(usage_norm),]
  spot_neighbours_df_wo_overlaps <- spot_neighbours_df_wo_overlaps[spot_neighbours_df_wo_overlaps$spot_of_interest %in% rownames(usage_norm),]

  program_var_names <- c(ls(), "program_var_names")
}

# Modify the radii to R0, 1, 2, and 3 ----
if (NEIGHBORHOOD_ANALYSIS_PROGRAMS || NEIGHBOURHOOD_ANALYSIS_GCAR_CELLS ||
    NEIGHBORHOOD_ANALYSIS_SUMMARY ||
    INSIDE_OUTSIDE_T_CELL || INSIDE_OUTSIDE_NICHE){
  sample_spot_radii <- list()

  for (sample in so_list) {
    sample_so <- eval(as.name(sample))

    if (length(sample_spot_radii) == 0){
      sample_spot_radii <- (sample_so@images[["slice1.024um"]]@boundaries[["centroids"]]@radius + 1)
    } else{
      sample_radii <- (sample_so@images[["slice1.024um"]]@boundaries[["centroids"]]@radius + 1)
      sample_spot_radii <- append(sample_spot_radii, sample_radii)
    }
  }

  spot_neighbours_df$r[spot_neighbours_df$r == 0] <- "R0"
  for (radii in sample_spot_radii){
    spot_neighbours_df$r[spot_neighbours_df$r == radii] <- "R1"
    spot_neighbours_df$r[spot_neighbours_df$r == (radii*2)] <- "R2"
    spot_neighbours_df$r[spot_neighbours_df$r == (radii*3)] <- "R3"
  }

  spot_neighbours_df_wo_overlaps$r[spot_neighbours_df_wo_overlaps$r == 0] <- "R0"
  for (radii in sample_spot_radii){
    spot_neighbours_df_wo_overlaps$r[spot_neighbours_df_wo_overlaps$r == radii] <- "R1"
    spot_neighbours_df_wo_overlaps$r[spot_neighbours_df_wo_overlaps$r == (radii*2)] <- "R2"
    spot_neighbours_df_wo_overlaps$r[spot_neighbours_df_wo_overlaps$r == (radii*3)] <- "R3"
  }
}

# 3.1 Neighbourhood of the GCAR+ / T-cell Program+ Cells ----
if (NEIGHBOURHOOD_ANALYSIS_GCAR_CELLS){
  # 1.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "K", RANK, "/neighbourhood_analysis/GCAR_T_Cell_Program_Cells/", sep = "")
  create_folder(save_path)

  usage_threshold <- GLOBAL_THRESHOLD

  # 1.1 Get the GCAR+ cells, the spot coordinates from the 2 samples, and their intersection ----
  GCAR_cell_so <- subset(x = combined_layers, subset = GCAR > 0, slot = "counts")
  GCAR_spots <- rownames(GCAR_cell_so@meta.data)

  # Get the spot names and merge them with the union_cell_coordinate dataframe
  GCAR_spots <- intersect(GCAR_spots, rownames(union_cell_coordinate))

  # Get the T-cell Program+ bins
  gep_pos_spots <- rownames(usage_norm[usage_norm[, "GEP5"] >= usage_threshold,])

  # Intersection of GCAR+ and T-cell Program+ Spots
  GCAR_spots <- intersect(GCAR_spots, gep_pos_spots)

  # 1.2 Subset the spot neighbors ----
  spot_neighbours_df_subset <- spot_neighbours_df_wo_overlaps[spot_neighbours_df_wo_overlaps$spot_of_interest %in% GCAR_spots,]

  # 1.3 Binarize our normalized usage ----
  binarized_gep <- usage_norm
  binarized_gep[binarized_gep >= 0.05] <- 1; binarized_gep[binarized_gep < 0.05] <- 0

  # 1.4 Make the spot_neighbour_programs_df ----
  # The final_df will have additional columns that are named on the unique program-radii combinations
  # Cbind the program_neighbors_df with an empty matrix w/ ncol = unique programs * unique radii
  program_radius_comb_df <- expand.grid(colnames(usage_norm), unique(spot_neighbours_df_subset$r))
  program_radius_comb_df$GEP_r_comb <- paste(program_radius_comb_df$Var1, program_radius_comb_df$Var2, sep = "_")
  temp_matrix <- matrix(0,
                        nrow = 1,
                        ncol = RANK*length(unique(spot_neighbours_df_subset$r)),
                        dimnames = list(NULL, program_radius_comb_df$GEP_r_comb))
  spot_neighbour_programs_df <- cbind(spot_neighbours_df_subset, temp_matrix)

  # 1.5 Get the binarized usage sum per program for all a spots neighbors ----
  # Iterate i over each row (spot-radius combination) of the program_neighbours_df
  for (i in seq_len(nrow(spot_neighbour_programs_df))){
    # Get the neighbors for the current row (spot-radius combination)
    neighbors <- unlist(spot_neighbour_programs_df$neighbours[i])

    # Get the binary gene expression programs for the neighbors
    neighbor_programs <- binarized_gep[rownames(binarized_gep) %in% neighbors, , drop = FALSE]

    # Compute the column sums for the neighbors
    neighbor_programs_sum <- colSums(neighbor_programs) / nrow(neighbor_programs)

    # Update spot_neighbour_programs_df with the calculated sums using column name matching
    colnames_to_paste <- paste(colnames(neighbor_programs), spot_neighbour_programs_df$r[i], sep = "_")
    spot_neighbour_programs_df[i, c(colnames_to_paste)] <- neighbor_programs_sum
  }

  # 1.6 Remove the spot-radii combs, thereby flattening the df ----
  spot_neighbour_programs_df_flat <- data.frame()
  for (soi in unique(spot_neighbour_programs_df$spot_of_interest)){
    soi_program_subset <- spot_neighbour_programs_df[spot_neighbour_programs_df$spot_of_interest == soi,]

    gep_sums <- soi_program_subset[, c(7:ncol(soi_program_subset))]
    gep_sums <- colSums(gep_sums)

    row_data <- c(soi_program_subset[1, c(1, 2, 4)], c(soi_program_subset$nr_neighbours), gep_sums)
    row_data <- as.data.frame(row_data)
    colnames(row_data)[4:7] <- c("nr_neigh_r0", "nr_neigh_r1", "nr_neigh_r2", "nr_neigh_r3")

    if (length(spot_neighbour_programs_df_flat) == 0){
      spot_neighbour_programs_df_flat <- row_data
    }
    else{
      spot_neighbour_programs_df_flat <- rbind(spot_neighbour_programs_df_flat, row_data)
    }
  }

  # 1.7 Modify the spot_neighbour_programs_df_flat for plotting ----
  # Add a tissue column and organize the columns
  spot_neighbour_programs_df_flat <- as.data.frame(spot_neighbour_programs_df_flat)
  spot_neighbour_programs_df_flat[, c(4:7)] <- as.numeric(unlist(spot_neighbour_programs_df_flat[, c(4:7)]))
  spot_neighbour_programs_df_flat$tissue <- spot_neighbour_programs_df_flat$sample
  spot_neighbour_programs_df_flat$tissue <- sub("bS1", "bS", spot_neighbour_programs_df_flat$tissue)
  spot_neighbour_programs_df_flat$tissue <- sub("bS2", "bS", spot_neighbour_programs_df_flat$tissue)
  spot_neighbour_programs_df_flat <- spot_neighbour_programs_df_flat[, c(1:7, ncol(spot_neighbour_programs_df_flat), 8:(ncol(spot_neighbour_programs_df_flat) - 1))]

  # Make the annotation and plotting df
  spot_neighbour_programs_df_flat_annot <- spot_neighbour_programs_df_flat[,1:8]
  spot_neighbour_programs_df_flat_plot <- t(spot_neighbour_programs_df_flat[, c(9:ncol(spot_neighbour_programs_df_flat))])
  colnames(spot_neighbour_programs_df_flat_plot) <- c(spot_neighbour_programs_df_flat$spot_of_interest)

  # 1.8 Create the annotations ----
  annotate_genes <- c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "GCAR")
  combined_count_table_subset <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[annotate_genes , spot_neighbour_programs_df_flat_annot$spot_of_interest])))
  spot_neighbour_programs_df_flat_annot <- cbind(spot_neighbour_programs_df_flat_annot, combined_count_table_subset)

  haT <- HeatmapAnnotation(sample = spot_neighbour_programs_df_flat_annot$sample,
                           nr_neigh_r0 = spot_neighbour_programs_df_flat_annot$nr_neigh_r0,
                           nr_neigh_r1 = spot_neighbour_programs_df_flat_annot$nr_neigh_r1,
                           nr_neigh_r2 = spot_neighbour_programs_df_flat_annot$nr_neigh_r2,
                           nr_neigh_r3 = spot_neighbour_programs_df_flat_annot$nr_neigh_r3,
                           CD3D = spot_neighbour_programs_df_flat_annot$CD3D,
                           CD3E = spot_neighbour_programs_df_flat_annot$CD3E,
                           CD3G = spot_neighbour_programs_df_flat_annot$CD3G,
                           CD4 = spot_neighbour_programs_df_flat_annot$CD4,
                           CD8A = spot_neighbour_programs_df_flat_annot$CD8A,
                           GCAR = spot_neighbour_programs_df_flat_annot$GCAR,
                           col = list(sample = sample_cols,
                                      nr_neigh_r0 = colorRamp2(c(range(c(0, spot_neighbour_programs_df_flat_annot$nr_neigh_r0)), range(spot_neighbour_programs_df_flat_annot$nr_neigh_r0) + 1), hcl_palette = "YlOrRd", reverse = TRUE),
                                      nr_neigh_r1 = colorRamp2(range(c(0, spot_neighbour_programs_df_flat_annot$nr_neigh_r1)), hcl_palette = "YlOrRd", reverse = TRUE),
                                      nr_neigh_r2 = colorRamp2(range(c(0, spot_neighbour_programs_df_flat_annot$nr_neigh_r2)), hcl_palette = "YlOrRd", reverse = TRUE),
                                      nr_neigh_r3 = colorRamp2(range(c(0, spot_neighbour_programs_df_flat_annot$nr_neigh_r3)), hcl_palette = "YlOrRd", reverse = TRUE),
                                      CD3D = colorRamp2(range(c(0, 1, spot_neighbour_programs_df_flat_annot$CD3D)), hcl_palette = "Reds", reverse = TRUE),
                                      CD3E = colorRamp2(range(c(0, 1, spot_neighbour_programs_df_flat_annot$CD3E)), hcl_palette = "Reds", reverse = TRUE),
                                      CD3G = colorRamp2(range(c(0, 1, spot_neighbour_programs_df_flat_annot$CD3G)), hcl_palette = "Reds", reverse = TRUE),
                                      CD4 = colorRamp2(range(c(0, 1, spot_neighbour_programs_df_flat_annot$CD4)), hcl_palette = "Reds", reverse = TRUE),
                                      CD8A = colorRamp2(range(c(0, 1, spot_neighbour_programs_df_flat_annot$CD8A)), hcl_palette = "Reds", reverse = TRUE),
                                      GCAR = colorRamp2(range(c(0, 1, spot_neighbour_programs_df_flat_annot$GCAR)), hcl_palette = "Reds", reverse = TRUE)),
                           which = "col")

  program_radius_comb_df$Var2 <- as.numeric(program_radius_comb_df$Var2)
  program_radius_comb_df$program <- gsub("GEP", "", program_radius_comb_df$Var1)
  program_radius_comb_df$program <- as.numeric(program_radius_comb_df$program)
  haL <- HeatmapAnnotation(radius = program_radius_comb_df$Var2, program = program_radius_comb_df$program,
                           col = list(radius = colorRamp2(range(program_radius_comb_df$Var2),
                                                          hcl_palette = "Blue-Yellow", reverse = TRUE),
                                      program = colorRamp2(range(program_radius_comb_df$program),
                                                           hcl_palette = "Inferno", reverse = TRUE)),
                           which = "row")

  col_spliter <- factor(spot_neighbour_programs_df_flat_annot$tissue, levels = c("pNA", "bS"))
  row_spliter <- factor(program_radius_comb_df$Var1)
  col_fun = colorRamp2(range(spot_neighbour_programs_df_flat_plot), hcl_palette = "Blues", reverse = TRUE)

  # Set the program annotation names
  rownames(spot_neighbour_programs_df_flat_plot) <- paste(program_annotations, c(rep(0, RANK), rep(75, RANK), rep(150, RANK), rep(225, RANK)), sep = "_")
  spot_neighbour_programs_df_flat_plot <- as.matrix(spot_neighbour_programs_df_flat_plot)

  # Determine the plot height
  plot_height <- (0.1969*nrow(spot_neighbour_programs_df_flat_plot)) + 1.3386 + (length(haT@anno_size) * 5 / 10)

  # 1.9 Plotting ----
  # 1.9.1 Plot the spot_neighbour_programs_df_flat ----
  ht_combined <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                         name = "Mean of Neighbour Program",
                         cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun,
                         show_row_names = TRUE, show_column_names = TRUE,
                         column_split = col_spliter, row_split = row_spliter)
  ht_combined_row_clust <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                                   name = "Mean of Neighbour Program",
                                   cluster_columns = TRUE, cluster_rows = TRUE, col = col_fun,
                                   show_row_names = TRUE, show_column_names = TRUE,
                                   column_split = col_spliter)
  ht_combined_col_clust <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                                   name = "Mean of Neighbour Program",
                                   cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun,
                                   show_row_names = TRUE, show_column_names = TRUE,
                                   row_split = row_spliter)
  ht_combined_clust <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                               name = "Mean of Neighbour Program",
                               cluster_columns = TRUE, cluster_rows = TRUE, col = col_fun,
                               show_row_names = TRUE, show_column_names = TRUE)
  pdf(paste(save_path, "GCAR_Cells_Mean_Neighbour_Programs_Binarized_0.05.pdf", sep = ""), width = 20, height = plot_height)
  draw(ht_combined, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_combined_row_clust, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_combined_col_clust, padding = unit(c(3, 3, 3, 3), "cm"))
  draw(ht_combined_clust, padding = unit(c(3, 3, 3, 3), "cm"))
  dev.off()

  # 1.9.2 Save the Rdata image, so plotting changes can be easily performed ----
  save.image(paste(save_path, "GCAR_Cell", "_neighbourhood_save.RData", sep = ""))
}

# 3.2 Neighbourhood of each program ----
if (NEIGHBORHOOD_ANALYSIS_PROGRAMS){
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  usage_threshold <- GLOBAL_THRESHOLD
  Binarized_threshold <- GLOBAL_THRESHOLD

  binarized_gep <- usage_norm
  binarized_gep[binarized_gep >= Binarized_threshold] <- 1; binarized_gep[binarized_gep < Binarized_threshold] <- 0

  # Set up the parallel backend with the number of available cores ----
  loaded_packages <- setdiff(loadedNamespaces(), c("base", "stats", "utils", "graphics", "grDevices", "methods", "datasets"))
  num_cores <- parallelly::availableCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  foreach(program = seq_len(ncol(usage_norm)), .packages = loaded_packages) %dopar% {
    print(program)

    # 2.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
    save_path <- paste(SAVE_FOLDER, "K", RANK, "/neighbourhood_analysis/Programs/GEP_", program, "/", sep = "")
    create_folder(save_path)

    # 2.1 Get the GEP+ cells, the spot coordinates from the 3 samples, and their intersection ----
    gep_pos_spots <- rownames(usage_norm[usage_norm[, program] >= usage_threshold,])

    # Get the spot names and merge them with the union_cell_coordinate dataframe
    gep_pos_spots <- intersect(gep_pos_spots, rownames(union_cell_coordinate))

    # 2.2 Subset the spot neighbors ----
    spot_neighbours_df_subset <- spot_neighbours_df_wo_overlaps[spot_neighbours_df_wo_overlaps$spot_of_interest %in% gep_pos_spots,]

    # 2.3 Make the spot_neighbour_programs_df ----
    # The final_df will have additional columns that are named on the unique program-radii combinations
    # Cbind the program_neighbors_df with an empty matrix w/ ncol = unique programs * unique radii
    program_radius_comb_df <- expand.grid(colnames(usage_norm), unique(spot_neighbours_df_subset$r))
    program_radius_comb_df$GEP_r_comb <- paste(program_radius_comb_df$Var1, program_radius_comb_df$Var2, sep = "_")
    temp_matrix <- matrix(0,
                          nrow = 1,
                          ncol = RANK*length(unique(spot_neighbours_df_subset$r)),
                          dimnames = list(NULL, program_radius_comb_df$GEP_r_comb))
    spot_neighbour_programs_df <- cbind(spot_neighbours_df_subset, temp_matrix)

    # 2.4 Get the binarized usage sum per program for all a spots neighbors ----
    # Iterate i over each row (spot-radius combination) of the program_neighbours_df
    for (i in seq_len(nrow(spot_neighbour_programs_df))){
      # Get the neighbors for the current row (spot-radius combination)
      neighbors <- unlist(spot_neighbour_programs_df$neighbours[i])

      # Get the binary gene expression programs for the neighbors
      neighbor_programs <- binarized_gep[rownames(binarized_gep) %in% neighbors, , drop = FALSE]

      # Compute the column sums for the neighbors
      neighbor_programs_sum <- colSums(neighbor_programs) / nrow(neighbor_programs)

      # Update spot_neighbour_programs_df with the calculated sums using column name matching
      colnames_to_paste <- paste(colnames(neighbor_programs), spot_neighbour_programs_df$r[i], sep = "_")
      spot_neighbour_programs_df[i, c(colnames_to_paste)] <- neighbor_programs_sum
    }

    # 2.5 Remove the spot-radii combs, thereby flattening the df ----
    spot_neighbour_programs_df_flat <- data.frame()
    for (soi in unique(spot_neighbour_programs_df$spot_of_interest)){
      soi_program_subset <- spot_neighbour_programs_df[spot_neighbour_programs_df$spot_of_interest == soi,]

      gep_sums <- soi_program_subset[, c(7:ncol(soi_program_subset))]
      gep_sums <- colSums(gep_sums)

      row_data <- c(soi_program_subset[1, c(1, 2, 4)], c(soi_program_subset$nr_neighbours), gep_sums)
      row_data <- as.data.frame(row_data)
      colnames(row_data)[4:7] <- c("nr_neigh_r0", "nr_neigh_r1", "nr_neigh_r2", "nr_neigh_r3")

      if (length(spot_neighbour_programs_df_flat) == 0){
        spot_neighbour_programs_df_flat <- row_data
      }
      else{
        spot_neighbour_programs_df_flat <- rbind(spot_neighbour_programs_df_flat, row_data)
      }
    }

    # 2.6 Modify the spot_neighbour_programs_df_flat for plotting ----
    # Add a tissue column and organize the columns
    spot_neighbour_programs_df_flat <- as.data.frame(spot_neighbour_programs_df_flat)
    spot_neighbour_programs_df_flat[, c(4:7)] <- as.numeric(unlist(spot_neighbour_programs_df_flat[, c(4:7)]))
    spot_neighbour_programs_df_flat$tissue <- spot_neighbour_programs_df_flat$sample
    spot_neighbour_programs_df_flat <- spot_neighbour_programs_df_flat[, c(1:7, ncol(spot_neighbour_programs_df_flat), 8:(ncol(spot_neighbour_programs_df_flat) - 1))]
    spot_neighbour_programs_df_flat[is.na(spot_neighbour_programs_df_flat)] <- 0

    # Make the annotation and plotting df
    spot_neighbour_programs_df_flat_annot <- spot_neighbour_programs_df_flat[,1:8]
    spot_neighbour_programs_df_flat_annot$GEP <- paste("GEP", program, sep = "_")
    spot_neighbour_programs_df_flat_plot <- t(spot_neighbour_programs_df_flat[, c(9:ncol(spot_neighbour_programs_df_flat))])
    colnames(spot_neighbour_programs_df_flat_plot) <- c(spot_neighbour_programs_df_flat$spot_of_interest)

    # 2.7 Plotting ----
    # 2.7.1 Create the annotations ----
    haT <- HeatmapAnnotation(sample = spot_neighbour_programs_df_flat_annot$sample,
                             nr_neigh_r0 = spot_neighbour_programs_df_flat_annot$nr_neigh_r0,
                             nr_neigh_r1 = spot_neighbour_programs_df_flat_annot$nr_neigh_r1,
                             nr_neigh_r2 = spot_neighbour_programs_df_flat_annot$nr_neigh_r2,
                             nr_neigh_r3 = spot_neighbour_programs_df_flat_annot$nr_neigh_r3,
                             col = list(sample = sample_cols,
                                        nr_neigh_r0 = colorRamp2(c(range(spot_neighbour_programs_df_flat_annot$nr_neigh_r0), range(spot_neighbour_programs_df_flat_annot$nr_neigh_r0)+1), hcl_palette = "YlOrRd", reverse = TRUE),
                                        nr_neigh_r1 = colorRamp2(range(spot_neighbour_programs_df_flat_annot$nr_neigh_r1), hcl_palette = "YlOrRd", reverse = TRUE),
                                        nr_neigh_r2 = colorRamp2(range(spot_neighbour_programs_df_flat_annot$nr_neigh_r2), hcl_palette = "YlOrRd", reverse = TRUE),
                                        nr_neigh_r3 = colorRamp2(range(spot_neighbour_programs_df_flat_annot$nr_neigh_r3), hcl_palette = "YlOrRd", reverse = TRUE)),
                             which = "col")

    program_radius_comb_df$Var2 <- as.numeric(program_radius_comb_df$Var2)
    program_radius_comb_df$program <- gsub("GEP", "", program_radius_comb_df$Var1)
    program_radius_comb_df$program <- as.numeric(program_radius_comb_df$program)
    haL <- HeatmapAnnotation(radius = program_radius_comb_df$Var2, program = program_radius_comb_df$program,
                             col = list(radius = colorRamp2(range(program_radius_comb_df$Var2),
                                                            hcl_palette = "Blue-Yellow", reverse = TRUE),
                                        program = colorRamp2(range(program_radius_comb_df$program),
                                                             hcl_palette = "Inferno", reverse = TRUE)),
                             which = "row")

    # Set the variables that the plots will be split by
    col_spliter <- factor(spot_neighbour_programs_df_flat_annot$tissue, levels = sort(unique(spot_neighbour_programs_df_flat_annot$tissue)))
    col_spliter_GEP <- as.factor(spot_neighbour_programs_df_flat_annot$GEP)
    row_spliter <- factor(program_radius_comb_df$Var1)
    col_fun = colorRamp2(range(spot_neighbour_programs_df_flat_plot), hcl_palette = "Blues", reverse = TRUE)

    # Set the program annotation names
    rownames(spot_neighbour_programs_df_flat_plot) <- paste(program_annotations, c(rep(0, ncol(usage_norm)), rep(75, ncol(usage_norm)), rep(150, ncol(usage_norm)), rep(225, ncol(usage_norm))), sep = "_")
    spot_neighbour_programs_df_flat_plot <- as.matrix(spot_neighbour_programs_df_flat_plot)

    # Determine the plot height
    plot_height <- (0.1969*nrow(spot_neighbour_programs_df_flat_plot)) + 1.3386 + (length(haT@anno_size) * 5 / 10)

    # 2.7.2 Plot but w/o the gene expression tracks ----
    ht_combined <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                           name = "Mean of Neighbour Program",
                           cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun,
                           show_row_names = TRUE, show_column_names = FALSE,
                           column_split = col_spliter, row_split = row_spliter)
    ht_combined_GEP_split <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                                     name = "Mean of Neighbour Program",
                                     cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun,
                                     show_row_names = TRUE, show_column_names = FALSE,
                                     column_split = col_spliter_GEP, row_split = row_spliter)
    ht_combined_row_clust <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                                     name = "Mean of Neighbour Program",
                                     cluster_columns = TRUE, cluster_rows = TRUE, col = col_fun,
                                     show_row_names = TRUE, show_column_names = FALSE,
                                     column_split = col_spliter)
    ht_combined_row_clust_GEP_split <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                                               name = "Mean of Neighbour Program",
                                               cluster_columns = TRUE, cluster_rows = TRUE, col = col_fun,
                                               show_row_names = TRUE, show_column_names = FALSE,
                                               column_split = col_spliter_GEP)
    ht_combined_col_clust <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                                     name = "Mean of Neighbour Program",
                                     cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun,
                                     show_row_names = TRUE, show_column_names = FALSE,
                                     row_split = row_spliter)
    ht_combined_clust <- Heatmap(spot_neighbour_programs_df_flat_plot, left_annotation = haL, top_annotation = haT,
                                 name = "Mean of Neighbour Program",
                                 cluster_columns = TRUE, cluster_rows = TRUE, col = col_fun,
                                 show_row_names = TRUE, show_column_names = FALSE)
    pdf(paste(save_path, "GEP", program, "_T", usage_threshold, "_Binarized_T", Binarized_threshold, "_Mean_Neighbour_Programs_No_Overlaps.pdf", sep = ""), width = 60, height = plot_height)
    draw(ht_combined, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_combined_GEP_split, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_combined_row_clust, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_combined_row_clust_GEP_split, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_combined_col_clust, padding = unit(c(3, 3, 3, 3), "cm"))
    draw(ht_combined_clust, padding = unit(c(3, 3, 3, 3), "cm"))
    dev.off()

    # 2.8 Save the Rdata image, so plotting changes can be easily performed ----
    save(list = ls(environment()), file = paste0(save_path, "GEP", program, "_neighbourhood_save.RData"))
  }
  stopCluster(cl)
}

# 3.3 Custom neighbourhood splitting and dot plot ----
if (NEIGHBORHOOD_ANALYSIS_SUMMARY){
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))

  for (program_selection in NEIGHBORHOOD_ANALYSIS_SUMMARY_PROGRAMS){
    # 3.1 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
    save_path <- paste(SAVE_FOLDER, "K", RANK, "/neighbourhood_analysis/Programs/GEP_",
                       program_selection, "/GEP", program_selection, sep = "")
    usage_threshold <- GLOBAL_THRESHOLD
    Binarized_threshold <- GLOBAL_THRESHOLD
    dendro_split <- 4

    load(paste(save_path, "_neighbourhood_save.RData", sep = ""))

    # 3.2 Make the heatmap, extract Kmeans column order, and then further differentiate large niches ----
    split <- data.frame(cutree(hclust(dist(t(spot_neighbour_programs_df_flat_plot))),
                               k = dendro_split),
                        col_spliter)
    ht_combined <- Heatmap(spot_neighbour_programs_df_flat_plot,
                           left_annotation = haL, top_annotation = haT,
                           name = "Mean of Neighbour Program",
                           cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun,
                           show_row_names = TRUE, show_column_names = FALSE,
                           column_split = split, row_split = row_spliter)
    ht_combined_drawn <- draw(ht_combined, padding = unit(c(3, 3, 3, 3), "cm"))

    # Plot the Kmeans split heatmap
    plot_height <- (0.1969*nrow(spot_neighbour_programs_df_flat_plot)) + 1.3386 + (length(haT@anno_size) * 5 / 10)
    pdf(paste0(save_path, "GEP5_T", usage_threshold, "_Binarized_T", Binarized_threshold, "_Dendro_", dendro_split, "_Split.pdf"), width = 40, height = plot_height)
    print(ht_combined_drawn)
    dev.off()

    heatmap_col_order <- column_order(ht_combined_drawn)

    # Filter the niches
    filtering_list <- list("4,SPSC_biopsy2" = c(23),
                           "4,SPSC_biopsy2" = c(24,26),
                           "2,SPSC_biopsy1" = c(13,14),
                           "2,SPSC_biopsy1" = c(26,27),
                           "1,SPSC_biopsy1" = c(24,26,27),
                           "1,SPSC_biopsy2" = c(10))

    temp <- spot_neighbour_programs_df_flat_plot[unique(unlist(filtering_list)),]

    for (filter in unique(names(filtering_list))){
      num_filters <- length(names(filtering_list)[names(filtering_list) == filter])
      heatmap_col_order_names <- colnames(spot_neighbour_programs_df_flat_plot)[unlist(heatmap_col_order[filter])]
      spot_neighbour_programs_df_flat_plot_subset <- spot_neighbour_programs_df_flat_plot[,colnames(spot_neighbour_programs_df_flat_plot) %in% heatmap_col_order_names]

      for (niche in seq(num_filters)){
        new_niche_name <- paste(filter, niche, sep = ",")
        niche_filter <- filtering_list[names(filtering_list) == filter][[niche]]

        if (length(niche_filter) == 1){
          temp <- spot_neighbour_programs_df_flat_plot_subset[niche_filter, , drop = F]
          temp <- temp[,colSums(temp) >= 1, drop = F]
        }
        if (length(niche_filter) == 2){
          temp <- spot_neighbour_programs_df_flat_plot_subset[niche_filter, , drop = F]
          temp <- temp[,colSums(temp) >= 1, drop = F]
        }
        if (length(niche_filter) == 3){
          temp <- spot_neighbour_programs_df_flat_plot_subset[niche_filter, , drop = F]
          temp <- temp[,colSums(temp) >= 0.75, drop = F]
        }
        if (length(niche_filter) == 4){
          temp <- spot_neighbour_programs_df_flat_plot_subset[niche_filter, , drop = F]
          temp <- temp[,colSums(temp) >= 2, drop = F]
        }

        valid_bins <- intersect(heatmap_col_order_names, colnames(temp))
        heatmap_col_order[[new_niche_name]] <- valid_bins
        heatmap_col_order[[new_niche_name]] <- which(colnames(spot_neighbour_programs_df_flat_plot) %in% heatmap_col_order[[new_niche_name]])
        heatmap_col_order_names <- heatmap_col_order_names[!heatmap_col_order_names %in% valid_bins]
      }

      heatmap_col_order[[filter]] <- heatmap_col_order_names
      heatmap_col_order[[filter]] <- which(colnames(spot_neighbour_programs_df_flat_plot) %in% heatmap_col_order[[filter]])
    }

    heatmap_col_order_final <- unlist(lapply(seq_along(heatmap_col_order), function(i) {
      # Get the name and values of the current list element
      name <- names(heatmap_col_order)[i]
      values <- heatmap_col_order[[i]]

      # Return the name repeated for each value
      rep(name, length(values))
    }))[order(unlist(heatmap_col_order))]

    # 3.3 Plot the final heatmap with the more granular niches ----
    ht_combined <- Heatmap(spot_neighbour_programs_df_flat_plot,
                           left_annotation = haL, top_annotation = haT,
                           name = "Mean of Neighbour Program",
                           cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun,
                           show_row_names = TRUE, show_column_names = FALSE,
                           column_split = heatmap_col_order_final, column_title_rot = 45,
                           row_split = row_spliter)

    plot_height <- (0.1969*nrow(spot_neighbour_programs_df_flat_plot)) + 1.3386 + (length(haT@anno_size) * 5 / 10)
    pdf(paste0(save_path, "GEP5_T", usage_threshold, "_Binarized_T", Binarized_threshold, "_Dendro_", dendro_split, "_Split_Custom.pdf"), width = 40, height = plot_height)
    draw(ht_combined, padding = unit(c(3, 3, 3, 3), "cm"))
    dev.off()

    # 3.4 Make a Dot plot with these niches ----
    # 3.4.1 Prepare the Seurat object with the niche info ----
    combined_layers_subset <- combined_layers
    combined_layers_subset@meta.data$spot <- rownames(combined_layers_subset@meta.data)
    combined_layers_subset@meta.data$sample <- paste(sapply(strsplit(combined_layers_subset@meta.data$spot, "_"), "[", 1))
    combined_layers_subset <- subset(combined_layers_subset,
                                     subset = spot %in% spot_neighbour_programs_df_flat$spot_of_interest)

    # Because the order of the rows to col transpose is the same, we can just set the col indices as rows
    setdiff(colnames(spot_neighbour_programs_df_flat_plot), spot_neighbour_programs_df_flat$spot_of_interest)
    spot_neighbour_programs_df_flat_mod <- spot_neighbour_programs_df_flat
    spot_neighbour_programs_df_flat_mod$clust <- "0"

    for (split_name in names(heatmap_col_order)){
      spot_neighbour_programs_df_flat_mod[heatmap_col_order[[split_name]], "clust"] <- split_name
    }
    combined_layers_subset@meta.data <- cbind(combined_layers_subset@meta.data,
                                              spot_neighbour_programs_df_flat_mod[9:ncol(spot_neighbour_programs_df_flat_mod)])
    cols_to_plot <- colnames(combined_layers_subset@meta.data)[7:(ncol(combined_layers_subset@meta.data)-1)]

    # 3.4.2 Create the dot plot ----
    pdf(paste(save_path, "Neighbourhood_dot_plot_Dendro_", dendro_split, "_Split_Custom.pdf", sep = ""), height = 14, width = 10)
    DotPlot(combined_layers_subset,
            features = cols_to_plot,
            group.by = "clust",
            scale = F) + coord_flip() + RotatedAxis()
    for (radii in unique(gsub(".*_", "", cols_to_plot))){
      cols_to_plot_subset <- cols_to_plot[grepl(paste("_", radii, sep = ""), cols_to_plot)]

      print(DotPlot(combined_layers_subset,
                    features = cols_to_plot_subset,
                    group.by = "clust",
                    scale = F) + coord_flip() + RotatedAxis())
    }
    dev.off()

    # 3.5 Save the spot id's per niche ----
    niche_spot_ids <- list()
    for (niche in names(heatmap_col_order)){
      spot_names <-  colnames(spot_neighbour_programs_df_flat_plot)
      spot_names_indicies <- heatmap_col_order[[niche]]
      spot_names <- spot_names[spot_names_indicies]

      niche_spot_ids[[niche]] <- spot_names
    }

    saveRDS(niche_spot_ids, file = paste(save_path, "Niche_spot_ids_split_custom_V2.rds", sep = ""))
    saveRDS(spot_neighbour_programs_df_flat_plot, file = paste(save_path, "Spot_neighbour_programs_df_flat_plot.rds", sep = ""))
  }
}

# 4. Gene Expression Analysis ----
# 4.1 Inside-Outside T-Cell Gene Expression Analysis ----
if (INSIDE_OUTSIDE_T_CELL){
  # 1.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "K", RANK, "/inside_outside/T_Cell/", sep = "")
  create_folder(paste(save_path, "diagnostic/", sep = ""))

  # 1.1 Get the T-cells, the spot coordinates from the 2 samples, and their intersection ----
  T_cell_so <- subset(x = combined_layers, subset = CD3D | CD3E | CD3G | CD8A | GCAR, slot = "counts")
  T_cell_spots <- rownames(T_cell_so@meta.data)

  # Get the spot names and merge them with the union_cell_coordinate dataframe
  T_cell_spots <- intersect(T_cell_spots, rownames(union_cell_coordinate))

  # 1.2 Subset the spot neighbours ----
  spot_neighbours_df_subset <- spot_neighbours_df[spot_neighbours_df$spot_of_interest %in% T_cell_spots,]

  # 1.3 Calculate the within and outside cells across radii ----
  combined_metadata_table_duplicate <- combined_metadata_table
  new_columns <- as.vector(unique(as.character(spot_neighbours_df_subset$r)))
  new_columns <- paste("dist", new_columns, sep = "_")
  combined_metadata_table_duplicate[, new_columns] <- "outside"

  # Set the neighbours of a spot to be "within", else make them "outside"
  for (radii in unique(spot_neighbours_df_subset$r)){
    # Get all the neighbours at a given radius
    radii_neighs <- unlist(spot_neighbours_df_subset[spot_neighbours_df_subset$r == radii, "neighbours"])

    # Set the column name to set cells to "within"
    radii_name <- paste("dist", radii, sep = "_")

    # Set the cells to "within"
    combined_metadata_table_duplicate[rownames(combined_metadata_table_duplicate) %in% radii_neighs, radii_name] <- "inside"
  }

  # 1.4 Merge the spots with specific gene counts ----
  count_table_subset <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[INSIDE_OUTSIDE_GENES$gene ,])))

  within_outside_df <- merge(combined_metadata_table_duplicate, count_table_subset, by = "row.names")
  colnames(within_outside_df)[1] <- "spot"

  within_outside_df <- merge(sample_metadata[, c("sample", "TREATMENT", "LOCATION")], within_outside_df,
                             by.x = "sample", by.y = "sample", all.x = TRUE)

  within_outside_df_spread <- within_outside_df %>%
    group_by(spot) %>%
    pivot_longer(c(dist_R0, dist_R1, dist_R2, dist_R3), names_to = "radius", values_to = "location") %>%
    pivot_longer(INSIDE_OUTSIDE_GENES$gene, names_to = "gene", values_to = "expression")

  # 1.5 Histogram of the fraction of gene+ bins in the inside vs outside bins ----
  histo_plot_data <- within_outside_df_spread %>%
    group_by(TREATMENT, radius, location, gene) %>%
    summarise(total_spots = n(), nr_gene_pos_bins = sum(expression > 0), gene_pos_bin_mean_exp = mean(expression[expression > 0])) %>%
    mutate(proportion_gene_pos_bins = nr_gene_pos_bins / total_spots,
           proportion_gene_pos_bin_mean_exp = gene_pos_bin_mean_exp / nrow(combined_metadata_table_duplicate)) %>%
    mutate(percent_gene_pos_bins = proportion_gene_pos_bins * 100,
           percent_gene_pos_bin_mean_exp = proportion_gene_pos_bin_mean_exp * 100) %>%
    ungroup()
  histo_plot_data[is.na(histo_plot_data)] <- 0
  histo_plot_data$split <- paste(histo_plot_data$TREATMENT, histo_plot_data$location)
  histo_plot_data$radius <- factor(histo_plot_data$radius, levels = c("dist_R0", "dist_R1", "dist_R2", "dist_R3"))

  # Modifications for plotting
  histo_plot_data$split <- sub("inside", "Inside", histo_plot_data$split)
  histo_plot_data$split <- sub("outside", "Outside", histo_plot_data$split)
  colnames(histo_plot_data)[12] <- "Split"

  histo_plot_data$radius <- sub("dist_R0", "R0", histo_plot_data$radius)
  histo_plot_data$radius <- sub("dist_R1", "R1", histo_plot_data$radius)
  histo_plot_data$radius <- sub("dist_R2", "R2", histo_plot_data$radius)
  histo_plot_data$radius <- sub("dist_R3", "R3", histo_plot_data$radius)
  histo_plot_data$radius <- factor(histo_plot_data$radius, levels = c("R0", "R1", "R2", "R3"))

  histo_plot_data$Split_V2 <- paste(histo_plot_data$location, histo_plot_data$radius, sep = "_")

  # 1.6 Run the stats ----
  save.image(paste(save_path, "T_Cell", "_IO_pre_stats_save.RData", sep = ""))

  stat_df <- histo_plot_data
  stat_df <- stat_df[!(stat_df$location == "inside" & stat_df$radius == "R3"),]
  stat_df <- stat_df[!(stat_df$location == "outside" & stat_df$radius == "R0"),]

  group_means <- stat_df %>%
    group_by(gene, Split) %>%
    summarise(group_mean = mean(proportion_gene_pos_bins), .groups = "drop")

  stats <- stat_df %>%
    group_by(gene, Split) %>%
    mutate(s = sd(proportion_gene_pos_bins)) %>%
    group_by(gene, Split) %>%
    filter(!any(s == 0)) %>%
    compare_means(proportion_gene_pos_bins ~ Split, .,
                  group.by = "gene", method = "t.test", var.equal = T)

  stats <- stats %>%
    left_join(group_means, by = c("gene", "group1" = "Split")) %>%
    rename(group_1_mean = group_mean) %>%
    left_join(group_means, by = c("gene", "group2" = "Split")) %>%
    rename(group_2_mean = group_mean)
  stats$group_mean_diff <- stats$group_1_mean - stats$group_2_mean

  # 1.7 Plotting ----
  # Get the range of the heatmap w/ p-values
  ht_col_range <- formatC(stats$p, format = "e", digits = 2)
  ht_col_range <- as.numeric(gsub(".*e[+|-]", "", ht_col_range))
  if (min(ht_col_range) != 0){
    if (min(ht_col_range) == max(ht_col_range)){
      col_fun <- colorRamp2(c(0, max(ht_col_range, na.rm = TRUE)+1), hcl_palette = "Blues", reverse = TRUE)
    }else{
      col_fun <- colorRamp2(c(0, max(ht_col_range, na.rm = TRUE)), hcl_palette = "Blues", reverse = TRUE)
    }
  }else{
    col_fun <- colorRamp2(range(ht_col_range, na.rm = TRUE), hcl_palette = "Blues", reverse = TRUE)
  }

  # Make the plots
  temp_df <- data.frame()
  write.csv(stats, paste(save_path, "stats_table.csv", sep = ""))
  write.csv(histo_plot_data, paste(save_path, "histogram_plot_proportions.csv", sep = ""))
  pdf(paste(save_path, "Inside_Outside_Gene_Expression_Histogram_T_Cell_Nominal_P.pdf", sep = ""), width = 40, height = 10)
  for (i in seq_len(nrow(INSIDE_OUTSIDE_GENES))){
    pair <- INSIDE_OUTSIDE_GENES$category[i]
    gene <- INSIDE_OUTSIDE_GENES$gene[i]

    # Subset the histogram data to the gene of interest
    temp_plot_data <- histo_plot_data[histo_plot_data$gene == gene,]

    # Create the histogram
    gp <- ggplot(temp_plot_data, aes(x = radius, y = proportion_gene_pos_bins, fill = Split_V2)) +
      geom_bar(stat = "identity", color = "black", position = "dodge", width = .8) + theme_classic2() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, colour="black"),
            axis.text.y=element_text(colour="black"), text = element_text(size = 20),
            legend.position = "bottom") +
      labs(x = "Distance From CD3D/E/G, CD8A, and GCAR Bin",
           y = "Proportion of Gene+ Spots\nof Total Spots Per Condition",
           title = paste(gene, pair, sep = ": "),
           subtitle = "Fraction of Adjacent and Distant Bins From CD3D/E/G, CD8A, and GCAR+ Cells") +
      scale_y_continuous(expand = c(0,0)) +
      facet_wrap( ~ TREATMENT + location, ncol = length(unique(temp_plot_data$Split)), nrow = 1) +
      scale_fill_manual(values = c("#2b465e", "#607486", "#95a2af", "#cad1d7",
                                   "#d2c691", "#ddd4ad", "#e9e2c8", "#f4f1e3")) +
      guides(fill = guide_legend(ncol = length(unique(temp_plot_data$Split)), byrow = TRUE))

    # Create the P-Value Heatmap in the same figure
    temp_stats_data <- as.data.frame(stats[stats$gene == gene, c(3, 4, 5)])
    temp_stats_data_rev <- temp_stats_data
    colnames(temp_stats_data_rev)[1:2] <- c("group2", "group1")
    temp_stats_data <- rbind(temp_stats_data, temp_stats_data_rev)
    setDT(temp_stats_data)

    temp_stats_data_spread <- dcast(temp_stats_data, group1 ~ group2, value.var = "p")
    temp_stats_data_spread <- as.data.frame(temp_stats_data_spread)
    rownames(temp_stats_data_spread) <- temp_stats_data_spread$group1
    temp_stats_data_spread <- temp_stats_data_spread[-1]
    temp_stats_data_spread[is.na(temp_stats_data_spread)] <- 0
    temp_stats_data_spread <- as.matrix(temp_stats_data_spread)
    temp_stats_data_spread <- formatC(temp_stats_data_spread, format = "e", digits = 2)

    temp_stats_data_spread <- gsub(".*e[+|-]", "", temp_stats_data_spread)

    temp_stats_data_spread <- as.data.frame(temp_stats_data_spread)
    stats_row_name_store <- rownames(temp_stats_data_spread)

    temp_stats_data_spread <- sapply(temp_stats_data_spread, as.numeric)
    rownames(temp_stats_data_spread) <- stats_row_name_store

    temp_stats_data_spread <- as.matrix(temp_stats_data_spread)
    diag(temp_stats_data_spread) <- NA

    data_cols <- colnames(temp_stats_data_spread)

    # Make the heatmap
    htp_name <- paste(gene, ": ", pair, ":\nAbsolute Value of The P-Value Scientific Notation Exponent\n[T.test]", sep = "")
    htp <- Heatmap(temp_stats_data_spread,
                   cluster_columns = FALSE, cluster_rows = FALSE,
                   show_row_names = TRUE, show_column_names = TRUE,
                   col = col_fun, na_col = "black")
    htp_grb <- grid.grabExpr(draw(htp, column_title=htp_name, padding = unit(c(3, 3, 3, 3), "cm")))

    grid.arrange(gp, htp_grb, ncol = 2)
  }
  dev.off()

  # 1.8 Save the Rdata image, so plotting changes can be easily performed ----
  save.image(paste(save_path, "T_Cell", "_neighbourhood_save.RData", sep = ""))
}

# 4.2 Inside-Outside Niche Gene Expression Analysis ----
if (INSIDE_OUTSIDE_NICHE){
  # 2.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste0(SAVE_FOLDER, "K", RANK, "/inside_outside/niche/")
  create_folder(save_path)
  create_folder(paste0(save_path, "diagnostic/"))

  # 2.1 Load the data ----
  var_save_root <- paste0(SAVE_FOLDER, "K", RANK, "/neighbourhood_analysis/Programs/GEP_5/")
  niche_spot_ids <- readRDS(file = paste0(var_save_root, "Niche_spot_ids_split_custom_V2.rds"))
  spot_neighbour_programs_df_flat_plot <- readRDS(file = paste0(var_save_root, "Spot_neighbour_programs_df_flat_plot.rds"))

  # 2.2 Inside outside analysis ----
  for (niche in names(niche_spot_ids)){
    # 2.2.1 Select the appropriate spots ----
    niche_name_underscore <- gsub(",", "_", niche)
    niche_save_path <- paste0(save_path, niche_name_underscore, "/")
    create_folder(niche_save_path)
    niche_spots <- niche_spot_ids[[niche]]

    # Get the spot names and merge them with the union_cell_coordinate dataframe
    niche_spots <- intersect(niche_spots, rownames(union_cell_coordinate))

    # 2.2.2 Subset the spot neighbours ----
    spot_neighbours_df_subset <- spot_neighbours_df[spot_neighbours_df$spot_of_interest %in% niche_spots,]

    # 2.2.3 Calculate the inside and outside cells across radii ----
    combined_metadata_table_duplicate <- combined_metadata_table
    new_columns <- as.vector(unique(as.character(spot_neighbours_df_subset$r)))
    new_columns <- paste("dist", new_columns, sep = "_")
    combined_metadata_table_duplicate[, new_columns] <- "outside"

    # Set the neighbours of a spot to be "inside", else make them "outside"
    for (radii in unique(spot_neighbours_df_subset$r)){
      # Get all the neighbours at a given radius
      radii_neighs <- unlist(spot_neighbours_df_subset[spot_neighbours_df_subset$r == radii, "neighbours"])

      # Set the column name to set cells to "inside"
      radii_name <- paste("dist", radii, sep = "_")

      # Set the cells to "inside"
      combined_metadata_table_duplicate[rownames(combined_metadata_table_duplicate) %in% radii_neighs, radii_name] <- "inside"
    }

    # Remove the other niche spots from this analysis
    niche_spots_to_remove <- unlist(niche_spot_ids[names(niche_spot_ids) != niche])
    combined_metadata_table_duplicate <- combined_metadata_table_duplicate[!rownames(combined_metadata_table_duplicate)
                                                                           %in% niche_spots_to_remove,]

    # 2.2.4 Make a linked dot plot to check that no inside converts to an outside across radii ----
    ohe <- data.frame(soi = rownames(combined_metadata_table_duplicate))
    rownames(ohe) <- ohe$soi
    ohe <- cbind(ohe, combined_metadata_table_duplicate[,5:8])
    ohe_spread <- ohe %>%
      pivot_longer(!soi, names_to = "radii", values_to = "location")
    ohe_spread$radii <- factor(ohe_spread$radii, levels = c("dist_R0", "dist_R1", "dist_R2", "dist_R3"))
    ohe_spread$location[ohe_spread$location == "inside"] <- 0
    ohe_spread$location[ohe_spread$location == "outside"] <- 1

    options(bitmapType = 'cairo')
    png(paste(save_path, "diagnostic/", "Inside_Outside_Diagnostic_", sub(",", "_", niche), "_subset.png", sep = ""), height = 1000, width = 1000)
    io_diag <- ggplot(ohe_spread, aes(x = radii, y = location)) +
      geom_line(aes(group = soi), color = "gray") + geom_point(shape = 16, color = "black", fill = "#69b3a2", size = 6) +
      theme_classic2() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, colour="black"),
            axis.text.y=element_text(colour="black"), text = element_text(size = 20),
            legend.position = "bottom") +
      labs(x = "Radius",
           y = paste(sub(",", "_", niche), "Spots:", sep = ""),
           title = paste("Number of Spots" , sub(",", "_", niche), sep = ": "),
           subtitle = paste("[Inside = 0, Outside = 1]"))
    print(io_diag)
    dev.off()

    # 2.2.5 Merge the spots with specific gene counts ----
    count_table_subset <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[INSIDE_OUTSIDE_GENES$gene ,])))

    within_outside_df <- merge(combined_metadata_table_duplicate, count_table_subset, by = "row.names")
    colnames(within_outside_df)[1] <- "spot"

    within_outside_df <- merge(sample_metadata[, c("sample", "TREATMENT", "LOCATION")], within_outside_df,
                               by.x = "sample", by.y = "sample", all.x = TRUE)

    # 2.2.6 Spread the data ----
    within_outside_df_spread <- within_outside_df %>%
      group_by(spot) %>%
      pivot_longer(c(dist_R0, dist_R1, dist_R2, dist_R3), names_to = "radius", values_to = "location") %>%
      pivot_longer(INSIDE_OUTSIDE_GENES$gene, names_to = "gene", values_to = "expression")

    # 2.2.7 Histogram of the fraction of PDL1+ bins in the inside vs outside bins ----
    histo_plot_data <- within_outside_df_spread %>%
      group_by(TREATMENT, radius, location, gene) %>%
      summarise(total_spots = n(), nr_gene_pos_bins = sum(expression > 0), gene_pos_bin_mean_exp = mean(expression[expression > 0])) %>%
      mutate(proportion_gene_pos_bins = nr_gene_pos_bins / total_spots,
             proportion_gene_pos_bin_mean_exp = gene_pos_bin_mean_exp / nrow(combined_metadata_table_duplicate)) %>%
      mutate(percent_gene_pos_bins = proportion_gene_pos_bins * 100,
             percent_gene_pos_bin_mean_exp = proportion_gene_pos_bin_mean_exp * 100) %>%
      ungroup()
    histo_plot_data[is.na(histo_plot_data)] <- 0
    histo_plot_data$split <- paste(histo_plot_data$TREATMENT, histo_plot_data$location)
    histo_plot_data$radius <- factor(histo_plot_data$radius, levels = c("dist_R0", "dist_R1", "dist_R2", "dist_R3"))

    # Modifications for plotting
    histo_plot_data$split <- sub("inside", "Inside", histo_plot_data$split)
    histo_plot_data$split <- sub("outside", "Outside", histo_plot_data$split)
    colnames(histo_plot_data)[12] <- "Split"

    histo_plot_data$radius <- sub("dist_R0", "R0", histo_plot_data$radius)
    histo_plot_data$radius <- sub("dist_R1", "R1", histo_plot_data$radius)
    histo_plot_data$radius <- sub("dist_R2", "R2", histo_plot_data$radius)
    histo_plot_data$radius <- sub("dist_R3", "R3", histo_plot_data$radius)
    histo_plot_data$radius <- factor(histo_plot_data$radius, levels = c("R0", "R1", "R2", "R3"))

    histo_plot_data$Split_V2 <- paste(histo_plot_data$location, histo_plot_data$radius, sep = "_")

    # 2.2.8 Run the stats ----
    stats <- histo_plot_data %>%
      group_by(gene, Split) %>%
      mutate(s = sd(proportion_gene_pos_bins)) %>%
      group_by(gene, Split) %>%
      filter(!any(s == 0)) %>%
      compare_means(proportion_gene_pos_bins ~ Split, .,
                    group.by = "gene", method = "t.test")

    # 2.2.9 Plotting ----
    # Get the range of the heatmap w/ p-values
    ht_col_range <- formatC(stats$p, format = "e", digits = 2)
    ht_col_range <- as.numeric(gsub(".*e[+|-]", "", ht_col_range))
    if (min(ht_col_range) != 0){
      if (min(ht_col_range) == max(ht_col_range)){
        col_fun <- colorRamp2(c(0, max(ht_col_range, na.rm = TRUE)+1), hcl_palette = "Blues", reverse = TRUE)
      }else{
        col_fun <- colorRamp2(c(0, max(ht_col_range, na.rm = TRUE)), hcl_palette = "Blues", reverse = TRUE)
      }
    }else{
      col_fun <- colorRamp2(range(ht_col_range, na.rm = TRUE), hcl_palette = "Blues", reverse = TRUE)
    }

    # Make the plots
    temp_df <- data.frame()
    write.csv(stats, paste(niche_save_path, "stats_table.csv", sep = ""))
    write.csv(histo_plot_data, paste(niche_save_path, "histogram_plot_proportions.csv", sep = ""))
    pdf(paste(niche_save_path, "Inside_Outside_Gene_Expression_Histogram_", niche_name_underscore, ".pdf", sep = ""), width = 30, height = 10)
    for (i in seq_len(nrow(INSIDE_OUTSIDE_GENES))){
      pair <- INSIDE_OUTSIDE_GENES$category[i]
      gene <- INSIDE_OUTSIDE_GENES$gene[i]

      if (gene %in% stats$gene){
        # Subset the histogram data to the gene of interest
        temp_plot_data <- histo_plot_data[histo_plot_data$gene == gene,]

        # Create the histogram
        gp <- ggplot(temp_plot_data, aes(x = radius, y = proportion_gene_pos_bins, fill = Split_V2)) +
          geom_bar(stat = "identity", color = "black", position = "dodge", width = .8) + theme_classic2() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, colour="black"),
                axis.text.y=element_text(colour="black"), text = element_text(size = 20),
                legend.position = "bottom") +
          labs(x = paste("Distance From GEP_5, Niche: ", niche, sep = ""),
               y = "Proportion of Gene+ Spots\nof Total Spots Per Condition",
               title = paste(gene, pair, sep = ": "),
               subtitle = paste("Fraction of Adjacent and Distant Bins GEP_5, Niche: ", niche, sep = "")) +
          scale_y_continuous(expand = c(0,0)) + facet_wrap( ~ TREATMENT + location, ncol = length(unique(temp_plot_data$Split)), nrow = 1) +
          scale_fill_manual(values = c("#2b465e", "#607486", "#95a2af", "#cad1d7",
                                       "#d2c691", "#ddd4ad", "#e9e2c8", "#f4f1e3")) +
          guides(fill = guide_legend(ncol = length(unique(temp_plot_data$Split)), byrow = TRUE))

        # Create the P-Value Heatmap in the same figure
        temp_stats_data <- as.data.frame(stats[stats$gene == gene, c(3, 4, 5)])
        temp_stats_data_rev <- temp_stats_data
        colnames(temp_stats_data_rev)[1:2] <- c("group2", "group1")
        temp_stats_data <- rbind(temp_stats_data, temp_stats_data_rev)
        setDT(temp_stats_data)

        temp_stats_data_spread <- dcast(temp_stats_data, group1 ~ group2, value.var = "p")
        temp_stats_data_spread <- as.data.frame(temp_stats_data_spread)
        rownames(temp_stats_data_spread) <- temp_stats_data_spread$group1
        temp_stats_data_spread <- temp_stats_data_spread[-1]
        temp_stats_data_spread[is.na(temp_stats_data_spread)] <- 0
        temp_stats_data_spread <- as.matrix(temp_stats_data_spread)
        temp_stats_data_spread <- formatC(temp_stats_data_spread, format = "e", digits = 2)

        temp_stats_data_spread <- gsub(".*e[+|-]", "", temp_stats_data_spread)

        temp_stats_data_spread <- as.data.frame(temp_stats_data_spread)
        stats_row_name_store <- rownames(temp_stats_data_spread)

        temp_stats_data_spread <- sapply(temp_stats_data_spread, as.numeric)
        rownames(temp_stats_data_spread) <- stats_row_name_store

        temp_stats_data_spread <- as.matrix(temp_stats_data_spread)
        diag(temp_stats_data_spread) <- NA

        data_cols <- colnames(temp_stats_data_spread)

        # Make the heatmap
        htp_name <- paste(gene, ": ", pair, ":\nAbsolute Value of The P-Value Scientific Notation Exponent\n[T.test]", sep = "")
        htp <- Heatmap(temp_stats_data_spread,
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       show_row_names = TRUE, show_column_names = TRUE,
                       col = col_fun, na_col = "black")
        htp_grb <- grid.grabExpr(draw(htp, column_title=htp_name, padding = unit(c(3, 3, 3, 3), "cm")))

        grid.arrange(gp, htp_grb, ncol = 2)
      }
    }
    dev.off()
  }

  # 2.3 Contingency table ----
  histo_plot_df <- data.frame()

  for (niche in names(niche_spot_ids)){
    # 2.3.1 Select the appropriate spots ----
    niche_name_underscore <- gsub(",", "_", niche)
    niche_save_path <- paste0(save_path, niche_name_underscore, "/")
    create_folder(niche_save_path)
    niche_spots <- niche_spot_ids[[niche]]

    # Get the spot names and merge them with the union_cell_coordinate dataframe
    niche_spots <- intersect(niche_spots, rownames(union_cell_coordinate))

    # 2.3.2 Subset the spot neighbours ----
    spot_neighbours_df_subset <- spot_neighbours_df[spot_neighbours_df$spot_of_interest %in% niche_spots,]

    # 2.3.3 Calculate the inside and outside cells across radii ----
    combined_metadata_table_duplicate <- combined_metadata_table
    new_columns <- as.vector(unique(as.character(spot_neighbours_df_subset$r)))
    new_columns <- paste("dist", new_columns, sep = "_")
    combined_metadata_table_duplicate[, new_columns] <- "outside"

    # Set the neighbours of a spot to be "inside", else make them "outside"
    for (radii in unique(spot_neighbours_df_subset$r)){
      # Get all the neighbours at a given radius
      radii_neighs <- unlist(spot_neighbours_df_subset[spot_neighbours_df_subset$r == radii, "neighbours"])

      # Set the column name to set cells to "inside"
      radii_name <- paste("dist", radii, sep = "_")

      # Set the cells to "inside"
      combined_metadata_table_duplicate[rownames(combined_metadata_table_duplicate) %in% radii_neighs, radii_name] <- "inside"
    }

    # Remove the other niche spots from this analysis
    niche_spots_to_remove <- unlist(niche_spot_ids[names(niche_spot_ids) != niche])
    combined_metadata_table_duplicate <- combined_metadata_table_duplicate[!rownames(combined_metadata_table_duplicate)
                                                                           %in% niche_spots_to_remove,]

    # 2.3.4 Merge the spots with specific gene counts ----
    count_table_subset <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[INSIDE_OUTSIDE_GENES$gene ,])))

    within_outside_df <- merge(combined_metadata_table_duplicate, count_table_subset, by = "row.names")
    colnames(within_outside_df)[1] <- "spot"

    within_outside_df <- merge(sample_metadata[, c("sample", "TREATMENT", "LOCATION")], within_outside_df,
                               by.x = "sample", by.y = "sample", all.x = TRUE)

    # 2.3.5 Spread the data ----
    within_outside_df_spread <- within_outside_df %>%
      group_by(spot) %>%
      pivot_longer(c(dist_R0, dist_R1, dist_R2, dist_R3), names_to = "radius", values_to = "location") %>%
      pivot_longer(INSIDE_OUTSIDE_GENES$gene, names_to = "gene", values_to = "expression")

    # 2.3.6 Histogram of the fraction of PDL1+ bins in the inside vs outside bins ----
    histo_plot_data <- within_outside_df_spread %>%
      group_by(TREATMENT, radius, location, gene) %>%
      summarise(total_spots = n(), nr_gene_pos_bins = sum(expression > 0), gene_pos_bin_mean_exp = mean(expression[expression > 0])) %>%
      mutate(proportion_gene_pos_bins = nr_gene_pos_bins / total_spots,
             proportion_gene_pos_bin_mean_exp = gene_pos_bin_mean_exp / nrow(combined_metadata_table_duplicate)) %>%
      mutate(percent_gene_pos_bins = proportion_gene_pos_bins * 100,
             percent_gene_pos_bin_mean_exp = proportion_gene_pos_bin_mean_exp * 100) %>%
      ungroup()
    histo_plot_data[is.na(histo_plot_data)] <- 0
    histo_plot_data$split <- paste(histo_plot_data$TREATMENT, histo_plot_data$location)
    histo_plot_data$radius <- factor(histo_plot_data$radius, levels = c("dist_R0", "dist_R1", "dist_R2", "dist_R3"))

    # Modifications for plotting
    histo_plot_data$split <- sub("inside", "Inside", histo_plot_data$split)
    histo_plot_data$split <- sub("outside", "Outside", histo_plot_data$split)
    colnames(histo_plot_data)[12] <- "Split"

    histo_plot_data$radius <- sub("dist_R0", "R0", histo_plot_data$radius)
    histo_plot_data$radius <- sub("dist_R1", "R1", histo_plot_data$radius)
    histo_plot_data$radius <- sub("dist_R2", "R2", histo_plot_data$radius)
    histo_plot_data$radius <- sub("dist_R3", "R3", histo_plot_data$radius)
    histo_plot_data$radius <- factor(histo_plot_data$radius, levels = c("R0", "R1", "R2", "R3"))

    histo_plot_data$Split_V2 <- paste(histo_plot_data$location, histo_plot_data$radius, sep = "_")

    histo_plot_data <- histo_plot_data[histo_plot_data$radius == "R0",]
    histo_plot_data$niche <- niche

    histo_plot_df <- rbind(histo_plot_df, histo_plot_data)
  }

  # 2.4 Data Manipulation ----
  histo_plot_df_final <- histo_plot_df[,c(14,3,4,6,5)]
  colnames(histo_plot_df_final)[4] <- "positive"
  histo_plot_df_final$negative <- histo_plot_df_final$total_spots - histo_plot_df_final$positive
  histo_plot_df_final <- histo_plot_df_final[,c(1:4,6,5)]
  histo_plot_df_final <- histo_plot_df_final[histo_plot_df_final$location == "inside",]
  histo_plot_df_final$niche <- gsub(",", "_", histo_plot_df_final$niche)
  write.table(histo_plot_df_final, paste0(save_path, "R0_Final_Table_Custom_Split.csv"),
              quote = FALSE, sep = ",", row.names = FALSE)

  # 2.5 Plotting ----
  col_fun <- colorRamp2(c(0, 10), hcl_palette = "Blues", reverse = TRUE)

  f_pw_df <- data.frame()
  pdf(paste(save_path, "R0_Gene_Expression_Fisher_Exact_Nominal_P.pdf", sep = ""), width = 15, height = 15)
  for (i in seq_len(nrow(INSIDE_OUTSIDE_GENES))){
    pair <- INSIDE_OUTSIDE_GENES$category[i]
    gene <- INSIDE_OUTSIDE_GENES$gene[i]

    # Subset the histogram data to the gene of interest ----
    histo_plot_df_subset <- histo_plot_df_final[histo_plot_df_final$gene == gene,]
    histo_plot_df_subset <- as.data.frame(histo_plot_df_subset)
    rownames(histo_plot_df_subset) <- histo_plot_df_subset$niche
    histo_plot_df_subset <- histo_plot_df_subset[,c(4,5)]

    temp <- fisher.test(histo_plot_df_subset, workspace = 2e9, simulate.p.value = T)
    temp_ph <- pairwise_fisher_test(histo_plot_df_subset, p.adjust.method = "fdr")

    # Create the P-Value Heatmap figure ----
    temp_stats_data <- temp_ph
    temp_stats_data_rev <- temp_stats_data
    colnames(temp_stats_data_rev)[1:2] <- c("group2", "group1")
    temp_stats_data <- rbind(temp_stats_data, temp_stats_data_rev)
    setDT(temp_stats_data)

    # Store the results in a dataframe ----
    temp_temp_ph <- temp_stats_data
    temp_temp_ph$gene <- gene
    f_pw_df <- rbind(f_pw_df, temp_temp_ph)

    temp_stats_data_spread <- dcast(temp_stats_data, group1 ~ group2, value.var = "p")
    temp_stats_data_spread <- as.data.frame(temp_stats_data_spread)
    rownames(temp_stats_data_spread) <- temp_stats_data_spread$group1
    temp_stats_data_spread <- temp_stats_data_spread[-1]
    temp_stats_data_spread[is.na(temp_stats_data_spread)] <- 0
    temp_stats_data_spread <- as.matrix(temp_stats_data_spread)
    temp_stats_data_spread <- formatC(temp_stats_data_spread, format = "e", digits = 2)

    temp_stats_data_spread <- gsub(".*e[+|-]", "", temp_stats_data_spread)

    temp_stats_data_spread <- as.data.frame(temp_stats_data_spread)
    stats_row_name_store <- rownames(temp_stats_data_spread)

    temp_stats_data_spread <- sapply(temp_stats_data_spread, as.numeric)
    rownames(temp_stats_data_spread) <- stats_row_name_store

    temp_stats_data_spread <- as.matrix(temp_stats_data_spread)
    diag(temp_stats_data_spread) <- NA

    # Make the heatmap ----
    if (min(range(temp_stats_data_spread, na.rm = TRUE)) == max(range(temp_stats_data_spread, na.rm = TRUE))){
      col_fun <- colorRamp2(c(min(temp_stats_data_spread, na.rm = TRUE), min(temp_stats_data_spread, na.rm = TRUE) + 1), hcl_palette = "Blues", reverse = TRUE)
    } else{col_fun <- colorRamp2(range(temp_stats_data_spread, na.rm = TRUE), hcl_palette = "Blues", reverse = TRUE)}

    htp_name <- paste(gene, "- ", pair, ":",
                      "\nAbsolute Value of The P-Value Scientific Notation Exponent [Capped at 10] \n[Post-Hoc Pairwise Fisher Exact P-Value]",
                      "\nFisher Exact P-Value: ", temp$p.value, sep = "")
    htp <- Heatmap(temp_stats_data_spread,
                   cluster_columns = FALSE, cluster_rows = FALSE,
                   show_row_names = TRUE, show_column_names = TRUE,
                   col = col_fun, na_col = "black")

    draw(htp, column_title = htp_name, padding = unit(c(3, 3, 3, 3), "cm"))
  }
  dev.off()
  write.table(f_pw_df, paste(save_path, "R0_Fischer_Pairwise_Stat_Table_Custom_Split.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
}

# 4.3 Differential ALDEx2 gene expression ----
if (DEG_ALDEX_ANALYSIS){
  # 3.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "K", RANK, "/differential_expression/aldex/", sep = "")
  create_folder(save_path)
  set.seed(123)

  usage_threshold <- 0.1

  # 3.1 Indicate the tissue and pre-/post-CAR infusion GEP+ bins----
  combined_layers_duplicate <- combined_layers
  combined_layers_duplicate@meta.data$pre_post[combined_layers_duplicate@meta.data$sample == "SPSC_biopsy1"] <- "Pre_CAR_Dose"
  combined_layers_duplicate@meta.data$pre_post[combined_layers_duplicate@meta.data$sample == "SPSC_biopsy2"] <- "Post_CAR_Dose"

  # 3.2 Create a list to store the results ----
  pre_post_deg_list <- vector("list", length = ncol(usage_norm)+2)
  names(pre_post_deg_list) <- c("all_bins_individual",
                                "all_bins_pseudobulk",
                                paste0("Program_", seq_len(ncol(usage_norm))))

  # 3.3 Function to generate pseudobulk's of the data ----
  make_random_pseudobulks <- function(count_mat, sample_col = "sample",
                                      cells_per_bulk = NA, samples_bulk = 7,
                                      equal_reads_per_bulk = NA){
    if (is.na(cells_per_bulk) & is.na(samples_bulk) & is.na(equal_reads_per_bulk)) {
      stop("One of cells_per_bulk, samples_bulk, or equal_reads_per_bulk must be set.")
    } else if (!is.na(cells_per_bulk) & !is.na(samples_bulk) & !is.na(equal_reads_per_bulk)) {
      stop("Set only one of cells_per_bulk, samples_bulk, equal_reads_per_bulk, but not two or more.")
    }

    # shuffle rows within each sample to avoid ordering bias
    count_mat <- count_mat %>%
      group_by(!!sym(sample_col)) %>%
      sample_frac() %>%
      ungroup()

    if (!is.na(cells_per_bulk)){
      count_mat <- count_mat %>%
        group_by(!!sym(sample_col)) %>%
        mutate(pseudobulk_id = paste0(!!sym(sample_col), "_",
                                      ceiling(row_number() / cells_per_bulk))) %>%
        ungroup()
    }

    if (!is.na(samples_bulk)){
      count_mat <- count_mat %>%
        group_by(!!sym(sample_col)) %>%
        mutate(pseudobulk_id = paste0(!!sym(sample_col), "_",
                                      ntile(row_number(), samples_bulk))) %>%
        ungroup()
    }

    if (!is.na(equal_reads_per_bulk)){
      count_mat <- count_mat %>%
        group_by(!!sym(sample_col)) %>%
        mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
        arrange(desc(total_counts)) %>%
        mutate(pseudobulk_id = paste0(!!sym(sample_col), "_",
                                      (row_number() - 1) %% samples_bulk + 1)) %>%
        ungroup()
    }

    # sum counts per pseudobulk
    bulked <- count_mat %>%
      group_by(pseudobulk_id) %>%
      summarise(across(where(is.numeric), sum), .groups = "drop")

    bulked <- as.data.frame(bulked)
    rownames(bulked) <- bulked$pseudobulk_id
    bulked$pseudobulk_id <- NULL

    # genes as rows, pseudobulks as columns
    pseudobulk_mat <- t(as.matrix(bulked))
    return(pseudobulk_mat)
  }

  make_count_balanced_pseudobulks <- function(count_mat, sample_col = "sample", target_counts = 1e6){
    # count_mat: data.frame with sample_col + numeric gene counts per row (spots)
    # target_counts: desired total counts per pseudobulk

    # ensure sample column is character
    count_mat[[sample_col]] <- as.character(count_mat[[sample_col]])

    pseudobulk_list <- list()

    for (s in unique(count_mat[[sample_col]])) {
      df <- count_mat[count_mat[[sample_col]] == s, ]

      # calculate total counts per spot
      gene_cols <- which(sapply(df, is.numeric))
      df$total_counts <- rowSums(df[, gene_cols[!grepl("usage", gene_cols)], drop = FALSE])

      # randomize spot order to avoid spatial/ordering bias
      df <- df[sample(nrow(df)), ]

      # cumulative binning: assign spots until target_counts is reached
      cum_counts <- cumsum(df$total_counts)
      bin_ids <- ceiling(cum_counts / target_counts)
      df$pseudobulk_id <- paste0(s, "_bulk", bin_ids)

      # sum counts per pseudobulk
      bulked <- df %>%
        group_by(pseudobulk_id) %>%
        summarise(across(all_of(colnames(df)[gene_cols]), sum),
                  n_spots = n(),
                  .groups = "drop")

      pseudobulk_list[[s]] <- bulked
    }

    # merge across samples
    bulked_all <- bind_rows(pseudobulk_list)

    # reshape into matrix (genes x pseudobulks)
    bulked_all <- as.data.frame(bulked_all)
    rownames(bulked_all) <- bulked_all$pseudobulk_id
    bulked_all$pseudobulk_id <- NULL
    pseudobulk_mat <- t(as.matrix(bulked_all))

    return(pseudobulk_mat)
  }

  # 3.4 Compute the DEG between the pre and post-CAR bins ----
  # 3.4.1 Compute the DEG between the pre and post-CAR individual bins ----
  count_mat <- as.data.frame(combined_layers_duplicate@assays$Spatial.024um$counts)

  # Had to run twice because the first removal means some more need to be removed
  count_mat_t <- t(count_mat)
  count_mat_subset <- cmultRepl(count_mat_t, label = 0, method = "CZM", suppress.print = TRUE)
  count_mat_subset <- t(count_mat_subset)
  subset_cols <- colnames(count_mat_subset)
  subset_rows <- rownames(count_mat_subset)
  count_mat <- count_mat[,colnames(count_mat) %in% subset_cols]
  count_mat <- count_mat[rownames(count_mat) %in% subset_rows,]

  count_mat_t <- t(count_mat)
  count_mat_subset <- cmultRepl(count_mat_t, label = 0, method = "CZM", suppress.print = TRUE)
  count_mat_subset <- t(count_mat_subset)
  subset_cols <- colnames(count_mat_subset)
  subset_rows <- rownames(count_mat_subset)
  count_mat <- count_mat[,colnames(count_mat) %in% subset_cols]
  count_mat <- count_mat[rownames(count_mat) %in% subset_rows,]

  count_mat <- as.data.frame(count_mat)
  conds <- gsub(paste0("(", pattern_regex, ").*"), "\\1", colnames(count_mat))
  conds <- ifelse(conds == "SPSC_biopsy1", "Pre_CAR_Dose", "Post_CAR_Dose")

  aldex_res <- aldex(count_mat, conds, mc.samples=128,
                     test="t", effect=TRUE, include.sample.summary=FALSE, denom="lvha",
                     verbose=FALSE, paired.test=FALSE, useMC = T)

  saveRDS(aldex_res, file = paste0(save_path, "aldex_res_all_bins_individual.rds"))

  group_counts <- table(conds)
  aldex_res$gene <- rownames(aldex_res)
  aldex_res$Pre_CAR_Dose <- group_counts[2]
  aldex_res$Post_CAR_Dose <- group_counts[1]
  pre_post_deg_list[["all_bins_individual"]] <- aldex_res

  # 3.4.2 Compute the DEG between the pre and post-CAR pseudobulked bins ----
  count_mat <- as.data.frame(t(combined_layers_duplicate@assays$Spatial.024um$counts))
  count_mat$pre_post <- c(combined_layers_duplicate$pre_post)

  pseudobulk_mat <- make_random_pseudobulks(count_mat, sample_col = "pre_post",
                                            cells_per_bulk = 20, samples_bulk = NA, equal_reads_per_bulk = NA)

  pseudobulk_mat_t <- t(pseudobulk_mat)
  pseudobulk_mat_subset <- cmultRepl(pseudobulk_mat_t, label = 0, method = "CZM", suppress.print = TRUE)
  pseudobulk_mat_subset <- t(pseudobulk_mat_subset)
  subset_cols <- colnames(pseudobulk_mat_subset)
  subset_rows <- rownames(pseudobulk_mat_subset)
  pseudobulk_mat <- pseudobulk_mat[,colnames(pseudobulk_mat) %in% subset_cols]
  pseudobulk_mat <- pseudobulk_mat[rownames(pseudobulk_mat) %in% subset_rows,]

  # Continue with the analysis
  pseudobulk_mat <- as.data.frame(pseudobulk_mat)
  conds <- gsub("_[^_]+$", "", colnames(pseudobulk_mat))

  aldex_res <- aldex(pseudobulk_mat, conds, mc.samples=128,
                     test="t", effect=TRUE, include.sample.summary=FALSE, denom="lvha",
                     verbose=FALSE, paired.test=FALSE, useMC = T)

  saveRDS(aldex_res, file = paste0(save_path, "aldex_res_all_bins_pseudobulk.rds"))

  group_counts <- table(conds)
  aldex_res$gene <- rownames(aldex_res)
  aldex_res$Pre_CAR_Dose <- group_counts[2]
  aldex_res$Post_CAR_Dose <- group_counts[1]
  pre_post_deg_list[["all_bins_pseudobulk"]] <- aldex_res

  # 3.5 Compute the DEG between the pre and post-CAR bins subset by program+ ----
  for (program in seq_len(ncol(usage_norm))){
    print(program)
    # Identify the GEP+ cells
    gep_pos_spots <- rownames(usage_norm[usage_norm[, program] >= usage_threshold,])

    if (length(gep_pos_spots) >= 120){
      # Subset and normalize
      subset_bins <- subset(combined_layers_duplicate, cells = gep_pos_spots)
      subset_bins_count_mat <- as.data.frame(t(subset_bins@assays$Spatial.024um$counts))
      subset_bins_count_mat$pre_post <- c(subset_bins$pre_post)
      subset_bins_count_mat$usage <- usage_norm[rownames(usage_norm) %in% rownames(subset_bins_count_mat), program]

      group_counts <- table(subset_bins@meta.data$pre_post)
      sufficient_groups <- all(c("Pre_CAR_Dose", "Post_CAR_Dose") %in% names(group_counts)) && all(group_counts[c("Pre_CAR_Dose", "Post_CAR_Dose")] >= 60)

      # Check if all pre post groups have >= 60 cells
      if (sufficient_groups){
        pseudobulk_mat <- make_random_pseudobulks(subset_bins_count_mat, sample_col = "pre_post", cells_per_bulk = 20, samples_bulk = NA, equal_reads_per_bulk = NA)

        # pseudobulk_mat <- make_count_balanced_pseudobulks(subset_bins_count_mat, sample_col = "pre_post")
        temp <- pseudobulk_mat[c((nrow(pseudobulk_mat)-1):nrow(pseudobulk_mat)),]
        pseudobulk_mat <- pseudobulk_mat[1:(nrow(pseudobulk_mat)-2),]

        # Aldex2 throws an error when an internal function call removes particular rows and columns
        # Here, these columns and rows are being removed
        pseudobulk_mat_t <- t(pseudobulk_mat)
        pseudobulk_mat_subset <- cmultRepl(pseudobulk_mat_t, label = 0, method = "CZM", suppress.print = TRUE)
        pseudobulk_mat_subset <- t(pseudobulk_mat_subset)
        subset_cols <- colnames(pseudobulk_mat_subset)
        subset_rows <- rownames(pseudobulk_mat_subset)
        pseudobulk_mat <- pseudobulk_mat[,colnames(pseudobulk_mat) %in% subset_cols]
        pseudobulk_mat <- pseudobulk_mat[rownames(pseudobulk_mat) %in% subset_rows,]

        # Continue with the analysis
        pseudobulk_mat <- as.data.frame(pseudobulk_mat)
        conds <- gsub("_[^_]+$", "", colnames(pseudobulk_mat))

        if (length(unique(conds))>= 2){
          if (table(conds)[[1]] >= 2 & table(conds)[[2]] >= 2){
            aldex_res <- aldex(pseudobulk_mat, conds, mc.samples=128,
                               test="t", effect=TRUE, include.sample.summary=FALSE, denom="lvha",
                               verbose=FALSE, paired.test=FALSE, useMC = T)

            aldex_res$gene <- rownames(aldex_res)
            aldex_res$Pre_CAR_Dose <- group_counts[2]
            aldex_res$Post_CAR_Dose <- group_counts[1]
            pre_post_deg_list[[paste0("Program_", program)]] <- aldex_res
          }
        }
      }
    }
  }

  # 3.6 Save the results ----
  rm(program, gep_pos_spots, subset_bins, sufficient_groups, tmp_deg, group_counts)

  deg_results <- list(pre_post = pre_post_deg_list)
  saveRDS(deg_results, file = paste0(save_path, "deg_results.rds"))

  # Save the pre- and post-CAR data
  pre_post_deg_list <- Filter(function(x) length(x) > 0, pre_post_deg_list)
  pre_post_deg_list <- bind_rows(pre_post_deg_list, .id = "column_label")
  rownames(pre_post_deg_list) <- NULL
  fwrite(pre_post_deg_list, file = paste0(save_path, "pre_post_deg_list.csv"), sep = ",")

  # 3.7 Function to generate a volcano plot ----
  plot_volcano <- function(deg_df, logfc_threshold = 0.25, pval_threshold = 0.05, caps = T, top_n_labels = 10,
                           title = "Volcano Plot", subtitle = "Subtitle Here"){
    # Prepare data ----
    # These next 3 lines of code are from the Aldex volcano plot code: https://github.com/ggloor/ALDEx_bioc/blob/0cffe7151bd00752719acf64fa62ced09482072d/R/plot.aldex.r#L11
    p.add <- min(deg_df$we.eBH[deg_df$we.eBH > 0])/10
    deg_df$p_val_adj <- deg_df$we.eBH + p.add

    colnames(deg_df)[5] <- "avg_log2FC"

    deg_df$significance <- "Not Significant"
    deg_df$significance[deg_df$p_val_adj <= pval_threshold & deg_df$avg_log2FC >= logfc_threshold] <- "Upregulated"
    deg_df$significance[deg_df$p_val_adj <= pval_threshold & deg_df$avg_log2FC <= -logfc_threshold] <- "Downregulated"
    deg_df$log_p <- -log10(deg_df$p_val_adj)

    # Set max limits to visualize the data better ----
    if (caps){
      deg_df$avg_log2FC[deg_df$avg_log2FC >= 5] <- 5
      deg_df$avg_log2FC[deg_df$avg_log2FC <= -5] <- -5
      deg_df$log_p[deg_df$log_p >= 10] <- 10
    }

    # Select top genes for labelling ----
    top_genes <- deg_df %>%
      filter(significance != "Not Significant") %>%
      arrange(p_val_adj) %>%
      slice_head(n = top_n_labels)

    # Plot ----
    ggplot(deg_df, aes(x = avg_log2FC, y = log_p)) +
      geom_point(aes(color = significance), alpha = 0.75, size = 1.8) +
      geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
                 linetype = "dashed", color = "darkgray", linewidth = 0.4) +
      geom_hline(yintercept = -log10(pval_threshold),
                 linetype = "dashed", color = "darkgray", linewidth = 0.4) +
      scale_color_manual(values = c("Upregulated" = "#D73027",
                                    "Downregulated" = "#4575B4",
                                    "Not Significant" = "gray70")) +
      ggrepel::geom_text_repel(data = top_genes,
                               aes(label = gene), size = 3,
                               max.overlaps = Inf,
                               box.padding = 0.5, point.padding = 0.3,
                               segment.color = "gray60") +
      theme_minimal(base_size = 14) +
      labs(title = title, subtitle = subtitle,
           x = "log2 Fold Change", y = expression(-log[10]~adjusted~italic(p)),
           color = "Significance") +
      theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            legend.position = "top",
            panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank())
  }

  # 6.9 Find the significant genes across all program DGE ----
  # 6.9.1 Check for gene significance ----
  pre_post_subset <- pre_post_deg_list[pre_post_deg_list$column_label == "all_bins_pseudobulk",]

  # These next 3 lines of code are from the Aldex volcano plot code: https://github.com/ggloor/ALDEx_bioc/blob/0cffe7151bd00752719acf64fa62ced09482072d/R/plot.aldex.r#L11
  p.add <- min(pre_post_subset$we.eBH[pre_post_subset$we.eBH > 0])/10
  pre_post_subset$p_val_adj <- pre_post_subset$we.eBH + p.add
  colnames(pre_post_subset)[5] <- "avg_log2FC"

  logfc_threshold = 1
  pval_threshold = 0.05
  pre_post_subset$significance <- "Not Significant"
  pre_post_subset$significance[pre_post_subset$p_val_adj <= pval_threshold & pre_post_subset$avg_log2FC >= logfc_threshold] <- "Upregulated"
  pre_post_subset$significance[pre_post_subset$p_val_adj <= pval_threshold & pre_post_subset$avg_log2FC <= -logfc_threshold] <- "Downregulated"
  pre_post_subset$log_p <- -log10(pre_post_subset$p_val_adj)
  pre_post_subset <- pre_post_subset[pre_post_subset$significance != "Not Significant",]

  # 3.9.2 Inverse gene rank line graph and heatmap plots ----
  spectra_df <- load_spectra(SPECTRA)
  gep_df <- ranked_spectra(spectra_df)

  gene_of_interest <- unique(pre_post_subset$gene)[unique(pre_post_subset$gene) %in% unique(as.vector(as.matrix(gep_df)))]
  gene_position_sig <- gene_position(gep_df, gene_of_interest)

  # 3.9.3 Inverse gene heatmap plot (custom, so we don't use the classic function) -----
  gep_gene_position_df <- gene_position_sig
  gep_gene_position_wide <- gep_gene_position_df
  gep_gene_position_wide <- gep_gene_position_wide %>%
    pivot_wider(names_from = col, values_from = "row_inverse", id_cols = gene_name)
  gep_gene_position_wide <- as.data.frame(gep_gene_position_wide)
  row.names(gep_gene_position_wide) <- gep_gene_position_wide$gene_name
  gep_gene_position_wide <- gep_gene_position_wide[,-1]
  gep_gene_position_wide <- as.matrix(gep_gene_position_wide)

  scaled_mat <- t(scale(t(gep_gene_position_wide)))

  colnames(scaled_mat) <- program_annotations
  label_colors <- c(brewer.pal(name = "Set3", n = 12), brewer.pal(name = "Set2", n = 3))
  names(label_colors) <- unique(program_annotations)
  haT <- HeatmapAnnotation(Label = program_annotations,
                           col = list(Label = label_colors),
                           which = "col")

  # This takes the gene names from scaled_mat, then uses the pre_post_subset to give each significance value
  # the name of the gene (name_sig), and then orders the name_sig by the gene order in scaled_mat
  gene_order <- rownames(scaled_mat)
  name_sig <- setNames(pre_post_subset$significance, pre_post_subset$gene)
  DE_gene_ordered <- name_sig[gene_order]
  de_cols <- c(Downregulated = "#6A4C93", Upregulated = "#F2B705")
  DE_gene_split <- factor(DE_gene_ordered, levels = c("Upregulated", "Downregulated"))

  haL <- HeatmapAnnotation(DE = DE_gene_ordered,
                           col = list(DE = de_cols),
                           which = "row")

  ht_clustered <- Heatmap(scaled_mat, top_annotation = haT,
                          left_annotation = haL, row_split = DE_gene_split,
                          cluster_columns = TRUE, cluster_rows = TRUE,
                          column_names_rot = 45,
                          name = "Scaled Inverse Gene Rank",
                          row_title = "Genes", column_title = "cNMF Programs")
  ht_col_clustered <- Heatmap(scaled_mat, top_annotation = haT,
                              left_annotation = haL, row_split = DE_gene_split,
                              cluster_columns = TRUE, cluster_rows = FALSE,
                              column_names_rot = 45,
                              name = "Scaled Inverse Gene Rank",
                              row_title = "Genes", column_title = "cNMF Programs")
  ht_row_clustered <- Heatmap(scaled_mat, top_annotation = haT,
                              left_annotation = haL, row_split = DE_gene_split,
                              cluster_columns = FALSE, cluster_rows = TRUE,
                              column_names_rot = 45,
                              name = "Scaled Inverse Gene Rank",
                              row_title = "Genes", column_title = "cNMF Programs")
  ht_unclustered <- Heatmap(scaled_mat, top_annotation = haT,
                            left_annotation = haL, row_split = DE_gene_split,
                            cluster_columns = FALSE, cluster_rows = FALSE,
                            column_names_rot = 45,
                            name = "Scaled Inverse Gene Rank",
                            row_title = "Genes", column_title = "cNMF Programs")

  pdf(paste(save_path, "Sig_gene_heatmap_pseudobulk.pdf", sep = ""), height = 80, width = 25)
  draw(ht_clustered, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "Scaled Inverse Gene Rank Of DGE")
  draw(ht_col_clustered, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "Scaled Inverse Gene Rank Of DGE")
  draw(ht_row_clustered, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "Scaled Inverse Gene Rank Of DGE")
  draw(ht_unclustered, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "Scaled Inverse Gene Rank Of DGE")
  dev.off()

  # 3.9.4 Gene expression to program usage correlation ----
  gene_panel <- data.frame(category = c(rep("Upregulated", length(names(DE_gene_ordered[DE_gene_ordered == "Upregulated"]))),
                                        rep("Downregulated", length(names(DE_gene_ordered[DE_gene_ordered == "Downregulated"])))),
                           gene = c(names(DE_gene_ordered[DE_gene_ordered == "Upregulated"]),
                                    names(DE_gene_ordered[DE_gene_ordered == "Downregulated"])))
  genes_in_data <- rownames(combined_layers@assays$Spatial.024um$counts)[rownames(combined_layers@assays$Spatial.024um$counts) %in% gene_panel$gene]
  combined_count_table <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[genes_in_data ,])))

  combined_count_table_subset <- combined_count_table[, colnames(combined_count_table) %in% gene_panel$gene]
  combined_count_table_subset <- combined_count_table_subset[rownames(combined_count_table_subset) %in% rownames(usage_norm),]

  corr_matrix <- cor(combined_count_table_subset, usage_norm)
  corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset, usage_norm)))
  corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

  gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
  corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
  corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

  colnames(corr_matrix) <- program_annotations
  colnames(corr_matrix_cosine) <- program_annotations

  # Make the left annotation and the split
  DE_gene_ordered <- gene_panel$category[match(gene_order, gene_panel$gene)]
  names(DE_gene_ordered) <- gene_order
  DE_gene_split <- factor(DE_gene_ordered, levels = c("Upregulated", "Downregulated"))
  de_cols <- c(Downregulated = "#6A4C93", Upregulated = "#F2B705")
  haL <- HeatmapAnnotation(DE = DE_gene_split,
                           col = list(DE = de_cols),
                           which = "row")
  split <- DE_gene_split

  ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_split = split, row_title_rot = 0, left_annotation = haL)
  ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = TRUE,
                    row_split = split, row_title_rot = 0, left_annotation = haL)
  ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = TRUE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0, left_annotation = haL)
  ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                     cluster_columns = TRUE, cluster_rows = TRUE,
                     row_split = split, row_title_rot = 0, left_annotation = haL)
  ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       row_split = split, row_title_rot = 0, left_annotation = haL)
  ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = TRUE,
                           row_split = split, row_title_rot = 0, left_annotation = haL)
  ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = TRUE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0, left_annotation = haL)
  ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                            cluster_columns = TRUE, cluster_rows = TRUE,
                            row_split = split, row_title_rot = 0, left_annotation = haL)

  p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
  p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
  p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
  p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

  pdf(paste(save_path, "Sig_genes_gene_expression_program_usage_correlation.pdf", sep = ""), width = 40, height = 80)
  print(plot_grid(p1, p2))
  print(plot_grid(p3, p4))
  print(plot_grid(p5, p6))
  print(plot_grid(p7, p8))
  dev.off()

  # 3.10 Generate a usage-expression correlation heatmap of DE genes that are also in the top N number of genes for X programs ----
  # Take the gep_df, rename the columns, subset to the top N genes, subset to the programs of interest, and then collapse to unique genes
  for (program in program_annotations){
    gep_df_subset <- gep_df
    colnames(gep_df_subset) <- program_annotations
    gep_df_subset <- gep_df_subset[c(1:250), program]

    top_genes <- unique(unlist(gep_df_subset))
    top_genes <- top_genes[top_genes %in% pre_post_subset$gene]

    if(!is.null(top_genes)){
      gene_panel <- data.frame(category = rep("Top 250 Genes of program also DE in Pre-Post Pseudobulk", length(top_genes)),
                               gene = top_genes)
      genes_in_data <- rownames(combined_layers@assays$Spatial.024um$counts)[rownames(combined_layers@assays$Spatial.024um$counts) %in% gene_panel$gene]
      combined_count_table <- as.data.frame(t(as.matrix(combined_layers@assays$Spatial.024um$counts[genes_in_data ,])))

      combined_count_table_subset <- combined_count_table[, colnames(combined_count_table) %in% gene_panel$gene]
      combined_count_table_subset <- combined_count_table_subset[rownames(combined_count_table_subset) %in% rownames(usage_norm),]

      corr_matrix <- cor(combined_count_table_subset, usage_norm)
      corr_matrix_cosine <- cosine(as.matrix(cbind(combined_count_table_subset, usage_norm)))
      corr_matrix_cosine <- corr_matrix_cosine[c(1:(nrow(corr_matrix))), c((nrow(corr_matrix)+1):ncol(corr_matrix_cosine))]

      gene_order <- gene_panel$gene[gene_panel$gene %in% rownames(corr_matrix)]
      corr_matrix <- corr_matrix[match(gene_order, rownames(corr_matrix)),]
      corr_matrix_cosine <- corr_matrix_cosine[match(gene_order, rownames(corr_matrix_cosine)),]

      colnames(corr_matrix) <- program_annotations
      colnames(corr_matrix_cosine) <- program_annotations

      # Make the left annotation and the split
      DE_gene_ordered <- setNames(pre_post_subset$significance, pre_post_subset$gene)
      DE_gene_ordered <- DE_gene_ordered[gene_order]
      DE_gene_split <- factor(DE_gene_ordered, levels = c("Upregulated", "Downregulated"))
      de_cols <- c(Downregulated = "#6A4C93", Upregulated = "#F2B705")
      haL <- HeatmapAnnotation(DE = DE_gene_split,
                               col = list(DE = de_cols),
                               which = "row")
      split <- DE_gene_split

      ht <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                    cluster_columns = FALSE, cluster_rows = FALSE,
                    row_split = split, row_title_rot = 0, left_annotation = haL)
      ht_row <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                        cluster_columns = FALSE, cluster_rows = TRUE,
                        row_split = split, row_title_rot = 0, left_annotation = haL)
      ht_col <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                        cluster_columns = TRUE, cluster_rows = FALSE,
                        row_split = split, row_title_rot = 0, left_annotation = haL)
      ht_both <- Heatmap(corr_matrix, name = "Gene Expression and Program Usage Pearson Correlation",
                         cluster_columns = TRUE, cluster_rows = TRUE,
                         row_split = split, row_title_rot = 0, left_annotation = haL)
      ht_cosine <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                           cluster_columns = FALSE, cluster_rows = FALSE,
                           row_split = split, row_title_rot = 0, left_annotation = haL)
      ht_cosine_row <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                               cluster_columns = FALSE, cluster_rows = TRUE,
                               row_split = split, row_title_rot = 0, left_annotation = haL)
      ht_cosine_col <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                               cluster_columns = TRUE, cluster_rows = FALSE,
                               row_split = split, row_title_rot = 0, left_annotation = haL)
      ht_cosine_both <- Heatmap(corr_matrix_cosine, name = "Gene Expression and Program Usage Cosine Correlation",
                                cluster_columns = TRUE, cluster_rows = TRUE,
                                row_split = split, row_title_rot = 0, left_annotation = haL)

      p1 <- grid.grabExpr(draw(ht, padding = unit(c(3, 3, 3, 3), "cm"))); p2 <- grid.grabExpr(draw(ht_cosine, padding = unit(c(3, 3, 3, 3), "cm")))
      p3 <- grid.grabExpr(draw(ht_row, padding = unit(c(3, 3, 3, 3), "cm"))); p4 <- grid.grabExpr(draw(ht_cosine_row, padding = unit(c(3, 3, 3, 3), "cm")))
      p5 <- grid.grabExpr(draw(ht_col, padding = unit(c(3, 3, 3, 3), "cm"))); p6 <- grid.grabExpr(draw(ht_cosine_col, padding = unit(c(3, 3, 3, 3), "cm")))
      p7 <- grid.grabExpr(draw(ht_both, padding = unit(c(3, 3, 3, 3), "cm"))); p8 <- grid.grabExpr(draw(ht_cosine_both, padding = unit(c(3, 3, 3, 3), "cm")))

      pdf(paste0(save_path, "T250_DE_", program, "_Intersecting_gene_expression_program_usage_correlation.pdf"), width = 40, height = 20)
      print(plot_grid(p1, p2))
      print(plot_grid(p3, p4))
      print(plot_grid(p5, p6))
      print(plot_grid(p7, p8))
      dev.off()
    }
  }
}

# 4.4 Program to Expression Associations ----
if (PROGRAM_EXPRESSION_ASSOCIATIONS){
  # 4.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "K", RANK, "/program_expression_association/", sep = "")
  create_folder(save_path)
  create_folder(paste0(save_path, "GCAR/"))
  create_folder(paste0(save_path, "Sample_Split/"))

  # 4.1 Modify the usage_norm df and subset the counts to genes of interest ----
  cutoff <- as.data.frame(usage_norm)
  cutoff$sample <- gsub(paste0("(", pattern_regex, ").*"), "\\1", rownames(cutoff))
  bins_above_cutoff <- lapply(seq_len(ncol(cutoff) - 1), function(col_idx) {
    rownames(cutoff)[cutoff[, col_idx] > GLOBAL_THRESHOLD]
  })

  genes_to_subset <- unique(c(GENE_EXPRESSION_LIST, as.vector(unlist(GENE_EXPRESSION_NORMALIZE))))
  combined_count_table_subset <- combined_layers@assays$Spatial.024um$counts[genes_to_subset, unique(unlist(bins_above_cutoff))]
  combined_count_table_subset <- as.data.frame(t(as.matrix(combined_count_table_subset)))
  combined_count_table_subset$sample <- gsub(paste0("(", pattern_regex, ").*"), "\\1", rownames(combined_count_table_subset))

  # 4.2 Check for an association between GCAR+ bins and some programs among all program bins ----
  program_2_list <- c(setNames(as.list(seq_along(colnames(usage_norm))),
                               as.character(seq_along(colnames(usage_norm)))),
                      list("9_11_12" = c(9, 11, 12), "11_12"   = c(11, 12)))

  gene <- "GCAR"
  ft_list <- vector(mode = "list", length = (ncol(usage_norm) * length(program_2_list)))
  ft_counter <- 1

  for (program_1 in seq_along(colnames(usage_norm))){
    create_folder(paste0(save_path, "GCAR/", program_1, "/"))

    for (program_2 in program_2_list){
      # 2.1 Subset to just 1) P05 (T_cell program) positive bins, 2) GCAR+ bins, and 3) All other programs ----
      program_2_collapse <- paste(program_2, collapse = "_")

      p1_pos_bins <- unlist(bins_above_cutoff[program_1])
      gene_pos_bins <- rownames(combined_count_table_subset[combined_count_table_subset[, gene] > 0,])
      p2_pos_bins <- unique(unlist(bins_above_cutoff[program_2]))

      # Need to subset all
      SPSC_bS2_bins <- rownames(combined_layers@meta.data[combined_layers@meta.data$sample == "SPSC_biopsy2",])
      p1_pos_bins <- intersect(p1_pos_bins, SPSC_bS2_bins)
      gene_pos_bins <- intersect(gene_pos_bins, SPSC_bS2_bins)
      p2_pos_bins <- intersect(p2_pos_bins, SPSC_bS2_bins)

      # 2.2 Calculate 1) Gene+ P1+ P2+, 2) Gene- P1+ P2+, 3) Gene+ P1+ P2-, and 4) Gene- P1+ P2- ----
      a <- length(intersect(p1_pos_bins, intersect(gene_pos_bins, p2_pos_bins))) # Gene+ P1+ P2+
      b <- length(setdiff(intersect(p1_pos_bins, p2_pos_bins), gene_pos_bins))  # Gene- P1+ P2+
      c <- length(setdiff(intersect(p1_pos_bins, gene_pos_bins), p2_pos_bins)) # Gene+ P1+ P2-
      d <- length(setdiff(p1_pos_bins, union(gene_pos_bins, p2_pos_bins))) # Gene- P1+ P2-

      if (sum(a, b, c, d) != 0){
        # 2.3 Make a contingency table and run stats ----
        contingency_table <- matrix(c(a, b, c, d),
                                    nrow = 2,
                                    byrow = T,
                                    dimnames = list(c("A+_B+", "A+_B-"),
                                                    c("A+_GCAR+", "A+_GCAR-")))

        ft_df <- pairwise_fisher_test(contingency_table, p.adjust.method = "fdr")

        ft_df$program_one <- program_1
        ft_df$program_two <- paste(program_2, collapse = "_")
        ft_df$GCAR_pos_A_pos_B_pos <- contingency_table["A+_B+", "A+_GCAR+"]
        ft_df$GCAR_pos_A_pos_B_neg <- contingency_table["A+_B-", "A+_GCAR+"]
        ft_df$A_pos_B_pos <- contingency_table["A+_B+", "A+_GCAR-"]
        ft_df$A_pos_B_neg <- contingency_table["A+_B-", "A+_GCAR-"]

        # 2.4 Plot the contingency table ----
        df_long <- contingency_table %>%
          as.data.frame() %>%
          rownames_to_column(var = "B_status") %>%
          pivot_longer(cols = -B_status, names_to = "GCAR_status", values_to = "Count") %>%
          group_by(GCAR_status) %>%
          mutate(GCAR_sum = sum(Count),
                 Proportion = Count / GCAR_sum) %>%
          ungroup()

        plot_title <- paste0(gene ,"+ vs P", program_2_collapse, " positivity among P", program_1 ,"+ bins")

        # Turn to proportions
        p1 <- ggplot(df_long, aes(x = GCAR_status, y = Count, fill = B_status)) +
          geom_bar(stat = "identity", position = "dodge") + theme_classic2() +
          labs(title = plot_title,
               subtitle = paste0("A = ", program_1, "; B = ", program_2_collapse, "; P-Value: ", ft_df$p),
               x = "GCAR status (A+ bins)", y = "Number of Bins of GCAR Status", fill = "B status") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
                axis.text.y = element_text(colour="black"),
                legend.position = "top", text = element_text(size = 12)) +
          theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) +
          scale_fill_manual(values = c("steelblue", "tomato"))
        p2 <- ggplot(df_long, aes(x = GCAR_status, y = Proportion, fill = B_status)) +
          geom_bar(stat = "identity", position = "dodge") + theme_classic2() +
          labs(title = plot_title,
               subtitle = paste0("A = ", program_1, "; B = ", program_2_collapse, "; P-Value: ", ft_df$p),
               x = "GCAR status (A+ bins)", y = "Proportion of Bins of GCAR Status", fill = "B status") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
                axis.text.y = element_text(colour="black"),
                legend.position = "top", text = element_text(size = 12)) +
          theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) +
          scale_fill_manual(values = c("steelblue", "tomato"))

        pdf(paste0(save_path, gene, "/", program_1, "/", gene, "_P", program_1, "_vs_P", program_2_collapse, ".pdf"))
        mosaicplot(contingency_table,
                   main = paste0("Mosaic Plot of: ", plot_title),
                   sub = paste0("A = ", program_1, "; B = ", program_2_collapse),
                   color = TRUE)
        print(p1)
        print(p2)
        dev.off()

        ft_df$empirical_pval_100000 <- NA
        if (sum(contingency_table[, "A+_GCAR+"]) != 0){
          # 2.5 Empirical P-value ----
          # Make a vector containing GCAR1- and P05+ bins
          gcar_neg_p1_pos <- setdiff(p1_pos_bins, gene_pos_bins)
          
          n_perm <- 100000L
          N <- length(gcar_neg_p1_pos)
          k <- sum(contingency_table[, "A+_GCAR+"])
          x_obs <- contingency_table["A+_B+", "A+_GCAR+"]
          is_p2 <- gcar_neg_p1_pos %in% p2_pos_bins
          idx <- replicate(n_perm, sample.int(N, k))
          p2_pos_rate <- colSums(matrix(is_p2[idx], nrow = k))
          
          ft_df$empirical_pval_100000 <- mean(p2_pos_rate == x_obs)
        }
        ft_list[[ft_counter]] <- as.data.frame(ft_df)

        ft_counter <-  ft_counter + 1
      }
    }
  }

  ft_df <- bind_rows(ft_list)

  fwrite(ft_df, file = paste0(save_path, gene, "/", gene, "_fisher_exact_test_results.csv"), sep = ",")

  # 4.3 Check for an association between P05+ bins and some programs between samples ----
  program_2_list <- c(setNames(as.list(seq_along(colnames(usage_norm))),
                               as.character(seq_along(colnames(usage_norm)))),
                      list("9_11_12" = c(9, 11, 12), "11_12"   = c(11, 12)))
  ft_list <- vector(mode = "list", length = (ncol(usage_norm) * length(program_2_list)))
  ft_counter <- 1

  for (program_1 in seq_along(colnames(usage_norm))){
    create_folder(paste0(save_path, "Sample_Split/", program_1, "/"))

    for (program_2 in program_2_list){
      # 2.1 Subset to just 1) P05 (T_cell program) positive bins, 2) GCAR+ bins, and 3) Vasculature unit positive bins (9, 11, 12) ----
      program_2_collapse <- paste(program_2, collapse = "_")

      p1_pos_bins <- unlist(bins_above_cutoff[program_1])
      SPSC_bS1_bins <- p1_pos_bins[grep("SPSC_biopsy1", p1_pos_bins)]
      SPSC_bS2_bins <- p1_pos_bins[grep("SPSC_biopsy2", p1_pos_bins)]
      p2_pos_bins <- unique(unlist(bins_above_cutoff[program_2]))

      # 2.2 Calculate 1) P1+, bS1+, P2+, 2) P1+, b2+, P2+, 3) P1+, b1+, P2-+, and 4) P1+, b2+, P2-+ ----
      a <- length(intersect(p1_pos_bins, intersect(SPSC_bS1_bins, p2_pos_bins))) # P1+, bS1+, P2+
      b <- length(intersect(p1_pos_bins, intersect(SPSC_bS2_bins, p2_pos_bins)))  # P1+, b2+, P2+
      c <- length(setdiff(intersect(p1_pos_bins, SPSC_bS1_bins), p2_pos_bins)) # P1+, b1+, P2-+
      d <- length(setdiff(intersect(p1_pos_bins, SPSC_bS2_bins), p2_pos_bins)) # P1+, b2+, P2-+

      # 2.3 Make a contingency table and run stats ----
      contingency_table <- matrix(c(a, b, c, d),
                                  nrow = 2,
                                  byrow = T,
                                  dimnames = list(c("A+_B+", "A+_B-"),
                                                  c("A+_bS1+", "A+_bS2+")))

      ft_df <- pairwise_fisher_test(contingency_table, p.adjust.method = "fdr")

      ft_df$program_one <- program_1
      ft_df$program_two <- program_2_collapse
      ft_df$bS1_A_pos_B_pos <- contingency_table["A+_B+", "A+_bS1+"]
      ft_df$bS1_A_pos_B_neg <- contingency_table["A+_B-", "A+_bS1+"]
      ft_df$bS2_A_pos_B_pos <- contingency_table["A+_B+", "A+_bS2+"]
      ft_df$bS2_A_pos_B_neg <- contingency_table["A+_B-", "A+_bS2+"]

      # 2.4 Plot the contingency table ----
      df_long <- contingency_table %>%
        as.data.frame() %>%
        rownames_to_column(var = "B_status") %>%
        pivot_longer(cols = -B_status, names_to = "Sample", values_to = "Count") %>%
        group_by(Sample) %>%
        mutate(Sample_Sum = sum(Count),
               Proportion = Count / Sample_Sum) %>%
        ungroup()

      plot_title <- paste0("Sample positivity vs P", program_2_collapse, " positivity among P", program_1 ,"+ bins")
      
      p1 <- ggplot(df_long, aes(x = Sample, y = Count, fill = B_status)) +
        geom_bar(stat = "identity", position = "dodge") + theme_classic2() +
        labs(title = plot_title,
             subtitle = paste0("A = ", program_1, "; B = ", program_2_collapse, "; P-Value: ", ft_df$p),
             x = "Sample (A+ bins)", y = "Number of Bins of Sample", fill = "P2 status") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
              axis.text.y = element_text(colour="black"),
              legend.position = "top", text = element_text(size = 12)) +
        theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) +
        scale_fill_manual(values = c("steelblue", "tomato"))
      p2 <- ggplot(df_long, aes(x = Sample, y = Proportion, fill = B_status)) +
        geom_bar(stat = "identity", position = "dodge") + theme_classic2() +
        labs(title = plot_title,
             subtitle = paste0("A = ", program_1, "; B = ", program_2_collapse, "; P-Value: ", ft_df$p),
             x = "Sample (A+ bins)", y = "Proportion of Bins of Sample", fill = "P2 status") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
              axis.text.y = element_text(colour="black"),
              legend.position = "top", text = element_text(size = 12)) +
        theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0)) +
        scale_fill_manual(values = c("steelblue", "tomato"))

      pdf(paste0(save_path, "Sample_Split/", program_1, "/Sample_Split_P", program_1, "_vs_P", paste(program_2, collapse = "_"), ".pdf"))
      mosaicplot(contingency_table,
                 main = paste0("Mosaic Plot of: ", plot_title),
                 sub = paste0("A = ", program_1, "; B = ", program_2_collapse),
                 color = TRUE)
      print(p1)
      print(p2)
      dev.off()
      
      ft_list[[ft_counter]] <- as.data.frame(ft_df)
    }
  }

  ft_df <- bind_rows(ft_list)

  fwrite(ft_df, file = paste0(save_path, "Sample_Split/", "Sample_Split_fisher_exact_test_results.csv"), sep = ",")
}

# 4.5 Proportion of GCAR+ bins among T-cell program+ bins ----
if (BIN_FRACTION){
  # 5.0 Clear existing non-function variables, set the save folders (create them if they doesn't exist) ----
  rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
  save_path <- paste(SAVE_FOLDER, "/K", RANK, "/bin_fraction/", sep = "")
  create_folder(save_path)
  
  # 5.1 Modify the usage_norm df and subset the counts to genes of interest ----
  cutoff <- as.data.frame(usage_norm)
  cutoff$sample <- gsub(paste0("(", pattern_regex, ").*"), "\\1", rownames(cutoff))
  bins_above_cutoff <- lapply(seq_len(ncol(cutoff) - 1), function(col_idx) {
    rownames(cutoff)[cutoff[, col_idx] > GLOBAL_THRESHOLD]
  })
  
  genes_to_subset <- unique(c(GENE_EXPRESSION_LIST, as.vector(unlist(GENE_EXPRESSION_NORMALIZE))))
  combined_count_table_subset <- combined_layers@assays$Spatial.024um$counts[genes_to_subset, unique(unlist(bins_above_cutoff))]
  combined_count_table_subset <- as.data.frame(t(as.matrix(combined_count_table_subset)))
  combined_count_table_subset$sample <- gsub(paste0("(", pattern_regex, ").*"), "\\1", rownames(combined_count_table_subset))
  
  # 5.2 Calculate GCAR+ proportion among program_1 positive bins ----
  gene <- "GCAR"
  program_1 <- 10
  p1_pos_bins <- unlist(bins_above_cutoff[[program_1]])
  
  # Subset only P1+ bins
  p1_bin_table <- combined_count_table_subset[p1_pos_bins, c(gene, "sample"), drop = FALSE]
  p1_bin_table$GCAR_pos <- p1_bin_table[[gene]] > 0
  
  # Calculate proportions per sample
  prop_tbl <- p1_bin_table %>%
    group_by(sample) %>%
    summarise(total_p1_bins = n(),
              gcar_pos_bins = sum(GCAR_pos),
              prop_gcar = gcar_pos_bins / total_p1_bins, 
              .groups = "drop") %>%
    filter(total_p1_bins > 0)
  
  # 5.3 Plot log-scaled bar plot ----
  p1 <- ggplot(prop_tbl, aes(x = sample, y = prop_gcar)) +
    geom_col(fill = "#2C7BB6") + theme_classic2() +
    labs(x = "Sample", y = "Proportion of GCAR+ bins",
         title = "Proportion of GCAR+ Bins of T-cell Program+ Bins",
         subtitle = paste0("Program = ", program_1, "; Usage Cutoff = ", GLOBAL_THRESHOLD)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour="black"),
          legend.position = "top", text = element_text(size = 12)) +
    theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_y_continuous(expand = c(0,0))
  
  fwrite(as.data.frame(prop_tbl), paste0(save_path, "GCAR_Bin_Fraction_Table.csv"), sep = ",")
  pdf(paste0(save_path, "GCAR_Bin_Fraction_Stacked.pdf"))
  print(p1)
  dev.off()
}
