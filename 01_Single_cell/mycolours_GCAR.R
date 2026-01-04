##############################
#### my colours           ####
#### - project: [ GCAR1 ] ####
##############################
# version: 1.14.1 (2026-01-04)
# by hyojinsong
today <- Sys.Date()
# ---

#### package ####
mypackages <- c(
  "RColorBrewer", "pals" #for colour palette
  ,"stringr" #for [str_pad] function
)

## [RColorBrewer] ##
# install.packages("RColorBrewer")

## [pals] ##
# install.packages("pals")

lapply(mypackages, library, character.only = TRUE)
# ---


#### GCAR ####
## #1 | longitudinal sample ##
# ---
## [Patient-01_Dose-01] - aka [SPS-Q] ##
## - note: created based on [Dark2] (n=7)
col_gcar_e <- "#666666" #dark_grey
col_gcar_h <- "#bf5b17" #brown
col_gcar_d08 <- "#a6d854" #light_green
col_gcar_d15 <- "#e7298a" #pink
col_gcar_d22 <- "#e6ab02" #yellow
col_gcar_d28 <- "#1b9e77" #dark_green
col_gcar_d42 <- "#386cb0" #blue
smpColor <- setNames(c(col_gcar_p,col_gcar_t,col_gcar_d08,col_gcar_d15,col_gcar_d22,col_gcar_d28,col_gcar_d42), c("Enriched_Apheresis","GCAR1_Product","D08","D15","D22","D28","D42")) #sample_new
# ---
## [Patient-01_Dose-02] - aka [SPS-Q2] ##
## - note: created based on [Dark2] + similar colours to differentiate between new samples (n=10*)
##         --> *[10 samples] = { 7 new samples from the same patient before/after 2nd infusion } + { 3 RERUN samples from the same patient before/after 1st infusion }
col_gcar_bx_pre <- "#8c6bb1" #light_purple
col_gcar_bx_post <- "#810f7c" #dark_purple

col_gcar_md7 <- "#666666" #dark_grey #minus_day7

col_gcar_d04 <- "#a6761d" #earthy_brown
col_gcar_d08 <- "#a6d854" #light_green
col_gcar_d11 <- "#1b9e77" #green
col_gcar_d14 <- "#80b1d3" #light_blue

col_gcar_p_re <- "#bf5b17" #brown #harvest
col_gcar_d15_re <- "#e7298a" #pink
col_gcar_d22_re <- "#e6ab02" #yellow

smpColor <- setNames(c(col_gcar_p_re,col_gcar_d15_re,col_gcar_d22_re, col_gcar_bx_pre,col_gcar_md7, col_gcar_d04,col_gcar_d08,col_gcar_d11,col_gcar_d14,col_gcar_bx_post), 
                     c("GCAR1product-RERUN","D15-RERUN","D22-RERUN", "Biopsy-PRE","DAY-7-20250210", "DAY04-20250221","DAY08-20250225","DAY11-20250228","DAY14-20250303","Biopsy-POST")) #nSmp=10
# ---


## #2 | expression ##
## - source: [https://xkcd.com/color/rgb/]
# pal_expr=sns.blend_palette(["lightgrey", sns.xkcd_rgb["scarlet"]], as_cmap=True)
col_pos <- "#be0119" #scarlet
col_neg <- "#d8dcd6" #lightgrey
exprColor <- setNames(c(col_pos,col_neg), c("pos","neg"))
# ---


## #3 | Seurat cell cluster ##
# nClst <- length(unique(sobj_flt_umap$seurat_clusters))
# if (nClst<=25) {
#     clstColor <- setNames(cols25(n = 25)[1:nClst], sort(unique(as.factor(sobj_flt_umap$seurat_clusters))))
# } else if (nClst>25) {
#     clstColor <- setNames(c(cols25(n = 25),kelly(n=22)[3:22])[1:nClst], sort(unique(as.factor(sobj_flt_umap$seurat_clusters))))
# } else if (nClst>45) {
#     clstColor <- setNames(c(cols25(n = 25),kelly(n=22))[1:nClst], sort(unique(as.factor(sobj_flt_umap$seurat_clusters))))
# }

# col_C_00   <- "#1F78C8"
# col_C_01   <- "#ff0000"
# col_C_02   <- "#33a02c"
# col_C_03   <- "#6A33C2"
# col_C_04   <- "#ff7f00"
# col_C_05   <- "#565656"
# col_C_06   <- "#FFD700"
# col_C_07   <- "#a6cee3"
# col_C_08   <- "#FB6496"
# col_C_09   <- "#b2df8a"
# col_C_10   <- "#CAB2D6"
# col_C_11   <- "#FDBF6F"
# col_C_12   <- "#999999"
# col_C_13   <- "#EEE685"
# col_C_14   <- "#C8308C"
# col_C_15   <- "#FF83FA"
# col_C_16   <- "#C814FA"
# col_C_17   <- "#0000FF"
# col_C_18   <- "#36648B"
# col_C_19   <- "#00E2E5"
# col_C_20   <- "#00FF00"
# col_C_21   <- "#778B00"
# col_C_22   <- "#BEBE00"
# col_C_23   <- "#8B3B00"
# col_C_24   <- "#A52A3C"
# col_C_25   <- "#F2F3F4"
# col_C_26   <- "#222222"
# col_C_27   <- "#F3C300"
# col_C_28   <- "#875692"
# col_C_29   <- "#F38400" #excluded from the analyses
# col_C_30   <- "#A1CAF1"
# col_C_31   <- "#BE0032"
# col_C_32   <- "#C2B280"
# col_C_33   <- "#848482"
# col_C_34   <- "#008856"
# col_C_35   <- "#E68FAC"
# col_C_36   <- "#0067A5"
# col_C_37   <- "#F99379"
# col_C_38   <- "#604E97"
# col_C_39   <- "#F6A600"
# col_C_40   <- "#B3446C"
# col_C_41   <- "#DCD300"
# col_C_42   <- "#882D17"
# col_C_43   <- "#8DB600" #excluded from the analyses
# col_C_44   <- "#654522"
# col_C_45   <- "#E25822"
# col_C_46   <- "#2B3D26"

nClst <- 47 #res=2.4
list_clst <- paste0("C_",str_pad(seq(0,(nClst-1),by=1),2,pad="0"))
clstColor <- setNames(c(cols25(n = 25),kelly(n=22))[1:nClst], list_clst) #up to 47 #from [C_00] to [C_46]
# ---


## #4 | cell class ##
list_CellClass_specific <- sort(unique(obs_c__exCC43_CC29$Specific_cell_class))
nCellClass_specific <- length(list_CellClass_specific) #23
ccColor21 <- setNames( cols25(n = 25)[1:nCellClass_specific], list_CellClass_specific ) #21
#                 ccColor21
# B-cells           #1F78C8
# CD4               #ff0000
# CD4_Treg          #33a02c
# CD4em             #6A33C2
# CD4naive          #ff7f00
# CD4naive_helper   #565656
# CD8e              #FFD700
# CD8e_td           #a6cee3
# CD8em             #FB6496
# CD8em_MAIT        #b2df8a
# CD8em_NK          #CAB2D6
# CD8em_td          #FDBF6F
# CD8em_td_MAIT     #999999
# DCs               #EEE685
# HSC               #C8308C
# Lymphocyte        #FF83FA
# MAIT              #C814FA
# Monocyte          #0000FF
# ncMonocyte        #36648B
# NK                #00E2E5
# NKT               #00FF00
# ---

## #5 | cell type ##
## - note: built based on this colour palette above: [ccColor21] -- union
fn_csv_ctColor49 <- "./list_ctColor49.csv"
csv_ctColor49 <- read.csv(file = fn_csv_ctColor49, header = T)
head(csv_ctColor49,3)
#           list_label_cellType ctColor49
# 1                      B_cell   #1F78C8
# 2                  Blood_cell   #3B00FB
# 3                         CD4   #ff0000
# 4                    CD4_Treg   #33a02c
# 5                       CD4cm   #66B0FF
# 6                CD4cytotoxic   #FA0087
# 7                       CD4em   #6A33C2
# 8                   CD4helper   #D85FF7
# 9                        CD4m   #1C7F93
# 10                   CD4naive   #ff7f00
# 11            CD4naive_helper   #565656
# 12                        CD8   #7ED7D1
# 13                      CD8cm   #B5EFB5
# 14               CD8cytotoxic   #822E1C
# 15                       CD8e   #FFD700
# 16                    CD8e_td   #a6cee3
# 17                      CD8em   #FB6496
# 18                 CD8em_MAIT   #b2df8a
# 19                   CD8em_NK   #CAB2D6
# 20                   CD8em_td   #FDBF6F
# 21              CD8em_td_MAIT   #999999
# 22                       CD8m   #BDCDFF
# 23                   CD8naive   #AAF400
# 24                        cDC   #782AB6
# 25                         DC   #EEE685
# 26           Endothelial_cell   #C075A6
# 27 Endothelial_cell_capillary   #F7E1A0
# 28            Epithelial_cell   #FC1CBF
# 29                 Fibroblast   #683B79
# 30                        HSC   #C8308C
# 31                 Lymphocyte   #FF83FA
# 32                 Macrophage   #1CBE4F
# 33                       MAIT   #C814FA
# 34                        mDC   #FBE426
# 35                   Monocyte   #0000FF
# 36                 ncMonocyte   #36648B
# 37                 Neutrophil   #B10DA1
# 38                    NK_cell   #00E2E5
# 39                        pDC   #85660D
# 40                   Pericyte   #1C8356
# 41                        SMC   #C4451C
# 42                     T_cell   #325A9B
# 43                  T_cell_ab   #F8A19F
# 44                  T_cell_gd   #AA0DFE
# 45                   T_cell_h   #DEA0FD
# 46                T_cell_MAIT   #2ED9FF
# 47                  T_cell_NK   #00FF00
# 48                       Treg   #90AD1C
# 49                Tumor_cells   #5A5156
ctColor49 <- setNames( csv_ctColor49$ctColor49, csv_ctColor49$list_label_cellType ) #49
# ---
