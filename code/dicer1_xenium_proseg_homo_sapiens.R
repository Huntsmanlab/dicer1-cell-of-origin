---
title: "DICER1_Xenium_proseg_homo_sapiens"
author: "Felix Kommoss"
date: "2024-09-16"
---

#################################################################################
## load libraries
library(Seurat)
library(tidyverse)
library(BPCells)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(SeuratWrappers)
library(monocle3)
library(reshape2)
library(MetBrewer)
library(RColorBrewer)
library(ggplot2)
library(viridis)
library(glmGamPoi)
library(future)
library(reshape2)
library(patchwork)
library(cowplot)
library(ggpubr)
library(yaml)
library(arrow)
library(ape)

#################################################################################
## set seed
set.seed(404)

#################################################################################
# output directory
output_dir <- "../output/plots/"

#################################################################################
## set themes
no_legend_theme <- theme(legend.position = "none")

colors <- c( "#C77CFF", "#00BB4E", "#00B0F6", "#00C1A3", "#FF61CC", 
             "#A3A500", "#D89000", "#FD7083", "#CE9500", "#39B606", 
             "#00B0E1", "#E28A00", "#D973FC", "#B2A000", "#00C091", 
             "#F27D53", "#FF6A98", "#B186FF", "#00BAE0", "#F265E8", 
             "#00B5ED", "#FF62BC", "#CD9600", "#00B92A", "#00C0B5", 
             "#EA8331", "#00A9FF", "#7099FF", "#00BE67", "#35A2FF", 
             "#E76BF3", "#91AA00", "#00BDD4", "#00C0B6", "#FF65AA", 
             "#FD6F88", "#61B200", "#C09B06", "#39B607", "#F8766D")
colors_secondary <- c( "#FC8D62", "#66C2A5", "#8DA0CB", "#E78AC3",
                       "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
                       "#fdb462", "#bc80bd")

cluster_colors <- list(
  "Endo" = "#F8766D",
  "Immune" = "#B186FF",
  "Fibro" = "#CE9500",
  "Epi" = "#FF6A98",
  "Uni Fibro" = "#00BAE0",
  "Ground" = "#00BB4E",
  "Ground 1" = "#00BB4E",
  "Ground 2" = "#00C0B6",
  "Progen" = "#C77CFF",
  "Diff Myo" = "#00C1A3",
  "Prolif" = "#FF61CC",
  "Prolif 1" = "#FF61CC",
  "Prolif 2" = "#F8766D",
  "TR Diff Myo" = "#7099FF",
  "Mural" = "#A3A500",
  "Unk" = "grey"
)

#################################################################################
## load data
# load Xenium object for human data
xenium.obj <- readRDS("../data/dicer1_xenium_homo_sapiens.rds")

#################################################################################
## QC
# remove cells with less than 3 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 3)
## remove cells with less than 3 genes
xenium.obj <- subset(xenium.obj, subset = nFeature_Xenium > 3)

# cell counts for xenium.obj
cell_count_tumour <- table(xenium.obj@meta.data$tumour)
print(cell_count_tumour)
cell_count_site <- table(xenium.obj@meta.data$site)
print(cell_count_site)
cell_count_class <- table(xenium.obj@meta.data$class)
print(cell_count_class)
cell_count_gender <- table(xenium.obj@meta.data$gender)
print(cell_count_gender)

#################################################################################
# SCTransform
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium") 

# RunPCA
xenium.obj <- RunPCA(xenium.obj, verbose = FALSE)
DimHeatmap(xenium.obj, dims = 1:20, cells = 500, balanced = TRUE)
# determine dims with ElbowPlot
ElbowPlot(xenium.obj, ndims = 50, reduction = "pca")
# RunUMAP
xenium.obj <- RunUMAP(xenium.obj, reduction = "pca", dims = 1:12)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:12)
# FindClusters using louvain algorithm
xenium.obj <- FindClusters(xenium.obj, resolution = 0.3)

#################################################################################
## plot cluster and tumour (ID)
# DimPlot cluster
(xenium.obj_umap_clusters = DimPlot(
  xenium.obj, label = T, repel = TRUE,reduction = "umap",
  shuffle = T, cols = colors, group.by = "seurat_clusters") +
  theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_hs_umap_clusters.svg", 
       plot = xenium.obj_umap_clusters, width = 5, height = 4, 
       device = "svg", path = output_dir)
# DimPlot tumour (ID)
(xenium.obj_umap_tumours = DimPlot(
  xenium.obj, label = T, repel = TRUE,reduction = "umap",
  shuffle = T, cols = colors, group.by = "tumour")+
    theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_hs_umap_tumours.svg", 
       plot = xenium.obj_umap_tumours, width = 5, height = 4, 
       device = "svg", path = output_dir)

## count cells in UMAP # 75,938
total_cells_all <- ncol(xenium.obj)
print(total_cells_all)

#################################################################################
## extract annotation for cells visualization in Xenium explorer
clusters_annotation = xenium.obj$SCT_snn_res.0.3
df_clusters = 
  data.frame(cell_id = names(clusters_annotation), group = clusters_annotation)
write.csv(x = df_clusters, file = "../output/tables/dicer1_xenium_hs_annotation_cluster.csv", 
          row.names = FALSE)

#################################################################################
## find markers for all 13 cluster
# FindAllMarkers
Idents(xenium.obj) <- xenium.obj@meta.data$seurat_clusters
markers_clusters <- FindAllMarkers(xenium.obj, only.pos = TRUE)
markers_clusters %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

## print markers for all 13 clusters
# write CSV
write_csv(x = markers_clusters, file = "dicer1_xenium_hs_markers_clusters.csv")

#################################################################################
## annotate cluster using broad cell type markers
# Path to YAML file for cell type markers
cell_type_markers <- read_yaml("../reference/cell_type_markers_homo_sapiens.yaml")

## plot cell type markers
# DotPlot
(xenium.obj_dotplot_clusters = DotPlot(xenium.obj, features = cell_type_markers) 
  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  + theme(strip.text.x = element_text(angle = 90, hjust = 0)))
ggsave(filename = "dicer1_xenium_hs_dotplot_clusters.svg", 
       plot = xenium.obj_dotplot_clusters, width = 10, height = 6, 
       device = "svg", path = output_dir)

#################################################################################
## batch and rename clusters based on broad cell type markers
# cluster names
xenium.obj_cell_type <- list(
  `0` = "Ground",    # Fibroblasts
  `1` = "Progen",    # Fibroblasts
  `2` = "Prolif",    # Fibroblasts
  `3` = "Prolif",    # Fibroblasts
  `4` = "Endo",    # Endothelium
  `5` = "TR Diff Myo",    # Fibroblasts
  `6` = "Immune",    # Macrophages
  `7` = "Epi",    # Epithelium
  `8` = "Mural", # Pericytes/Smooth msucle cells (Mural cells)
  `9` = "Immune",   # B lymphocytes
  `10` = "Fibro",   # Fibroblasts
  `11` = "Unk",   # Fibroblasts
  `12` = "Immune", # T lymphocytes/NK
  `13` = "Diff Myo"   # Skeletal muscle (Myocytes)
)

xenium.obj@meta.data$cell_type <- 
  factor(x = 
           xenium.obj_cell_type[as.character(xenium.obj$seurat_clusters)],
         levels = c("Progen", "Ground", "Prolif", "TR Diff Myo", "Diff Myo", 
                    "Uni Fibro", "Fibro", "Mural", "Endo", "Epi", "Immune", "Unk"))


#################################################################################
## plot cell type
# DimPlot cell type
Idents(xenium.obj) <- xenium.obj@meta.data$cell_type
(xenium.obj_umap_cell_type = 
   DimPlot(xenium.obj, label = T, label.box = F, repel = TRUE, 
           reduction = "umap", shuffle = T, cols =
             unlist(cluster_colors[levels(xenium.obj@active.ident)]),
           group.by = "cell_type")+
    theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_hs_umap_cell_type.svg", 
       plot = xenium.obj_umap_cell_type, width = 5.8, height = 4, 
       device = "svg", path = output_dir)

#################################################################################
## extract annotation for cells visualization in Xenium explorer
clusters_annotation_cell_type = xenium.obj$cell_type
df_cell_type = 
  data.frame(cell_id = names(clusters_annotation_cell_type), 
             group = clusters_annotation_cell_type)
write.csv(x = df_cell_type, file = "../output/tables/dicer1_xenium_hs_annotation_cell_type.csv", 
          row.names = FALSE)

#################################################################################
## find markers for all cell types
# FindAllMarkers
Idents(xenium.obj) <- xenium.obj@meta.data$cell_type
markers_cell_types <- FindAllMarkers(xenium.obj, only.pos = TRUE)
markers_cell_types %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

## print markers for all 13 clusters
# write CSV
write_csv(x = markers_cell_types, file = "dicer1_xenium_hs_obj_cell_type.csv")


#################################################################################
#################################################################################
## subset to mesenchymal lineage
xenium.obj.filter <- 
  subset(xenium.obj, idents = c("Progen", "Ground", "Prolif", "TR Diff Myo", 
                                "Diff Myo", "Fibro"))

# SCTransform for normalization
xenium.obj.filter <- SCTransform(xenium.obj.filter, assay = "Xenium") 
#(, vars.to.regress = "cell_area") could use to regress out cell area
# run PCA
xenium.obj.filter <- RunPCA(xenium.obj.filter, verbose = FALSE)
DimHeatmap(xenium.obj.filter, dims = 1:20, cells = 500, balanced = TRUE)
# determine dims
ElbowPlot(xenium.obj.filter, ndims = 50, reduction = "pca")
# run UMAP
xenium.obj.filter <- RunUMAP(xenium.obj.filter, reduction = "pca", dims = 1:10)
xenium.obj.filter <- FindNeighbors(xenium.obj.filter, reduction = "pca", dims = 1:10)
# find clusters using louvain algorithm
xenium.obj.filter <- FindClusters(xenium.obj.filter, resolution = 0.3)

#################################################################################
## plot cluster
(xenium.obj.filter_umap_cluster = DimPlot(
  xenium.obj.filter, label = T, repel = TRUE, reduction = "umap",
  shuffle = T, cols = colors, group.by = "seurat_clusters") +
  theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_hs_filter_umap_cluster.svg", 
       plot = xenium.obj.filter_umap_cluster, width = 4, height = 4, 
       device = "svg", path = output_dir)

## count cells in UMAP # 57,724
total_cells_filter <- ncol(xenium.obj.filter)
print(total_cells_filter)

# Count cells in each cluster
batch_counts_filter_all <- table(xenium.obj.filter@meta.data$seurat_clusters)
print(batch_counts_filter_all)

#################################################################################
## extract annotation for cells visualization in Xenium explorer
annotations_filter_clusters = xenium.obj.filter$SCT_snn_res.0.2
df_filter = data.frame(cell_id = names(annotations_filter_clusters), 
                       group = annotations_filter_clusters)
write.csv(x = df_filter, file = "../output/tables/dicer1_xenium_hs_filter_annotation_clusters.csv", 
          row.names = FALSE)

#################################################################################
## find markers for 11 clusters
# FindAllMarkers
Idents(xenium.obj.filter) <- xenium.obj.filter@meta.data$seurat_clusters
markers.filter_clusters <- FindAllMarkers(xenium.obj.filter, only.pos = TRUE)
markers.filter_clusters %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# top 5 markers
markers.filter_clusters %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top5.filter_cluster

## plot top 5 markers 
# DoHeatmap (downsampled)
maxcells = min(table(Idents(xenium.obj.filter)))
xenium.obj.filter_top5_genes_clusters = 
  DoHeatmap(subset(xenium.obj.filter, downsample = maxcells),
            features = top5.filter_cluster$gene, raster = TRUE, angle = 45, 
            draw.line = F, label = T) + 
  viridis::scale_fill_viridis(option = "G" ) 
ggsave(filename = "dicer1_xenium_hs_filter_top5_genes_clusters.svg", 
       plot = xenium.obj.filter_top5_genes_clusters, width = 15, height = 14, 
       device = "svg", path = output_dir)

## print markers for all 11 clusters
# write CSV
write_csv(x = markers.filter_clusters, file = "dicer1_xenium_hs_markers_filter_clusters.csv")

#################################################################################
## annotate cluster using broad cell type markers
# Path to YAML file for cell type markers
cell_type_markers_mesenchymal <- read_yaml("../reference/cell_type_markers_homo_sapiens_mesenchymal.yaml")

## plot cell type markers
# DotPlot
(xenium.obj.filter_dotplot_cluster = 
    DotPlot(xenium.obj.filter, features = cell_type_markers_mesenchymal) 
  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  + theme(strip.text.x = element_text(angle = 90, hjust = 0)))
ggsave(filename = "dicer1_xenium_hs_filter_dotplot_cluster.svg", 
       plot = xenium.obj.filter_dotplot_cluster, width = 6, height = 6, 
       device = "svg", path = output_dir)


#################################################################################
## batch and rename FB.clusters based on shared markers
# cluster names
cluster.names.filter_cell_type <- list(
  `0` = "Ground",              # Ground tumor cells
  `1` = "Progen",              # DPT+ Progenitor tumor cells
  `2` = "Progen",       # Proliferative tumor cells
  `3` = "TR Diff Myo",   # Transiting-Differentiated tumor cells
  `4` = "Ground",       # Proliferative tumor cells
  `5` = "Prolif",       # Proliferative tumor cells
  `6` = "Prolif",    # C7+ universal fibroblasts
  `7` = "Progen",          # CXCL10+ fibroblasts
  `8` = "Fibro",
  `9` = "Progen",
  `10` = "Diff Myo"
)
xenium.obj.filter@meta.data$cell_type <- 
  factor(x = 
           cluster.names.filter_cell_type[as.character(xenium.obj.filter$seurat_clusters)],
         levels = c("Progen", "Ground", "Prolif","TR Diff Myo",
                    "Diff Myo", "Fibro", "Mural"))

#################################################################################
## plot cell types
Idents(xenium.obj.filter) <- xenium.obj.filter@meta.data$cell_type
(xenium.obj.filter_clusters_umap_cell_type = 
   DimPlot(xenium.obj.filter, label = T, label.box = F, repel = TRUE,
           reduction = "umap", shuffle = T, cols =
           unlist(cluster_colors[levels(xenium.obj.filter@active.ident)])))
ggsave(filename = "dicer1_xenium_hs_filter_clusters_umap_cell_type.svg", 
       plot = xenium.obj.filter_clusters_umap_cell_type, width = 5.2, height = 4, 
       device = "svg", path = output_dir)

## count cells in UMAP # 52,624
total_cells_filter <- ncol(xenium.obj.filter)
print(total_cells_filter)

# Count cells in each cluster
batch_counts_filter <- table(xenium.obj.filter@meta.data$cell_type)
print(batch_counts_filter)

#################################################################################
## extract annotation for cells visualization in Xenium explorer
annotation_filter_cell_type = xenium.obj.filter$cell_type
df_filter_cell_type = data.frame(cell_id = names(annotation_filter_cell_type), 
                         group = annotation_filter_cell_type)
write.csv(x = df_filter_cell_type, file = "../output/tables/dicer1_xenium_hs_filter_annotation_cell_type.csv", 
          row.names = FALSE)

#################################################################################
## plot key genes genes
(xenium.obj.filter_featureplot = 
   FeaturePlot(xenium.obj.filter, order = T,
               features = 
                 c("PDGFRA", "DPT", "C7", "DES", "MYBPC1", "MKI67")))
ggsave(filename = "dicer1_xenium_hs_filter_featureplot.svg", 
       plot = dicer1_xenium_hs_filter_featureplot, width = 6, height = 8, 
       device = "svg", path = output_dir)


#################################################################################
## find markers for all 6 cell states
# FindAllMarkers
Idents(xenium.obj.filter) <- xenium.obj.filter@meta.data$cell_type
markers.filter_cell_type <- FindAllMarkers(xenium.obj.filter, only.pos = TRUE)
markers.filter_cell_type %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# top 5 markers
markers.filter_cell_type %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15.filter_cell_type

## plot top 15 markers 
# DoHeatmap (downsampled)
maxcells = min(table(Idents(xenium.obj.filter)))
xenium.obj.filter_top15_genes_cell_type = 
  DoHeatmap(subset(xenium.obj.filter, downsample = maxcells),
            features = top15.filter_cell_type$gene, raster = TRUE, angle = 45, 
            draw.line = F, label = T, group.colors = 
              unlist(cluster_colors[levels(xenium.obj.filter@active.ident)])) + 
  viridis::scale_fill_viridis(option = "G" ) 
ggsave(filename = "dicer1_xenium_hs_filter_top15_genes_cell_type.svg", 
       plot = xenium.obj.filter_top5_genes_cell_type, width = 15, height = 14, 
       device = "svg", path = output_dir)

## print markers for all 6 cell states
# write CSV
write_csv(x = markers.filter_cell_type, file = "dicer1_xenium_hs_filter_cell_type.csv")


#################################################################################
#################################################################################
## Pseudotime analyses using tumor clusters
## create file for monocle3
xenium.obj.filter[["UMAP"]] <- xenium.obj.filter[["umap"]]
cds <- 
  as.cell_data_set(xenium.obj.filter, group.by="cell_type")

#################################################################################
# cluster cells
cds <- 
  cluster_cells(cds, reduction_method = c("UMAP"))
# learn the trajectory graph
cds <- learn_graph(cds, use_partition = T, close_loop = F, 
                   learn_graph_control =  
                     list(minimal_branch_len = 30, nn.k = NULL))
# order cells
cds <- order_cells(cds)

# plot_cells with pseudotime
cds_trajecotry_pseudotime = 
  plot_cells(cds, 
             color_cells_by = 
               "pseudotime",
             label_groups_by_cluster = FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             show_trajectory_graph = TRUE,
             group_label_size = 0,
             alpha = 0.5,
             ) +
  theme(legend.position = "right")
ggsave(filename = "dicer1_xenium_hs_cds_trajectory_pseudotime.svg", 
       plot = cds_trajecotry_pseudotime, width = 5, height = 4, 
       device = "svg", path = output_dir)

# plot_cells with clusters
cds_trajecotry_cluster = 
  plot_cells(cds, color_cells_by = 
               "cell_type",
             label_groups_by_cluster = FALSE,
             label_cell_groups=FALSE,
             label_leaves=T,
             label_branch_points=FALSE,
             show_trajectory_graph = TRUE,
             group_label_size = 0,
             alpha = 1)
ggsave(filename = "dicer1_xenium_hs_cds_trajectory_cluster.svg", 
       plot = cds_trajecotry_cluster, width = 5.2, height = 4, 
       device = "svg", path = output_dir)



#################################################################################
#################################################################################
## compare clusters to signatures of RMS
# create the gene list for "Progenitor_RMS"
progenitor_RMS_genes <- c(
  "CD44", "COL3A1", "COL6A2", "COL5A2", "COL6A3", "COL1A1", "POSTN", "COL6A1", "COL5A1", 
  "LAMA4", "DCN", "COL11A1", "COL1A2", "LRP1", "LAMC1", "SCN7A", "CCDC102B", "COL4A1", 
  "VCAN", "COL5A3", "COL4A2", "NAV3", "SULF1", "GSN", "EPS8", "FAM129A", "LAMB1", "EGFR", 
  "EBF1", "ABI3BP", "COL15A1", "NR4A1", "MEOX2", "FRMD4A", "FBN1", "PPAP2B", "BICC1", "TSHZ2", 
  "COL12A1", "AHNAK", "EMP1", "ANXA1", "LOXL2", "RRBP1", "FN1", "TFPI", "MMP2", "LUM", 
  "SERPINH1", "SERPINE2", "RHOC", "LGALS1", "EVA1B", "MGP", "IGFBP6", "SPARC", "EDNRA", "TIMP1", 
  "CMYA5", "ELN", "NRP1", "CPQ", "GALNT10", "ADD3", "FNDC3B", "DHRS3", "SERPINE1", "PDGFRB", 
  "ZFHX4", "EPHA3", "FBLN1", "CACNB4", "SP100", "UTRN", "DLC1", "DUSP1", "FBN2", "FOSB", "CALD1", 
  "FOS", "MMP16", "AKAP12", "PRICKLE1", "SLIT2", "MIR4435-1HG", "ENPP2", "ADAMTS9", "THBS1", 
  "PLOD2", "NEAT1", "CTGF", "COLEC12", "EGR1", "SAT1", "UBC", "KLF6", "PRRX1", "CREB3L1", 
  "FIBIN", "SFRP2", "PTN", "MAP1B", "ALDH1A3", "PPP1R14A", "PDGFRL", "MATN2", "ZFP36", "ID3", 
  "CRABP2", "NR4A2", "OGN", "PCOLCE", "THY1", "S100A4", "TNMD", "HIC1", "SELM", "PLAC9", 
  "PTGS2", "NOV", "TAGLN2", "ENC1", "THBS4", "APOE", "ATP1B1", "TGFB1I1", "ITM2C", "C11orf96", 
  "LGALS3", "CCDC80", "HES4", "ID1", "PAG1", "NEFM", "NOTCH3", "SOCS3", "CDC42EP5", "TGFBI", 
  "CFH", "MTUS1", "IGFBP7", "SPARCL1", "CHODL", "NUPR1", "HTRA1", "IGFBP5", "FSTL1", "ANXA2", 
  "FTH1", "SH3BGRL3", "DAB2", "TSC22D1", "BST2", "HLA-A", "ISLR", "SQSTM1", "B2M", "RRAS", 
  "DLK1", "LRRC17", "ERRFI1", "EMP3", "S100A6", "CAV1", "MT2A", "APLP2", "CLIC1", "CCL2", "TXNIP", 
  "S100A10", "FTL"
)

# filter out missing genes
missing_genes <- setdiff(progenitor_RMS_genes, rownames(xenium.obj.filter))
present_genes <- intersect(progenitor_RMS_genes, rownames(xenium.obj.filter))

# Print how many genes are missing or present
cat("Number of missing genes:", length(missing_genes), "\n")
cat("Number of present genes:", length(present_genes), "\n")

# proceed with filtering (make sure genes are present in the dataset)
progenitor_RMS_genes_filtered <- present_genes

# create the gene list for "Differentiated_RMS"
differentiated_RMS_genes <- c(
  "MYL4", "MACF1", "ARPP21", "CHRNA1", "RBM24", "CYB5R1", "TNNT3", "ERBB3", "FNDC5", "BLCAP", 
  "DCLK1", "SRPK3", "RYR1", "ACTC1", "CHRND", "TNNI1", "ZNF106", "RASSF4", "TTN", "COBL", 
  "MYH3", "DES", "MYOG", "TNNT2", "NES", "STAC3", "MYL1", "MYLPF", "ACTN2", "FILIP1", "NEB", 
  "GPC1", "ATP2B1", "BIN1", "MEF2C", "SYNPO2L", "SHD", "HSPB3", "KLHL41", "TNNC2", "FLNC", 
  "TNNC1", "TNNI2", "TPM1", "SMPX", "DOK7", "CKM", "KREMEN2", "38231", "CTD-2545M3.8", "DLG2", 
  "LRRN1", "UNC45B", "NEXN", "PDE4DIP", "MYOZ2", "SEPW1", "ACTA1", "HIPK3", "COL25A1", "LMO7", 
  "UBE2E3", "A2M", "PALLD", "SERPINI1", "DYSF", "CKB", "CDH15", "SPTBN1", "CASZ1", "CAP2", 
  "CHRNB1", "NCALD", "MTSS1", "SETD7", "PDLIM3", "ITGA6", "CASQ2", "CHRNG", "SMYD1", "PDLIM5", 
  "RGR", "PLS3", "RP1-302G2.5", "VGLL2", "NEU4", "IL17B", "SCRIB", "MAP1A", "TSPAN33", "INPPL1", 
  "SIRT2", "PRKAR1A", "HSPB1", "FITM1", "PRDX6", "RGS16", "PDLIM7", "S100A4", "BTG1", "LMOD3", 
  "CRYAB", "TUBB6", "KLF5", "APOBEC2", "PHLDA1", "ZBTB18", "METRN", "GLRX", "C1orf105", "CSRP3", 
  "MYOD1", "ENO3", "MYBPH", "NUAK1", "MAP3K7CL", "EHBP1L1", "NGFR", "TPM2", "COX6A2", "HRC", 
  "SHISA2", "MYL6B", "LMNA", "NNMT", "S100A13", "TNNT1", "VASH2", "RTN4", "FAM178A", "GPC4", 
  "LDB3", "LAMA2", "FRAS1", "CACNA2D1", "SYNPO2", "SASH1", "ST6GALNAC5", "MYOM1", "INPP4B", "HIVEP2", 
  "FAM65B", "RP11-1L12.3", "TRDN", "TMEM38A", "FGF7", "CCDC141", "AKAP6", "HES6", "RNF217", "CADM2", 
  "RP11-14N7.2", "SHISA9", "AC007970.1", "SHROOM3", "CACNA1S", "AFAP1L1", "NCAM1", "CLSTN2", 
  "HIST1H2AC", "NCOA1", "DMPK", "ANTXR2", "MYO10", "RUNX1", "TRIM55", "DST", "PTGFR", "SORBS2", 
  "ACVR2A", "SGCD", "TMEM108", "PKIA", "MIR503HG", "DEPTOR", "MYPN", "LIMA1", "SVIL", "MYH8", 
  "TTLL7", "FGF13", "WWTR1", "JPH2", "MEG3", "CAV3", "ANK2", "ITGB6", "LZTS1", "RBM20", "DMD", 
  "MAP2K1", "ASB5", "SNTB1", "DAPK2", "ARHGAP24", "ANKRD1", "OBSCN", "LINC00882", "NPNT", "SEMA6B", 
  "PRUNE2", "MAMDC2", "MYOM3", "FBXO32", "MYH7B", "MYO18B", "CTNNA3", "VAV3", "TP63", "REEP1", 
  "SEMA6A", "SRL", "ALPK2", "TEAD4", "SCN3B"
)

# filter out missing genes
missing_genes <- setdiff(differentiated_RMS_genes, rownames(xenium.obj.filter))
present_genes <- intersect(differentiated_RMS_genes, rownames(xenium.obj.filter))

# print how many genes are missing or present
cat("Number of missing genes:", length(missing_genes), "\n")
cat("Number of present genes:", length(present_genes), "\n")

# proceed only if some genes are present
differentiated_RMS_genes_filtered <- present_genes

# create the gene list for "Proliferative_RMS"
proliferative_RMS_genes <- c(
  "ASPM", "TOP2A", "MIS18BP1", "SMC4", "MKI67", "CENPF", "CENPE", "NUCKS1", "TUBA1B", "SYNE2", 
  "ATAD2", "HELLS", "DHFR", "CENPK", "BUB3", "ARHGAP11A", "TPX2", "CCNB1", "HMGB3", "UBE2S", 
  "ARL6IP1", "CDC25B", "CDCA3", "KPNA2", "RNF26", "TACC3", "BIRC5", "PLK1", "UBE2C", "CCNB2", 
  "PRC1", "CDK1", "PIF1", "KIF11", "CKAP2", "CEP55", "MZT1", "RACGAP1", "KNSTRN", "TUBA1C", 
  "NCAPD2", "SAPCD2", "DBF4", "AURKA", "KIF14", "BUB1B", "DLGAP5", "HJURP", "AURKB", "CKAP2L", 
  "FAM83D", "KIF20B", "CDCA2", "HMMR", "KIF4A", "CCNF", "TUBB4B", "KIF2C", "KIF20A", "RAD21", 
  "CKAP5", "CENPA", "KIF22", "CDC20", "NDC80", "CCNA2", "TTK", "GTSE1", "CDKN3", "BUB1", "NUF2", 
  "NCAPG", "PSRC1", "PRR11", "NEK2", "TROAP", "HMGB2", "UBE2T", "CKS2", "ECT2", "LBR", "KIF23", 
  "DEPDC1", "CDCA8", "PTMA", "DIAPH3", "EZH2", "TUBB", "NNAT", "SPC25", "GPC3", "GAPDH", "FGFR4", 
  "HIST1H4C", "ACTG1", "C21orf58", "ARHGAP11B", "HMGA2", "HSP90AA1", "ANP32E", "HMGN2", "HNRNPA2B1", 
  "ANLN", "FOXM1", "DTYMK", "PTMS", "CDC25C", "KIFC1", "HP1BP3", "MAD2L1", "TMPO", "CKS1B", "GAS2L3", 
  "SPDL1", "HDGF", "TRIM59", "SUN2", "CALM2", "PTTG1", "STMN1", "RANGAP1", "CCDC88A", "HMGB1", "NDE1", 
  "GPSM2", "G2E3", "LMNB1", "DNMT1", "MCM7", "CASC5", "RRM2", "FEN1", "E2F1", "MXD3", "CENPU", "DUT", 
  "KIF18A", "MYBL2", "RRM1", "MCM3", "CLSPN", "SGOL1", "REEP4", "HIST1H3D", "PCNA", "MCM6", "FBXO5", 
  "SGOL2", "TYMS", "FAM111B", "ZWINT", "KIAA0101", "H2AFX", "GINS2", "FAM64A", "NUSAP1", "CIT"
)

# filter out missing genes
missing_genes <- setdiff(proliferative_RMS_genes, rownames(xenium.obj.filter))
present_genes <- intersect(proliferative_RMS_genes, rownames(xenium.obj.filter))

# print how many genes are missing or present
cat("Number of missing genes:", length(missing_genes), "\n")
cat("Number of present genes:", length(present_genes), "\n")

# proceed only if some genes are present
proliferative_RMS_genes_filtered <- present_genes

# add module score
xenium.obj.filter <- AddModuleScore(
  object = xenium.obj.filter,
  features = list(progenitor_RMS_genes_filtered),  # Add the gene signature
  name = "Progenitor_RMS",  # Name the new metadata column
  ctrl = 15  # Disable random background sampling
)
xenium.obj.filter <- AddModuleScore(
  object = xenium.obj.filter,
  features = list(proliferative_RMS_genes_filtered),  # Add the gene signature
  name = "Proliferative_RMS",  # Name the new metadata column
  ctrl = 15  # Disable random background sampling
)
xenium.obj.filter <- AddModuleScore(
  object = xenium.obj.filter,
  features = list(differentiated_RMS_genes_filtered),  # Add the gene signature
  name = "Differentiated_RMS",  # Name the new metadata column
  ctrl = 15  # Disable random background sampling
)
# visualize the module score using VlnPlot
VlnPlot(xenium.obj.filter, features, pt.size = 0.0)

# visualize the scores on a UMAP or PCA plot
signature_prog = FeaturePlot(xenium.obj.filter, order = T, features = "Progenitor_RMS1") +
  scale_color_gradient2(low = "#f1a340" , mid = "#fee0b6", high = "#542788")
# Step 6: Optionally, visualize the scores on a UMAP or PCA plot
signature_prolif =FeaturePlot(xenium.obj.filter, order = T, features = "Proliferative_RMS1") +
  scale_color_gradient2(low = "#f1a340" , mid = "#fee0b6", high = "#542788")
# Step 6: Optionally, visualize the scores on a UMAP or PCA plot
signature_diff =FeaturePlot(xenium.obj.filter, order = T, features = "Differentiated_RMS1") +
  scale_color_gradient2(low = "#f1a340" , mid = "#fee0b6", high = "#542788")

ggsave(filename = "dicer1_xenium_hs_signature_prog.png", 
       plot = signature_prog, width = 5, height = 4, 
       device = "png", path = output_dir)
ggsave(filename = "dicer1_xenium_hs_signature_prolif.png", 
       plot = signature_prolif, width = 5, height = 4, 
       device = "png", path = output_dir)
ggsave(filename = "dicer1_xenium_hs_signature_diff.png", 
       plot = signature_diff, width = 5, height = 4, 
       device = "png", path = output_dir)

#################################################################################
#################################################################################
## calculate cell numbers per cluster for all 16 tumors

## analysis of all cell numbers per type
# cell numbers
xenium.obj.filter@meta.data %>%
  group_by(tumour) %>%
  summarise(n = n())
xenium.obj.filter_clusters_sample = xenium.obj.filter@meta.data %>%
  group_by(tumour, cell_type) %>%
  summarise(n = n())
xenium.obj.filter_clusters_sample$tumour <- factor(xenium.obj.filter_clusters_sample$tumour)
xenium.obj.filter_clusters_sample <- xenium.obj.filter_clusters_sample %>%
  group_by(tumour) %>%
  mutate(percentage = n / sum(n) * 100)

# plotting of cell numbers per cluster
tumour_order <- c(43, 44, 45, 47, 5, 8, 9, 12, 15, 46, 48, 49, 7, 11, 4, 18)
(xenium.obj.filter_clusters_sample_bar_plot <- ggplot(xenium.obj.filter_clusters_sample, 
                                                  aes(x = factor(tumour, levels = tumour_order), 
                                                      y = percentage, 
                                                      fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = colors) +
    labs(x = "Tumour", y = "Percentage", title = "Stacked Bar Plot of Tumours") +
    theme_minimal())
ggsave(filename = "dicer1_xenium_hs_filter_clusters_sample_bar_plot.svg", 
       plot = xenium.obj.filter_clusters_sample_bar_plot , width = 10, height = 10, 
       device = "svg", path = output_dir)