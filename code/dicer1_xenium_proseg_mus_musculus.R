---
title: "Dicer1_Xenium_proseg_mus_musculus"
author: "Felix Kommoss"
date: "2024-10-23"
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

#################################################################################
## set seed
set.seed(404)

#################################################################################
# output directory
output_dir <- "../output/plots/"

#################################################################################
## load data 
# readRDS
xenium.obj <- readRDS("../data/dicer1_xenium_mus_musculus.rds")

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
  "Prox Tub" = "#B2A000",
  "Loop Henle" = "#00C091",
  "Endo" = "#F8766D",
  "Immune" = "#B186FF",
  "Fibro" = "#CE9500",
  "Collect Duct" = "#F27D53",
  "Trans Epi" = "#FF6A98",
  "Uni Fibro" = "#00BAE0",
  "Podo" = "#D973FC",
  "Ground" = "#00BB4E",
  "Progen" = "#C77CFF",
  "Diff Myo" = "#00C1A3",
  "Prolif" = "#FF61CC",
  "Diff Fibro" = "#A3A500",
  "TR Diff Myo" = "#7099FF",
  "TR Diff Chondro" = "#00A9FF",
  "Mural" = "#39B606",
  "Uni Fibro (Pi16)" = "#B186FF",
  "Uni Fibro (Col15a1)" = "#FD7083",
  "Cycle" = "#00BAE0",
  "Mac Densa" = "#39B607",
  "Adipo" = "#C09B06",
  "Prolif 1" = "#FF61CC",
  "Prolif 2" = "#FF65AA"
)

colors_syngenic <- c( "#C77CFF", "#FF61CC", "#E78AC3", "#00A9FF",
                      "#F8766D", "#00C0B6", "#C09B06")

#################################################################################
## QC
# remove cells with less than 3 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 3)
# remove cells with less than 3 genes
xenium.obj <- subset(xenium.obj, subset = nFeature_Xenium > 3)

# cell counts for xenium.obj
cell_count_core <- table(xenium.obj@meta.data$core)
print(cell_count_core)
cell_count_Lesion_ID <- table(xenium.obj@meta.data$Lesion_ID)
print(cell_count_Lesion_ID)
cell_count_Sample_type <- table(xenium.obj@meta.data$Sample_type)
print(cell_count_Sample_type)
cell_count_Gender <- table(xenium.obj@meta.data$Gender)
print(cell_count_Gender)


#################################################################################
#################################################################################
# Control 
#################################################################################
## subset xenium.obj to samples with control type genotype
# subset
xenium.obj.wt <- xenium.obj [, xenium.obj$Inclusion_control != "no"]
# SCTransform
xenium.obj.wt <- SCTransform(xenium.obj.wt, assay = "Xenium")
# RunPCA
xenium.obj.wt <- RunPCA(xenium.obj.wt, verbose = FALSE)
DimHeatmap(xenium.obj.wt, dims = 1:20, cells = 500, balanced = TRUE)
# determine dims with ElbowPlot
ElbowPlot(xenium.obj.wt, ndims = 50, reduction = "pca")
# RunUMAP
xenium.obj.wt <- RunUMAP(xenium.obj.wt, reduction = "pca", dims = 1:30)
xenium.obj.wt <- FindNeighbors(xenium.obj.wt, reduction = "pca", dims = 1:30)
# FindClusters using louvain algorithm
xenium.obj.wt <- FindClusters(xenium.obj.wt, resolution = 0.6)

#################################################################################
## plot clusters and lesion_ID
# DimPLot clusters
(xenium.obj.wt_umap_clusters <- 
   DimPlot(xenium.obj.wt, label = T, repel = TRUE,
           reduction = "umap", shuffle = T, cols = colors,
           group.by = "seurat_clusters") +
   theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_mm_wt_umap_clusters.svg", 
       plot = xenium.obj.wt_umap_clusters, width = 5, height = 4, 
       device = "svg", path = output_dir)
# DimPLot lesion_ID
(xenium.obj.wt_umap_lesion_ID <- 
    DimPlot(xenium.obj.wt, label = F, repel = TRUE,reduction = "umap", 
            shuffle = T, cols = colors_secondary, group.by = "Lesion_ID") +
    theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_mm_wt_umap_lesion_id.svg", 
       plot = xenium.obj.wt_umap_lesion_ID, width = 5.5, height = 4, 
       device = "svg", path = output_dir)

## count cells in DimPlot 
# n = 15414
total_cells_wt <- ncol(xenium.obj.wt)
print(total_cells_wt)

#################################################################################
## extract annotation for cells visualization in Xenium explorer
# write CSV
clusters_wt_annotation <- xenium.obj.wt$SCT_snn_res.0.6
df_wt_clusters <- data.frame(cell_id = names(clusters_wt_annotation), 
                             group = clusters_wt_annotation)
write.csv(x = df_wt_clusters , file = "../output/tables/dicer1_xenium_mm_wt_annotation_clusters.csv", 
          row.names = FALSE)

#################################################################################
## find markers for all 22 cluster
# FindAllMarkers
Idents(xenium.obj.wt) <- xenium.obj.wt@meta.data$seurat_clusters
markers_wt_clusters <- FindAllMarkers(xenium.obj.wt, only.pos = TRUE)
markers_wt_clusters %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

## print markers for all 22 clusters
# write CSV
write_csv(x = markers_wt_clusters, file = "../output/tables/dicer1_xenium_mm_wt_markers_clusters.csv")

#################################################################################
## annotate cluster using broad cell type markers
# Path to YAML file for cell type markers
cell_type_markers_mus_musculus_kidney <- 
  read_yaml("../reference/cell_type_markers_mus_musculus_kidney.yaml")

## plot cell type markers
# DotPlot
(xenium.obj.wt_dotplot_clusters <- 
    DotPlot(xenium.obj.wt, features = cell_type_markers_mus_musculus_kidney) 
  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  + theme(strip.text.x = element_text(angle = 90, hjust = 0)))
ggsave(filename = "dicer1_xenium_mm_wt_dotplot_clusters.svg", 
       plot = xenium.obj.wt_dotplot_clusters, width = 10, height = 8, 
       device = "svg", path = output_dir)

#################################################################################
## batch and rename clusters based on broad cell type markers
# cluster names
xenium.obj.wt_cell_type <- list(
  "0" = "Prox Tub",       # Proximal tubule
  "1" = "Prox Tub",       # Proximal tubule
  "2" = "Loop Henle",     # Loop of Henle
  "3" = "Loop Henle",     # Loop of Henle
  "4" = "Endo",           # Endothelium
  "5" = "Immune",         # Immune cells
  "6" = "Fibro",          # Fibroblasts
  "7" = "Prox Tub",       # Proximal tubule
  "8" = "Prox Tub",       # Proximal tubule
  "9" = "Loop Henle",     # Loop of Henle
  "10" = "Mural",         # Mural cells
  "11" = "Loop Henle",    # Loop of Henle
  "12" = "Fibro",         # Fibroblasts
  "13" = "Collect Duct",  # Collecting duct
  "14" = "Trans Epi",     # Transitional epithelium
  "15" = "Collect Duct",  # Collecting duct
  "16" = "Loop Henle",    # Loop of Henle
  "17" = "Endo",          # Endothelium
  "18" = "Loop Henle",    # Loop of Henle
  "19" = "Prox Tub",      # Proximal tubule
  "20" = "Uni Fibro",     # Universal fibroblasts
  "21" = "Podo"           # Podocytes
)

xenium.obj.wt@meta.data$cell_type <- 
  factor(x = 
           xenium.obj.wt_cell_type[as.character(xenium.obj.wt$seurat_clusters)],
         levels = c("Uni Fibro", "Fibro", "Mural", "Endo", "Podo", "Prox Tub", 
                    "Loop Henle", "Collect Duct", "Trans Epi", "Immune"))

#################################################################################
## plot cell type
# DimPlot cell types
(xenium.obj.wt_umap_cell_type <- DimPlot(
  xenium.obj.wt, label = T, repel = TRUE,reduction = "umap",
  shuffle = T, cols = cluster_colors, group.by = "cell_type") +
   theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_mm_wt_umap_cell_type.svg", 
       plot = xenium.obj.wt_umap_cell_type, width = 5, height = 4, 
       device = "svg", path = output_dir)

#################################################################################
## extract annotation for cells visualization in Xenium explorer
# write CSV
cell_types_wt_annotation <- xenium.obj.wt$cell_type
df_wt_batch <- data.frame(cell_id = names(cell_types_wt_annotation), 
                          group = cell_types_wt_annotation)
write.csv(x = df_wt_batch, file = "../output/tables/dicer1_xenium_mm_wt_annotation_cell_type.csv", 
          row.names = FALSE)

#################################################################################
## find markers for cell types
# FindAllMarkers
Idents(xenium.obj.wt) <- xenium.obj.wt@meta.data$cell_type
markers_wt_cell_type <- FindAllMarkers(xenium.obj.wt, only.pos = TRUE)
markers_wt_cell_type %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# top 5 markers
markers_wt_cell_type %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5_wt_cell_type

## plot top 5 markers
# DoHeatmap (downsampled)
maxcells.wt <- min(table(Idents(xenium.obj.wt)))
xenium.obj.wt_top5_genes_cell_type <- 
  DoHeatmap(subset(xenium.obj.wt, downsample = maxcells.wt), 
            features = top5_wt_cell_type$gene, raster = T, angle = 45, 
            draw.line = F, label = T, group.colors = 
              unlist(cluster_colors[levels(xenium.obj.wt@active.ident)])) + 
  viridis::scale_fill_viridis(option = "G") 
ggsave(filename = "dicer1_xenium_mm_wt_top5_genes_cell_type.svg", 
       plot = xenium.obj.wt_top5_genes_cell_type, width = 15, height = 7, 
       device = "svg", path = output_dir)

## print markers for cell types
# write CSV
write_csv(x = markers_wt_cell_type, file = "../output/tables/dicer1_xenium_mm_wt_markers_cell_type.csv")


#################################################################################
#################################################################################
#HDT 
#################################################################################
## subset xenium.obj to samples with mutant genotype
# subset
xenium.obj.mut <- xenium.obj [, xenium.obj$Inclusion_mutant_primary != "no"]
# SCTransform
xenium.obj.mut <- SCTransform(xenium.obj.mut, assay = "Xenium")
# RunPCA
xenium.obj.mut <- RunPCA(xenium.obj.mut, verbose = FALSE)
DimHeatmap(xenium.obj.mut, dims = 1:20, cells = 500, balanced = TRUE)
# determine dims with ElbowPlot
ElbowPlot(xenium.obj.mut, ndims = 50, reduction = "pca")
# RunUMAP
xenium.obj.mut <- RunUMAP(xenium.obj.mut, reduction = "pca", dims = 1:30)
xenium.obj.mut <- FindNeighbors(xenium.obj.mut, reduction = "pca", dims = 1:30)
# FindClusters using louvain algorithm
xenium.obj.mut <- FindClusters(xenium.obj.mut, resolution = 0.6)

FeaturePlot(xenium.obj.mut, features = c("nFeature_Xenium"))

#################################################################################
## exclude cluster 8
## low feature count detected and absence of positive markers above background expression
# subset
xenium.obj.mut <- xenium.obj.mut [, xenium.obj.mut$seurat_clusters != 8]
# SCTransform
xenium.obj.mut <- SCTransform(xenium.obj.mut, assay = "Xenium")
# RunPCA
xenium.obj.mut <- RunPCA(xenium.obj.mut, verbose = FALSE)
DimHeatmap(xenium.obj.mut, dims = 1:20, cells = 500, balanced = TRUE)
# determine dims with ElbowPlot
ElbowPlot(xenium.obj.mut, ndims = 50, reduction = "pca")
# RunUMAP
xenium.obj.mut <- RunUMAP(xenium.obj.mut, reduction = "pca", dims = 1:30)
xenium.obj.mut <- FindNeighbors(xenium.obj.mut, reduction = "pca", dims = 1:30)
# FindClusters using louvain algorithm
xenium.obj.mut <- FindClusters(xenium.obj.mut, resolution = 0.6)

#################################################################################
## plot clusters, lesion_ID and sample type
# DimPlot clusters
(xenium.obj.mut_umap_clusters <- DimPlot(
  xenium.obj.mut, label = T, repel = TRUE, reduction = "umap",
  shuffle = T, cols = colors, group.by = "seurat_clusters") +
   theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_mm_mut_umap_clusters.svg", 
       plot = xenium.obj.mut_umap_clusters, width = 5.5, height = 4, 
       device = "svg", path = output_dir)
# DimPlot lesion_ID
(xenium.obj.mut_umap_lesion_ID <- DimPlot(
  xenium.obj.mut, label = F, repel = TRUE, reduction = "umap",
  shuffle = T, cols = colors, group.by = "Lesion_ID") + 
    guides(colour = guide_legend(ncol = 2)) +
    theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_mm_mut_umap_lesion_id.svg", 
       plot = dicer1_xenium_mm_mut_umap_lesion_ID, width = 9, height = 4, 
       device = "svg", path = output_dir)
# DimPlot histotype
(xenium.obj.mut_umap_sample_type <- DimPlot(
  xenium.obj.mut, label = F, repel = TRUE, reduction = "umap",
  shuffle = T, cols = colors_secondary, group.by = "Histotype") +
    theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_mm_mut_umap_histotype.svg", 
       plot = dicer1_xenium_mm_mut_umap_sample_type, width = 7.5, height = 4, 
       device = "svg", path = output_dir)

## count cells in DimPlot
# n = 250153
total_cells_mut <- ncol(xenium.obj.mut)
print(total_cells_mut)

#################################################################################
## extract annotation for cells visualization in Xenium explorer
# write CSV
clusters_mut_annotation <- xenium.obj.mut$SCT_snn_res.0.6
df_mut_clusters <- data.frame(cell_id = names(clusters_mut_annotation), 
                              group = clusters_mut_annotation)
write.csv(x = df_mut_clusters, file = "../output/tables/dicer1_xenium_mm_mut_annotation_clusters.csv", 
          row.names = FALSE)

#################################################################################
## find markers for all 27 cluster
# FindAllMarkers
Idents(xenium.obj.mut) <- xenium.obj.mut$seurat_clusters
markers_mut_clusters <- FindAllMarkers(xenium.obj.mut, only.pos = TRUE)
markers_mut_clusters %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

## print markers for all 25 clusters
# write CSV
write_csv(x = markers_mut_clusters, file = "../output/tables/dicer1_xenium_mm_mut_markers_clusters.csv")

#################################################################################
## annotate cluster using broad cell type markers
# Path to YAML file for cell type markers
cell_type_markers_mus_musculus_kidney <- 
  read_yaml("../reference/cell_type_markers_mus_musculus_kidney.yaml")

## plot cell type markers
# DotPlot
Idents(xenium.obj.mut) <- xenium.obj.mut$seurat_clusters
(xenium.obj.mut_dotplot_all <- 
    DotPlot(xenium.obj.mut, features = cell_type_markers_mus_musculus_kidney) 
  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  + theme(strip.text.x = element_text(angle = 90, hjust = 0)))
ggsave(filename = "dicer1_xenium_mm_mut_dotplot_all.svg", 
       plot = xenium.obj.mut_dotplot_all, width = 11, height = 10, 
       device = "svg", path = output_dir)

#################################################################################
## batch and rename clusters based on broad cell type markers
# cluster names
xenium.obj.mut_cell_type <- list(
  `0` = "Progen",        # Progenitor
  `1` = "Ground",        # Ground
  `2` = "Ground",        # Ground
  `3` = "Immune",        # Immune cells
  `4` = "Prox Tub",      # Proximal tubule
  `5` = "Prox Tub",      # Proximal tubule
  `6` = "Fibro",         # Fibroblasts
  `7` = "Endo",          # Endothelium
  `8` = "Prolif",        # Proliferative
  `9` = "Diff Myo",      # Differentiated myogenic
  `10` = "Trans Epi",    # Transitional epithelium
  `11` = "Endo",         # Endothelium
  `12` = "Loop Henle",   # Loop of Henle
  `13` = "Diff Myo",     # Differentiated myogenic
  `14` = "Loop Henle",   # Loop of Henle
  `15` = "Diff Myo",     # Differentiated myogenic
  `16` = "Prox Tub",     # Proximal tubule
  `17` = "Progen",       # Progenitor
  `18` = "Immune",       # Immune cells
  `19` = "Collect Duct", # Collecting duct
  `20` = "Ground",       # Ground
  `21` = "Immune",       # Immune cells
  `22` = "TR Diff Myo",  # Transiting differentiated myogenic
  `23` = "Collect Duct", # Collecting duct
  `24` = "Mural",        # Mural cells
  `25` = "Podo",         # Podocytes
  `26` = "Diff Myo"      # Differentiated myogenic
)

xenium.obj.mut@meta.data$cell_type <- 
  factor(x = 
           xenium.obj.mut_cell_type[as.character(xenium.obj.mut$seurat_clusters)],
         levels = c("Progen", "Ground", "Prolif", "TR Diff Myo", "Diff Myo", 
                    "Fibro", "Mural", "Endo", "Podo", "Prox Tub", "Loop Henle", 
                    "Collect Duct", "Trans Epi", "Immune"))

# count cells in cell types
xenium.obj.mut_cell_count_cell_type <- table(xenium.obj.mut@meta.data$cell_type)
print(xenium.obj.mut_cell_count_cell_type)

#################################################################################
## plot cell types
# DimPLot
Idents(xenium.obj.mut) <- xenium.obj.mut$cell_type
(xenium.obj.mut_umap_cell_type <- DimPlot(
  xenium.obj.mut, label = T, repel = TRUE, reduction = "umap",
  shuffle = T, group.by = "cell_type", cols = 
    unlist(cluster_colors[levels(xenium.obj.mut@active.ident)])) +
    theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_mm_mut_umap_cell_type.svg", 
       plot = xenium.obj.mut_umap_cell_type, width = 5.7, height = 4, 
       device = "svg", path = output_dir)

#################################################################################
## subset to MSC lineage
# subset
xenium.obj.mut.subset <- 
  subset(xenium.obj.mut, idents = c("Progen", "Ground", "Prolif", "TR Diff Myo",
                                    "Diff Myo", "Fibro", "Mural"))

#################################################################################
## plot key genes
# FeaturePlot
(xenium.obj.mut.subset_feature_plot = 
   FeaturePlot(xenium.obj.mut.subset, 
               order = T, 
               features = c("Mfap4", "Itga8", "Myoz2", "Cenpf", "Cfh", "Myh11"), 
               keep.scale = "all"))
ggsave(filename = "dicer1_xenium_mm_mut_subset_feature_plot.svg", 
       plot = xenium.obj.mut.subset_feature_plot, width = 6.66, 
       height = 8, device = "svg", path = output_dir)

# VlnPlot
(xenium.obj.mut.subset_Vln_plot = 
    VlnPlot(xenium.obj.mut.subset, 
            pt.size = 0, cols = colors_secondary, stack = TRUE, sort = F, 
            flip = TRUE, features = 
              c("Pi16", "Mfap4", "Itga8", "Cfh", "Myh11", "Des", "Car3", "Myoz2", 
                "Cenpf", "Sox9")))
ggsave(filename = "dicer1_xenium_mm_mut_subset_vln_plot.svg", 
       plot = xenium.obj.mut.subset_Vln_plot, width = 9, 
       height = 6, device = "svg", path = output_dir)

#################################################################################
## extract annotation for cells visualization in Xenium explorer
# write CSV
cell_type_mut_annotation <- xenium.obj.mut$cell_type
df_mut_cell_type <- data.frame(cell_id = names(cell_type_mut_annotation), 
                               group = cell_type_mut_annotation)
write.csv(x = df_mut_cell_type, file = "../output/tables/dicer1_xenium_mm_mut_annotation_cell_type.csv", 
          row.names = FALSE)

#################################################################################
## find markers for cell types
# FindAllMarkers
Idents(xenium.obj.mut) <- xenium.obj.mut$cell_type
markers_mut_cell_type <- FindAllMarkers(xenium.obj.mut, only.pos = TRUE)
markers_mut_cell_type %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

## print markers for cell types
# write CSV
write_csv(x = markers_mut_cell_type, file = "../output/tables/dicer1_xenium_mm_mut_markers_cell_types.csv")

#################################################################################
#################################################################################
## Pseudotime analysis using monocle3 on xenium.obj.mut.subset

#################################################################################
## create file for monocle3
# cds
xenium.obj.mut.subset[["UMAP"]] <- xenium.obj.mut.subset[["umap"]]
cds <- 
  as.cell_data_set(xenium.obj.mut.subset, group.by="cell_type")

#################################################################################
## cluster cells
# cluster_cells
cds <- 
  cluster_cells(cds, reduction_method = c("UMAP"))

#################################################################################
## learn the trajectory graph
#learn_graph
cds <- learn_graph(cds, use_partition = T, close_loop = F, 
                   learn_graph_control =  list(minimal_branch_len = 40, nn.k = NULL))

# order_cells
cds <- order_cells(cds)

#################################################################################
## plot annotations with trajectory graph
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
             alpha = 0.5) +
  theme(legend.position = "right")
ggsave(filename = "dicer1_xenium_mm_cds_trajectory_pseudotime.svg", 
       plot = cds_trajecotry_pseudotime, width = 5, height = 4, 
       device = "svg", path = output_dir)

# plot_cells with clusters
cds_trajecotry_cluster = 
  plot_cells(cds, color_cells_by = 
               "cell_type",
             label_groups_by_cluster = FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             show_trajectory_graph = TRUE,
             group_label_size = 0,
             alpha = 0.5)
ggsave(filename = "dicer1_xenium_mm_cds_trajectory_cluster.svg", 
       plot = cds_trajecotry_cluster, width = 6, height = 4, 
       device = "svg", path = output_dir)