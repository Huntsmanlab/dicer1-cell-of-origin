### CONTROLS
# Script to normalize data via Seurat's integration, then cluster to identify cell types and find markers

# Required packages
library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
library(yaml)

# For reproducibility
set.seed(1337)

# Read in our previously QC'd data
metadata <- read.table(file = "../metadata/dicer1_scrna_sample_metadata.csv", sep = ",", header = TRUE)
sce_qc <- readRDS(file = "../output/scrna/dicer1_scrna_sce_qc.rds")

# We want to use gene symbols
rownames(sce_qc) <- rowData(sce_qc)$Symbol

# Uniquely bar code cells for Seurat by combining bar code with sample ID
cell_ids <- paste(sce_qc$Sample, sce_qc$Barcode, sep = ".")
colnames(sce_qc) <- cell_ids

### UNFILTERED
# Subset to the control samples
keep_cells <- sce_qc$genotype == "Control"
sce_qc_filter_sample <- sce_qc[, keep_cells]

# Convert to a Seurat object
sce_seurat_sample <- CreateSeuratObject(counts = counts(sce_qc_filter_sample), names.field = 1, names.delim = "/")

# Add our metadata to the Seurat object
seurat_metadata_sample = list("Sample" = sce_qc_filter_sample$Sample,
                              "Barcode" = sce_qc_filter_sample$Barcode,
                              "Genotype" = sce_qc_filter_sample$genotype,
                              "Time" = sce_qc_filter_sample$time,
                              "Tomato" = sce_qc_filter_sample$subsets_tomato_sum,
                              "mito_percent" = sce_qc_filter_sample$subsets_mito_percent,
                              "Name" = sce_qc_filter_sample$name,
                              "Delta" = sce_qc_filter_sample$delta_days,
                              "ribo_percent" = sce_qc_filter_sample$subsets_ribo_percent)
sce_seurat_sample <- AddMetaData(object = sce_seurat_sample, metadata = seurat_metadata_sample)

# Split into one Seurat object per sample and normalize
seurat_list_sample <- SplitObject(sce_seurat_sample, split.by = "Sample")

# Normalize each sample with Seurat's method
seurat_list_sample <- lapply(seurat_list_sample, function(seurat_obj) NormalizeData(seurat_obj, verbose = FALSE))

# Identify variable features (genes) per sample; 2,000 is roughly in line with the top 20% of HVGs
seurat_list_sample <- lapply(seurat_list_sample, function(seurat_obj) FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000))

# Select variable features by their rank across samples, prioritizing conserved features
features_sample <- SelectIntegrationFeatures(object.list = seurat_list_sample, nfeatures = 2000)

# Do per-sample PCA so we can use rPCA to integrate
seurat_list_sample <- lapply(seurat_list_sample, function(seurat_obj) ScaleData(seurat_obj, features = features_sample))
seurat_list_sample <- lapply(seurat_list_sample, function(seurat_obj) RunPCA(seurat_obj, features = features_sample))

# Find anchors across pairwise combinations of our datasets: pairs of similar/conserved cells
anchors_sample <- FindIntegrationAnchors(object.list = seurat_list_sample, anchor.features = features_sample, reduction = "rpca")

# Create the integrated dataset, using the anchors we found; default features used are those we used to find anchors
sce_integrated_sample <- IntegrateData(anchorset = anchors_sample)

# Center and scale the integrated values for improved performance with dim. reduction and clustering
sce_integrated_sample <- ScaleData(sce_integrated_sample, do.scale = TRUE, do.center = TRUE, assay = "integrated")

# Now we can start looking at some dim. reduction and clustering of the integrated dataset
sce_integrated_sample <- RunPCA(sce_integrated_sample, npcs = 50, assay = "integrated", seed.use = 1337)
sce_integrated_sample <- RunTSNE(sce_integrated_sample, dims = 1:30, assay = "integrated", dim.embed = 3, seed.use = 300, perplexity = 20)
sce_integrated_sample <- RunUMAP(sce_integrated_sample, reduction = "pca", dims = 1:30, assay = "integrated", seed.use = 1337)
sce_integrated_sample <- FindNeighbors(sce_integrated_sample, reduction = "pca", dims = 1:30, assay = "integrated")
sce_integrated_sample <- FindClusters(sce_integrated_sample, resolution = 0.6, random.seed = 1337)

# Visualize clusters and set some plotting defaults
DefaultAssay(object = sce_integrated_sample) <- "RNA"
sce_integrated_sample <- ScaleData(object = sce_integrated_sample)

global_theme <- theme(plot.title = element_blank(),
  axis.line = element_line(size = 0.3))

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
                       "#FDB462", "#BC80BD")

plot_clusters <- DimPlot(object = sce_integrated_sample, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE) +
  global_theme +
  scale_colour_manual(values = colors)

plot_samples <- DimPlot(object = sce_integrated_sample, reduction = "umap", group.by = "Name", shuffle = TRUE) +
  global_theme +
  scale_colour_manual(values = colors_secondary)

plot_time <- DimPlot(object = sce_integrated_sample, reduction = "umap", group.by = "Time", shuffle = TRUE) +
  global_theme +
  scale_colour_manual(values = colors_secondary)

# Visualize genes of interest
dot_plot_genes <- read_yaml(file = "../reference/cell_type_markers_mus_musculus_kidney.yaml")

plot_cell_dot <- DotPlot(object = sce_integrated_sample, features = dot_plot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

# Find marker genes per cluster
markers_cluster <- FindAllMarkers(object = sce_integrated_sample, only.pos = TRUE)

top_markers_cluster <- markers_cluster %>%
  group_by(cluster) %>%
  slice_head(n = 3)

# Plot broad cell types derived from key marker genes
cluster_labels <- c("0" = "Prox Tub",
                    "1" = "Prox Tub",
                    "2" = "Prox Tub",
                    "3" = "Immune",
                    "4" = "Prox Tub",
                    "5" = "Endo",
                    "6" = "Prox Tub",
                    "7" = "Fibro",
                    "8" = "Uni Fibro",
                    "9" = "Endo",
                    "10" = "Cycle",
                    "11" = "Podo",
                    "12" = "Immune",
                    "13" = "Collect Duct",
                    "14" = "Mac Densa",
                    "15" = "Cycle",
                    "16" = "Loop Henle",
                    "17" = "Prox Tub",
                    "18" = "Endo",
                    "19" = "Fibro",
                    "20" = "Endo",
                    "21" = "Trans Epi",
                    "22" = "Mural",
                    "23" = "Endo")

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
  "Mural" = "#39B606",
  "Uni Fibro (Pi16)" = "#FD7083",
  "Fibro Perturb" = "#EA8331",
  "Uni Fibro (Col15a1)" = "#B186FF",
  "Cycle" = "#00BAE0",
  "Mac Densa" = "#39B607",
  "Peri" = "#C09B06",
  "SM Vasc" = "#91AA00",
  "Unk" = "black"
)

sce_integrated_sample$cell_type <- factor(x = as.character(cluster_labels[as.character(sce_integrated_sample$seurat_clusters)]),
                                          levels = c("Uni Fibro",
                                                     "Fibro",
                                                     "Mural",
                                                     "Endo",
                                                     "Podo",
                                                     "Prox Tub",
                                                     "Loop Henle",
                                                     "Mac Densa",
                                                     "Collect Duct",
                                                     "Trans Epi",
                                                     "Immune",
                                                     "Cycle"))

plot_cell_type <- DimPlot(object = sce_integrated_sample, reduction = "umap", label = TRUE, group.by = "cell_type", label.size = 3) +
  global_theme +
  scale_colour_manual(values = cluster_colors)

downsample_clusters <- min(table(sce_integrated_sample$seurat_clusters))
plot_heatmap <- DoHeatmap(object = subset(sce_integrated_sample, downsample = downsample_clusters),
                          features = top_markers_cluster$gene,
                          raster = TRUE,
                          angle = 45,
                          draw.lines = FALSE,
                          label = TRUE,
                          group.colors = colors) +
  viridis::scale_fill_viridis(option = "G")

Idents(sce_integrated_sample) <- sce_integrated_sample$cell_type

markers_batch <- FindAllMarkers(object = sce_integrated_sample, only.pos = TRUE)

top_markers_batch <- markers_batch %>%
  group_by(cluster) %>%
  slice_head(n = 5)

downsample_batch <- min(table(sce_integrated_sample$cell_type))
plot_heatmap_batch <- DoHeatmap(object = subset(sce_integrated_sample, downsample = downsample_clusters),
                          features = top_markers_batch$gene,
                          raster = TRUE,
                          angle = 45,
                          draw.lines = FALSE,
                          label = TRUE,
                          group.colors =  unlist(cluster_colors[levels(sce_integrated_sample$cell_type)])) +
  viridis::scale_fill_viridis(option = "G")

plot_cell_dot_batch <- DotPlot(object = sce_integrated_sample, features = dot_plot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

# Save the unfiltered, integrated Seurat object
saveRDS(object = sce_integrated_sample, file = "../output/scrna/dicer1_scrna_unfilter_control.rds")

# Write the table(s) of markers out
write.csv(x = markers_cluster, file = "../output/tables/dicer1_scrna_unfilter_control_cluster_markers_table.csv", row.names = FALSE)
write.csv(x = markers_batch, file = "../output/tables/dicer1_scrna_unfilter_control_batch_markers_table.csv", row.names = FALSE)

# Save plots for the unfiltered data
out_dir <- "../output/plots/"

ggsave(filename = "dicer1_scrna_unfilter_control_clusters.png", plot = plot_clusters, path = out_dir, width = 5.5, height = 4)
ggsave(filename = "dicer1_scrna_unfilter_control_clusters.svg", plot = plot_clusters, path = out_dir, width = 5.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_control_samples.png", plot = plot_samples, path = out_dir, width = 6.5, height = 4)
ggsave(filename = "dicer1_scrna_unfilter_control_samples.svg", plot = plot_samples, path = out_dir, width = 6.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_control_cell_type.png", plot = plot_cell_type, path = out_dir, width = 6, height = 4)
ggsave(filename = "dicer1_scrna_unfilter_control_cell_type.svg", plot = plot_cell_type, path = out_dir, width = 6, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_control_time.png", plot = plot_time, path = out_dir, width = 5, height = 4)
ggsave(filename = "dicer1_scrna_unfilter_control_time.svg", plot = plot_time, path = out_dir, width = 5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_control_dot_plot.png", plot = plot_cell_dot, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_unfilter_control_dot_plot.svg", plot = plot_cell_dot, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_control_dot_plot_batch.png", plot = plot_cell_dot_batch, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_unfilter_control_dot_plot_batch.svg", plot = plot_cell_dot_batch, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_control_heatmap_clusters.png", plot = plot_heatmap, path = out_dir, width = 15, height = 9)
ggsave(filename = "dicer1_scrna_unfilter_control_heatmap_clusters.svg", plot = plot_heatmap, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_control_heatmap_batch.png", plot = plot_heatmap_batch, path = out_dir, width = 15, height = 9)
ggsave(filename = "dicer1_scrna_unfilter_control_heatmap_batch.svg", plot = plot_heatmap_batch, path = out_dir, width = 15, height = 9, device = "svg")

### FILTERED
# Filter for the mesenchymal Hic1 lineage only, which is enriched in tdTomato+ cells
mes_cells <- sce_integrated_sample$cell_type %in% c("Uni Fibro", "Fibro", "Mural")
sce_integrated_sample_filter <- sce_integrated_sample[, mes_cells]

seurat_list_sample_filter <- SplitObject(sce_integrated_sample_filter, split.by = "Sample")
seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) NormalizeData(seurat_obj, verbose = FALSE))
seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000))

features_sample_filter <- SelectIntegrationFeatures(object.list = seurat_list_sample_filter, nfeatures = 2000)

seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) ScaleData(seurat_obj, features = features_sample_filter))
seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) RunPCA(seurat_obj, features = features_sample_filter))

anchors_sample_filter <- FindIntegrationAnchors(object.list = seurat_list_sample_filter, anchor.features = features_sample_filter, reduction = "rpca")

sce_integrated_sample_filter <- IntegrateData(anchorset = anchors_sample_filter, k.weight = 80)
sce_integrated_sample_filter <- ScaleData(sce_integrated_sample_filter, do.scale = TRUE, do.center = TRUE, assay = "integrated")
sce_integrated_sample_filter <- RunPCA(sce_integrated_sample_filter, npcs = 13, assay = "integrated", seed.use = 1337)
sce_integrated_sample_filter <- RunTSNE(sce_integrated_sample_filter, dims = 1:13, assay = "integrated", dim.embed = 3, seed.use = 300, perplexity = 20)
sce_integrated_sample_filter <- RunUMAP(sce_integrated_sample_filter, reduction = "pca", dims = 1:13, assay = "integrated", seed.use = 1337)
sce_integrated_sample_filter <- FindNeighbors(sce_integrated_sample_filter, reduction = "pca", dims = 1:13, assay = "integrated")
sce_integrated_sample_filter <- FindClusters(sce_integrated_sample_filter, resolution = 0.2, random.seed = 1337)

# Visualize clusters
DefaultAssay(object = sce_integrated_sample_filter) <- "RNA"
sce_integrated_sample_filter <- ScaleData(object = sce_integrated_sample_filter)

plot_clusters_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE) +
  global_theme +
  scale_colour_manual(values = colors)

plot_samples_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", group.by = "Name", shuffle = TRUE) +
  global_theme +
  scale_colour_manual(values = colors_secondary)

plot_time_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", group.by = "Time", shuffle = TRUE) +
  global_theme +
  scale_colour_manual(values = colors_secondary)

# Find marker genes
markers_filter <- FindAllMarkers(object = sce_integrated_sample_filter, only.pos = TRUE)

top_markers_filter <- markers_filter %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), desc(pct.1)) %>%
  slice_head(n = 15)

# Visualize by cell type
cluster_labels_filter <- c("0" = "Fibro",
                           "1" = "Uni Fibro (Pi16)",
                           "2" = "Cycle",
                           "3" = "Mural",
                           "4" = "Uni Fibro (Col15a1)",
                           "5" = "Fibro")

sce_integrated_sample_filter$cell_type <- factor(x = as.character(cluster_labels_filter[as.character(sce_integrated_sample_filter$seurat_clusters)]),
                                          levels = c("Uni Fibro (Pi16)",
                                                     "Uni Fibro (Col15a1)",
                                                     "Fibro",
                                                     "Mural",
                                                     "Cycle"))

plot_cell_type_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", label = TRUE, group.by = "cell_type", label.size = 3) +
  global_theme +
  scale_colour_manual(values = cluster_colors)

dot_plot_genes_mes <- read_yaml(file = "../reference/cell_type_markers_mus_musculus_mesenchymal.yaml")
plot_cell_dot_filter <- DotPlot(object = sce_integrated_sample_filter, features = dot_plot_genes_mes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

plot_features_filter <- FeaturePlot(sce_integrated_sample_filter, features = c("tdtomato", "Hic1", "Pdgfra", "Myh11", "Dpt", "Pi16", "Col15a1", "Mfap4"),
                                    ncol = 2,
                                    order = TRUE,
                                    keep.scale = "feature")

# Find marker genes for the broad labels
Idents(sce_integrated_sample_filter) <- sce_integrated_sample_filter$cell_type
markers_batch_filter <- FindAllMarkers(object = sce_integrated_sample_filter, only.pos = TRUE)

top_markers_batch_filter <- markers_batch_filter %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), desc(pct.1)) %>%
  slice_head(n = 50)

downsample_batch_filter <- min(table(sce_integrated_sample_filter$cell_type))
plot_filter_heatmap <- DoHeatmap(object = subset(sce_integrated_sample_filter, downsample = downsample_batch_filter),
                                 features = top_markers_batch_filter$gene,
                                 raster = TRUE,
                                 angle = 45,
                                 draw.lines = FALSE,
                                 label = TRUE,
                                 group.colors =  unlist(cluster_colors[levels(sce_integrated_sample_filter$cell_type)])) +
  viridis::scale_fill_viridis(option = "G")

plot_cell_dot_batch_filter <- DotPlot(object = sce_integrated_sample_filter, features = dot_plot_genes_mes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

# Save the final, filtered and integrated Seurat object
saveRDS(object = sce_integrated_sample_filter, file = "../output/scrna/dicer1_scrna_filter_control.rds")

# Write the table(s) of markers out
write.csv(x = markers_filter, file = "../output/tables/dicer1_scrna_filter_control_cluster_markers_table.csv", row.names = FALSE)
write.csv(x = markers_batch_filter, file = "../output/tables/dicer1_scrna_filter_control_batch_markers_table.csv", row.names = FALSE)

# Save plots for the filtered data
out_dir <- "../output/plots/"

ggsave(filename = "dicer1_scrna_filter_control_clusters.png", plot = plot_clusters_filter, path = out_dir, width = 5.5, height = 4)
ggsave(filename = "dicer1_scrna_filter_control_clusters.svg", plot = plot_clusters_filter, path = out_dir, width = 5.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_samples.png", plot = plot_samples_filter, path = out_dir, width = 6.5, height = 4)
ggsave(filename = "dicer1_scrna_filter_control_samples.svg", plot = plot_samples_filter, path = out_dir, width = 6.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_cell_type.png", plot = plot_cell_type_filter, path = out_dir, width = 6, height = 4)
ggsave(filename = "dicer1_scrna_filter_control_cell_type.svg", plot = plot_cell_type_filter, path = out_dir, width = 6, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_time.png", plot = plot_time_filter, path = out_dir, width = 5, height = 4)
ggsave(filename = "dicer1_scrna_filter_control_time.svg", plot = plot_time_filter, path = out_dir, width = 5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_dot_plot.png", plot = plot_cell_dot_filter, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_filter_control_dot_plot.svg", plot = plot_cell_dot_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_dot_plot_batch.png", plot = plot_cell_dot_batch_filter, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_filter_control_dot_plot_batch.svg", plot = plot_cell_dot_batch_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_heatmap_batch.png", plot = plot_filter_heatmap, path = out_dir, width = 15, height = 9)
ggsave(filename = "dicer1_scrna_filter_control_heatmap_batch.svg", plot = plot_filter_heatmap, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_features.png", plot = plot_features_filter, path = out_dir, width = 12, height = 20)
ggsave(filename = "dicer1_scrna_filter_control_features.svg", plot = plot_features_filter, path = out_dir, width = 12, height = 20, device = "svg")
