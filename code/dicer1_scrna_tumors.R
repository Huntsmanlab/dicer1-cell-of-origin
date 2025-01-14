### MUTANTS
# Script to normalize data via Seurat's integration, then cluster to identify cell types and find markers

# Required packages
library(dplyr)
library(Seurat)
library(scater)
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

### Tumor samples, no filtering by expressed genes, Dicer1 CH
keep_samples <- sce_qc$genotype == "Dicer1 CH"
sce_qc_filter_geno <- sce_qc[, keep_samples]

# Convert to a Seurat object
sce_seurat_geno <- CreateSeuratObject(counts = counts(sce_qc_filter_geno), names.field = 1, names.delim = "/")

# Add our metadata to the Seurat object
seurat_metadata_geno = list("Sample" = sce_qc_filter_geno$Sample,
                            "Barcode" = sce_qc_filter_geno$Barcode,
                            "Genotype" = sce_qc_filter_geno$genotype,
                            "Time" = sce_qc_filter_geno$time,
                            "Tomato" = sce_qc_filter_geno$subsets_tomato_sum,
                            "mito_percent" = sce_qc_filter_geno$subsets_mito_percent,
                            "Name" = sce_qc_filter_geno$name,
                            "Delta" = sce_qc_filter_geno$delta_days,
                            "ribo_percent" = sce_qc_filter_geno$subsets_ribo_percent)
sce_seurat_geno <- AddMetaData(object = sce_seurat_geno, metadata = seurat_metadata_geno)

# Split into one Seurat object per sample and normalize
seurat_list_geno <- SplitObject(sce_seurat_geno, split.by = "Sample")

# Normalize each sample with Seurat's method
seurat_list_geno <- lapply(seurat_list_geno, function(seurat_obj) NormalizeData(seurat_obj, verbose = FALSE))

# Identify variable features (genes) per sample; 2,000 is roughly in line with the top 20% of HVGs
seurat_list_geno <- lapply(seurat_list_geno, function(seurat_obj) FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000))

# Select variable features by their rank across samples, prioritizing conserved features
features_geno <- SelectIntegrationFeatures(object.list = seurat_list_geno, nfeatures = 2000)

# Do per-sample PCA so we can use rPCA to integrate
seurat_list_geno <- lapply(seurat_list_geno, function(seurat_obj) ScaleData(seurat_obj, features = features_geno))
seurat_list_geno <- lapply(seurat_list_geno, function(seurat_obj) RunPCA(seurat_obj, features = features_geno))

# Find anchors across pairwise combinations of our datasets: pairs of similar/conserved cells
anchors_geno <- FindIntegrationAnchors(object.list = seurat_list_geno, anchor.features = features_geno, reduction = "rpca")

# Create the integrated dataset, using the anchors we found; default features used are those we used to find anchors
sce_integrated_geno <- IntegrateData(anchorset = anchors_geno)

# Center and scale the integrated values for improved performance with dim. reduction and clustering
sce_integrated_geno <- ScaleData(sce_integrated_geno, do.scale = TRUE, do.center = TRUE, assay = "integrated")

# Now we can start looking at some dim. reduction and clustering of the integrated dataset
sce_integrated_geno <- RunPCA(sce_integrated_geno, npcs = 50, assay = "integrated", seed.use = 1337)
sce_integrated_geno <- RunTSNE(sce_integrated_geno, dims = 1:30, assay = "integrated", dim.embed = 3, seed.use = 300, perplexity = 20)
sce_integrated_geno <- RunUMAP(sce_integrated_geno, reduction = "pca", dims = 1:30, assay = "integrated", seed.use = 1337)
sce_integrated_geno <- FindNeighbors(sce_integrated_geno, reduction = "pca", dims = 1:30, assay = "integrated")
sce_integrated_geno <- FindClusters(sce_integrated_geno, resolution = 0.6, random.seed = 1337)

# Visualize clusters and set some plotting defaults
DefaultAssay(object = sce_integrated_geno) <- "RNA"
sce_integrated_geno <- ScaleData(object = sce_integrated_geno)

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

plot_clusters <- DimPlot(object = sce_integrated_geno, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE) +
  global_theme +
  scale_colour_manual(values = colors)

plot_samples <- DimPlot(object = sce_integrated_geno, reduction = "umap", group.by = "Name", shuffle = TRUE) +
  global_theme +
  scale_colour_manual(values = colors_secondary)

plot_time <- DimPlot(object = sce_integrated_geno, reduction = "umap", group.by = "Time", shuffle = TRUE) +
  global_theme +
  scale_colour_manual(values = colors_secondary)

# Visualize genes of interest
dot_plot_genes <- read_yaml(file = "../reference/cell_type_markers_mus_musculus_kidney.yaml")
dot_plot_genes_mes <- read_yaml(file = "../reference/cell_type_markers_mus_musculus_mesenchymal.yaml")

plot_cell_dot_kidney <- DotPlot(object = sce_integrated_geno, features = dot_plot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

plot_cell_dot_mes <- DotPlot(object = sce_integrated_geno, features = dot_plot_genes_mes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

# Find marker genes
markers_cluster <- FindAllMarkers(object = sce_integrated_geno, only.pos = TRUE)

top_markers_cluster <- markers_cluster %>%
  group_by(cluster) %>%
  slice_head(n = 3)

# Plot broad cell types derived from key marker genes
cluster_labels <- c("0" = "Trans Epi",
                    "1" = "Immune",
                    "2" = "Immune",
                    "3" = "Fibro",
                    "4" = "Trans Epi",
                    "5" = "Fibro",
                    "6" = "Immune",
                    "7" = "Trans Epi",
                    "8" = "Immune",
                    "9" = "TR Diff Myo",
                    "10" = "Trans Epi",
                    "11" = "Trans Epi",
                    "12" = "Immune",
                    "13" = "Immune",
                    "14" = "Immune",
                    "15" = "Immune",
                    "16" = "Podo",
                    "17" = "Uni Fibro",
                    "18" = "Immune",
                    "19" = "Fibro",
                    "20" = "Diff Myo",
                    "21" = "Cycle",
                    "22" = "Mural",
                    "23" = "Immune",
                    "24" = "Fibro")

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

sce_integrated_geno$cell_type <- factor(x = as.character(cluster_labels[as.character(sce_integrated_geno$seurat_clusters)]),
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
                                                   "Cycle",
                                                   "TR Diff Myo",
                                                   "Diff Myo"))

plot_cell_type <- DimPlot(object = sce_integrated_geno, reduction = "umap", label = TRUE, group.by = "cell_type", label.size = 3) +
  global_theme +
  scale_colour_manual(values = cluster_colors)

downsample_clusters <- min(table(sce_integrated_geno$seurat_clusters))
plot_heatmap <- DoHeatmap(object = subset(sce_integrated_geno, downsample = downsample_clusters),
                          features = top_markers_cluster$gene,
                          raster = TRUE,
                          angle = 45,
                          draw.lines = FALSE,
                          label = TRUE,
                          group.colors = colors) +
  viridis::scale_fill_viridis(option = "G")

Idents(sce_integrated_geno) <- sce_integrated_geno$cell_type

markers_batch <- FindAllMarkers(object = sce_integrated_geno, only.pos = TRUE)

top_markers_batch <- markers_batch %>%
  group_by(cluster) %>%
  slice_head(n = 5)

downsample_batch <- min(table(sce_integrated_geno$cell_type))
plot_heatmap_batch <- DoHeatmap(object = subset(sce_integrated_geno, downsample = downsample_clusters),
                                features = top_markers_batch$gene,
                                raster = TRUE,
                                angle = 45,
                                draw.lines = FALSE,
                                label = TRUE,
                                group.colors =  unlist(cluster_colors[levels(sce_integrated_geno$cell_type)])) +
  viridis::scale_fill_viridis(option = "G")

plot_cell_dot_batch_kidney <- DotPlot(object = sce_integrated_geno, features = dot_plot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

plot_cell_dot_batch_mes <- DotPlot(object = sce_integrated_geno, features = dot_plot_genes_mes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

# Save the unfiltered, integrated Seurat object
saveRDS(object = sce_integrated_geno, file = "../output/scrna/dicer1_scrna_unfilter_tumor.rds")

# Write the table(s) of markers out
write.csv(x = markers_cluster, file = "../output/tables/dicer1_scrna_unfilter_tumor_cluster_markers_table.csv", row.names = FALSE)
write.csv(x = markers_batch, file = "../output/tables/dicer1_scrna_unfilter_tumor_batch_markers_table.csv", row.names = FALSE)

# Save plots for the unfiltered data
out_dir <- "../output/plots/"

ggsave(filename = "dicer1_scrna_unfilter_tumor_clusters.png", plot = plot_clusters, path = out_dir, width = 5.5, height = 4)
ggsave(filename = "dicer1_scrna_unfilter_tumor_clusters.svg", plot = plot_clusters, path = out_dir, width = 5.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_tumor_samples.png", plot = plot_samples, path = out_dir, width = 5.5, height = 4)
ggsave(filename = "dicer1_scrna_unfilter_tumor_samples.svg", plot = plot_samples, path = out_dir, width = 5.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_tumor_cell_type.png", plot = plot_cell_type, path = out_dir, width = 6, height = 4)
ggsave(filename = "dicer1_scrna_unfilter_tumor_cell_type.svg", plot = plot_cell_type, path = out_dir, width = 6, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_tumor_time.png", plot = plot_time, path = out_dir, width = 5, height = 4)
ggsave(filename = "dicer1_scrna_unfilter_tumor_time.svg", plot = plot_time, path = out_dir, width = 5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_tumor_dot_plot_kidney.png", plot = plot_cell_dot_kidney, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_unfilter_tumor_dot_plot_kidney.svg", plot = plot_cell_dot_kidney, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_tumor_dot_plot_mes.png", plot = plot_cell_dot_mes, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_unfilter_tumor_dot_plot_mes.svg", plot = plot_cell_dot_mes, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_tumor_dot_plot_batch_kidney.png", plot = plot_cell_dot_batch_kidney, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_unfilter_tumor_dot_plot_batch_kidney.svg", plot = plot_cell_dot_batch_kidney, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_tumor_dot_plot_batch_mes.png", plot = plot_cell_dot_batch_mes, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_unfilter_tumor_dot_plot_batch_mes.svg", plot = plot_cell_dot_batch_mes, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_tumor_heatmap_clusters.png", plot = plot_heatmap, path = out_dir, width = 15, height = 9)
ggsave(filename = "dicer1_scrna_unfilter_tumor_heatmap_clusters.svg", plot = plot_heatmap, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_unfilter_tumor_heatmap_batch.png", plot = plot_heatmap_batch, path = out_dir, width = 15, height = 9)
ggsave(filename = "dicer1_scrna_unfilter_tumor_heatmap_batch.svg", plot = plot_heatmap_batch, path = out_dir, width = 15, height = 9, device = "svg")

### FILTERED
# Filter for the Mesenchymal/no_1mo_3mo Hic1 lineage only, and ensure we are only looking at tdTomato+ cells
mes_cells <- sce_integrated_geno$cell_type %in% c("Uni Fibro",
                                                  "Fibro",
                                                  "Cycle",
                                                  "TR Diff Myo",
                                                  "Diff Myo",
                                                  "Mural")

# Keep tdTomato+ cells, and discard Ptprc+ (Immune, n = 113), Epcam+ (Epithelial, n = 102) cells, and Krt19+ (Transitional Epithelieum, n = 376);
# Otherwise we will observe small amounts of residual immmune/epi contamination, which are likely doublets with fibro that can be seen near these groups in the unfiltered plot
keep_cells <- sce_integrated_geno$Tomato > 0

counts <- GetAssayData(object = sce_integrated_geno, assay = "RNA", slot = "counts")
ptprc_cells <- counts["Ptprc", ] > 0
epcam_cells <- counts["Epcam", ] > 0
trans_epi_cells <- counts["Krt19", ] > 0
discard_cells <- ptprc_cells | epcam_cells | trans_epi_cells

sce_integrated_sample_filter <- sce_integrated_geno[, mes_cells & keep_cells & !discard_cells]

# We will only keep samples that have 50 or more tdTomato+ cells remaining, as integration is poor with very few cells
sce_integrated_sample_filter <- sce_integrated_sample_filter[, sce_integrated_sample_filter$Sample %in%
                                                               names(table(sce_integrated_sample_filter$Sample))[table(sce_integrated_sample_filter$Sample) > 50]]

seurat_list_sample_filter <- SplitObject(sce_integrated_sample_filter, split.by = "Sample")
seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) NormalizeData(seurat_obj, verbose = FALSE))
seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000))

features_sample_filter <- SelectIntegrationFeatures(object.list = seurat_list_sample_filter, nfeatures = 2000)

seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) ScaleData(seurat_obj, features = features_sample_filter))
seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) RunPCA(seurat_obj, features = features_sample_filter))

anchors_sample_filter <- FindIntegrationAnchors(object.list = seurat_list_sample_filter, anchor.features = features_sample_filter, reduction = "rpca")

sce_integrated_sample_filter <- IntegrateData(anchorset = anchors_sample_filter, k.weight = 25)
sce_integrated_sample_filter <- ScaleData(sce_integrated_sample_filter, do.scale = TRUE, do.center = TRUE, assay = "integrated")
sce_integrated_sample_filter <- RunPCA(sce_integrated_sample_filter, npcs = 50, assay = "integrated", seed.use = 1337)

# Fewer PCs are required to discriminate clusters, as expected given that we have far fewer cell types after filtering
plot_elbow_filter <- ElbowPlot(sce_integrated_sample_filter, ndims = 50)

sce_integrated_sample_filter <- RunTSNE(sce_integrated_sample_filter, dims = 1:20, assay = "integrated", dim.embed = 3, seed.use = 300, perplexity = 20)
sce_integrated_sample_filter <- RunUMAP(sce_integrated_sample_filter, reduction = "pca", dims = 1:20, assay = "integrated", seed.use = 1337)
sce_integrated_sample_filter <- FindNeighbors(sce_integrated_sample_filter, reduction = "pca", dims = 1:20, assay = "integrated")
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

plot_cell_dot_kidney_filter <- DotPlot(object = sce_integrated_sample_filter, features = dot_plot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

plot_cell_dot_mes_filter <- DotPlot(object = sce_integrated_sample_filter, features = dot_plot_genes_mes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

# Find marker genes
markers_cluster_filter <- FindAllMarkers(object = sce_integrated_sample_filter, only.pos = TRUE)

top_markers_cluster_filter <- markers_cluster_filter %>%
  group_by(cluster) %>%
  slice_head(n = 30)

downsample_cluster_filter <- min(table(sce_integrated_sample_filter$seurat_clusters))
plot_heatmap_cluster_filter <- DoHeatmap(object = subset(sce_integrated_sample_filter, downsample = downsample_cluster_filter),
                                         features = top_markers_cluster_filter$gene,
                                         raster = TRUE,
                                         angle = 45,
                                         draw.lines = FALSE,
                                         label = TRUE) +
  viridis::scale_fill_viridis(option = "G")

# Plot broad cell types derived from key marker genes
cluster_labels_filter <- c("0" = "Progen",
                           "1" = "Ground",
                           "2" = "TR Diff Myo",
                           "3" = "Ground",
                           "4" = "Prolif",
                           "5" = "Diff Myo",
                           "6" = "Progen",
                           "7" = "Peri",
                           "8" = "SM Vasc",
                           "9" = "Fibro Perturb")

sce_integrated_sample_filter$cell_type <- factor(x = as.character(cluster_labels_filter[as.character(sce_integrated_sample_filter$seurat_clusters)]),
                                                 levels = c("Progen",
                                                            "Ground",
                                                            "Prolif",
                                                            "TR Diff Myo",
                                                            "Diff Myo",
                                                            "Fibro Perturb",
                                                            "Peri",
                                                            "SM Vasc"))

plot_cell_type_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", label = TRUE, group.by = "cell_type", label.size = 3) +
  global_theme +
  scale_colour_manual(values = cluster_colors)

plot_cell_type_time_split_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", label = FALSE, group.by = "cell_type", split.by = "Time", ncol = 2) +
  global_theme +
  scale_colour_manual(values = cluster_colors)

Idents(sce_integrated_sample_filter) <- sce_integrated_sample_filter$cell_type

markers_batch_filter <- FindAllMarkers(object = sce_integrated_sample_filter, only.pos = TRUE)
all_de_batch_filter <- FindAllMarkers(object = sce_integrated_sample_filter, only.pos = FALSE)

top_markers_batch_filter <- markers_batch_filter %>%
  group_by(cluster) %>%
  slice_head(n = 35)

downsample_batch_filter <- min(table(sce_integrated_sample_filter$cell_type))
plot_heatmap_batch_filter <- DoHeatmap(object = subset(sce_integrated_sample_filter, downsample = downsample_batch_filter),
                                       features = top_markers_batch_filter$gene,
                                       raster = TRUE,
                                       angle = 45,
                                       draw.lines = FALSE,
                                       label = TRUE,
                                       group.colors =  unlist(cluster_colors[levels(sce_integrated_sample_filter$cell_type)])) +
  viridis::scale_fill_viridis(option = "G")

plot_cell_dot_batch_kidney_filter <- DotPlot(object = sce_integrated_sample_filter, features = dot_plot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

plot_cell_dot_batch_mes_filter <- DotPlot(object = sce_integrated_sample_filter, features = dot_plot_genes_mes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

plot_features_filter <- FeaturePlot(object = sce_integrated_sample_filter,
                                    features = c("tdtomato", "Hic1", "Pdgfra", "Dpt", "Pi16", "Col15a1"),
                                    ncol = 3,
                                    order = TRUE)

plot_features_xenium_filter <- FeaturePlot(object = sce_integrated_sample_filter,
                                           features = c("Pi16", "Lum", "Col15a1", "Mfap4", "Tagln", "Myh11", "Cenpf", "Mki67", "Car3", "Myoz2", "Itga8", "Rgs5"),
                                           ncol = 3,
                                           order = TRUE)

# Save the filtered, integrated Seurat object
saveRDS(object = sce_integrated_sample_filter, file = "../output/scrna/dicer1_scrna_filter_tumor.rds")

# Write the table(s) of markers out
write.csv(x = markers_cluster_filter, file = "../output/tables/dicer1_scrna_filter_tumor_cluster_markers_table.csv", row.names = FALSE)
write.csv(x = markers_batch_filter, file = "../output/tables/dicer1_scrna_filter_tumor_batch_markers_table.csv", row.names = FALSE)
write.csv(x = all_de_batch_filter, file = "../output/tables/dicer1_scrna_filter_tumor_batch_all_de_table.csv", row.names = FALSE)

# Save plots for the filtered data
out_dir <- "../output/plots/"

ggsave(filename = "dicer1_scrna_filter_tumor_clusters.png", plot = plot_clusters_filter, path = out_dir, width = 5.5, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_clusters.svg", plot = plot_clusters_filter, path = out_dir, width = 5.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_samples.png", plot = plot_samples_filter, path = out_dir, width = 6.5, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_samples.svg", plot = plot_samples_filter, path = out_dir, width = 6.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_cell_type.png", plot = plot_cell_type_filter, path = out_dir, width = 6, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_cell_type.svg", plot = plot_cell_type_filter, path = out_dir, width = 6, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_time.png", plot = plot_time_filter, path = out_dir, width = 5, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_time.svg", plot = plot_time_filter, path = out_dir, width = 5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_kidney.png", plot = plot_cell_dot_kidney_filter, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_kidney.svg", plot = plot_cell_dot_kidney_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_mes.png", plot = plot_cell_dot_mes_filter, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_mes.svg", plot = plot_cell_dot_mes_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_batch_kidney.png", plot = plot_cell_dot_batch_kidney_filter, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_batch_kidney.svg", plot = plot_cell_dot_batch_kidney_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_batch_mes.png", plot = plot_cell_dot_batch_mes_filter, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_batch_mes.svg", plot = plot_cell_dot_batch_mes_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_heatmap_clusters.png", plot = plot_heatmap_cluster_filter, path = out_dir, width = 15, height = 9)
ggsave(filename = "dicer1_scrna_filter_tumor_heatmap_clusters.svg", plot = plot_heatmap_cluster_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_heatmap_batch.png", plot = plot_heatmap_batch_filter, path = out_dir, width = 15, height = 9)
ggsave(filename = "dicer1_scrna_filter_tumor_heatmap_batch.svg", plot = plot_heatmap_batch_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_features.png", plot = plot_features_filter, path = out_dir, width = 18, height = 10)
ggsave(filename = "dicer1_scrna_filter_tumor_features.svg", plot = plot_features_filter, path = out_dir, width = 18, height = 10, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_xenium_features.png", plot = plot_features_xenium_filter, path = out_dir, width = 18, height = 20)
ggsave(filename = "dicer1_scrna_filter_tumor_xenium_features.svg", plot = plot_features_xenium_filter, path = out_dir, width = 18, height = 20, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_cell_type_time_split.png", plot = plot_cell_type_time_split_filter, path = out_dir, width = 10, height = 8)
ggsave(filename = "dicer1_scrna_filter_tumor_cell_type_time_split.svg", plot = plot_cell_type_time_split_filter, path = out_dir, width = 10, height = 8, device = "svg")
