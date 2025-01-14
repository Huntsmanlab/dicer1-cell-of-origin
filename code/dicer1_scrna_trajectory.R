### Script to perform trajectory analysis on the tumor samples

# Required packages
library(Seurat)
library(monocle3)
library(dplyr)
library(biomaRt)
library(clusterProfiler)
library(ggplot2)

# Seed for reproducibility
set.seed(seed = 1337)

# Read in the integrated data from Seurat
seurat_integrated <- readRDS(file = "../output/scrna/dicer1_scrna_filter_tumor.rds")

# Cluster colours remain the same
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

# Convert to the CellDataSet class used by Monocle
gene_annotation <- as.data.frame(rownames(seurat_integrated@assays$RNA@data))
colnames(gene_annotation) = "gene_short_name"
rownames(gene_annotation) = gene_annotation$gene_short_name

cell_metadata <- seurat_integrated@meta.data

exp_matrix <- GetAssayData(object = seurat_integrated, assay = "RNA", slot = "data")

cds_integrated <- new_cell_data_set(exp_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# Dummy paritions to prevent monocle errors; these cannot be used, hence we set use_partition to FALSE when calling learn_graph
recreate_partition <- c(rep(1, length(cds_integrated@colData@rownames)))
names(recreate_partition) <- cds_integrated@colData@rownames
recreate_partition <- as.factor(recreate_partition)
cds_integrated@clusters@listData[["UMAP"]][["partitions"]] <- recreate_partition

clusters <- seurat_integrated$cell_type
names(clusters) <- seurat_integrated@assays[["RNA"]]@data@Dimnames[[2]]
cds_integrated@clusters@listData[["UMAP"]][["clusters"]] <- clusters
cds_integrated@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
cds_integrated@int_colData@listData$reducedDims@listData[["UMAP"]] <-seurat_integrated@reductions[["umap"]]@cell.embeddings

# Learn graph and infer pseudotime
cds_integrated <- learn_graph(cds = cds_integrated,
                              close_loop = TRUE,
                              use_partition = FALSE,
                              learn_graph_control = list("minimal_branch_len" = 15))

plot_cell_type <- plot_cells(cds = cds_integrated, reduction_method = "UMAP", color_cells_by = "cell_type", label_cell_groups = FALSE)

# Infer pseudo-time positions; choose start node as node nearest the Dpt+ region of the Progen cluster
cds_integrated <- order_cells(cds = cds_integrated)
plot_pseudotime <- plot_cells(cds = cds_integrated, color_cells_by = "pseudotime") +
  theme( legend.position = "none")
plot_pseudotime_no_graph <- plot_cells(cds = cds_integrated, color_cells_by = "pseudotime", show_trajectory_graph = FALSE)
plot_cell_type_graph <- plot_cells(cds = cds_integrated, color_cells_by = "cell_type", label_cell_groups = FALSE) +
  scale_colour_manual(name = "Cell type", values = cluster_colors) +
  theme( legend.position = "none")

###  Subset to our branches of interest so we can show key genes over each trajectory
# Trajectory to diff myo; choose terminal nodes in TR Diff Myo and Diff Myo as endpoints
cds_to_muscle <- choose_graph_segments(cds = cds_integrated, clear_cds = FALSE)
cds_to_muscle <- order_cells(cds = cds_to_muscle)

# Plot key genes over the trajectory
muscle_genes <- c("Dpt", "Mfap4", "Mki67", "Pax7", "Myog")
cds_to_muscle_genes <- cds_to_muscle[rowData(cds_to_muscle)$gene_short_name %in% muscle_genes, ]

plot_genes_rhabdo <- plot_genes_in_pseudotime(cds_subset = cds_to_muscle_genes, panel_order = muscle_genes)
plot_genes_clusters_rhabdo <- plot_genes_in_pseudotime(cds_subset = cds_to_muscle_genes, color_cells_by = "cell_type", panel_order = muscle_genes) +
  scale_colour_manual(values = cluster_colors, name = "Cell type") +
  theme(legend.position = "none")

# Trajectory to ground; choose terminal nodes in ground as endpoints
cds_to_ground <- choose_graph_segments(cds = cds_integrated, clear_cds = FALSE)
cds_to_ground <- order_cells(cds = cds_to_ground)

# Plot key genes over the trajectory
ground_genes <- c("Mfap4", "Itga8")
cds_to_ground_genes <- cds_to_ground[rowData(cds_to_ground)$gene_short_name %in% ground_genes, ]

plot_genes_ground <- plot_genes_in_pseudotime(cds_subset = cds_to_ground_genes)
plot_genes_clusters_ground <- plot_genes_in_pseudotime(cds_subset = cds_to_ground_genes, color_cells_by = "cell_type", panel_order = ground_genes) +
  scale_colour_manual(values = cluster_colors, name = "Cell type") +
  theme(legend.position = "none")

# Save some plots
out_dir <- "../output/plots/"

ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_pseudotime.png", plot = plot_pseudotime, path = out_dir, width = 5.5, height = 3.5)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_pseudotime.svg", plot = plot_pseudotime, path = out_dir, width = 5.5, height = 3.5, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_pseudotime_no_graph.png", plot = plot_pseudotime_no_graph, path = out_dir, width = 5.5, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_pseudotime_no_graph.svg", plot = plot_pseudotime_no_graph, path = out_dir, width = 5.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_cell_type_graph.png", plot = plot_cell_type_graph, path = out_dir, width = 5.5, height = 3.5)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_cell_type_graph.svg", plot = plot_cell_type_graph, path = out_dir, width = 5.5, height = 3.5, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_myo_diff_genes_time.png", plot = plot_genes_rhabdo, path = out_dir, width = 10, height = 12)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_myo_diff_genes_time.svg", plot = plot_genes_rhabdo, path = out_dir, width = 10, height = 12, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_myo_diff_genes_cell_type.png", plot = plot_genes_clusters_rhabdo, path = out_dir, width = 4, height = 5)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_myo_diff_genes_cell_type.svg", plot = plot_genes_clusters_rhabdo, path = out_dir, width = 4, height = 5, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_ground_genes_time.png", plot = plot_genes_ground, path = out_dir, width = 10, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_ground_genes_time.svg", plot = plot_genes_ground, path = out_dir, width = 10, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_ground_genes_cell_type.png", plot = plot_genes_clusters_ground, path = out_dir, width = 4, height = 2)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_ground_genes_cell_type.svg", plot = plot_genes_clusters_ground, path = out_dir, width = 4, height = 2, device = "svg")
