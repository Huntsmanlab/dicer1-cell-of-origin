# Script to get a broad-strokes view of the data, during and after QC

# Required packages
library(DropletUtils)
library(scuttle)
library(scran)
library(scater)
library(cowplot)
library(ggplot2)

# Read in the 10X data and metadata
metadata <- read.csv(file = "../metadata/dicer1_scrna_sample_metadata.csv")
sample <- metadata$sample
genotype <- metadata$genotype
time <- metadata$time_months
batch <- metadata$batch
name <- gsub(pattern = "/", replacement = "_", x = metadata$name)
delta_days <- metadata$delta_days

sample_dirs <- paste0("../data/", sample)
sce <- read10xCounts(samples = sample_dirs, sample.names = sample)

# Compute cells per sample to build the metadata per cell
n_cells <- unlist(purrr::map(sample, function(x) sum(sce$Sample == x))) 
genotype_per_cell <- rep(genotype, n_cells)
sce$genotype <- genotype_per_cell
time_per_cell <- rep(time, n_cells)
sce$time <- time_per_cell
batch_per_cell <- rep(batch, n_cells)
sce$batch <- batch_per_cell
name_per_cell <- rep(name, n_cells)
sce$name <- name_per_cell
delta_days_per_cell <- rep(delta_days, n_cells)
sce$delta_days <- delta_days_per_cell

# Key genes for QC
mito_genes <- grepl("^mt", rowData(sce)$Symbol)
ribo_genes <- grepl("^Rp[sl][[:digit:]]", rowData(sce)$Symbol)
tomato_genes <- grepl("tomato", rowData(sce)$Symbol)

# Add QC data to the colData of the sce; sums of the relevant counts and respective percentages will be computed for our gene subsets
feature_controls <- list(mito = rownames(sce)[mito_genes], 
  ribo = rownames(sce)[ribo_genes], 
  tomato = rownames(sce)[tomato_genes])

sce <- addPerCellQC(x = sce, subsets = feature_controls)

# Visualize the distributions of various features before QC
qc_data <- as.data.frame(colData(sce))

umi_plot <- ggplot(data = qc_data, aes(x = sum)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  facet_wrap(facets = vars(name)) +
  xlab(label = "UMIs per cell") +
  ylab(label = "Density")

gene_plot <- ggplot(data = qc_data, aes(x = detected)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  facet_wrap(facets = vars(name)) +
  xlab(label = "Genes per cell") +
  ylab(label = "Density")

genes_vs_umis <- ggplot(data = qc_data, aes(x = sum, y = detected, color = subsets_mito_percent)) +
  geom_point(size = 1) +
  scale_colour_gradient(low = "gray", high = "black") +
  stat_smooth(method = lm) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(facets = vars(name)) +
  xlab(label = "UMIs per cell") +
  ylab(label = "Genes per cell") +
  geom_hline(yintercept = 500, linetype = "dashed") +
  geom_vline(xintercept = 1000, linetype = "dashed")

mito_hist <- ggplot(data = qc_data, aes(x = subsets_mito_percent)) +
  geom_histogram() +
  xlab(label = "Percentage of mitochondrial reads per cell") +
  ylab(label = "Count") +
  facet_wrap(facets = vars(name)) +
  geom_vline(xintercept = 25, linetype = "dashed")

ribo_hist <- ggplot(data = qc_data, aes(x = subsets_ribo_percent)) +
  geom_histogram() +
  xlab(label = "Percentage of ribosomal reads per cell") +
  ylab(label = "Count") +
  facet_wrap(facets = vars(name)) +
  geom_vline(xintercept = 30, linetype = "dashed")

tomato_percent <- sum(sce$subsets_tomato_detected) / length(sce$subsets_tomato_detected)

# Drop cells based on percentage of mitochondrial reads
mito_drop <- sce$subsets_mito_percent > 25

# Drop cells with fewer than 1,000 UMIs or fewer than 500 genes detected
umi_drop <- sce$sum < 1000
gene_drop <- sce$detected < 500

sce_qc <- sce[, !(mito_drop | umi_drop | gene_drop)]

# Normalize by library size
sce_qc <- logNormCounts(sce_qc)

# Get most variable genes, all together and incorporating batch
dec <- modelGeneVar(sce_qc)
hvg <- getTopHVGs(dec, prop = 0.2)

dec_batch <- modelGeneVar(sce_qc, block = sce_qc$batch)
hvg_batch <- getTopHVGs(dec_batch, prop = 0.2)

# Visualize overall spread of samples via UMAP
set.seed(1337)

# Euclidean distance is inaccurate for compositional data but since we have log-normalized by library size, we may use it
sce_qc <- runUMAP(sce_qc, subset_row = hvg_batch, exprs_values = "logcounts", ncomponents = 2, min_dist = 0.75, n_neighbors = 15, metric = "euclidean")

# Visualize the distribution of various features after QC
post_qc_data <- as.data.frame(colData(sce_qc))

umi_plot_qc <- ggplot(data = post_qc_data, aes(x = sum)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  facet_wrap(facets = vars(name)) +
  xlab(label = "UMIs per cell") +
  ylab(label = "Density")

gene_plot_qc <- ggplot(data = post_qc_data, aes(x = detected)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  facet_wrap(facets = vars(name)) +
  xlab(label = "Genes per cell") +
  ylab(label = "Density")

genes_vs_umis_qc <- ggplot(data = post_qc_data, aes(x = sum, y = detected, color = subsets_mito_percent)) +
  geom_point(size = 1) +
  scale_colour_gradient(low = "gray", high = "black") +
  stat_smooth(method = lm) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(facets = vars(name)) +
  xlab(label = "UMIs per cell") +
  ylab(label = "Genes per cell") +
  geom_hline(yintercept = 500, linetype = "dashed") +
  geom_vline(xintercept = 1000, linetype = "dashed")

mito_hist_qc <- ggplot(data = post_qc_data, aes(x = subsets_mito_percent)) +
  geom_histogram() +
  xlab(label = "Percentage of mitochondrial reads per cell") +
  ylab(label = "Count") +
  facet_wrap(facets = vars(name)) +
  geom_vline(xintercept = 25, linetype = "dashed")

ribo_hist_qc <- ggplot(data = post_qc_data, aes(x = subsets_ribo_percent)) +
  geom_histogram() +
  xlab(label = "Percentage of ribosomal reads per cell") +
  ylab(label = "Count") +
  facet_wrap(facets = vars(name)) +
  geom_vline(xintercept = 30, linetype = "dashed")

tomato_percent_qc <- round(sum(sce_qc$subsets_tomato_detected) / length(sce_qc$subsets_tomato_detected) * 100, digits = 1)

umap_time <- plotUMAP(sce_qc, shape_by = "genotype", colour_by = "time", text_size = 12)
umap_tomato <- plotUMAP(sce_qc, colour_by = "subsets_tomato_detected", text_size = 12) +
  ggplot2::ggtitle(label = paste0(tomato_percent_qc, "% of cells positive for tdTomato after QC"))
umap_sample <- plotUMAP(sce_qc, shape_by = "genotype", colour_by = "name", size_by = "detected", text_size = 12)
umap_genotype <- plotUMAP(sce_qc, colour_by = "genotype", text_size = 12)
plot_grid(umap_time, umap_tomato)

# Check the amount of variance in the normalized counts that is accounted for by various variables
variables_plot <- plotExplanatoryVariables(sce_qc, exprs_values = "logcounts", subset_row = hvg_batch, variables = c("sum", "detected", "Sample",
  "time", "genotype", "subsets_mito_percent", "subsets_ribo_percent", "batch"), theme_size = 12)

# Write one SCE per sample, if desired
#for (s in sample) {
#  saveRDS(object = sce_qc[, sce_qc$Sample == s], file = paste0("../Output/QC/", s, "_sce_qc.rds"))
#}

# Save the QC'd SCE for further processing
saveRDS(object = sce_qc, file = "../output/scrna/dicer1_scrna_sce_qc.rds")
