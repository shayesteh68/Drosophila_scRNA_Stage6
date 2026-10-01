# ==============================================================================
# Pipeline: Drosophila Embryogenesis Single-Cell Transcriptomics (Stage 6)
# Script 02: Normalization, HVG Identification, PCA, UMAP, and Graph Clustering
# Dataset: GEO GSE95025 (Karaiskos et al., Science 2017)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# 1. Define Paths
project_dir <- "/mnt/f/Drosophila_scRNA_Muller"
input_rds <- file.path(project_dir, "data", "01_seurat_qc_filtered.rds")
fig_dir <- file.path(project_dir, "results", "figures")
data_out_dir <- file.path(project_dir, "data")

message("==> Step 1: Loading QC-filtered Seurat object...")
seurat_obj <- readRDS(input_rds)

# 2. Normalization
message("==> Step 2: Normalizing expression counts...")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

# 3. Identify Highly Variable Features (HVGs)
message("==> Step 3: Finding Highly Variable Features (HVGs)...")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 1000, verbose = FALSE)

top10_genes <- head(VariableFeatures(seurat_obj), 10)
message(sprintf("Top 10 HVGs identified: %s", paste(top10_genes, collapse = ", ")))

# Plot Variable Features
hvg_plot1 <- VariableFeaturePlot(seurat_obj)
hvg_plot2 <- LabelPoints(plot = hvg_plot1, points = top10_genes, repel = TRUE, max.overlaps = 20)
ggsave(filename = file.path(fig_dir, "02_variable_genes_plot.png"), plot = hvg_plot2, width = 9, height = 5, dpi = 300)

# 4. Scaling
message("==> Step 4: Scaling data across all features...")
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes, verbose = FALSE)

# 5. Principal Component Analysis (PCA)
message("==> Step 5: Running PCA...")
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), verbose = FALSE)

# Save PCA Diagnostics
pca_dim_plot <- DimPlot(seurat_obj, reduction = "pca", pt.size = 1.2) + ggtitle("PCA - Drosophila Embryo (Stage 6)")
ggsave(filename = file.path(fig_dir, "03_pca_dimension_plot.png"), plot = pca_dim_plot, width = 7, height = 5, dpi = 300)

elbow_plot <- ElbowPlot(seurat_obj, ndims = 30) + ggtitle("Elbow Plot - Dimensionality Assessment")
ggsave(filename = file.path(fig_dir, "04_pca_elbow_plot.png"), plot = elbow_plot, width = 6, height = 4, dpi = 300)

# 6. Graph-based Clustering and Non-linear Dimensionality Reduction (UMAP)
# Using top 15 PCs (captures major spatial axes of embryo blastoderm)
dims_to_use <- 1:15
message(sprintf("==> Step 6: Constructing SNN Graph and Running UMAP (Dims: 1 to %d)...", max(dims_to_use)))

seurat_obj <- FindNeighbors(seurat_obj, dims = dims_to_use, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.6, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = dims_to_use, verbose = FALSE)

# 7. Generate UMAP Visualizations
message("==> Step 7: Saving UMAP Visualizations...")
umap_plot <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  pt.size = 1.2
) +
  ggtitle("Single-Cell Transcriptomic Map: Drosophila Stage 6 Embryo") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

ggsave(filename = file.path(fig_dir, "05_umap_clusters.png"), plot = umap_plot, width = 8, height = 6, dpi = 300)

# 8. Save Clustered Object
saveRDS(seurat_obj, file = file.path(data_out_dir, "02_seurat_dimreduced.rds"))
message(sprintf("==> Step 8: Clustering complete. Identified %d clusters across %d cells.",
                length(unique(Idents(seurat_obj))), ncol(seurat_obj)))
message("==> Processed Seurat object saved to: data/02_seurat_dimreduced.rds")
