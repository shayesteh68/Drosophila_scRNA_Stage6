# ==============================================================================
# Pipeline: Drosophila Embryogenesis Single-Cell Transcriptomics (Stage 6)
# Script 03: Cluster Marker Identification & Morphogenetic/Cytoskeletal Profiling
# Focus: Src42A, Adherens Junctions, and Cytoskeletal Remodeling (Muller Lab)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# 1. Define Paths
project_dir <- "/mnt/f/Drosophila_scRNA_Muller"
input_rds <- file.path(project_dir, "data", "02_seurat_dimreduced.rds")
fig_dir <- file.path(project_dir, "results", "figures")
tbl_dir <- file.path(project_dir, "results", "tables")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(tbl_dir)) dir.create(tbl_dir, recursive = TRUE)

message("==> Step 1: Loading dimreduced Seurat object...")
seurat_obj <- readRDS(input_rds)

# 2. Identify All Cluster Markers
message("==> Step 2: Finding differentially expressed markers for all 7 clusters...")
cluster_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)

# Save cluster markers to CSV
markers_out_file <- file.path(tbl_dir, "all_cluster_markers_stage6.csv")
write.csv(cluster_markers, file = markers_out_file, row.names = FALSE)
message(sprintf("Found %d marker records. Saved to: %s", nrow(cluster_markers), markers_out_file))

# Extract Top 10 Markers per Cluster
top10_per_cluster <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(top10_per_cluster, file = file.path(tbl_dir, "top10_cluster_markers.csv"), row.names = FALSE)

# 3. Muller Lab Target Genes: Cytoskeletal, Polarity, and Junctional Regulators
# Checking which target genes are present in the dataset matrix
potential_targets <- c("Src42A", "arm", "zip", "shot", "twi", "sna", "zen", "sim", "rho", "bkg", "Cka", "shg")
available_targets <- intersect(potential_targets, rownames(seurat_obj))

message(sprintf("Target genes identified in dataset: %s", paste(available_targets, collapse = ", ")))

# 4. Visualization: Feature Plots (UMAP localization)
message("==> Step 3: Generating Spatial Feature Plots for Muller Lab Targets...")
feature_plot <- FeaturePlot(
  seurat_obj,
  features = available_targets[1:min(6, length(available_targets))],
  cols = c("lightgrey", "#1E88E5", "#D81B60"),
  ncol = 3,
  pt.size = 1.0
) & theme(
  plot.title = element_text(size = 12, face = "bold.italic"),
  axis.text = element_blank(),
  axis.ticks = element_blank()
)

ggsave(filename = file.path(fig_dir, "06_target_genes_featureplot.png"), plot = feature_plot, width = 12, height = 7, dpi = 300)

# 5. Visualization: Dot Plot of Morphogenetic Targets
message("==> Step 4: Generating Target Dot Plot across Embryonic Clusters...")
dot_plot <- DotPlot(
  seurat_obj,
  features = available_targets,
  cols = c("lightgrey", "#0D47A1")
) +
  RotatedAxis() +
  ggtitle("Expression of Cytoskeletal & Patterning Regulators across Clusters") +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave(filename = file.path(fig_dir, "07_morphogenetic_targets_dotplot.png"), plot = dot_plot, width = 10, height = 5, dpi = 300)

# 6. Visualization: Heatmap of Top Markers
message("==> Step 5: Generating Marker Heatmap...")
top5_per_cluster <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

heatmap_plot <- DoHeatmap(
  seurat_obj,
  features = unique(top5_per_cluster$gene),
  size = 3.5,
  angle = 45
) + NoLegend() + ggtitle("Top 5 Marker Genes per Embryonic Cluster")

ggsave(filename = file.path(fig_dir, "08_top_markers_heatmap.png"), plot = heatmap_plot, width = 11, height = 8, dpi = 300)

message("==> Step 6: Marker profiling and Muller targets analysis successfully completed!")
