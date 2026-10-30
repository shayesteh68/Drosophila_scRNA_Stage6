# ==============================================================================
# Pipeline: Drosophila Embryogenesis Single-Cell Transcriptomics (Stage 6)
# Script 04: Cell Type / Spatial Domain Annotation & Final Publication Figures
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# 1. Paths
script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
project_dir <- if (length(script_file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_file_arg[1])), ".."))
} else {
  normalizePath(getwd())
}
input_rds <- file.path(project_dir, "data", "02_seurat_dimreduced.rds")
output_rds <- file.path(project_dir, "data", "03_seurat_annotated.rds")
fig_dir <- file.path(project_dir, "results", "figures")
tbl_dir <- file.path(project_dir, "results", "tables")

message("==> Step 1: Loading dimreduced Seurat object...")
seurat_obj <- readRDS(input_rds)

# 2. Assign Biological Annotations based on Top Stage 6 Markers
cluster_annotations <- c(
  "0" = "Mesoderm",
  "1" = "Neuroectoderm",
  "2" = "Ectoderm, subtype unresolved",
  "3" = "Anterior / Head",
  "4" = "Amnioserosa / Dorsal",
  "5" = "Mixed lineage, unresolved"
)

annotation_status_map <- c(
  "0" = "provisional_marker_supported",
  "1" = "provisional_marker_supported",
  "2" = "lineage_unresolved",
  "3" = "provisional_marker_supported",
  "4" = "provisional_marker_supported",
  "5" = "mixed_marker_state_unresolved"
)

annotation_evidence <- c(
  "0" = "twi, sna, and tin are enriched; zfh1 and stumps are top markers",
  "1" = "SoxN, ths, and sog are enriched; ths and SoxN are top markers",
  "2" = "SoxN and sog are detected, but cluster-specific lineage markers are inconclusive",
  "3" = "toy and tll are enriched; both are among top cluster markers",
  "4" = "zen and zen2 are strongly enriched and are top cluster markers",
  "5" = "sog, SoxN, twi, and sna co-occur; mixed lineage signature remains unresolved"
)

observed_clusters <- unique(as.character(Idents(seurat_obj)))
if (!setequal(observed_clusters, names(cluster_annotations))) {
  stop("Current cluster IDs differ from the marker-reviewed map. Review current markers and update cluster_annotations.")
}
if (!setequal(names(annotation_status_map), names(cluster_annotations)) ||
    !setequal(names(annotation_evidence), names(cluster_annotations))) {
  stop("Annotation labels, statuses, and evidence must cover the same cluster IDs.")
}
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
message("==> Step 2: Annotating embryonic cell populations...")
cluster_id <- as.character(Idents(seurat_obj))
seurat_obj$cluster_id <- cluster_id
seurat_obj <- RenameIdents(seurat_obj, cluster_annotations)
seurat_obj$cell_type <- Idents(seurat_obj)
seurat_obj$annotation_status <- unname(annotation_status_map[cluster_id])
warning("Labels are provisional. Clusters 2 and 5 remain unresolved; review current marker evidence before interpretation.")

cluster_sizes <- as.data.frame(table(cluster_id))
colnames(cluster_sizes) <- c("cluster_id", "n_cells")
annotation_evidence_table <- data.frame(
  cluster_id = names(cluster_annotations),
  cell_type = unname(cluster_annotations),
  annotation_status = unname(annotation_status_map[names(cluster_annotations)]),
  marker_evidence = unname(annotation_evidence[names(cluster_annotations)]),
  n_cells = cluster_sizes$n_cells[match(names(cluster_annotations), cluster_sizes$cluster_id)],
  stringsAsFactors = FALSE
)
write.csv(annotation_evidence_table, file.path(tbl_dir, "cluster_annotation_evidence.csv"), row.names = FALSE)

# 3. Custom Color Palette for Embryonic Germ Layers
embryo_cols <- c(
  "Mesoderm" = "#E64B35",
  "Neuroectoderm" = "#4DBBD5",
  "Ectoderm, subtype unresolved" = "#00A087",
  "Anterior / Head" = "#3C5488",
  "Amnioserosa / Dorsal" = "#8491B4",
  "Mixed lineage, unresolved" = "#F39B7F"
)
if (!setequal(names(embryo_cols), unname(cluster_annotations))) {
  stop("The annotation palette does not match the cluster annotation labels.")
}

# 4. Provisional Annotated UMAP Plot
message("==> Step 3: Generating annotated UMAP...")
umap_annotated <- DimPlot(
  seurat_obj,
  reduction = "umap",
  cols = embryo_cols,
  pt.size = 1.0,
  label = FALSE
) +
  ggtitle("Drosophila Stage 6: Provisional Marker Annotations") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 9)
  )

ggsave(file.path(fig_dir, "09_annotated_embryo_umap.png"), plot = umap_annotated, width = 11, height = 6.5, dpi = 300)

# 5. Cell Type Proportion Chart
message("==> Step 4: Calculating and Plotting Cell Type Proportions...")
prop_df <- as.data.frame(table(seurat_obj$cell_type))
colnames(prop_df) <- c("CellType", "Count")
prop_df$Percentage <- (prop_df$Count / sum(prop_df$Count)) * 100

prop_plot <- ggplot(prop_df, aes(x = reorder(CellType, Percentage), y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = embryo_cols) +
  theme_minimal() +
  labs(
    title = paste0("Embryonic Domain Proportions (", format(ncol(seurat_obj), big.mark = ","), " cells)"),
    x = "",
    y = "Percentage of Total Cells (%)"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

ggsave(file.path(fig_dir, "10_cell_type_proportions.png"), plot = prop_plot, width = 9, height = 5, dpi = 300)
write.csv(prop_df, file = file.path(tbl_dir, "embryonic_cell_proportions.csv"), row.names = FALSE)

# 6. Save Annotated RDS Object
message("==> Step 5: Saving final annotated Seurat RDS object...")
saveRDS(seurat_obj, file = output_rds)

message("==> Step 6: Cluster annotation and final visualization completed successfully!")
