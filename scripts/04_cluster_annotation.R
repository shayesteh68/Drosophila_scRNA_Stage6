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
project_dir <- "/mnt/f/Drosophila_scRNA_Muller"
input_rds <- file.path(project_dir, "data", "02_seurat_dimreduced.rds")
output_rds <- file.path(project_dir, "data", "03_seurat_annotated.rds")
fig_dir <- file.path(project_dir, "results", "figures")
tbl_dir <- file.path(project_dir, "results", "tables")

message("==> Step 1: Loading dimreduced Seurat object...")
seurat_obj <- readRDS(input_rds)

# 2. Assign Biological Annotations based on Top Stage 6 Markers
cluster_annotations <- c(
  "0" = "Mesoderm (twi+/sna+)",
  "1" = "Neurogenic Ectoderm (vnd+/SoxN+)",
  "2" = "Lateral Trunk Ectoderm",
  "3" = "Anterior / Head (toy+/tll+)",
  "4" = "Lateral Ectoderm (sog+)",
  "5" = "Amnioserosa / Dorsal (zen+/zen2+)",
  "6" = "Mesectoderm / Midline (sim+)"
)

message("==> Step 2: Annotating embryonic cell populations...")
seurat_obj <- RenameIdents(seurat_obj, cluster_annotations)
seurat_obj$cell_type <- Idents(seurat_obj)

# 3. Custom Color Palette for Embryonic Germ Layers
embryo_cols <- c(
  "Mesoderm (twi+/sna+)"            = "#E64B35", # Red / Ventral
  "Neurogenic Ectoderm (vnd+/SoxN+)"= "#4DBBD5", # Blue / Ventral-Lateral
  "Lateral Trunk Ectoderm"          = "#00A087", # Teal / Lateral
  "Anterior / Head (toy+/tll+)"     = "#3C5488", # Dark Blue / Anterior
  "Lateral Ectoderm (sog+)"         = "#F39B7F", # Coral / BMP gradient
  "Amnioserosa / Dorsal (zen+/zen2+)"= "#8491B4", # Lavender / Dorsal-most
  "Mesectoderm / Midline (sim+)"    = "#91D1C2"  # Light Teal / Midline
)

# 4. Publication-Ready UMAP Plot
message("==> Step 3: Generating Annotated Publication UMAP...")
umap_annotated <- DimPlot(
  seurat_obj,
  reduction = "umap",
  cols = embryo_cols,
  pt.size = 1.0,
  label = TRUE,
  label.size = 3.5,
  repel = TRUE
) +
  ggtitle("Drosophila Blastoderm (Stage 6) - Spatial Germ Layers") +
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
    title = "Embryonic Domain Proportions (1,297 High-Quality Cells)",
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
