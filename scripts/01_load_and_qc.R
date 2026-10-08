# ==============================================================================
# Pipeline: Drosophila Embryogenesis Single-Cell Transcriptomics (Stage 6)
# Script 01: Data Ingestion, Quality Control, and Initial Filtering
# Dataset: GEO GSE95025 (Karaiskos et al., Science 2017)
# Target: Cytoskeletal & Morphogenetic Dynamics (Muller Lab Collaboration)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# 1. Define Paths
project_dir <- "."
raw_data_path <- file.path(project_dir, "raw_data", "GSE95025_high_quality_cells_digital_expression.txt.gz")
fig_dir <- file.path(project_dir, "results", "figures")
data_out_dir <- file.path(project_dir, "data")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(data_out_dir)) dir.create(data_out_dir, recursive = TRUE)

message("==> Step 1: Loading raw count matrix...")
# Read gz compressed table natively using gzfile
con <- gzfile(raw_data_path, "rt")
counts_raw <- read.table(con, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
close(con)

counts_matrix <- as.matrix(counts_raw)
sparse_counts <- as(counts_matrix, "CsparseMatrix")

message(sprintf("Raw matrix loaded: %d genes across %d cells.", nrow(sparse_counts), ncol(sparse_counts)))

# 2. Create Seurat Object
message("==> Step 2: Creating Seurat Object...")
seurat_obj <- CreateSeuratObject(
  counts = sparse_counts,
  project = "Drosophila_Stage6_Karaiskos",
  min.cells = 3,
  min.features = 200
)

# 3. Calculate Mitochondrial and Ribosomal QC Metrics
# In Drosophila (FlyBase annotation), mitochondrial genes start with 'mt:'
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt:")
# Ribosomal protein genes in Drosophila typically start with RpL or RpS
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[LS]")

# 4. Generate Pre-filter QC Plots
message("==> Step 3: Generating QC Violin Plots...")
qc_vln <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  ncol = 4,
  pt.size = 0.1
) & theme(plot.title = element_text(size = 11, face = "bold"))

ggsave(filename = file.path(fig_dir, "01_qc_metrics_prefilter.png"), plot = qc_vln, width = 12, height = 4, dpi = 300)

# QC Scatter Plots
scatter1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
scatter2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
qc_scatter <- scatter1 + scatter2

ggsave(filename = file.path(fig_dir, "01_qc_scatter_prefilter.png"), plot = qc_scatter, width = 10, height = 4, dpi = 300)

# 5. Apply Stringent Filtering
# Retain high-confidence embryonic cells
message("==> Step 4: Filtering low quality cells...")
seurat_filtered <- subset(
  seurat_obj,
  subset = nFeature_RNA >= 500 & nFeature_RNA <= 6000 & percent.mt < 5
)

message(sprintf("Cells before filtering: %d | Cells after filtering: %d", ncol(seurat_obj), ncol(seurat_filtered)))

# 6. Save Processed Object
saveRDS(seurat_filtered, file = file.path(data_out_dir, "01_seurat_qc_filtered.rds"))
message("==> Step 5: Processed Seurat object successfully saved to: data/01_seurat_qc_filtered.rds")
