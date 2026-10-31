# ==============================================================================
# Pipeline: Drosophila Embryogenesis Single-Cell Transcriptomics (Stage 6)
# Script 01: Data Ingestion, Quality Control, and Initial Filtering
# Dataset: GEO GSE95025 (Karaiskos et al., Science 2017)
# Scope: Independent expression analysis of a public Stage 6 dataset
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

# 1. Define Paths
script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
project_dir <- if (length(script_file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_file_arg[1])), ".."))
} else {
  normalizePath(getwd())
}
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
if (!is.numeric(counts_matrix) || anyNA(counts_matrix) ||
    any(!is.finite(counts_matrix)) || any(counts_matrix < 0) ||
    any(counts_matrix != floor(counts_matrix)) || anyDuplicated(colnames(counts_matrix))) {
  stop("The input must contain non-negative integer counts and unique cell names.")
}
sparse_counts <- as(counts_matrix, "CsparseMatrix")
rm(counts_raw, counts_matrix)

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
# In Drosophila (FlyBase annotation), mitochondrial genes start with 'mt:'.
mt_features <- grep("^mt:", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
mito_qc_available <- length(mt_features) > 0L
if (mito_qc_available) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = mt_features)
} else {
  seurat_obj[["percent.mt"]] <- rep(NA_real_, ncol(seurat_obj))
  warning("No mitochondrial features were found; percent.mt is unavailable and the mitochondrial filter will be skipped.")
}
# Ribosomal protein genes in Drosophila typically start with RpL or RpS
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[LS]")

# 4. Generate Pre-filter QC Plots
message("==> Step 3: Generating QC Violin Plots...")
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.ribo")
if (mito_qc_available) qc_features <- c(qc_features, "percent.mt")
qc_vln <- VlnPlot(
  seurat_obj,
  features = qc_features,
  ncol = length(qc_features),
  pt.size = 0.1
) & theme(plot.title = element_text(size = 11, face = "bold"))

ggsave(filename = file.path(fig_dir, "01_qc_metrics_prefilter.png"), plot = qc_vln, width = 12, height = 4, dpi = 300)

# QC Scatter Plots
scatter2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
if (mito_qc_available) {
  scatter1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  qc_scatter <- scatter1 + scatter2
} else {
  qc_scatter <- scatter2
}

ggsave(filename = file.path(fig_dir, "01_qc_scatter_prefilter.png"), plot = qc_scatter, width = 10, height = 4, dpi = 300)

# 5. Apply Stringent Filtering
# Retain high-confidence embryonic cells
message("==> Step 4: Filtering low quality cells...")
cell_keep <- seurat_obj$nFeature_RNA >= 500 & seurat_obj$nFeature_RNA <= 6000
if (mito_qc_available) {
  cell_keep <- cell_keep & seurat_obj$percent.mt < 5
} else {
  message("Mitochondrial filtering skipped because this matrix has no mitochondrial genes.")
}
seurat_filtered <- subset(seurat_obj, cells = colnames(seurat_obj)[which(cell_keep)])

message(sprintf("Cells before filtering: %d | Cells after filtering: %d", ncol(seurat_obj), ncol(seurat_filtered)))

# 6. Save Processed Object
saveRDS(seurat_filtered, file = file.path(data_out_dir, "01_seurat_qc_filtered.rds"))
message("==> Step 5: Processed Seurat object successfully saved to: data/01_seurat_qc_filtered.rds")
