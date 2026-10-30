# Drosophila Melanogaster Stage 6 Embryogenesis scRNA-seq Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R: 4.4.3](https://img.shields.io/badge/R-4.4.3-blue.svg)](https://www.r-project.org/)
[![Seurat: v5.5.0](https://img.shields.io/badge/Seurat-v5.5.0-brightgreen.svg)](https://satijalab.org/seurat/)

A single-cell RNA sequencing (scRNA-seq) workflow analyzing *Drosophila melanogaster* embryonic development at **Stage 6**. Implemented using **Seurat v5**, this project covers quality control, log normalization, dimensionality reduction, unsupervised clustering, differential expression, and lineage annotation with stage-specific developmental markers.

---

## Scientific Background

*Drosophila melanogaster* embryonic Stage 6 marks the onset of gastrulation directly following cellular blastoderm formation. Single-cell transcriptomic profiling at this developmental stage enables high-resolution dissection of:
- Primary germ layers: mesoderm, anterior/posterior midgut invaginations, and neuroectoderm.
- Spatial morphogenetic regulators (e.g., *twist*, *snail*, *single-minded*).
- Lineage commitment dynamics across distinct embryonic territories.

---

## Repository Structure

```text
Drosophila_scRNA_Stage6/
+-- .gitignore
+-- LICENSE
+-- README.md
+-- environment.yml
+-- scripts/
|   +-- 00_download_data.R
|   +-- 01_load_and_qc.R
|   +-- 02_normalization_dimreduction.R
|   +-- 03_marker_and_target_analysis.R
|   +-- 04_cluster_annotation.R
+-- raw_data/ (downloaded input; ignored by Git)
+-- results/ (historical snapshots; regenerate before citing)
    +-- figures/
    +-- tables/
```

---

## Analysis Workflow

The analysis is modularized into four consecutive steps located in `scripts/`:

### 1. Data Ingestion & Quality Control (`scripts/01_load_and_qc.R`)
- Loads count matrices and constructs Seurat objects.
- Calculates detected features (`nFeature_RNA`), total transcripts (`nCount_RNA`), and ribosomal percentage. It calculates mitochondrial percentage only when mitochondrial genes are present in the input.
- Applies feature-count filters and applies the mitochondrial-percentage filter only when that metric is available. The current GSE95025 matrix contains no genes matching the Drosophila mitochondrial prefix, so `percent.mt` is marked unavailable and the mitochondrial filter is skipped with a warning. No dedicated doublet classifier is implemented.
- Outputs: `01_qc_metrics_prefilter.png` and `01_qc_scatter_prefilter.png`.

### 2. Normalization & Dimensionality Reduction (`scripts/02_normalization_dimreduction.R`)
- Normalizes expression data and identifies top variable features.
- Runs Principal Component Analysis (PCA) and generates elbow plots to determine dimensionality.
- Performs non-linear manifold approximation via UMAP and computes graph-based Louvain clustering.
- Outputs: `02_variable_genes_plot.png`, `03_pca_dimension_plot.png`, `04_pca_elbow_plot.png`, and `05_umap_clusters.png`.

### 3. Marker Discovery & Target Analysis (`scripts/03_marker_and_target_analysis.R`)
- Discovers within-cluster marker genes with Seurat's Wilcoxon rank-sum test via `FindAllMarkers`. These cell-level marker p-values are exploratory and do not estimate replicate-level or condition effects.
- Profiles expression patterns of key morphogenetic and lineage target genes.
- Outputs: `06_target_genes_featureplot.png`, `07_morphogenetic_targets_dotplot.png`, `08_top_markers_heatmap.png`, and marker summary tables (`all_cluster_markers_stage6.csv`, `top10_cluster_markers.csv`).

### 4. Embryonic Cluster Annotation (`scripts/04_cluster_annotation.R`)
- Requires the exact current cluster IDs to match the marker-reviewed map. The audited environment produces clusters 0 through 5; cluster 2 has an unresolved ectodermal subtype and cluster 5 has a mixed marker state, so both remain explicitly unresolved.
- Computes cell-type population frequencies across clusters.
- Outputs: `09_annotated_embryo_umap.png`, `10_cell_type_proportions.png`, `embryonic_cell_proportions.csv`, and `cluster_annotation_evidence.csv`.

---

## Installation & Environment Setup

### 1. Clone the Repository
```bash
git clone https://github.com/shayesteh68/Drosophila_scRNA_Stage6.git
cd Drosophila_scRNA_Stage6
```

### 2. Create and Activate the Conda Environment
```bash
conda env create -f environment.yml
conda activate drosophila-scrna-env
```

For an exact Linux x86_64 environment used in the audit, create the environment from the lock file instead:

```bash
conda create -n drosophila-scrna-env --file environment-linux-64.lock.txt
conda activate drosophila-scrna-env
```

---

## Execution Guide

All analysis scripts are executed from the repository root:

```bash
Rscript scripts/00_download_data.R
Rscript scripts/01_load_and_qc.R
Rscript scripts/02_normalization_dimreduction.R
Rscript scripts/03_marker_and_target_analysis.R
Rscript scripts/04_cluster_annotation.R
```

---

## Key Results & Figures

| Figure | Description |
| :--- | :--- |
| **Annotated Embryo UMAP** | `results/figures/09_annotated_embryo_umap.png` |
| **Target Genes FeaturePlot** | `results/figures/06_target_genes_featureplot.png` |
| **Top Markers Heatmap** | `results/figures/08_top_markers_heatmap.png` |
| **Cell Type Proportions** | `results/figures/10_cell_type_proportions.png` |

---

## Dependencies

- **R** 4.4.3
- **Seurat** 5.5.0
- **ggplot2** 4.0.3
- **Matrix**
- **dplyr**
- **patchwork**
- **cowplot**
- **scales**
- **readr**

R, Seurat, and ggplot2 are pinned in environment.yml. Seurat 5.5.0 is used with ggplot2 4.0.3; earlier Seurat releases can fail in plotting functions under ggplot2 4.x ([upstream compatibility issue](https://github.com/satijalab/seurat/issues/10188)). Other R dependencies may change over time when the YAML is resolved. Use the lock file to reproduce the exact audited Linux x86_64 package set.

The workflow reads one expression matrix and has no sample-level replicate design. Cluster marker p-values support exploratory cell annotation only; they are not independent-replicate evidence.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Author

**Narges Shayesteh**  
Bioinformatics & Computational Biology
