# Drosophila Melanogaster Stage 6 Embryogenesis scRNA-seq Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R: >=4.3.0](https://img.shields.io/badge/R-%3E%3D4.3.0-blue.svg)](https://www.r-project.org/)
[![Seurat: v5](https://img.shields.io/badge/Seurat-v5.0.1-brightgreen.svg)](https://satijalab.org/seurat/)

An end-to-end, fully reproducible single-cell RNA sequencing (scRNA-seq) workflow analyzing *Drosophila melanogaster* embryonic development at **Stage 6**. Implemented using **Seurat v5**, this project covers quality control, variance stabilization, dimensionality reduction, unsupervised clustering, differential expression, and lineage annotation with stage-specific developmental markers.

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
├── .gitignore
├── LICENSE
├── README.md
├── environment.yml
├── scripts/
│   ├── 01_load_and_qc.R
│   ├── 02_normalization_dimreduction.R
│   ├── 03_marker_and_target_analysis.R
│   └── 04_cluster_annotation.R
└── results/
    ├── figures/
    │   ├── 01_qc_metrics_prefilter.png
    │   ├── 01_qc_scatter_prefilter.png
    │   ├── 02_variable_genes_plot.png
    │   ├── 03_pca_dimension_plot.png
    │   ├── 04_pca_elbow_plot.png
    │   ├── 05_umap_clusters.png
    │   ├── 06_target_genes_featureplot.png
    │   ├── 07_morphogenetic_targets_dotplot.png
    │   ├── 08_top_markers_heatmap.png
    │   ├── 09_annotated_embryo_umap.png
    │   └── 10_cell_type_proportions.png
    └── tables/
        ├── all_cluster_markers_stage6.csv
        ├── embryonic_cell_proportions.csv
        └── top10_cluster_markers.csv
```

---

## Analysis Workflow

The analysis is modularized into four consecutive steps located in `scripts/`:

### 1. Data Ingestion & Quality Control (`scripts/01_load_and_qc.R`)
- Loads count matrices and constructs Seurat objects.
- Calculates quality metrics: detected features (`nFeature_RNA`), total transcripts (`nCount_RNA`), and mitochondrial percentage (`percent.mt`).
- Applies QC filters to remove low-quality droplets and potential doublets.
- Outputs: `01_qc_metrics_prefilter.png` and `01_qc_scatter_prefilter.png`.

### 2. Normalization & Dimensionality Reduction (`scripts/02_normalization_dimreduction.R`)
- Normalizes expression data and identifies top variable features.
- Runs Principal Component Analysis (PCA) and generates elbow plots to determine dimensionality.
- Performs non-linear manifold approximation via UMAP and computes graph-based Louvain clustering.
- Outputs: `02_variable_genes_plot.png`, `03_pca_dimension_plot.png`, `04_pca_elbow_plot.png`, and `05_umap_clusters.png`.

### 3. Marker Discovery & Target Analysis (`scripts/03_marker_and_target_analysis.R`)
- Performs differential gene expression (Wilcoxon rank-sum test) via `FindAllMarkers`.
- Profiles expression patterns of key morphogenetic and lineage target genes.
- Outputs: `06_target_genes_featureplot.png`, `07_morphogenetic_targets_dotplot.png`, `08_top_markers_heatmap.png`, and marker summary tables (`all_cluster_markers_stage6.csv`, `top10_cluster_markers.csv`).

### 4. Embryonic Cluster Annotation (`scripts/04_cluster_annotation.R`)
- Assigns definitive developmental cell-type identities using canonical Stage 6 markers.
- Computes cell-type population frequencies across clusters.
- Outputs: `09_annotated_embryo_umap.png`, `10_cell_type_proportions.png`, and `embryonic_cell_proportions.csv`.

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

---

## Execution Guide

All analysis scripts are executed from the repository root:

```bash
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

- **R** >= 4.3.0
- **Seurat** >= 5.0.1
- **Matrix**
- **dplyr**
- **ggplot2**
- **patchwork**
- **cowplot**
- **scales**
- **readr**

All packages and version pins are defined in `environment.yml`.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Author

**Narges Shayesteh**  
Bioinformatics & Computational Biology
