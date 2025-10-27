# Single-cell RNA-seq Analysis Pipeline (Python)

This repository contains Python scripts for single-cell RNA sequencing data analysis, migrated from R/Seurat to Python/scanpy.

## Overview

The pipeline consists of two main scripts:

1. **`cellbender_qc_annotation.py`** - Data preprocessing, quality control, and cell type annotation
2. **`limma_voom_gsea.py`** - Differential expression analysis and pathway enrichment

## Installation

We use `uv` for fast Python package management. If you don't have `uv` installed:

```bash
# Install uv (macOS/Linux)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Or with pip
pip install uv
```

### Set up the environment:

```bash
# Create virtual environment and install dependencies
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -e .
```

### Alternative: Use uv to run scripts directly (recommended):

```bash
# uv will automatically manage the virtual environment
uv run python cellbender_qc_annotation.py
uv run python limma_voom_gsea.py
```

**Why use uv?**

- 🚀 **10-100x faster** than pip for package resolution and installation
- 🔒 **Automatic dependency resolution** and conflict detection
- 📦 **Built-in virtual environment management** - no need to manually activate/deactivate
- 🔄 **Reproducible builds** with lockfile generation

## Usage

### Step 1: Data Preprocessing and Annotation

```bash
# Run with uv (recommended - automatically manages environment)
uv run python cellbender_qc_annotation.py
```

This script:

- Loads CellBender-processed H5 files from multiple samples
- Performs quality control filtering
- Detects and removes doublets using Scrublet
- Normalizes data and performs dimensionality reduction
- Clusters cells and annotates cell types
- Saves results to `annotated_cellbender_data.h5ad`

### Step 2: Differential Expression and GSEA

```bash
# Run with uv (recommended - automatically manages environment)
uv run python limma_voom_gsea.py
```

This script:

- Loads the annotated data from Step 1
- Calculates pathway module scores
- Creates pseudobulk samples for differential expression
- Performs statistical testing across conditions
- Runs gene set enrichment analysis (GSEA)
- Saves results to CSV files

## Data Structure

The scripts expect data in the following structure:

```
base/
├── D25-2675/
│   └── D25-2675_processed_feature_bc_matrix_filtered.h5
├── D25-2676/
│   └── D25-2676_processed_feature_bc_matrix_filtered.h5
└── ...
```

## Key Features

### Quality Control

- Mitochondrial gene percentage calculation
- Ribosomal gene percentage calculation
- Doublet detection with Scrublet
- Cell and gene filtering based on expression thresholds

### Cell Type Annotation

- Automatic annotation based on marker genes
- Support for major brain cell types:
  - Excitatory neurons (ExN)
  - Inhibitory neurons (InN_SST, InN_VIP)
  - Astrocytes (AST)
  - Oligodendrocytes (ODC)
  - Oligodendrocyte precursors (OPC)
  - Endothelial cells (Endo)
  - Pericytes (Peri)
  - Vascular leptomeningeal cells (VLMC)

### Differential Expression

- Pseudobulk aggregation by sample and cell type
- Multiple contrast testing:
  - Treatment effects within genotypes
  - Genotype effects within conditions
- False discovery rate correction

### Pathway Analysis

- Module score calculation for pathway activity
- Gene Set Enrichment Analysis (GSEA)
- Support for multiple gene set databases:
  - MSigDB Hallmark pathways
  - KEGG pathways
  - GO Biological Process
  - Custom pathway definitions

## Output Files

- `annotated_cellbender_data.h5ad` - Annotated single-cell data
- `differential_expression_results.csv` - DE analysis results
- `gsea_results.csv` - Pathway enrichment results

## Customization

### Modify file paths:

Edit the `base_path` variable in `1_cellbender_qc_annotation.py`:

```python
base_path = "/path/to/your/data/"
```

### Adjust filtering parameters:

Modify QC thresholds in the `filter_cells_and_genes()` function:

```python
adata = adata[adata.obs.n_genes_by_counts < 8000, :]  # Max genes per cell
adata = adata[adata.obs.percent_mt < 10, :]           # Max MT percentage
```

### Add custom gene sets:

Edit the `astrocyte_sets` dictionary in `2_limma_voom_gsea.py`:

```python
astrocyte_sets = {
    'YOUR_PATHWAY': ['Gene1', 'Gene2', 'Gene3'],
    # ...
}
```

## Notes

- The scripts use simplified statistical methods compared to the original R implementation
- For production use, consider implementing more sophisticated DE methods (e.g., using diffxpy or similar packages)
- GSEA implementation is simplified; for comprehensive analysis, consider using dedicated tools like GSEApy with proper gene set databases
