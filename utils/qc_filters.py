#!/usr/bin/env python3
"""
Quality control filter parameters for single-cell RNA-seq analysis

This file centralizes all QC thresholds used in the pipeline.
Modify these values to adjust filtering stringency.
"""

# Cell-level filters
CELL_FILTERS = {
    "min_genes": 200,  # Minimum genes detected per cell (increased from 200)
    "max_genes": 8000,  # Maximum genes detected per cell (reduced from 8000)
    "min_counts": 1000,  # Minimum total counts per cell (increased from 1000)
    "max_counts": 50000,  # Maximum total counts per cell (reduced from 50000)
    "max_mt_pct": 10,  # Maximum mitochondrial gene percentage (reduced from 10)
    "max_ribo_pct": None,  # Maximum ribosomal gene percentage (None = no filter)
}

# Gene-level filters
GENE_FILTERS = {
    "min_cells": 10,  # Minimum cells expressing a gene
}

# Doublet detection parameters
# NOTE: Improve these parameters. Understand what's going on here.
DOUBLET_PARAMS = {
    "expected_doublet_rate": 0.1,  # Expected doublet rate (6%)
    "min_counts": 2,  # Minimum counts for Scrublet
    "min_cells": 3,  # Minimum cells for Scrublet
    "min_gene_variability_pctl": 85,  # Gene variability percentile
    "n_prin_comps": 30,  # Number of principal components
}

"""
# Conceptual workflow:
1. Filter genes/cells (min_counts=2, min_cells=3)
2. Select highly variable genes (top 15%)
3. Simulate artificial doublets by adding pairs of cells
4. Reduce dimensions to 30 PCs
5. Calculate doublet scores (KNN-based)
6. Set threshold based on expected_doublet_rate (6%)
7. Flag cells above threshold as doublets
"""

# Sample-specific adaptive filtering (optional)
ADAPTIVE_FILTERING = {
    "use_adaptive": False,  # Whether to use IQR-based adaptive thresholds
    "iqr_multiplier": 3,  # Number of IQRs for outlier detection
}

# Mitochondrial and ribosomal gene patterns
GENE_PATTERNS = {
    "mt_pattern": "mt-",  # Mouse mitochondrial genes (use "MT-" for human)
    "ribo_pattern": r"^Rp[sl]",  # Ribosomal protein genes
}


# Filtering summary messages
def get_filter_summary():
    """Return a formatted summary of current filter settings"""
    summary = [
        "=== QC Filter Settings ===",
        "\nCell-level filters:",
        f"  - Genes per cell: {CELL_FILTERS['min_genes']} - {CELL_FILTERS['max_genes']}",
        f"  - Counts per cell: {CELL_FILTERS['min_counts']} - {CELL_FILTERS['max_counts']}",
        f"  - Max mitochondrial %: {CELL_FILTERS['max_mt_pct']}%",
    ]

    if CELL_FILTERS["max_ribo_pct"]:
        summary.append(f"  - Max ribosomal %: {CELL_FILTERS['max_ribo_pct']}%")

    summary.extend(
        [
            "\nGene-level filters:",
            f"  - Min cells expressing: {GENE_FILTERS['min_cells']}",
            "\nDoublet detection:",
            f"  - Expected rate: {DOUBLET_PARAMS['expected_doublet_rate']*100}%",
        ]
    )

    return "\n".join(summary)


# Validation function
def validate_filters():
    """Validate that filter parameters make sense"""
    errors = []

    # Check min/max relationships
    if CELL_FILTERS["min_genes"] >= CELL_FILTERS["max_genes"]:
        errors.append("min_genes must be less than max_genes")

    if CELL_FILTERS["min_counts"] >= CELL_FILTERS["max_counts"]:
        errors.append("min_counts must be less than max_counts")

    # Check percentage bounds
    if not 0 <= CELL_FILTERS["max_mt_pct"] <= 100:
        errors.append("max_mt_pct must be between 0 and 100")

    if CELL_FILTERS["max_ribo_pct"] and not 0 <= CELL_FILTERS["max_ribo_pct"] <= 100:
        errors.append("max_ribo_pct must be between 0 and 100")

    # Check doublet parameters
    if not 0 < DOUBLET_PARAMS["expected_doublet_rate"] < 1:
        errors.append("expected_doublet_rate must be between 0 and 1")

    if errors:
        raise ValueError("Filter validation failed:\n" + "\n".join(errors))

    return True


# Run validation on import
validate_filters()
