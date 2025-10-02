#!/usr/bin/env python3
"""
Quality control utilities for single-cell RNA-seq analysis
Handles QC metrics calculation, doublet detection, and filtering
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scrublet as scr


def calculate_qc_metrics(adata):
    """Calculate QC metrics
    
    Args:
        adata: AnnData object
        
    Returns:
        AnnData object with QC metrics added
    """
    print("Calculating QC metrics...")

    # Mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    # Ribosomal genes
    adata.var["ribo"] = adata.var_names.str.match(r"^Rp[sl]")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, percent_top=None, log1p=False, inplace=True, var_type="genes"
    )

    # Add mitochondrial and ribosomal percentages
    adata.obs["percent_mt"] = (
        adata[:, adata.var["mt"]].X.sum(axis=1).A1 / adata.obs["total_counts"]
    ) * 100

    adata.obs["percent_ribo"] = (
        adata[:, adata.var["ribo"]].X.sum(axis=1).A1 / adata.obs["total_counts"]
    ) * 100

    return adata


def plot_qc_metrics(adata):
    """Plot QC metrics
    
    Args:
        adata: AnnData object with QC metrics
    """
    print("Plotting QC metrics...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Violin plots
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts"],
        jitter=0.4,
        multi_panel=True,
        ax=axes[0],
    )

    sc.pl.violin(
        adata, ["percent_mt", "percent_ribo"], jitter=0.4, multi_panel=True, ax=axes[1]
    )

    plt.tight_layout()
    plt.show()

    # Scatter plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    sc.pl.scatter(adata, x="total_counts", y="percent_mt", ax=axes[0])
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", ax=axes[1])

    plt.tight_layout()
    plt.show()


def detect_doublets_scrublet(adata, sample_col="orig.ident"):
    """Detect doublets using Scrublet for each sample
    
    Args:
        adata: AnnData object
        sample_col: Column name for sample identification
        
    Returns:
        AnnData object with doublet predictions added
    """
    print("Detecting doublets with Scrublet...")

    doublet_scores = []
    predicted_doublets = []

    for sample in adata.obs[sample_col].unique():
        print(f"Processing sample: {sample}")

        # Subset to sample
        mask = adata.obs[sample_col] == sample
        X_sample = adata[mask].X.copy()

        # Run Scrublet
        scrub = scr.Scrublet(X_sample, expected_doublet_rate=0.06)
        doublet_scores_sample, predicted_doublets_sample = scrub.scrub_doublets(
            min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30
        )

        doublet_scores.extend(doublet_scores_sample)
        predicted_doublets.extend(predicted_doublets_sample)

    # Add to adata
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    # Summary
    doublet_summary = (
        adata.obs.groupby("orig.ident")
        .agg({"predicted_doublet": ["count", "sum"]})
        .round(2)
    )
    doublet_summary.columns = ["n_cells", "n_doublets"]
    doublet_summary["pct_doublets"] = (
        doublet_summary["n_doublets"] / doublet_summary["n_cells"] * 100
    ).round(2)

    print("Doublet summary:")
    print(doublet_summary)

    return adata


def filter_cells_and_genes(adata, min_genes=200, max_genes=8000, max_mt_pct=10):
    """Apply QC filtering
    
    Args:
        adata: AnnData object
        min_genes: Minimum genes per cell
        max_genes: Maximum genes per cell
        max_mt_pct: Maximum mitochondrial percentage
        
    Returns:
        Filtered AnnData object
    """
    print("Applying QC filters...")

    print(f"Starting with {adata.n_obs} cells and {adata.n_vars} genes")

    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)  # Filter cells with too few genes

    # Filter genes expressed in at least 10 cells
    sc.pp.filter_genes(adata, min_cells=10)

    # Filter cells based on QC metrics
    adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    adata = adata[adata.obs.percent_mt < max_mt_pct, :]

    # Remove predicted doublets
    adata = adata[~adata.obs.predicted_doublet, :]

    print(f"After filtering: {adata.n_obs} cells and {adata.n_vars} genes")

    return adata
