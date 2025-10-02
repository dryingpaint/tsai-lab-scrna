#!/usr/bin/env python3
"""
Data loading utilities for single-cell RNA-seq analysis
Handles CellBender H5 file loading and merging
"""

import pandas as pd
import h5py
from scipy import sparse
import anndata
from pathlib import Path


def load_cellbender_h5(file_path):
    """Load CellBender processed h5 file

    Args:
        file_path: Path to the CellBender H5 file

    Returns:
        AnnData object with loaded data
    """
    with h5py.File(file_path, "r") as f:
        # Get the matrix data
        matrix = f["matrix"]
        features = f["matrix"]["features"]
        barcodes = f["matrix"]["barcodes"]
        data = f["matrix"]["data"]
        indices = f["matrix"]["indices"]
        indptr = f["matrix"]["indptr"]
        shape = f["matrix"]["shape"]

        # Read the actual values and convert shape properly
        data_vals = data[:]
        indices_vals = indices[:]
        indptr_vals = indptr[:]
        shape_vals = tuple(shape[:])  # Convert to tuple and read values

        # Create sparse matrix (CellBender stores as genes x cells, so we reverse)
        X = sparse.csc_matrix((data_vals, indices_vals, indptr_vals), shape=shape_vals)

        # Get feature names and barcodes
        gene_names = [x.decode("utf-8") for x in features["name"][:]]
        gene_ids = [x.decode("utf-8") for x in features["id"][:]]
        cell_barcodes = [x.decode("utf-8") for x in barcodes[:]]

        # Create AnnData object (transpose if needed to get cells x genes)
        if X.shape[0] == len(gene_names) and X.shape[1] == len(cell_barcodes):
            # Matrix is genes x cells, transpose to cells x genes
            adata = anndata.AnnData(X.T.tocsr())
        else:
            # Matrix is already cells x genes
            adata = anndata.AnnData(X.tocsr())

        adata.var_names = gene_names
        adata.var["gene_ids"] = gene_ids
        adata.obs_names = cell_barcodes

        # Make gene names unique
        adata.var_names_make_unique()

    return adata


def load_and_merge_cellbender_data(base_path, sample_names, custom_name):
    """Load and merge multiple CellBender h5 files

    Args:
        base_path: Base directory path
        sample_names: List of sample names
        custom_name: Custom filename suffix

    Returns:
        Merged AnnData object
    """
    print("Loading CellBender data...")

    adatas = []
    for sample in sample_names:
        file_path = Path(base_path) / sample / f"{sample}{custom_name}"
        print(f"Loading {file_path}")

        adata = load_cellbender_h5(file_path)
        adata.obs["sample"] = sample
        adata.obs["orig.ident"] = sample

        # Add sample prefix to cell barcodes
        adata.obs_names = [f"{sample}_{barcode}" for barcode in adata.obs_names]

        adatas.append(adata)

    # Ensure all samples have the same gene names and make them unique
    for adata in adatas:
        adata.var_names_make_unique()

    # Concatenate all samples
    adata_merged = anndata.concat(adatas, join="outer", fill_value=0)
    adata_merged.var_names_make_unique()

    return adata_merged


def add_metadata(adata, sample_ids):
    """Add experimental metadata to AnnData object

    Args:
        adata: AnnData object
        sample_ids: List of sample identifiers

    Returns:
        AnnData object with added metadata
    """
    print("Adding metadata...")

    # Create metadata mapping
    metadata_df = pd.DataFrame(
        {
            "orig.ident": sample_ids,
            "Genotype": ["E3", "E4", "E3", "E4"] * 4,
            "Stimulation": ["Ctrl"] * 8 + ["GENUS"] * 8,
            "Sex": ["M", "M", "F", "F"] * 4,
        }
    )

    # Map metadata to observations
    for col in ["Genotype", "Stimulation", "Sex"]:
        adata.obs[col] = adata.obs["orig.ident"].map(
            metadata_df.set_index("orig.ident")[col]
        )

    return adata
