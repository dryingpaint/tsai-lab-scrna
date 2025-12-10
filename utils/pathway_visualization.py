"""
Visualization functions for pathway enrichment analysis across cell types.

This module provides functions to create comprehensive visualizations of
pathway enrichment results across multiple cell types and contrasts.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _recover_pathway_names(
    results_df: pd.DataFrame,
    gene_sets: Optional[Dict[str, Dict[str, List[str]]]] = None,
) -> pd.DataFrame:
    """
    Recover pathway names for results that have numeric pathway identifiers.
    
    This function checks if pathway names are numeric and attempts to map them
    back to actual pathway names from gene_sets.
    
    Parameters:
    -----------
    results_df : pd.DataFrame
        GSEA summary dataframe with potentially numeric pathway identifiers
    gene_sets : Optional[Dict[str, Dict[str, List[str]]]]
        The gene sets dictionary used for GSEA (keys are collection names, values are pathway dicts)
        If None, will attempt to detect numeric pathway names but cannot recover them
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with pathway names recovered (if possible)
    """
    if results_df.empty or "pathway" not in results_df.columns:
        return results_df
    
    recovered_df = results_df.copy()
    
    # Check if pathway column contains numeric values
    pathway_sample = recovered_df["pathway"].iloc[0] if len(recovered_df) > 0 else None
    
    # Check if pathways are numeric (integers or numeric strings)
    is_numeric = False
    if pathway_sample is not None:
        try:
            int(pathway_sample)
            is_numeric = True
        except (ValueError, TypeError):
            # Check if it's a numeric string representation
            if isinstance(pathway_sample, str) and pathway_sample.isdigit():
                is_numeric = True
    
    if not is_numeric:
        # Pathway names are already strings, no recovery needed
        return recovered_df
    
    if gene_sets is None:
        print("⚠️  Warning: Pathway names appear to be numeric, but gene_sets not provided.")
        print("   Cannot recover pathway names. Please provide gene_sets parameter or")
        print("   use recover_pathway_names() function from the notebook before plotting.")
        return recovered_df
    
    # Group by collection to map pathway indices to names
    for collection_name in recovered_df["collection"].unique():
        collection_mask = recovered_df["collection"] == collection_name
        
        if collection_name in gene_sets:
            pathway_list = list(gene_sets[collection_name].keys())
            
            # Map numeric pathway identifiers to actual names
            def map_pathway_name(row):
                try:
                    pathway_idx = int(row["pathway"])
                    if 0 <= pathway_idx < len(pathway_list):
                        return pathway_list[pathway_idx]
                    else:
                        return f"{collection_name}_pathway_{pathway_idx}"
                except (ValueError, TypeError):
                    # If pathway is already a string, return as-is
                    return row["pathway"]
            
            recovered_df.loc[collection_mask, "pathway"] = recovered_df[collection_mask].apply(
                map_pathway_name, axis=1
            )
        else:
            print(f"⚠️  Warning: Collection '{collection_name}' not found in gene_sets. Cannot recover pathway names.")
    
    return recovered_df


def plot_pathways_across_cell_types(
    results_df: pd.DataFrame,
    fdr_threshold: float = 0.1,
    max_pathways: int = 20,
    aggregation_method: str = "max_nes",
    collection_filter: Optional[str] = None,
    contrast_filter: Optional[str] = None,
    figsize: tuple = (14, 10),
    output_path: Optional[Path] = None,
    gene_sets: Optional[Dict[str, Dict[str, List[str]]]] = None,
) -> plt.Figure:
    """
    Create a comprehensive visualization of pathways across all cell types.
    
    Creates a heatmap showing pathways (rows) x cell types (columns) with
    color intensity representing Normalized Enrichment Score (NES).
    
    Parameters:
    -----------
    results_df : pd.DataFrame
        GSEA summary dataframe with columns: cell_type, contrast, pathway, nes, fdr, etc.
    fdr_threshold : float
        FDR threshold for filtering pathways
    max_pathways : int
        Maximum number of pathways to display
    aggregation_method : str
        How to aggregate pathways across contrasts: "max_nes", "mean_nes", or "min_fdr"
    collection_filter : Optional[str]
        Filter to specific collection (e.g., "Hallmark", "GO_Biological_Process")
    contrast_filter : Optional[str]
        Filter to specific contrast (e.g., "E3_GENUS_vs_Ctrl")
    figsize : tuple
        Figure size (width, height)
    output_path : Optional[Path]
        Path to save the figure
    gene_sets : Optional[Dict[str, Dict[str, List[str]]]]
        Gene sets dictionary used for GSEA. If provided, will automatically recover
        pathway names if they are numeric identifiers.
    
    Returns:
    --------
    plt.Figure
        The matplotlib figure object
    """
    if results_df.empty:
        print("No enrichment results to plot.")
        return None
    
    # Recover pathway names if needed
    results_df = _recover_pathway_names(results_df, gene_sets)
    
    # Filter data
    filtered_df = results_df[results_df["fdr"] <= fdr_threshold].copy()
    
    if collection_filter:
        filtered_df = filtered_df[filtered_df["collection"] == collection_filter]
    
    if contrast_filter:
        filtered_df = filtered_df[filtered_df["contrast"] == contrast_filter]
    
    if filtered_df.empty:
        print(f"No pathways pass filters (FDR ≤ {fdr_threshold})")
        return None
    
    # Aggregate pathways across contrasts for each cell type
    if aggregation_method == "max_nes":
        # For each pathway-cell_type combination, take the maximum absolute NES
        agg_df = (
            filtered_df.groupby(["pathway", "cell_type"])
            .agg({"nes": lambda x: x.loc[x.abs().idxmax()], "fdr": "min"})
            .reset_index()
        )
    elif aggregation_method == "mean_nes":
        agg_df = (
            filtered_df.groupby(["pathway", "cell_type"])
            .agg({"nes": "mean", "fdr": "min"})
            .reset_index()
        )
    elif aggregation_method == "min_fdr":
        # Take the NES value corresponding to the minimum FDR
        agg_df = (
            filtered_df.groupby(["pathway", "cell_type"])
            .apply(lambda x: x.loc[x["fdr"].idxmin()])
            .reset_index(drop=True)
        )
    else:
        raise ValueError(f"Unknown aggregation_method: {aggregation_method}")
    
    # Select top pathways based on maximum absolute NES across all cell types
    pathway_scores = agg_df.groupby("pathway")["nes"].apply(lambda x: x.abs().max()).sort_values(ascending=False)
    top_pathways = pathway_scores.head(max_pathways).index.tolist()
    
    # Create pivot table for heatmap: pathways x cell_types
    pivot_data = agg_df[agg_df["pathway"].isin(top_pathways)].pivot(
        index="pathway", columns="cell_type", values="nes"
    )
    
    # Sort pathways by their maximum absolute NES
    pivot_data["max_abs_nes"] = pivot_data.abs().max(axis=1)
    pivot_data = pivot_data.sort_values("max_abs_nes", ascending=False).drop(columns="max_abs_nes")
    
    # Create the figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap
    vmax = pivot_data.abs().max().max()
    vmin = -vmax
    
    im = ax.imshow(pivot_data.values, aspect="auto", cmap="RdBu_r", vmin=vmin, vmax=vmax)
    
    # Set ticks and labels
    ax.set_xticks(range(len(pivot_data.columns)))
    ax.set_xticklabels(pivot_data.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(pivot_data.index)))
    ax.set_yticklabels(pivot_data.index)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label="Normalized Enrichment Score (NES)")
    
    # Add text annotations for significant values
    for i in range(len(pivot_data.index)):
        for j in range(len(pivot_data.columns)):
            value = pivot_data.iloc[i, j]
            if not pd.isna(value):
                # Only show text if value is significant
                if abs(value) > 1.0:
                    text_color = "white" if abs(value) > vmax * 0.6 else "black"
                    ax.text(j, i, f"{value:.2f}", ha="center", va="center", 
                           color=text_color, fontsize=8, fontweight="bold")
    
    ax.set_xlabel("Cell Type", fontsize=12, fontweight="bold")
    ax.set_ylabel("Pathway", fontsize=12, fontweight="bold")
    
    title = f"Pathway Enrichment Across Cell Types (FDR ≤ {fdr_threshold})"
    if collection_filter:
        title += f" - {collection_filter}"
    if contrast_filter:
        title += f" - {contrast_filter}"
    ax.set_title(title, fontsize=14, fontweight="bold", pad=20)
    
    plt.tight_layout()
    
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Saved figure to {output_path}")
    
    return fig


def plot_pathways_by_cell_type_grid(
    results_df: pd.DataFrame,
    fdr_threshold: float = 0.1,
    max_pathways_per_cell: int = 10,
    collection_filter: Optional[str] = None,
    contrast_filter: Optional[str] = None,
    figsize_per_subplot: tuple = (8, 6),
    max_cols: int = 3,
    output_path: Optional[Path] = None,
    gene_sets: Optional[Dict[str, Dict[str, List[str]]]] = None,
) -> plt.Figure:
    """
    Create a grid of bar plots showing top pathways for each cell type.
    
    Parameters:
    -----------
    results_df : pd.DataFrame
        GSEA summary dataframe
    fdr_threshold : float
        FDR threshold for filtering pathways
    max_pathways_per_cell : int
        Maximum number of pathways to show per cell type
    collection_filter : Optional[str]
        Filter to specific collection
    contrast_filter : Optional[str]
        Filter to specific contrast
    figsize_per_subplot : tuple
        Size of each subplot
    max_cols : int
        Maximum number of columns in grid
    output_path : Optional[Path]
        Path to save the figure
    gene_sets : Optional[Dict[str, Dict[str, List[str]]]]
        Gene sets dictionary used for GSEA. If provided, will automatically recover
        pathway names if they are numeric identifiers.
    
    Returns:
    --------
    plt.Figure
        The matplotlib figure object
    """
    if results_df.empty:
        print("No enrichment results to plot.")
        return None
    
    # Recover pathway names if needed
    results_df = _recover_pathway_names(results_df, gene_sets)
    
    # Filter data
    filtered_df = results_df[results_df["fdr"] <= fdr_threshold].copy()
    
    if collection_filter:
        filtered_df = filtered_df[filtered_df["collection"] == collection_filter]
    
    if contrast_filter:
        filtered_df = filtered_df[filtered_df["contrast"] == contrast_filter]
    
    if filtered_df.empty:
        print(f"No pathways pass filters (FDR ≤ {fdr_threshold})")
        return None
    
    # Get unique cell types
    cell_types = sorted(filtered_df["cell_type"].unique())
    n_cells = len(cell_types)
    
    if n_cells == 0:
        print("No cell types found.")
        return None
    
    n_cols = min(max_cols, n_cells)
    n_rows = int(np.ceil(n_cells / n_cols))
    
    fig_width = figsize_per_subplot[0] * n_cols
    fig_height = figsize_per_subplot[1] * n_rows
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
    
    if n_cells == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # Plot each cell type
    for idx, cell_type in enumerate(cell_types):
        ax = axes[idx]
        
        cell_data = filtered_df[filtered_df["cell_type"] == cell_type].copy()
        
        # Aggregate across contrasts if needed
        if not contrast_filter and len(cell_data["contrast"].unique()) > 1:
            # Take pathway with maximum absolute NES
            cell_data = (
                cell_data.groupby("pathway")
                .apply(lambda x: x.loc[x["nes"].abs().idxmax()])
                .reset_index(drop=True)
            )
        
        # Get top pathways
        top_up = cell_data[cell_data["nes"] > 0].sort_values("nes", ascending=False).head(max_pathways_per_cell // 2)
        top_down = cell_data[cell_data["nes"] < 0].sort_values("nes", ascending=True).head(max_pathways_per_cell // 2)
        top_pathways = pd.concat([top_up, top_down]).sort_values("nes")
        
        if top_pathways.empty:
            ax.text(0.5, 0.5, f"No pathways\nFDR ≤ {fdr_threshold}", 
                   ha="center", va="center", transform=ax.transAxes)
            ax.set_title(cell_type, fontsize=10, fontweight="bold")
            ax.axis("off")
            continue
        
        # Create bar plot
        pathway_names = top_pathways["pathway"].astype(str)
        colors = top_pathways["nes"].apply(lambda x: "#d7301f" if x > 0 else "#225ea8")
        
        bars = ax.barh(range(len(pathway_names)), top_pathways["nes"], color=colors)
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
        ax.set_yticks(range(len(pathway_names)))
        
        # Truncate long pathway names for display
        display_names = [name[:60] + "..." if len(name) > 60 else name for name in pathway_names]
        ax.set_yticklabels(display_names, fontsize=8)
        ax.set_xlabel("NES", fontsize=9)
        ax.set_title(cell_type, fontsize=10, fontweight="bold")
        ax.grid(axis="x", alpha=0.3, linestyle="--")
    
    # Hide unused subplots
    for idx in range(n_cells, len(axes)):
        axes[idx].axis("off")
    
    title = f"Top Pathways by Cell Type (FDR ≤ {fdr_threshold})"
    if collection_filter:
        title += f" - {collection_filter}"
    if contrast_filter:
        title += f" - {contrast_filter}"
    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.995)
    
    plt.tight_layout()
    
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Saved figure to {output_path}")
    
    return fig
