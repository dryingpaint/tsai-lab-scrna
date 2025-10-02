#!/usr/bin/env python3
"""
Cell type annotation utilities for single-cell RNA-seq analysis
Handles marker gene analysis and cell type assignment
"""

import scanpy as sc
import matplotlib.pyplot as plt


def get_marker_genes():
    """Get dictionary of marker genes for different cell types
    
    Returns:
        Dictionary mapping cell types to marker genes
    """
    marker_genes = {
        "Neuron": ["Snap25", "Rbfox3", "Syp"],
        "Excit": ["Slc17a7", "Camk2a", "Satb2"],
        "Inhib": ["Gad1", "Gad2", "Slc6a1"],
        "Astro": ["Slc1a2", "Slc1a3", "Aqp4", "Aldh1l1", "Gfap"],
        "Oligo": ["Plp1", "Mog", "Mobp", "Mbp"],
        "OPC": ["Pdgfra", "Cspg4"],
        "Micro": ["P2ry12", "Tmem119", "Cx3cr1", "Csf1r", "Sall1", "Aif1"],
        "Endo": ["Pecam1", "Kdr", "Flt1", "Klf2", "Slco1a4"],
        "Peri": ["Pdgfrb", "Rgs5", "Kcnj8", "Abcc9"],
        "VLMC": ["Col1a1", "Col1a2", "Lum", "Dcn"],
    }
    return marker_genes


def get_cluster_to_celltype_mapping():
    """Get mapping from cluster numbers to cell types
    
    Returns:
        Dictionary mapping cluster IDs to cell type names
    """
    cluster_to_celltype = {
        "0": "ExN",
        "1": "ExN",
        "2": "ExN",
        "3": "ExN",
        "4": "InN_SST",
        "5": "AST",
        "6": "ExN",
        "7": "ExN",
        "8": "InN_VIP",
        "9": "ExN",
        "10": "ExN",
        "11": "VLMC",
        "12": "Peri",
        "13": "ODC",
        "14": "Endo",
        "15": "OPC",
        "16": "ExN",
    }
    return cluster_to_celltype


def plot_marker_genes(adata):
    """Plot marker genes across clusters
    
    Args:
        adata: AnnData object with clustering results
    """
    marker_genes = get_marker_genes()
    
    # Plot marker genes
    available_markers = []
    for cell_type, genes in marker_genes.items():
        available = [g for g in genes if g in adata.var_names]
        available_markers.extend(available)

    if available_markers:
        sc.pl.dotplot(
            adata,
            available_markers,
            groupby="leiden",
            standard_scale="var",
            figsize=(15, 8),
        )
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.show()


def annotate_cell_types(adata):
    """Annotate cell types based on marker genes and clustering
    
    Args:
        adata: AnnData object with clustering results
        
    Returns:
        AnnData object with cell type annotations
    """
    print("Annotating cell types...")

    # Plot marker genes first
    plot_marker_genes(adata)

    # Get cluster to cell type mapping
    cluster_to_celltype = get_cluster_to_celltype_mapping()

    # Map clusters to cell types
    adata.obs["celltype"] = adata.obs["leiden"].map(cluster_to_celltype)

    # Plot annotated cell types
    sc.pl.umap(
        adata, color="celltype", legend_loc="on data", title="Cell type annotation"
    )
    plt.show()

    return adata


def plot_cell_type_summary(adata):
    """Plot summary of cell types across samples
    
    Args:
        adata: AnnData object with cell type annotations
    """
    # Cell type counts
    celltype_counts = adata.obs.groupby(['orig.ident', 'celltype']).size().unstack(fill_value=0)
    
    # Plot stacked bar chart
    fig, ax = plt.subplots(figsize=(12, 6))
    celltype_counts.plot(kind='bar', stacked=True, ax=ax)
    plt.title('Cell type distribution across samples')
    plt.xlabel('Sample')
    plt.ylabel('Number of cells')
    plt.xticks(rotation=45, ha='right')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()
    
    # Print summary table
    print("\nCell type summary:")
    print(adata.obs['celltype'].value_counts().sort_index())
