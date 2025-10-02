#!/usr/bin/env python3
"""
Pathway analysis utilities for single-cell RNA-seq analysis
Handles gene set loading, module scores, and GSEA
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from statsmodels.stats.multitest import multipletests


def load_gene_sets():
    """Load gene sets for pathway analysis
    
    Returns:
        Dictionary of gene sets organized by collection
    """
    print("Loading gene sets...")

    # Get gene sets from Enrichr/MSigDB
    gene_sets = {}

    # Hallmark pathways
    try:
        hallmark = gp.get_library(name="MSigDB_Hallmark_2020", organism="Mouse")
        gene_sets["Hallmark"] = hallmark
    except:
        print("Could not load Hallmark gene sets")

    # KEGG pathways
    try:
        kegg = gp.get_library(name="KEGG_2021_Mouse", organism="Mouse")
        gene_sets["KEGG"] = kegg
    except:
        print("Could not load KEGG gene sets")

    # GO Biological Process
    try:
        gobp = gp.get_library(name="GO_Biological_Process_2021", organism="Mouse")
        gene_sets["GO_BP"] = gobp
    except:
        print("Could not load GO BP gene sets")

    # Custom astrocyte-related gene sets
    astrocyte_sets = {
        "ASTROCYTE_ACTIVATION": [
            "Aqp4",
            "Gfap",
            "Vegfa",
            "Edn1",
            "Ptgs2",
            "Pla2g4a",
            "Kcnj10",
            "Slc1a2",
            "Aldh1l1",
        ],
        "CIRCADIAN_CLOCK": [
            "Clock",
            "Bmal1",
            "Per1",
            "Per2",
            "Cry1",
            "Cry2",
            "Nr1d1",
            "Dbp",
            "Tef",
        ],
        "VASOMOTION": [
            "Aqp4",
            "Pla2g",
            "Ptgs1",
            "Cyp2c",
            "Slc16a3",
            "Adora1",
            "Ip3r2",
            "Gfap",
            "Pgt",
            "Kcnj10",
            "Ptgs2",
            "Pla2g4a",
            "Edn1",
            "Vegfa",
        ],
    }
    gene_sets["Custom"] = astrocyte_sets

    return gene_sets


def calculate_module_scores(adata, gene_sets, cell_types=None):
    """Calculate module scores for pathway activity
    
    Args:
        adata: AnnData object
        gene_sets: Dictionary of gene sets
        cell_types: List of cell types to focus on
        
    Returns:
        AnnData object with module scores added
    """
    print("Calculating module scores...")

    if cell_types is None:
        cell_types = ["AST", "ExN", "InN_SST", "InN_VIP"]

    # Calculate scores for custom gene sets
    for pathway_name, genes in gene_sets["Custom"].items():
        # Filter genes present in dataset
        genes_present = [g for g in genes if g in adata.var_names]

        if len(genes_present) > 0:
            print(f"Calculating score for {pathway_name} ({len(genes_present)} genes)")
            sc.tl.score_genes(adata, genes_present, score_name=f"{pathway_name}_score")

    return adata


def plot_module_scores(adata, score_names, groupby="celltype"):
    """Plot module scores across conditions
    
    Args:
        adata: AnnData object with module scores
        score_names: List of score column names
        groupby: Column to group by for plotting
    """
    print("Plotting module scores...")

    fig, axes = plt.subplots(len(score_names), 1, figsize=(10, 4 * len(score_names)))
    if len(score_names) == 1:
        axes = [axes]

    for i, score_name in enumerate(score_names):
        if score_name in adata.obs.columns:
            sc.pl.violin(adata, score_name, groupby=groupby, ax=axes[i], show=False)
            axes[i].set_title(f"{score_name} by {groupby}")

    plt.tight_layout()
    plt.show()


def run_gsea_analysis(de_results, gene_sets):
    """Run GSEA analysis on differential expression results
    
    Args:
        de_results: DataFrame with differential expression results
        gene_sets: Dictionary of gene sets
        
    Returns:
        DataFrame with GSEA results
    """
    print("Running GSEA analysis...")

    gsea_results = []

    for cell_type in de_results["cell_type"].unique():
        print(f"GSEA for {cell_type}")
        ct_data = de_results[de_results["cell_type"] == cell_type]

        for contrast in ct_data["contrast"].unique():
            print(f"  {contrast}")
            contrast_data = ct_data[ct_data["contrast"] == contrast]

            # Create ranked gene list
            contrast_data = contrast_data.dropna(subset=["logFC", "P.Value"])
            if len(contrast_data) < 15:
                continue

            # Rank by signed -log10(p-value)
            contrast_data["rank"] = contrast_data["logFC"] * -np.log10(
                contrast_data["P.Value"]
            )

            # Remove duplicates, keep gene with highest absolute rank
            contrast_data = contrast_data.loc[
                contrast_data.groupby("gene")["rank"].apply(lambda x: x.abs().idxmax())
            ]

            gene_list = contrast_data.set_index("gene")["rank"].sort_values(
                ascending=False
            )

            # Run GSEA for each gene set collection
            for collection_name, pathways in gene_sets.items():
                if collection_name == "Custom":
                    for pathway_name, pathway_genes in pathways.items():
                        try:
                            # Simple GSEA implementation
                            overlap_genes = set(pathway_genes) & set(gene_list.index)
                            if len(overlap_genes) < 5:
                                continue

                            # Get ranks for pathway genes
                            pathway_ranks = gene_list[list(overlap_genes)]

                            # Calculate enrichment score (simplified)
                            es = pathway_ranks.mean()

                            # Simple p-value estimation
                            n_permutations = 1000
                            null_scores = []
                            for _ in range(n_permutations):
                                random_genes = np.random.choice(
                                    gene_list.index, len(overlap_genes), replace=False
                                )
                                null_scores.append(gene_list[random_genes].mean())

                            null_scores = np.array(null_scores)
                            if es >= 0:
                                pval = (null_scores >= es).sum() / n_permutations
                            else:
                                pval = (null_scores <= es).sum() / n_permutations

                            gsea_results.append(
                                {
                                    "pathway": pathway_name,
                                    "NES": es,
                                    "pval": max(pval, 1 / n_permutations),
                                    "padj": pval,  # Will adjust later
                                    "cell_type": cell_type,
                                    "contrast": contrast,
                                    "collection": collection_name,
                                    "size": len(overlap_genes),
                                }
                            )
                        except Exception as e:
                            print(f"Error in GSEA for {pathway_name}: {e}")
                            continue

    if gsea_results:
        gsea_df = pd.DataFrame(gsea_results)
        # Adjust p-values
        gsea_df["padj"] = multipletests(gsea_df["pval"], method="fdr_bh")[1]
        return gsea_df
    else:
        return pd.DataFrame()


def plot_gsea_results(gsea_results, cell_type, contrast, top_n=10):
    """Plot top GSEA results for a specific cell type and contrast
    
    Args:
        gsea_results: DataFrame with GSEA results
        cell_type: Cell type to plot
        contrast: Contrast to plot
        top_n: Number of top pathways to show
    """

    if gsea_results.empty:
        print("No GSEA results to plot")
        return

    # Filter to specific cell type and contrast
    subset = gsea_results[
        (gsea_results["cell_type"] == cell_type)
        & (gsea_results["contrast"] == contrast)
        & (gsea_results["padj"] < 0.1)
    ].copy()

    if subset.empty:
        print(f"No significant pathways for {cell_type} - {contrast}")
        return

    # Get top pathways by absolute NES
    subset = subset.nlargest(top_n, "NES")

    # Plot
    plt.figure(figsize=(10, 6))
    colors = ["red" if x > 0 else "blue" for x in subset["NES"]]
    plt.barh(range(len(subset)), subset["NES"], color=colors, alpha=0.7)
    plt.yticks(range(len(subset)), subset["pathway"])
    plt.xlabel("Normalized Enrichment Score (NES)")
    plt.title(f"Top GSEA results: {cell_type} - {contrast}")
    plt.axvline(x=0, color="black", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()
