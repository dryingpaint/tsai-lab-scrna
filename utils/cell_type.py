import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from utils.processing import choose_leiden_resolution

def plot_cell_type_summary(adata, save_dir=None):
    """Plot summary of cell types across samples

    Args:
        adata: AnnData object with cell type annotations
        save_dir: Directory to save plots (optional). If provided, plots are saved without display.
    """
    # Cell type counts
    celltype_counts = (
        adata.obs.groupby(["orig.ident", "celltype"]).size().unstack(fill_value=0)
    )

    # Plot stacked bar chart
    fig, ax = plt.subplots(figsize=(12, 6))
    celltype_counts.plot(kind="bar", stacked=True, ax=ax)
    plt.title("Cell type distribution across samples")
    plt.xlabel("Sample")
    plt.ylabel("Number of cells")
    plt.xticks(rotation=45, ha="right")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    if save_dir:
        fig.savefig(
            save_dir / "celltype_distribution.png", dpi=300, bbox_inches="tight"
        )
        print(f"  Saved: {save_dir}/celltype_distribution.png")
        plt.close(fig)
    else:
        plt.show()

    # Print summary table
    print("\nCell type summary:")
    print(adata.obs["celltype"].value_counts().sort_index())


def compute_top_markers_per_cluster(
    adata,
    groupby="leiden",
    method="wilcoxon",
    n_top=30,
    pval_adj_cutoff=None,
    save_dir=None,
    plot=False,
):
    """Compute top marker genes per cluster using differential expression.

    Args:
        adata: AnnData object with clustering results.
        groupby: Column in adata.obs to group by (default: "leiden").
        method: DE method passed to scanpy (e.g., "wilcoxon", "t-test").
        n_top: Number of top genes to rank per group.
        pval_adj_cutoff: Optional adjusted p-value cutoff to filter results.
        save_dir: Optional Path to save a CSV summary and optional plots.
        plot: If True, create a rank_genes_groups plot (saved if save_dir provided).

    Returns:
        Pandas DataFrame with ranked markers across all groups.
    """
    if groupby not in adata.obs:
        raise KeyError(f"Groupby key '{groupby}' not found in adata.obs")

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        n_genes=int(n_top),
        pts=True,
    )

    markers_df = sc.get.rank_genes_groups_df(adata, None)
    # Ensure consistent columns exist across scanpy versions
    if "pvals_adj" in markers_df.columns and pval_adj_cutoff is not None:
        markers_df = markers_df[markers_df["pvals_adj"] <= float(pval_adj_cutoff)]

    if save_dir is not None:
        out_csv = save_dir / "top_markers_by_cluster.csv"
        markers_df.to_csv(out_csv, index=False)
        print(f"  Saved: {out_csv}")

    if plot:
        sc.pl.rank_genes_groups(adata, n_genes=min(n_top, 20), sharey=False, show=False)
        if save_dir is not None:
            out_png = save_dir / "top_markers_ranked.png"
            plt.savefig(out_png, dpi=300, bbox_inches="tight")
            print(f"  Saved: {out_png}")
            plt.close()
        else:
            plt.show()

    return markers_df


def compare_top_markers_to_expected(
    adata,
    markers_df=None,
    groupby="leiden",
    top_n=10,
    panels=None,
    save_dir=None,
    plot=True,
):
    """Compare top DE genes per cluster with expected marker panels.

    Builds overlap metrics between each cluster's top-N DE genes and each
    expected marker gene panel.

    Args:
        adata: AnnData with DE results in .uns["rank_genes_groups"] or provide markers_df.
        markers_df: Optional DataFrame from sc.get.rank_genes_groups_df(adata, None).
        groupby: Cluster column used for DE (default: "leiden").
        top_n: Number of top genes per cluster to evaluate.
        panels: Optional dict mapping panel name -> list of genes. Defaults to MARKER_GENES.
        save_dir: Optional Path to write CSVs and heatmap.
        plot: If True, save a heatmap of precision (overlap/top_n).

    Returns:
        A tuple of (long_df, precision_matrix) where long_df is a list of dicts
        usable as rows for a DataFrame, and precision_matrix is a dict of
        {group: {panel: precision}}.
    """
    if panels is None:
        panels = MARKER_GENES

    # Make sets for faster overlap
    panel_to_genes = {k: set(v) for k, v in panels.items()}

    # Acquire ranked genes per group
    if markers_df is None:
        markers_df = sc.get.rank_genes_groups_df(adata, None)

    # Determine sort key preference
    sort_key = (
        "scores"
        if "scores" in markers_df.columns
        else ("logfoldchanges" if "logfoldchanges" in markers_df.columns else None)
    )
    if sort_key is None:
        raise KeyError("markers_df must contain 'scores' or 'logfoldchanges' column")

    # Build top-N gene sets per group
    var_names = set(adata.var_names)
    groups = sorted(markers_df["group"].unique())
    group_to_top = {}
    for g in groups:
        sub = markers_df[markers_df["group"] == g]
        sub = sub.sort_values(sort_key, ascending=False).head(int(top_n))
        genes = [nm for nm in sub["names"].tolist() if nm in var_names]
        group_to_top[g] = set(genes)

    # Compute overlaps
    long_rows = []
    precision_matrix = {g: {} for g in groups}
    for g in groups:
        top_set = group_to_top[g]
        for panel_name, panel_genes in panel_to_genes.items():
            overlap = len(top_set & panel_genes)
            precision = overlap / max(1, len(top_set))
            recall = overlap / max(1, len(panel_genes))
            denom = len(top_set | panel_genes)
            jaccard = overlap / max(1, denom)
            long_rows.append(
                {
                    "group": g,
                    "panel": panel_name,
                    "overlap": overlap,
                    "top_n": len(top_set),
                    "panel_size": len(panel_genes),
                    "precision": precision,
                    "recall": recall,
                    "jaccard": jaccard,
                }
            )
            precision_matrix[g][panel_name] = precision

    # Optional outputs
    if save_dir is not None:
        try:
            import pandas as _pd  # Local import to avoid hard dependency elsewhere

            long_df = _pd.DataFrame(long_rows)
            out_csv = save_dir / "expected_marker_overlap_long.csv"
            long_df.to_csv(out_csv, index=False)
            print(f"  Saved: {out_csv}")

            # Pivot to matrix for heatmap
            mat = long_df.pivot(index="group", columns="panel", values="precision")
            out_mat = save_dir / "expected_marker_overlap_matrix_precision.csv"
            mat.to_csv(out_mat)
            print(f"  Saved: {out_mat}")

            if plot:
                fig, ax = plt.subplots(
                    figsize=(
                        max(6, len(mat.columns) * 0.6),
                        max(4, len(mat.index) * 0.4),
                    )
                )
                im = ax.imshow(
                    mat.values, aspect="auto", cmap="viridis", vmin=0, vmax=1
                )
                ax.set_xticks(range(len(mat.columns)))
                ax.set_xticklabels(mat.columns, rotation=45, ha="right")
                ax.set_yticks(range(len(mat.index)))
                ax.set_yticklabels(mat.index)
                ax.set_title("Precision: overlap of top-N vs expected markers")
                cbar = fig.colorbar(im, ax=ax)
                cbar.set_label("precision (overlap/top_n)")
                plt.tight_layout()
                out_png = save_dir / "expected_marker_overlap_heatmap.png"
                plt.savefig(out_png, dpi=300, bbox_inches="tight")
                print(f"  Saved: {out_png}")
                plt.close(fig)
        except Exception as e:
            print(f"Warning: could not write overlap CSVs/plot: {e}")

    return long_rows, precision_matrix
