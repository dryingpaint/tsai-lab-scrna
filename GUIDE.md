## Quickstart

Run the end‑to‑end pipeline (uses uv to manage the env):

```bash
uv run python cellbender_qc_annotation.py
```

Outputs are written to `plots/` and `annotated_cellbender_data.h5ad`.

## What the pipeline does

1. Load CellBender H5 files and merge samples
2. Compute QC metrics and plot QC summaries
3. Detect doublets (Scrublet) on QC‑filtered cells and transfer labels
4. Filter cells/genes using centralized thresholds
5. Normalize, log1p, select HVGs, scale
6. PCA → kNN graph → UMAP
7. Leiden clustering with automatic resolution selection
8. UMAP plots (clusters and metadata)
9. Marker genes per cluster (csv + plots)
10. Cell type annotation and summary plots
11. Save annotated `.h5ad`

## Where to configure things

- Data paths and samples: `cellbender_qc_annotation.py`
  - `base_path`, `sample_names`, `custom_name`
- QC thresholds and doublet params: `utils/qc_filters.py`
  - `CELL_FILTERS`, `GENE_FILTERS`, `DOUBLET_PARAMS`
- Dimensionality reduction and clustering: `utils/processing.py`
  - `run_pca_umap_clustering(..., n_pcs=15, auto_resolution=True, resolution_grid=None, min_cluster_size=20)`
- Annotation behavior: `utils/annotation.py`
  - used via `annotate_cell_types(adata, label_mode, margin, cluster_agg)` in `cellbender_qc_annotation.py`

## Key outputs (plots/)

- `pca_elbow_plot.png`: choose PCs
- `leiden_resolution_sweep.csv`: per‑resolution metrics
- `clustree_leiden_labels.csv`: wide table of labels for clustree
- `leiden_sweep_diagnostics.png`: silhouette and #clusters vs resolution
- `umap_embeddings.png`: UMAP colored by clusters/metadata
- `marker_genes_dotplot.png`, `top_markers_by_cluster.csv`: marker support
- QC plots: `qc_violin_plots.png`, `qc_scatter_plots.png`

## How parameter selection works (and how to adjust)

### QC thresholds (`utils/qc_filters.py`)

- **min_genes / max_genes**: start with 200–8000; inspect `qc_violin_plots.png` and adjust to trim outliers without removing valid cells.
- **percent_mt (max_mt_pct)**: 5–15 is typical; increase slightly if tissue has high MT content.
- **min_counts / max_counts**: use scatter plots to set reasonable bounds based on library depth.
- **max_ribo_pct**: often left `None`; set if ribosomal content is a concern.

### Doublets (`utils/qc_filters.py` + script)

- **expected_doublet_rate**: platform‑dependent (e.g., 0.06 for 10x high‑throughput; start with 0.06–0.1). We currently use `0.1`.
- The script applies a `manual_threshold=0.35` in doublet detection; lower it for stricter removal if you see residual doublet clusters on UMAP.

### PCs and neighbors (`utils/processing.py`)

- **n_pcs**: default 15. Use `pca_elbow_plot.png`; pick at or just past the elbow (often 20–40 for complex tissues).
- **kNN neighbors**: set to 10 in code. Increase (e.g., 15–30) for smoother manifolds; decrease for finer local structure.

### Leiden clustering (automatic selection)

- The pipeline sweeps resolutions 0.2–2.0 (step 0.1), computes silhouette on PCA and the fraction of small clusters, and picks a stable resolution.
- Selection rule: maximize silhouette; among near‑ties (≤0.02), minimize small‑cluster fraction, then fewer clusters, then lower resolution.
- Chosen value is stored in `adata.uns["leiden_optimal_resolution"]`; `adata.obs["leiden"]` is set to the chosen labels.
- Adjust via `run_pca_umap_clustering(auto_resolution=True, resolution_grid=..., min_cluster_size=...)`.
- To force a resolution: call with `auto_resolution=False, resolution=0.6`.

### Annotation (`cellbender_qc_annotation.py`)

- `label_mode`: `"cell"` (default) for per‑cell labels or `"cluster"` to assign cluster‑level labels.
- `margin`: confidence margin for label assignment; increase slightly to be more conservative.
- `cluster_agg`: summary statistic used when labeling clusters (e.g., `"median"`).

## Using clustree (optional, in R)

You can visualize how clusters split/merge across resolutions with `clustree` using the exported CSV:

```r
library(clustree)
labels <- read.csv("plots/clustree_leiden_labels.csv")
clustree(labels, prefix = "leiden_")
```

Pick the smallest resolution before unstable branching and an explosion of tiny clusters; confirm with marker genes and silhouette.

## Typical tuning workflow

1. Inspect QC plots → adjust `CELL_FILTERS` if needed, re‑run.
2. Check `pca_elbow_plot.png` → set `n_pcs` accordingly.
3. Review `leiden_sweep_diagnostics.png` and UMAPs → accept auto resolution or tweak `min_cluster_size`/grid.
4. Validate clusters with `top_markers_by_cluster.csv` and dotplots → iterate if splits look unbiological.
5. Proceed to DE/GSEA with `limma_voom_gsea.py` when satisfied.
