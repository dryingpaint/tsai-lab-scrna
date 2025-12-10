
## Pipeline Overview

1. Quality Control (Notebook `cell_annotation_colab/1_setup_qc_filtering.ipynb`)
   
   a. Load and merge CellBender outputs per sample, then attach external metadata columns. The loader prints per-sample cell counts plus a dataset summary so you can spot obvious outliers before QC.
   
   b. Set QC parameters once up-front: pick `default`, `stringent`, or `permissive` thresholds for `min_genes`, `max_genes`, `min_counts`, `max_counts`, and `max_mt_pct`, then override individual values as needed for your tissue/chemistry. Re-run Stage 2 whenever these change so the dashed red guide-lines in the violin/scatter plots update.

```263:303:cell_annotation_colab/1_setup_qc_filtering.ipynb
QC_PRESETS = {
    'default': {'min_genes': 200, 'max_genes': 8000, 'min_counts': 1000,
                'max_counts': 50000, 'max_mt_pct': 10, 'max_ribo_pct': None},
    'stringent': {'min_genes': 500, 'max_genes': 6000, 'min_counts': 1500,
                 'max_counts': 40000, 'max_mt_pct': 5, 'max_ribo_pct': None},
    'permissive': {'min_genes': 100, 'max_genes': 10000, 'min_counts': 500,
                  'max_counts': 60000, 'max_mt_pct': 15, 'max_ribo_pct': None}
}
CELL_FILTERS = QC_PRESETS[SELECTED_PRESET]
```

   c. Compute QC metrics with `calculate_qc_metrics` and visualize violin + scatter plots via `plot_qc_metrics`. Use those plots to place thresholds: keep the main violin peak, trim obvious empty droplets/doublets, and set `max_mt_pct` just above the main mitochondrial peak. The scatter of counts vs. MT% highlights dying cells; counts vs. genes shows doublet shoulders.

```1160:1184:cell_annotation_colab/1_setup_qc_filtering.ipynb
### 💡 Interpreting QC Plots
1. Genes per cell: <200 → empty droplets, >8000 → doublets; move `min_genes`/`max_genes`.
2. Total counts: align with genes, trim extreme outliers via count bounds.
3. Mitochondrial %: >10–20% denotes stress; set `max_mt_pct` per tissue.
4. Scatter plots: inspect counts-vs-genes and counts-vs-MT% for aberrant populations.
```

   d. Detect doublets sample-by-sample with Scrublet. Inspect the per-sample histograms the notebook produces; lower `manual_threshold` (0.25–0.30) when distributions overlap or batch doublet rates fall below expectations, raise it (0.40–0.50) when you are over-calling.

```1579:1609:cell_annotation_colab/1_setup_qc_filtering.ipynb
### 💡 Tuning Doublet Detection
- Expected rates: 10x v2 ≈4–6%, v3 ≈6–8%, high-throughput ≈8–10%.
- Manual threshold: 0.35 default; ↓ for stricter removal, ↑ for permissive runs.
- Always process samples separately; re-check UMAP for remaining doublet clusters.
```

   e. Apply filters in the prescribed order so downstream stats (retention tables, per-sample bar charts) make sense. The code removes low-gene cells first, then rare genes, high-gene cells, high-MT cells, count outliers, ribosomal outliers, and finally doublets; the notebook prints removal counts for each filter so you can judge which threshold bites hardest. Violin QC, doublet histograms, and retention plots (cells per sample, genes per cell) are the key visuals to revisit when tweaking.

```1678:1752:cell_annotation_colab/1_setup_qc_filtering.ipynb
def filter_cells_and_genes(...):
    """Apply QC filtering in the correct order
    1. min_genes → 2. min_cells → 3. max_genes → 4. max_mt_pct →
    5. count filters → 6. ribo% → 7. doublets LAST"""
    ...
```

   f. Save `outputs/qc_filtered_data.h5ad` plus QC plots (`plots/qc_violin_plots.png`, `plots/qc_scatter_plots.png`) and summary tables to feed Notebook 2.

2. Clustering & Marker Discovery (Notebook `cell_annotation_colab/2_clustering_markers.ipynb`)
   
   a. Load the QC-filtered AnnData, confirm `percent_mt`, `doublet_score`, and `orig.ident` exist, and set `N_PCS`, `N_NEIGHBORS`, and `CLUSTERING_PARAMS['resolution']` (the default run used 33 PCs, 10 neighbors, resolution 0.8).
   
   b. Stage 5 normalization rescales per-cell depth to 10k, log1p transforms, finds HVGs, saves the HVG diagnostic plot, and scales to unit variance while preserving raw counts in `adata.raw`. Use the HVG plot to confirm you have enough variable genes (expect a few thousand).

```218:266:cell_annotation_colab/2_clustering_markers.ipynb
def normalize_and_scale(adata):
    print("NORMALIZATION AND SCALING")
    [1] save raw counts → [2] `sc.pp.normalize_total` (target_sum=1e4) →
    [3] `sc.pp.log1p` → [4] `sc.pp.highly_variable_genes` (plots saved) →
    [5] subset to HVGs and `sc.pp.scale` (max_value=10).
```

   c. Stage 6 performs PCA, kNN graph construction, UMAP, and Leiden clustering. Inspect the PCA elbow plot to decide how many PCs to keep, then either run the built-in resolution sweep (silhouette vs. cluster count plot plus CSV) or adjust manually using the notebook’s decision tree for UMAP structure, cluster counts, and cluster sizes. Typical plots: PCA variance curve, UMAP colored by cluster/sample, resolution-sweep diagnostics, and cluster-size bar chart.

```561:643:cell_annotation_colab/2_clustering_markers.ipynb
## 🔍 How to Determine Optimal Leiden Resolution
- Automatic sweep (`choose_leiden_resolution`) optimizes silhouette, small-cluster fraction, and cluster count, writing CSV + diagnostics plot.
- Manual guidance: start 0.6–0.8, inspect UMAP, adjust ±0.2–0.3 if clusters are merged, tiny, or overlapping; keep clusters >50–100 cells.
```

   d. Stage 7 runs `sc.tl.rank_genes_groups` per cluster, exports CSV/plots, and provides a marker-validation checklist (e.g., neuron markers Snap25/Slc17a7, astro markers Gfap/Aqp4). The dotplot and cluster heatmaps are what you use to verify marker specificity and to name clusters before annotation.

   e. Save `outputs/clustered_data.h5ad` together with PCA/UMAP/resolution plots under `plots/notebook2/`.

3. Annotation & Export (Notebook `cell_annotation_colab/3_annotation_export.ipynb`)
   
   a. Configure marker panels spanning major neuron/glia/vascular classes plus subtype markers, and set annotation knobs like `label_mode`, `margin`, `cluster_agg`, and optional neuron re-clustering parameters.

```153:238:cell_annotation_colab/3_annotation_export.ipynb
MARKER_GENES = {
    "Neuron": ["Snap25","Rbfox3","Syp"],
    "Excit": ["Slc17a7","Camk2a","Satb2"],
    "Inhib": ["Gad1","Gad2","Slc6a1",...],
    ...
    "SMC": ["Acta2","Myh11","Tagln"],
}
ANNOTATION_PARAMS = {'margin': 0.05}
RECLUSTERING_PARAMS enable auto-resolution neuron sub-clustering if desired.
```

   b. Stage 8 scores both major cell-type and subtype modules, assigns the highest-scoring label per cell or per cluster, optionally re-clusters excitatory/inhibitory populations with their own HVGs/PCs, and visualizes markers on UMAP.

   c. Use the annotation QC dashboard—UMAP colored by cell type, composition heatmap, stacked bar of per-sample cell types—to evaluate whether ≥80 % of cells receive confident labels. Adjust the confidence margin, extend `MARKER_GENES`, or fall back to cluster-level labeling when >20 % of cells remain “Unlabeled”; if clusters mix multiple types, revisit Notebook 1 thresholds or Notebook 2 resolution. The decision tree also tells you when to raise/lower doublet thresholds or clustering resolution based on heatmap patterns.

```976:1070:cell_annotation_colab/3_annotation_export.ipynb
### 🎛️ Parameter Tuning Guide: Annotation Quality
- If >20% cells are “Unlabeled”, lower `ANNOTATION_PARAMS['margin']` or add tissue-specific markers.
- Mixed clusters: inspect doublet scores; lower `DOUBLET_PARAMS['manual_threshold']` or raise Leiden resolution.
- Entire clusters with one type: ensure marker specificity or reduce resolution; decision tree summarizes fixes.
```

   d. Stage 9 exports `outputs/annotated_data.h5ad`, `cell_metadata.csv`, `celltype_counts.csv`, `analysis_summary.csv`, and plots (`marker_genes_dotplot.png`, `celltype_umap.png`, `composition_heatmap.png`, `celltype_distribution.png`). These are the hand-off inputs for DE analysis.

4. Differential Expression (Notebook `cell_annotation_colab/4_differential_gene.ipynb`)
   
   a. Load `annotated_data.h5ad`, verify required columns, then compute cluster purity: clusters with ≥50 % of one cell type are retained, mixed clusters are documented and dropped. This protects pseudobulk DE from mislabeled or doublet-rich clusters; review the printed purity table and adjust `PURITY_THRESHOLD` if needed.

```120:226:cell_annotation_colab/4_differential_gene.ipynb
# CLUSTER PURITY ANALYSIS
cluster_composition = pd.crosstab(leiden, celltype, normalize='index')
PURITY_THRESHOLD = 0.50
pure_clusters = cluster_purity[cluster_purity['purity'] >= PURITY_THRESHOLD]
adata = adata[adata.obs['leiden'].isin(pure_clusters)].copy()
print retention stats and save `cluster_purity.csv`.
```

   b. Prepare metadata (`Genotype`, `Stimulation`, combined contrasts), create pseudobulk counts per cell type × sample, and check that each contrast has enough replicates; the notebook prints per-cell-type condition counts.

   c. Run `utils.differential_expression.run_de_for_celltype`, which supports PyDESeq2. For each cell type/contrast, the notebook writes full DE tables plus visualization-ready summaries, then generates heatmaps (top genes per condition) and volcano plots (per contrast). Use these plots to decide fold-change and FDR cutoffs to carry forward.

   d. Export everything under `outputs/differential_expression_results/` (full tables, summary stats, top genes, pseudobulk matrices) so downstream notebooks—or R/DESeq2 workflows—can consume the results.

5. Pathway Enrichment (Notebook `cell_annotation_colab/5_pathway_enrichment.ipynb`)
   
   a. Turn each `(cell_type, contrast)` DE table into a preranked list using the signed score `log2FC × –log10(p-value)`; duplicates keep the highest absolute score so the ranking stays monotonic.

```132:199:cell_annotation_colab/5_pathway_enrichment.ipynb
def compute_rank_vector(df, config=RankingConfig(logfc_col="logFC", pval_col="P.Value")):
    rank_score = logFC * safe_negative_log10(p-value);
    deduplicate by abs(score); enforce strict ordering for GSEA.
```

   b. Fetch curated Hallmark, KEGG (mouse), and GO BP libraries through `gseapy.get_library`, optionally filtering to specific cell types or contrasts.

   c. Run `run_prerank_for_group` for every group: 1,000 permutations, gene-set size 15–500, multi-collection support. The notebook saves raw GSEApy result tables, a merged summary CSV, and per-group visualizations under `cell_annotation_colab/outputs/pathway_enrichment/`.

```21:34:cell_annotation_colab/5_pathway_enrichment.ipynb
## Overview & references
1. Assemble ranked lists from `differential_expression_results.csv`.
2. Load Hallmark/KEGG/GO libraries via `gseapy.get_library`.
3. Run `gseapy.prerank` (1,000 perms) per `(cell_type, contrast)`.
4. Save raw + summary outputs under `outputs/pathway_enrichment/`.
5. Plot the top pathways for QC and reporting.
```

   d. Use the `plot_top_pathways` helper (example included) to render bar charts of normalized enrichment scores for any cell type/contrast pair, ensuring the pathways you discuss align with the ranked statistics.

### What to watch as you rerun notebooks
- Rerun Notebook 1 from the parameter cell whenever you change QC thresholds; it regenerates violin/scatter plots and doublet histograms so you can iterate.
- After adjusting `N_PCS`, `N_NEIGHBORS`, or `resolution`, rerun Stage 6 in Notebook 2 so PCA/UMAP plots and the resolution sweep regenerate; use the provided heuristics to keep cluster counts reasonable.
- Tune `ANNOTATION_PARAMS['margin']`/`MARKER_GENES` in Notebook 3 whenever >20 % cells remain unlabeled or when composition heatmaps look mixed.
- In Notebook 4, re-check cluster purity if you change annotations or resolutions—the pseudobulk DE only sees clusters passing the purity filter.
- Notebook 5 assumes the DE CSV schema; regenerate DE results first if you add new contrasts or rename columns.