#!/usr/bin/env python3
"""01_scrna_analysis.py: Scanpy pipeline on 10x PBMC 3k.

Steps: load 10x matrix -> QC metrics + violins -> filter cells/genes ->
normalize + log1p -> highly variable genes -> regress out technical
covariates -> scale -> PCA -> neighbors -> Leiden -> UMAP ->
marker genes (Wilcoxon) -> cell-type annotation from canonical markers.

Run from the repo root: ../venv/bin/python scripts/01_scrna_analysis.py
Outputs: results/qc_metrics.csv, results/cluster_markers.csv,
results/annotation.csv, figures/00_qc_violins.png, 01_umap_leiden.png,
02_umap_celltype.png, 03_marker_dotplot.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=120, facecolor="white")
np.random.seed(0)

# ---------------- load ----------------
adata = sc.read_10x_mtx(
    "data/pbmc3k/filtered_gene_bc_matrices/hg19/",
    var_names="gene_symbols",
    cache=True,
)
adata.var_names_make_unique()
print(f"loaded: {adata.n_obs} cells x {adata.n_vars} genes")

# ---------------- QC ----------------
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

qc = adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
qc.to_csv("results/qc_metrics.csv")
print("QC summary (before filtering):")
print(qc.describe().round(1).to_string())

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
             jitter=0.4, multi_panel=True, ax=axes, show=False)
fig.tight_layout()
fig.savefig("figures/00_qc_violins.png", bbox_inches="tight")
plt.close(fig)

# ---------------- filter ----------------
print(f"before filter: {adata.n_obs} cells")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
print(f"after filter (min_genes=200, min_cells=3, MT%<5): {adata.n_obs} cells x {adata.n_vars} genes")

# ---------------- normalize / HVG / scale ----------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # keep full log-normalized matrix for marker plots
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print(f"highly variable genes: {int(adata.var['highly_variable'].sum())}")
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)

# ---------------- PCA / neighbors / Leiden / UMAP ----------------
sc.tl.pca(adata, svd_solver="arpack", random_state=0)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, random_state=0)
sc.tl.leiden(adata, resolution=0.5, random_state=0)
sc.tl.umap(adata, random_state=0)
n_clusters = adata.obs["leiden"].nunique()
print(f"Leiden clusters (res=0.5): {n_clusters}")

fig = sc.pl.umap(adata, color="leiden", legend_loc="on data",
                 title="PBMC 3k: Leiden clusters", show=False, return_fig=True)
fig.savefig("figures/01_umap_leiden.png", bbox_inches="tight")
plt.close(fig)

# ---------------- marker genes ----------------
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv("results/cluster_markers.csv", index=False)
print("top markers per cluster:")
for cl in sorted(adata.obs["leiden"].unique(), key=int):
    top = markers[markers["group"] == cl].head(5)
    print(f"  cluster {cl}: " + ", ".join(
        f"{r['names']} ({r['scores']:.1f})" for _, r in top.iterrows()))

# ---------------- annotation from canonical markers ----------------
# Canonical immune markers; a cluster is labeled by the marker set whose
# genes are most enriched among its top-10 Wilcoxon markers.
CANONICAL = {
    "T cells": ["CD3D", "CD3E", "LTB", "IL32", "TRAC", "IL7R", "CCR7", "CD8A", "CD8B"],
    "CD14+ Monocytes": ["CD14", "LYZ", "TYROBP", "AIF1", "LST1", "S100A8", "S100A9", "CST3"],
    "B cells": ["MS4A1", "CD79A", "CD79B"],
    "NK cells": ["GNLY", "NKG7", "GZMA", "CTSW"],
    "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7"],
    "Dendritic cells": ["FCER1A", "HLA-DRA", "HLA-DPB1", "CD74"],
    "Megakaryocytes": ["PPBP", "PF4"],
}
top10 = {cl: set(markers[markers["group"] == cl].head(10)["names"]) for cl in markers["group"].unique()}
rows = []
for cl in sorted(top10, key=int):
    best, best_hits = "Unassigned", []
    for ctype, genes in CANONICAL.items():
        hits = [g for g in genes if g in top10[cl]]
        if len(hits) > len(best_hits):
            best, best_hits = ctype, hits
    rows.append({"cluster": cl, "cell_type": best,
                 "supporting_markers": ";".join(best_hits),
                 "n_cells": int((adata.obs["leiden"] == cl).sum())})
annot = pd.DataFrame(rows)
annot.to_csv("results/annotation.csv", index=False)
print(annot.to_string(index=False))

adata.obs["cell_type"] = adata.obs["leiden"].map(
    dict(zip(annot["cluster"], annot["cell_type"]))).astype("category")

fig = sc.pl.umap(adata, color="cell_type", title="PBMC 3k: annotated cell types",
                 show=False, return_fig=True)
fig.savefig("figures/02_umap_celltype.png", bbox_inches="tight")
plt.close(fig)

# dot plot of one representative marker per annotated type
rep = {"T cells": "CD3D", "CD14+ Monocytes": "CD14", "B cells": "MS4A1",
       "NK cells": "NKG7", "FCGR3A+ Monocytes": "FCGR3A",
       "Dendritic cells": "FCER1A", "Megakaryocytes": "PPBP"}
genes = [g for t, g in rep.items() if t in set(annot["cell_type"])]
adata_plot = adata.raw.to_adata()
adata_plot.obs["cell_type"] = adata.obs["cell_type"].values
dp = sc.pl.dotplot(adata_plot, genes, groupby="cell_type", show=False, return_fig=True)
dp.savefig("figures/03_marker_dotplot.png", bbox_inches="tight")

print("DONE: scRNA-seq pipeline complete.")
