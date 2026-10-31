#!/usr/bin/env python3
"""02_spatial_analysis.py: Squidpy pipeline on 10x Visium human lymph node.

Steps: load visium_hne_adata (cached under data/visium/) -> per-spot QC ->
normalize + log1p -> HVG -> PCA -> neighbors -> Leiden -> spatial maps of
clusters and of B/T-cell marker genes -> Moran's I spatial autocorrelation.

Run from the repo root: ../venv/bin/python scripts/02_spatial_analysis.py
Outputs: results/visium_qc.csv, results/visium_markers.csv,
results/moransI_top.csv, figures/04_visium_clusters.png,
figures/05_visium_genes.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

sc.settings.verbosity = 1
np.random.seed(0)

# ---------------- load ----------------
adata = sq.datasets.visium_hne_adata(path="data/visium/")
print(f"loaded: {adata.n_obs} spots x {adata.n_vars} genes")
print("obs columns:", list(adata.obs.columns))

# ---------------- QC ----------------
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].to_csv("results/visium_qc.csv")
print("QC summary:")
print(adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].describe().round(1).to_string())

sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_genes(adata, min_cells=10)
print(f"after filter: {adata.n_obs} spots x {adata.n_vars} genes")

# ---------------- normalize / cluster ----------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
print(f"HVGs: {int(adata.var['highly_variable'].sum())}")
sc.pp.pca(adata, random_state=0)
sc.pp.neighbors(adata, random_state=0)
sc.tl.leiden(adata, resolution=0.5, random_state=0, key_added="leiden")
adata.uns.pop("leiden_colors", None)  # drop stale palette shipped with the dataset
print(f"Leiden clusters: {adata.obs['leiden'].nunique()}")

# ---------------- spatial plots (matplotlib over the H&E image) ----------------
# Rendered manually: squidpy loads coordinates/image, matplotlib draws.
lib_id = list(adata.uns["spatial"].keys())[0]
sp = adata.uns["spatial"][lib_id]
img = sp["images"]["hires"]
scalef = sp["scalefactors"]["tissue_hires_scalef"]
coords = adata.obsm["spatial"] * scalef

fig, ax = plt.subplots(figsize=(7, 7))
ax.imshow(img)
cats = adata.obs["leiden"].astype(str)
for cl in sorted(cats.unique(), key=lambda x: int(x)):
    m = (cats == cl).values
    ax.scatter(coords[m, 0], coords[m, 1], s=14, label=f"cluster {cl}", alpha=0.9)
ax.legend(markerscale=2, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
ax.set_title("Visium mouse brain: Leiden clusters")
ax.axis("off")
fig.savefig("figures/04_visium_clusters.png", bbox_inches="tight")
plt.close(fig)

# Two contrasting brain markers in space (mouse symbols; chosen from the data)
for g in ["Mbp", "Gad1"]:
    assert g in adata.var_names, f"{g} missing from Visium genes"
expr = {g: adata[:, g].X.toarray().ravel() for g in ["Mbp", "Gad1"]}
fig, axes = plt.subplots(1, 2, figsize=(13, 6))
for axi, (g, title) in zip(axes, [("Mbp", "Mbp (myelin, white matter)"), ("Gad1", "Gad1 (inhibitory interneurons)")] ):
    axi.imshow(img)
    sc_ = axi.scatter(coords[:, 0], coords[:, 1], s=14, c=expr[g], cmap="viridis", alpha=0.9)
    plt.colorbar(sc_, ax=axi, shrink=0.7, label="log1p expression")
    axi.set_title(title)
    axi.axis("off")
fig.savefig("figures/05_visium_genes.png", bbox_inches="tight")
plt.close(fig)

# ---------------- markers per spatial cluster ----------------
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
vm = sc.get.rank_genes_groups_df(adata, group=None)
vm.to_csv("results/visium_markers.csv", index=False)
print("top spatial-cluster markers:")
for cl in sorted(adata.obs["leiden"].unique(), key=int):
    top = vm[vm["group"] == cl].head(3)
    print(f"  cluster {cl}: " + ", ".join(top["names"]))

# ---------------- Moran's I spatial autocorrelation ----------------
sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)
hvg = adata.var_names[adata.var["highly_variable"]].tolist()
sq.gr.spatial_autocorr(adata, mode="moran", genes=hvg)
mi = adata.uns["moranI"][["I", "pval_norm"]].copy()
mi = mi.rename(columns={"I": "moransI", "pval_norm": "pval"})
mi["gene"] = mi.index
mi = mi.sort_values("moransI", ascending=False).reset_index(drop=True)
mi[["gene", "moransI", "pval"]].head(20).to_csv("results/moransI_top.csv", index=False)
print("top 5 spatially autocorrelated genes (Moran's I):")
print(mi[["gene", "moransI", "pval"]].head(5).round(4).to_string(index=False))

print("DONE: spatial pipeline complete.")
