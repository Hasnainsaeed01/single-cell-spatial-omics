#!/usr/bin/env bash
# Download the two public 10x Genomics demo datasets used in this project.
# Run from the repo root: bash data/download.sh
set -euo pipefail
mkdir -p data/pbmc3k data/visium

echo "--- PBMC 3k filtered feature-barcode matrix (10x Genomics) ---"
if [ ! -f data/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz ]; then
  curl -L --fail -o data/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz \
    "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
fi
tar -xzf data/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/pbmc3k

echo "--- Visium spatial dataset (via squidpy.datasets, cached under data/visium) ---"
echo "Downloaded programmatically by scripts/02_spatial_analysis.py using"
echo "squidpy.datasets.visium_hne_adata(); the raw files are cached by scanpy"
echo "under ./data/visium on first run. See data/METADATA for details."
echo "Done."
