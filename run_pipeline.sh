#!/usr/bin/env bash
set -euo pipefail

echo "[1/4] Fetching MPXV reference CDS and complete genomes from NCBI..."
python3 scripts/01_fetch_mpox_data.py

echo "[2/4] Building local BLAST+ database and extracting orthologous CDS..."
python3 scripts/02_extract_genes_blast.py

echo "[3/4] Validating open reading frames and counting reference SNPs..."
python3 scripts/03_validate_orfs.py

echo "[4/4] Generating conservation heatmap and SNP divergence plots..."
python3 scripts/04_plot_gene_conservation.py

echo "[PIPELINE COMPLETE] All outputs saved to data/extracted_cds/ and results/"
