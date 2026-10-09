#!/usr/bin/env bash
# =============================================================================
# SYNC DFT RESULTS & FIGURES DARI HPC KE PC LOKAL
# =============================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HPC_HOST="hpc-ub"
HPC_BASE="/media/user/uid1083/dft_gpaw_graphene_tpms_testing"

echo "🔄 Menarik hasil DFT & figures dari HPC ($HPC_HOST)..."

# Tarik results (kecualikan file .gpw besar berukuran gigabyte)
rsync -avz --progress --exclude '*.gpw' \
    "$HPC_HOST:$HPC_BASE/results/" "$SCRIPT_DIR/results/"

# Tarik figures/gambar hasil plot
rsync -avz --progress \
    "$HPC_HOST:$HPC_BASE/figures/" "$SCRIPT_DIR/figures/"

echo "✅ Sinkronisasi berhasil! File terbaru tersimpan di dft/results/ dan dft/figures/"
