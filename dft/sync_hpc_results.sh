#!/usr/bin/env bash
# Auto-sync results from HPC to Local machine
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOTE_DIR="hpc-ub:/media/user/uid1083/dft_gpaw_graphene_tpms_testing/results/"
LOCAL_DIR="${SCRIPT_DIR}/results/"
LOG_FILE="${SCRIPT_DIR}/sync.log"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Syncing DFT results from HPC..."
rsync -avz --update "$REMOTE_DIR" "$LOCAL_DIR"
cp -f "${LOCAL_DIR}/results.json" "${SCRIPT_DIR}/results.json" 2>/dev/null || true

# Otomatis perbarui file Excel jika ada script export
if [ -f "${SCRIPT_DIR}/export_to_excel.py" ]; then
    python3 "${SCRIPT_DIR}/export_to_excel.py"
fi
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Sync complete!"
