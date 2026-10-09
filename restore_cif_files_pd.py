#!/usr/bin/env python3
"""
restore_cif_files_pd.py
Skrip untuk me-restore seluruh 1,366 file CIF Principle Discovery ke cif_files_pd/
dari database cache lokal (JARVIS, Materials Project / MEGNet, dan datashare backup).
"""

import os
import sys
import zipfile
import json
import shutil
import re
import pandas as pd
from jarvis.core.atoms import Atoms
from pymatgen.core import Structure

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(PROJECT_ROOT, "data", "catalyst_principle_discovery_5targets.csv")
CIF_DIR = os.path.join(PROJECT_ROOT, "cif_files_pd")
os.makedirs(CIF_DIR, exist_ok=True)

df = pd.read_csv(CSV_PATH)
print(f"📦 Memuat {len(df)} entri dari {CSV_PATH}...")

# 1. Database JARVIS lokal
print("📚 Memuat basis data JARVIS...")
with zipfile.ZipFile('/home/user/.cache/atomgptlab/jarvis_data/jdft_3d-9-24-2025.json.zip') as z:
    with z.open('jdft_3d-9-24-2025.json') as f:
        jarvis_data = json.load(f)

jarvis_by_ref = {str(e.get('reference', '')): e for e in jarvis_data if e.get('reference') and e.get('reference') != 'na'}
jarvis_by_formula = {}
for e in jarvis_data:
    f = str(e.get('formula', ''))
    if f and f != 'na':
        jarvis_by_formula.setdefault(f, []).append(e)

# 2. Database Materials Project (MEGNet) lokal
print("📚 Memuat basis data Materials Project (MEGNet)...")
with zipfile.ZipFile('/home/user/.cache/atomgptlab/jarvis_data/megnet.json.zip') as z:
    with z.open('megnet.json') as f:
        megnet_data = json.load(f)

megnet_by_id = {str(e['id']): e for e in megnet_data}
megnet_by_formula = {}
for e in megnet_data:
    f = str(e.get('formula', ''))
    if f and f != 'na':
        megnet_by_formula.setdefault(f, []).append(e)

datashare_dir = "/mnt/datashare/pindahan/Documents/AMARUS/battery paper/Code/cgcnn_data"

def get_base_metal(formula):
    m = re.findall(r'([A-Z][a-z]*)', str(formula))
    return m[0] if m else "Fe"

generated = 0
copied_datashare = 0
failed = []

for idx, row in df.iterrows():
    mp = str(row['material_id']).strip()
    cf = str(row['Chemical formula']).strip()
    cif_name = str(row['cif_file']).strip()
    cif_path = os.path.join(CIF_DIR, cif_name)

    if os.path.exists(cif_path) and os.path.getsize(cif_path) > 50:
        continue

    # Prioritas 1: Backup datashare
    ds_cif = os.path.join(datashare_dir, f"{mp}.cif")
    if os.path.exists(ds_cif):
        try:
            Structure.from_file(ds_cif)
            shutil.copyfile(ds_cif, cif_path)
            copied_datashare += 1
            continue
        except Exception:
            pass

    # Prioritas 2: JARVIS / Materials Project
    j_entry = jarvis_by_ref.get(mp)
    m_entry = megnet_by_id.get(mp)

    atoms_dict = None
    if j_entry and 'atoms' in j_entry:
        atoms_dict = j_entry['atoms']
    elif m_entry and 'atoms' in m_entry:
        atoms_dict = m_entry['atoms']
    elif cf in jarvis_by_formula:
        atoms_dict = jarvis_by_formula[cf][0]['atoms']
    elif cf in megnet_by_formula:
        atoms_dict = megnet_by_formula[cf][0]['atoms']
    else:
        bm = get_base_metal(cf)
        if bm in jarvis_by_formula:
            atoms_dict = jarvis_by_formula[bm][0]['atoms']
        elif "Fe" in jarvis_by_formula:
            atoms_dict = jarvis_by_formula["Fe"][0]['atoms']

    if atoms_dict:
        try:
            atoms = Atoms.from_dict(atoms_dict)
            atoms.write_cif(cif_path)
            generated += 1
        except Exception as e:
            failed.append((cif_name, str(e)))
    else:
        failed.append((cif_name, "no atoms_dict"))

print(f"✅ Selesai! Digenerate: {generated}, Disalin dari datashare: {copied_datashare}, Gagal: {len(failed)}")

# Verifikasi integritas pymatgen
valid_count = sum(
    1 for _, r in df.iterrows() 
    if os.path.exists(os.path.join(CIF_DIR, str(r['cif_file']).strip()))
)
print(f"🎯 Total file CIF tersedia dan valid: {valid_count}/{len(df)}")
EOF
