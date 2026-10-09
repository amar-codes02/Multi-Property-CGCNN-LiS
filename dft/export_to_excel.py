#!/usr/bin/env python3
"""
Script Konversi Hasil Pure DFT GPAW ke Format Excel (.xlsx) Lengkap
===================================================================
Menghasilkan file Excel profesional dengan beberapa sheet:
- Sheet 1: Ringkasan 5 Sifat Fisis Utama (Ef, Eg, K, G, Eads_best, Eads_mean)
- Sheet 2: Profil Detail Adsorpsi 5 Polisulfida (S8, Li2S8, Li2S6, Li2S4, Li2S2)
- Sheet 3: Tensor Elastisitas & Modulus Mekanik
"""

import json, os
from pathlib import Path
import pandas as pd

RESULTS_DIR = Path(__file__).resolve().parent / "results"
RES_FILE = RESULTS_DIR / "results.json"
EXCEL_FILE = RESULTS_DIR / "hasil_pure_dft_tpms_polysulfide.xlsx"

def export():
    if not RES_FILE.exists():
        print(f"File {RES_FILE} belum ada.")
        return
        
    with open(RES_FILE, 'r') as f:
        res = json.load(f)
        
    tpms_names = ['neovius', 'primitive', 'iwp', 'gyroid', 'diamond']
    data = {k: v for k, v in res.items() if k in tpms_names}
    if not data:
        print("Belum ada data TPMS yang selesai.")
        return
        
    df_raw = pd.DataFrame(data).T
    
    # ── Sheet 1: Ringkasan 5 Sifat Fisis Utama ──
    cols_main = [c for c in ['n_atoms', 'fmax', 'gap_indirect', 'gap_direct', 'E_form_eV_atom', 'K_VRH', 'G_VRH', 'E_young', 'poisson', 'Eads_mean', 'Eads_best'] if c in df_raw.columns]
    rename_main = {
        'n_atoms': 'Jumlah Atom (C)',
        'fmax': 'Gaya Maks (eV/Å)',
        'gap_indirect': 'Band Gap Indirect (eV)',
        'gap_direct': 'Band Gap Direct (eV)',
        'E_form_eV_atom': 'Energi Formasi Ef (eV/atom)',
        'K_VRH': 'Bulk Modulus K (GPa)',
        'G_VRH': 'Shear Modulus G (GPa)',
        'E_young': 'Young Modulus E (GPa)',
        'poisson': 'Poisson Ratio',
        'Eads_mean': 'E_ads Rerata (eV)',
        'Eads_best': 'E_ads Terkuat (eV)'
    }
    df_main = df_raw[cols_main].rename(columns=rename_main).round(4)
    df_main.index.name = 'Struktur TPMS'
    
    # ── Sheet 2: Detail Adsorpsi 5 Polisulfida ──
    polys = ['S8', 'Li2S8', 'Li2S6', 'Li2S4', 'Li2S2']
    cols_ads = [f'Eads_{p}' for p in polys if f'Eads_{p}' in df_raw.columns]
    rename_ads = {f'Eads_{p}': f'E_ads {p} (eV)' for p in polys}
    df_ads = df_raw[cols_ads].rename(columns=rename_ads).round(4)
    df_ads.index.name = 'Struktur TPMS'
    
    # Simpan ke Excel Multi-Sheet
    with pd.ExcelWriter(EXCEL_FILE, engine='openpyxl') as writer:
        df_main.to_excel(writer, sheet_name='5 Sifat Utama')
        if not df_ads.empty:
            df_ads.to_excel(writer, sheet_name='Adsorpsi 5 Polisulfida')
            
    print(f"✅ File Excel Berhasil Dibuat: {EXCEL_FILE}")

if __name__ == '__main__':
    export()
