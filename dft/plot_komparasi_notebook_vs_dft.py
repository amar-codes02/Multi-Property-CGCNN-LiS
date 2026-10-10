#!/usr/bin/env python3
"""
VISUALISASI KOMPARASI: NOTEBOOK (Thesis_Principle_Discovery.ipynb) VS PURE DFT HPC
==================================================================================
Membuat visualisasi multi-panel kualitas publikasi membandingkan 5 sifat fisis utama.
"""
import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#94a3b8'
plt.rcParams['axes.linewidth'] = 1.2

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / 'results'
FIG_DIR = BASE_DIR / 'figures'
FIG_DIR.mkdir(parents=True, exist_ok=True)

df_comp = pd.read_csv(RESULTS_DIR / 'komparasi_notebook_vs_dft.csv')

# Urutkan berdasarkan stabilitas termodinamika (Ef terendah ke tertinggi)
df_comp = df_comp.sort_values(by='Ef_DFT [eV/atom]').reset_index(drop=True)

topos = df_comp['Topology'].tolist()
x = np.arange(len(topos))
width = 0.35

palette_topos = {
    'Diamond': '#2563eb',
    'Gyroid': '#f97316',
    'IWP': '#059669',
    'Neovius': '#dc2626',
    'Primitive': '#7c3aed'
}

fig, axes = plt.subplots(2, 2, figsize=(16, 12), dpi=300)

# ── Panel 1: Formation Energy (Ef) ──────────────────────────────────────────
ax1 = axes[0, 0]
rects1_nb = ax1.bar(x - width/2, df_comp['Ef_Notebook [eV/atom]'], width,
                    label='Notebook (Benchmark / ML)', color='#64748b', edgecolor='#1e293b', alpha=0.9, zorder=3)
rects1_dft = ax1.bar(x + width/2, df_comp['Ef_DFT [eV/atom]'], width,
                     label='Pure DFT GPAW (HPC)', color='#2563eb', edgecolor='#1e293b', alpha=0.9, zorder=3)

ax1.set_title('(a) Formation Energy ($E_f$) [Konsistensi Sangat Tinggi]', fontsize=13, fontweight='bold', pad=12)
ax1.set_ylabel('$E_f$ (eV/atom)', fontsize=11, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(topos, fontsize=11, fontweight='bold')
ax1.set_ylim(0, 4.0)
ax1.grid(axis='y', linestyle='--', alpha=0.5, zorder=0)
ax1.legend(frameon=True, facecolor='white', edgecolor='#cbd5e1', fontsize=10, loc='upper left')

for rect in rects1_nb:
    h = rect.get_height()
    ax1.annotate(f'{h:.3f}', (rect.get_x() + rect.get_width() / 2, h),
                 xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=9, color='#334155')

for rect in rects1_dft:
    h = rect.get_height()
    ax1.annotate(f'{h:.3f}', (rect.get_x() + rect.get_width() / 2, h),
                 xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=9, fontweight='bold', color='#1d4ed8')

# ── Panel 2: Band Gap (Eg) ──────────────────────────────────────────────────
ax2 = axes[0, 1]
rects2_nb = ax2.bar(x - width/2, df_comp['Eg_Notebook [eV]'], width,
                    label='Notebook (Asumsi Semimetal = 0 eV)', color='#94a3b8', edgecolor='#1e293b', alpha=0.9, zorder=3)
rects2_dft = ax2.bar(x + width/2, df_comp['Eg_DFT [eV]'], width,
                     label='Pure DFT (Kelengkungan TPMS Quasi-Semimetal)', color='#059669', edgecolor='#1e293b', alpha=0.9, zorder=3)

ax2.set_title('(b) Electronic Band Gap ($E_g$)', fontsize=13, fontweight='bold', pad=12)
ax2.set_ylabel('Band Gap ($E_g$, eV)', fontsize=11, fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels(topos, fontsize=11, fontweight='bold')
ax2.set_ylim(0, 0.22)
ax2.grid(axis='y', linestyle='--', alpha=0.5, zorder=0)
ax2.legend(frameon=True, facecolor='white', edgecolor='#cbd5e1', fontsize=10, loc='upper left')

for rect in rects2_dft:
    h = rect.get_height()
    ax2.annotate(f'{h:.3f} eV', (rect.get_x() + rect.get_width() / 2, h),
                 xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=9.5, fontweight='bold', color='#047857')

# ── Panel 3: Adsorption Energy Affinity ──────────────────────────────────────
ax3 = axes[1, 0]
rects3_nb = ax3.bar(x - width/2, df_comp['Eads_Notebook [eV]'], width,
                    label='Notebook (Model Intrinsik)', color='#f59e0b', edgecolor='#1e293b', alpha=0.9, zorder=3)
# Magnitudo ikatan terkuat DFT (|Eads_best| untuk Li2S2)
rects3_dft = ax3.bar(x + width/2, df_comp['Eads_DFT_Best [eV]'].abs(), width,
                     label='Pure DFT (|Eads, best| Li2S2)', color='#dc2626', edgecolor='#1e293b', alpha=0.9, zorder=3)

ax3.set_title(r'(c) Kekuatan Penjeratan Polisulfida ($|E_{ads}|$, eV)', fontsize=13, fontweight='bold', pad=12)
ax3.set_ylabel('Energi Ikatan Polisulfida (eV)', fontsize=11, fontweight='bold')
ax3.set_xticks(x)
ax3.set_xticklabels(topos, fontsize=11, fontweight='bold')
ax3.set_ylim(0, 5.2)
ax3.grid(axis='y', linestyle='--', alpha=0.5, zorder=0)
ax3.legend(frameon=True, facecolor='white', edgecolor='#cbd5e1', fontsize=10, loc='upper left')

for rect in rects3_nb:
    h = rect.get_height()
    ax3.annotate(f'{h:.2f}', (rect.get_x() + rect.get_width() / 2, h),
                 xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=9, color='#b45309')

for rect in rects3_dft:
    h = rect.get_height()
    ax3.annotate(f'{h:.2f} eV', (rect.get_x() + rect.get_width() / 2, h),
                 xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=9.5, fontweight='bold', color='#b91c1c')

# ── Panel 4: Parity Plot Formation Energy Ef ─────────────────────────────────
ax4 = axes[1, 1]
x_ef = df_comp['Ef_Notebook [eV/atom]']
y_ef = df_comp['Ef_DFT [eV/atom]']

for i, topo in enumerate(topos):
    ax4.scatter(x_ef[i], y_ef[i], s=140, color=palette_topos.get(topo, '#2563eb'),
                edgecolors='#0f172a', linewidth=1.5, zorder=5, label=f"{topo}")
    ax4.annotate(f" {topo} (Δ={y_ef[i]-x_ef[i]:+.3f})", (x_ef[i], y_ef[i]),
                 xytext=(5, -2), textcoords='offset points', fontsize=9.5, fontweight='bold')

min_val = 1.9
max_val = 3.5
ax4.plot([min_val, max_val], [min_val, max_val], 'k--', lw=1.5, alpha=0.7, zorder=2, label='Ideal Parity (y=x)')

# Hitung R2 dan MAE
r2 = np.corrcoef(x_ef, y_ef)[0, 1]**2
mae = np.mean(np.abs(y_ef - x_ef))

ax4.set_title(f'(d) Validasi Paritas $E_f$ ($R^2 = {r2:.4f}$, MAE = {mae:.3f} eV/atom)',
              fontsize=13, fontweight='bold', pad=12)
ax4.set_xlabel('Notebook $E_f$ (eV/atom)', fontsize=11, fontweight='bold')
ax4.set_ylabel('Pure DFT $E_f$ (eV/atom)', fontsize=11, fontweight='bold')
ax4.set_xlim(min_val, max_val)
ax4.set_ylim(min_val, max_val)
ax4.grid(True, linestyle='--', alpha=0.5, zorder=0)
ax4.legend(frameon=True, facecolor='white', edgecolor='#cbd5e1', fontsize=9.5, loc='upper left')

plt.tight_layout(pad=2.5)
out_png = FIG_DIR / 'komparasi_notebook_vs_dft.png'
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.close()

print(f"✅ Gambar perbandingan tersimpan ke: {out_png}")
