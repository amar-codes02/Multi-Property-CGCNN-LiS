#!/usr/bin/env python3
"""
MASTER PURE DFT GPAW RUNNER — 5 GRAPHENE TPMS & 5 POLYSULFIDES (S8 -> Li2S2)
[VERSI SINGLE-POINT DIRECT SCF — TANPA RELAKSASI GEOMETRI]
=============================================================================
Struktur TPMS : Neovius (188 C), Primitive (200 C), IWP (228 C),
                Gyroid (244 C), Diamond (332 C)
Adsorbat      : S8, Li2S8, Li2S6, Li2S4, Li2S2
Target        : 1. Ground-State SCF  → E_total, fmax
                2. Band Gap & DOS    → gap_indirect, gap_direct
                3. Formation Energy  → E_form_eV_atom (vs 2D Graphene)
                4. Adsorption Energy → Eads per polysulfide
                5. Bulk Modulus (K) & Shear Modulus (G) via EOS + shear strain

Output Layout (relatif terhadap lokasi skrip ini):
    <script_dir>/results/  ← semua data numerik (.json, .dat, .cif, .xlsx, .csv)
    <script_dir>/figures/  ← semua grafik PNG (DOS, adsorpsi)
"""

import os
import json
import time
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

from ase.io import read, write
from ase.build import graphene
from ase.units import GPa
from ase.dft.bandgap import bandgap
from ase.dft.dos import DOS
from gpaw import GPAW, PW, FermiDirac

warnings.filterwarnings('ignore')

# =============================================================================
# DETEKSI DIREKTORI (otomatis: lokal & HPC)
# =============================================================================
SCRIPT_DIR  = Path(__file__).resolve().parent   # folder tempat skrip ini berada
RESULTS_DIR = SCRIPT_DIR / 'results'            # → dft/results/
FIG_DIR     = SCRIPT_DIR / 'figures'            # → dft/figures/
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Cari folder CIF pristine TPMS
# Prioritas: repo/graphene_tpms/ → symlink di HPC → path absolut HPC
_cif_cands = [
    SCRIPT_DIR.parent / 'graphene_tpms',          # repo lokal: cgcnn_data_multiproperty/graphene_tpms
    SCRIPT_DIR / 'graphene_tpms',                  # fallback: dft/graphene_tpms (jika ada symlink)
    Path('/media/user/uid1083/graphene_tpms'),      # path absolut HPC
    Path('/home/user/Amarus/cgcnn_data_multiproperty/graphene_tpms'),
]
CIF_DIR = next((p for p in _cif_cands if (p / 'graphene_sheet_neovius.cif').exists()), None)
if CIF_DIR is None:
    raise FileNotFoundError(
        "Folder graphene_tpms tidak ditemukan!\n"
        "Pastikan folder 'graphene_tpms' ada di root repo (satu level di atas dft/).\n"
        f"Script ini berada di: {SCRIPT_DIR}"
    )

# Cari folder CIF adsorbat polysulfida
_ads_cands = [
    SCRIPT_DIR.parent / 'cifs_graphene_tpms_adsorbate',
    Path('/media/user/uid1083/cifs_graphene_tpms_adsorbate'),
    Path('/home/user/Amarus/cgcnn_data_multiproperty/cifs_graphene_tpms_adsorbate'),
    Path('../cifs_graphene_tpms_adsorbate'),
    Path('cifs_graphene_tpms_adsorbate'),
]
ADS_DIR = next((p for p in _ads_cands if (p / 'graphene_neovius_S8.cif').exists()), None)
if ADS_DIR is None:
    raise FileNotFoundError("Folder cifs_graphene_tpms_adsorbate tidak ditemukan!")

print(f"📂 CIF_DIR  : {CIF_DIR}")
print(f"📂 ADS_DIR  : {ADS_DIR}")
print(f"📂 RESULTS  : {RESULTS_DIR}")
print(f"📂 FIGURES  : {FIG_DIR}")

# =============================================================================
# PARAMETER GLOBAL
# =============================================================================
POLYSULFIDES = ['S8', 'Li2S8', 'Li2S6', 'Li2S4', 'Li2S2']

STRUCTS = {
    'neovius':   'graphene_sheet_neovius.cif',
    'primitive': 'graphene_sheet_primitive.cif',
    'iwp':       'graphene_sheet_iwp.cif',
    'gyroid':    'graphene_sheet_gyroid.cif',
    'diamond':   'graphene_sheet_diamond.cif',
}
RUN = ['neovius', 'primitive', 'iwp', 'gyroid', 'diamond']

# Mode DFT: 'lcao' (cepat, ~20-60 menit/struktur) atau 'pw' (akurat, berjam-jam)
DFT_MODE     = os.environ.get('DFT_MODE', 'lcao').lower()
XC           = 'PBE'
ECUT         = 400          # eV, hanya dipakai di mode 'pw'
ELASTIC_DELTA = 0.02        # strain magnitude — 0.02 stabil untuk LCAO
FORCE        = False         # True = hitung ulang meski sudah ada

# =============================================================================
# MANAJEMEN FILE HASIL (results.json)
# =============================================================================
RES_FILE = RESULTS_DIR / 'results.json'
RES: dict = json.loads(RES_FILE.read_text()) if RES_FILE.exists() else {}

def save_res() -> None:
    """Simpan state RES ke results/results.json secara atomik."""
    tmp = RES_FILE.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(RES, indent=2))
    tmp.replace(RES_FILE)

def rec(name: str) -> dict:
    """Ambil/buat sub-dict hasil untuk struktur `name`."""
    return RES.setdefault(name, {})

# =============================================================================
# PEMBUATAN KALKULATOR GPAW
# =============================================================================
def make_calc(atoms, log_path: Path, kpts=None, molecule: bool = False):
    """
    Buat kalkulator GPAW.
    - kpts=None  → Gamma-point (1,1,1): cukup untuk sel TPMS besar >12Å
    - molecule=True → smearing lebih ketat (FermiDirac 0.01 eV)
    """
    if kpts is None:
        kpts = (1, 1, 1)

    common = dict(
        xc=XC,
        kpts=kpts,
        occupations=FermiDirac(0.01 if molecule else 0.05),
        txt=str(log_path),
    )
    if DFT_MODE == 'lcao':
        return GPAW(mode='lcao', basis='dzp',
                    convergence={'energy': 1e-4, 'density': 1e-3},
                    **common)
    else:
        return GPAW(mode=PW(ECUT),
                    convergence={'energy': 1e-5, 'density': 1e-4},
                    **common)

def load_cif(name: str):
    """Baca CIF pristine TPMS."""
    return read(str(CIF_DIR / STRUCTS[name]))

# =============================================================================
# LANGKAH 1: GROUND-STATE SCF
# =============================================================================
def scf_ground_state(name: str) -> None:
    r   = rec(name)
    gpw = RESULTS_DIR / f'{name}.gpw'

    if gpw.exists() and not FORCE and 'E_total' in r:
        print(f"✨ [{name}] SCF sudah ada — E = {r['E_total']:.4f} eV. Dilewati.")
        return

    atoms = load_cif(name)
    print(f"\n{'='*70}")
    print(f"🚀 LANGKAH 1: SCF {name.upper()} ({len(atoms)} atom C) | Mode: {DFT_MODE.upper()}")
    print(f"   Cell: {atoms.cell.lengths().round(3)} Å")
    print(f"{'='*70}")

    t0          = time.time()
    atoms.calc  = make_calc(atoms, RESULTS_DIR / f'{name}_scf.txt')
    E_tot       = float(atoms.get_potential_energy())
    F           = atoms.get_forces()
    fmax        = float(np.linalg.norm(F, axis=1).max())
    frms        = float(np.sqrt((F**2).mean()))

    # Stress (hanya tersedia di PW; LCAO tidak mendukung → fallback 0)
    try:
        stress_GPa = [float(x) for x in atoms.get_stress() / GPa]
        stress_max = float(max(abs(s) for s in stress_GPa))
    except Exception:
        stress_GPa = [0.0] * 6
        stress_max = 0.0

    elapsed = time.time() - t0
    r.update(
        n_atoms=len(atoms),
        cell_A=float(atoms.cell.lengths()[0]),
        E_total=E_tot,
        fmax=fmax,
        frms=frms,
        stress_GPa=stress_GPa,
        stress_max_GPa=stress_max,
        scf_time_s=round(elapsed, 1),
    )
    atoms.calc.write(str(gpw), mode='all')
    save_res()
    print(f"✅ [{name}] SCF selesai ({elapsed:.0f} s) | E = {E_tot:.4f} eV | fmax = {fmax:.3f} eV/Å")

# =============================================================================
# LANGKAH 2: BAND GAP & DOS
# =============================================================================
def gap_dos(name: str) -> None:
    r   = rec(name)
    gpw = RESULTS_DIR / f'{name}.gpw'

    if not gpw.exists():
        print(f"⚠️  [{name}] File .gpw belum ada — jalankan SCF lebih dulu.")
        return
    if 'gap_indirect' in r and not FORCE:
        print(f"✨ [{name}] Band gap & DOS sudah ada. Dilewati.")
        return

    print(f"⚡ [{name}] Menghitung Band Gap & DOS...")
    calc = GPAW(str(gpw), txt=None)
    ef   = float(calc.get_fermi_level())
    g,  _, _ = bandgap(calc, direct=False, efermi=ef)
    gd, _, _ = bandgap(calc, direct=True,  efermi=ef)

    # Simpan data DOS numerik → results/
    try:
        dos  = DOS(calc, width=0.1, npts=2000)
        e_   = dos.get_energies() - ef
        d_   = dos.get_dos()
        dat_path = RESULTS_DIR / f'{name}_dos.dat'
        np.savetxt(str(dat_path), np.c_[e_, d_], header='E-EF(eV)  DOS(states/eV)')

        # Plot DOS individu → figures/
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.fill_between(e_, d_, alpha=0.25, color='#1f77b4')
        ax.plot(e_, d_, color='#1f77b4', lw=1.5, label=f'Graphene {name.capitalize()}')
        ax.axvline(0, color='red', linestyle='--', lw=1.2, alpha=0.8, label=r'$E_F$')
        ax.set_xlim(-5, 5)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(r'$E - E_F$ (eV)', fontweight='bold')
        ax.set_ylabel('DOS (states/eV)', fontweight='bold')
        ax.set_title(f'Electronic DOS — Graphene {name.capitalize()} TPMS', fontweight='bold')
        ax.grid(True, linestyle=':', alpha=0.5)
        ax.legend()
        plt.tight_layout()
        plt.savefig(str(FIG_DIR / f'dos_{name}.png'), dpi=300)
        plt.close()
        print(f"   📊 DOS plot tersimpan → figures/dos_{name}.png")
    except Exception as exc:
        print(f"   ⚠️  DOS gagal [{name}]: {exc}")

    r.update(E_fermi=ef, gap_indirect=float(g), gap_direct=float(gd))
    save_res()
    print(f"✅ [{name}] E_F = {ef:.3f} eV | gap = {g:.4f} eV (indirect) | {gd:.4f} eV (direct)")

# =============================================================================
# LANGKAH 3: ENERGI FORMASI (vs flat graphene)
# =============================================================================
def graphene_ref() -> float:
    """Hitung atau ambil cache energi referensi 2D graphene flat."""
    if 'graphene_ref' in RES and not FORCE:
        return float(RES['graphene_ref']['E_per_atom'])

    print("📏 Menghitung referensi energi 2D Graphene (flat)...")
    g       = graphene(formula='C2', a=2.46, vacuum=8.0)
    g.pbc   = True
    g.calc  = make_calc(g, RESULTS_DIR / 'graphene_ref.txt', kpts=(6, 6, 1))
    E_atom  = float(g.get_potential_energy() / len(g))
    RES['graphene_ref'] = {'E_per_atom': E_atom}
    save_res()
    print(f"✅ E(graphene) = {E_atom:.4f} eV/atom")
    return E_atom

def calc_formation_energy(name: str) -> None:
    r    = rec(name)
    if 'E_form_eV_atom' in r and not FORCE:
        print(f"✨ [{name}] Energi formasi sudah ada. Dilewati.")
        return
    if 'E_total' not in r:
        print(f"⚠️  [{name}] E_total belum ada — lewati energi formasi.")
        return

    E_gr = graphene_ref()
    E_f  = float(r['E_total'] / r['n_atoms'] - E_gr)
    r['E_form_eV_atom'] = E_f
    save_res()
    print(f"✅ [{name}] E_form = {E_f:.4f} eV/atom")

# =============================================================================
# LANGKAH 4: ADSORPSI 5 POLISULFIDA
# =============================================================================
def free_molecule(mol_name: str):
    """
    Hitung atau ambil energi fasa gas polysulfida.
    Struktur diambil dari ADS_DIR (neovius complex), lalu atom non-C diekstrak.
    """
    tag      = f'gas_{mol_name}'
    xyz_path = RESULTS_DIR / f'{tag}.xyz'

    if tag in RES and xyz_path.exists() and not FORCE:
        return read(str(xyz_path)), float(RES[tag]['E_free'])

    ref_cif  = ADS_DIR / f'graphene_neovius_{mol_name}.cif'
    if not ref_cif.exists():
        raise FileNotFoundError(f"Tidak ada CIF referensi gas: {ref_cif}")

    cplx     = read(str(ref_cif))
    syms     = cplx.get_chemical_symbols()
    mol      = cplx[[i for i, s in enumerate(syms) if s != 'C']]
    mol.center(vacuum=10.0)
    mol.pbc  = True

    print(f"   🧪 [Gas SCF] {mol_name} ({mol.get_chemical_formula()})...")
    mol.calc  = make_calc(mol, RESULTS_DIR / f'{tag}_scf.txt', molecule=True)
    E_free    = float(mol.get_potential_energy())
    mol.calc  = None
    write(str(xyz_path), mol)
    RES[tag]  = {'E_free': E_free, 'formula': mol.get_chemical_formula()}
    save_res()
    print(f"   ✅ Gas {mol_name}: E = {E_free:.4f} eV")
    return mol, E_free

def adsorb_all_polysulfides(name: str) -> None:
    r       = rec(name)
    E_host  = r.get('E_total')
    if E_host is None:
        print(f"⚠️  [{name}] E_total host belum ada — lewati adsorpsi.")
        return

    print(f"\n🧲 [{name}] Adsorpsi 5 Polisulfida [Single-Point SCF]")
    eads_all = {}

    for mol in POLYSULFIDES:
        key     = f'Eads_{mol}'
        cif_out = RESULTS_DIR / f'{name}_{mol}_adsorbed.cif'

        # Skip jika sudah ada
        if key in r and cif_out.exists() and not FORCE:
            print(f"   ✨ [{name}+{mol}] sudah ada → Eads = {r[key]:.3f} eV. Dilewati.")
            eads_all[mol] = r[key]
            continue

        cplx_cif = ADS_DIR / f'graphene_{name}_{mol}.cif'
        if not cplx_cif.exists():
            print(f"   ⚠️  CIF kompleks tidak ada: {cplx_cif.name}")
            continue

        try:
            _, E_gas = free_molecule(mol)
            cplx     = read(str(cplx_cif))
            t0       = time.time()
            cplx.calc = make_calc(cplx, RESULTS_DIR / f'{name}_{mol}_scf.txt')
            E_cplx   = float(cplx.get_potential_energy())
            E_ads    = float(E_cplx - E_host - E_gas)

            r[key]          = E_ads
            eads_all[mol]   = E_ads
            write(str(cif_out), cplx)
            save_res()
            print(f"   ✅ [{name}+{mol}] Eads = {E_ads:.3f} eV ({time.time()-t0:.0f} s)")
        except Exception as exc:
            print(f"   ❌ [{name}+{mol}] Error: {exc}")

    # Hitung statistik Eads dari semua 5 polisulfida yang tersedia
    all_vals = [r[f'Eads_{m}'] for m in POLYSULFIDES if f'Eads_{m}' in r]
    if all_vals:
        r['Eads_mean'] = round(float(np.mean(all_vals)), 4)
        r['Eads_best'] = round(float(min(all_vals)), 4)
        save_res()

# =============================================================================
# LANGKAH 5: BULK MODULUS (K) & SHEAR MODULUS (G)
# =============================================================================
def calc_elastic(name: str) -> None:
    r = rec(name)
    if 'K_VRH' in r and not FORCE:
        print(f"✨ [{name}] Modulus elastisitas sudah ada. Dilewati.")
        return
    if 'E_total' not in r:
        print(f"⚠️  [{name}] E_total belum ada — lewati elastisitas.")
        return

    base = load_cif(name)
    v0   = float(base.get_volume())
    e0   = float(r['E_total'])

    # ─── METODE LCAO: EOS (Birch-Murnaghan) + Shear Strain Energy ───
    if DFT_MODE == 'lcao':
        print(f"\n🔩 [{name}] Menghitung K (EOS) & G (Shear Strain) — LCAO...")

        # ── Bulk Modulus K via 5-point EOS ──
        # Skala volume: 0.98, 0.99, 1.00, 1.01, 1.02
        scale_factors = [s**(1/3) for s in [0.98, 0.99, 1.00, 1.01, 1.02]]
        vols, energies = [], []
        for sf in scale_factors:
            v_s = float(v0 * sf**3)
            vols.append(v_s)
            if abs(sf - 1.0) < 1e-5:
                energies.append(e0)
            else:
                at = base.copy()
                at.set_cell(at.cell * sf, scale_atoms=True)
                at.calc = make_calc(at, RESULTS_DIR / f'{name}_eos_sf{sf**3:.3f}.txt')
                energies.append(float(at.get_potential_energy()))

        # K = V0 * d²E/dV² (konversi eV/Å³ → GPa: × 160.2176)
        coeffs   = np.polyfit(vols, energies, 2)
        K_gpa    = float(v0 * 2.0 * coeffs[0] * 160.2176)

        # ── Shear Modulus G via Pure Shear Strain ──
        gammas   = [-0.015, -0.008, 0.0, 0.008, 0.015]
        e_shears = []
        for g in gammas:
            if abs(g) < 1e-8:
                e_shears.append(0.0)
            else:
                at = base.copy()
                sm = np.array([[1.0,  g,   0.0],
                               [g,    1.0, 0.0],
                               [0.0,  0.0, 1.0 / (1.0 - g**2)]])
                at.set_cell(base.cell.array @ sm.T, scale_atoms=True)
                at.calc = make_calc(at, RESULTS_DIR / f'{name}_shear_g{g:+.3f}.txt')
                e_shears.append(float(at.get_potential_energy()) - e0)

        poly_g = np.polyfit(np.array(gammas)**2, e_shears, 1)
        G_gpa  = float(poly_g[0] / (0.5 * v0) * 160.2176)

        denom   = 3 * K_gpa + G_gpa
        E_young = float(9 * K_gpa * G_gpa / denom) if denom > 0 else 0.0
        poisson = float((3 * K_gpa - 2 * G_gpa) / (2 * denom)) if denom > 0 else 0.0

        r.update(K_VRH=round(K_gpa, 2), G_VRH=round(G_gpa, 2),
                 E_young=round(E_young, 2), poisson=round(poisson, 4))
        save_res()
        print(f"✅ [{name}] K = {K_gpa:.1f} GPa | G = {G_gpa:.1f} GPa | E = {E_young:.1f} GPa | ν = {poisson:.3f}")
        return

    # ─── METODE PW: Tensor Elastisitas C_ij via Stress ───
    print(f"\n🔩 [{name}] Menghitung Tensor Elastisitas C_ij — PW...")

    def strain_mat(j, d):
        e = np.zeros((3, 3))
        if j < 3:
            e[j, j] = d
        else:
            a, b = {3: (1, 2), 4: (0, 2), 5: (0, 1)}[j]
            e[a, b] = e[b, a] = d / 2
        return e

    C = np.zeros((6, 6))
    for j in range(6):
        stresses = []
        for sgn in (+1, -1):
            at = base.copy()
            at.set_cell(at.cell.array @ (np.eye(3) + strain_mat(j, sgn * ELASTIC_DELTA)).T,
                        scale_atoms=True)
            at.calc = make_calc(at, RESULTS_DIR / f'{name}_el{j}_{sgn:+d}.txt')
            stresses.append(at.get_stress() / GPa)
        C[:, j] = (stresses[0] - stresses[1]) / (2 * ELASTIC_DELTA)

    C = 0.5 * (C + C.T)
    S = np.linalg.inv(C)
    Kv = (C[0,0]+C[1,1]+C[2,2] + 2*(C[0,1]+C[0,2]+C[1,2])) / 9
    Gv = (C[0,0]+C[1,1]+C[2,2] - (C[0,1]+C[0,2]+C[1,2]) + 3*(C[3,3]+C[4,4]+C[5,5])) / 15
    Kr = 1.0 / (S[0,0]+S[1,1]+S[2,2] + 2*(S[0,1]+S[0,2]+S[1,2]))
    Gr = 15.0 / (4*(S[0,0]+S[1,1]+S[2,2]) - 4*(S[0,1]+S[0,2]+S[1,2]) + 3*(S[3,3]+S[4,4]+S[5,5]))
    K  = float((Kv + Kr) / 2)
    G  = float((Gv + Gr) / 2)

    denom   = 3 * K + G
    E_young = float(9 * K * G / denom) if denom > 0 else 0.0
    poisson = float((3 * K - 2 * G) / (2 * denom)) if denom > 0 else 0.0

    r.update(K_VRH=round(K, 2), G_VRH=round(G, 2),
             E_young=round(E_young, 2), poisson=round(poisson, 4),
             C_min_eig=round(float(np.linalg.eigvalsh(C).min()), 2),
             C_GPa=C.round(2).tolist())
    save_res()
    np.savetxt(str(RESULTS_DIR / f'{name}_Cij_GPa.dat'), C, fmt='%9.2f')
    print(f"✅ [{name}] K = {K:.1f} GPa | G = {G:.1f} GPa | E = {E_young:.1f} GPa | ν = {poisson:.3f}")

# =============================================================================
# LANGKAH 6: RINGKASAN, EXCEL, & GRAFIK
# =============================================================================
def summarize_and_plot() -> None:
    """Buat ringkasan CSV, laporan Excel multi-sheet, dan grafik DOS + adsorpsi."""

    # ── Bangun DataFrame ringkasan ──
    tpms_data = {k: v for k, v in RES.items() if k in STRUCTS}
    if not tpms_data:
        print("⚠️  Belum ada data TPMS yang selesai.")
        return

    df = pd.DataFrame(tpms_data).T
    cols_phys = [c for c in ['n_atoms', 'fmax', 'gap_indirect', 'gap_direct',
                              'E_form_eV_atom', 'K_VRH', 'G_VRH', 'E_young', 'poisson']
                 if c in df.columns]
    cols_ads  = [c for c in [f'Eads_{m}' for m in POLYSULFIDES] + ['Eads_mean', 'Eads_best']
                 if c in df.columns]

    avail_cols = [c for c in cols_phys + cols_ads if c in df.columns]
    summary    = df[avail_cols].apply(pd.to_numeric, errors='coerce').round(4)
    summary.index.name = 'struktur_tpms'

    # ── Simpan CSV ringkasan → results/ ──
    csv_path = RESULTS_DIR / 'summary.csv'
    summary.to_csv(csv_path)
    print(f"\n📋 Ringkasan CSV tersimpan → results/summary.csv")
    print(summary.to_string())

    # ── Ekspor Excel multi-sheet → results/ ──
    try:
        excel_path  = RESULTS_DIR / 'hasil_pure_dft_tpms_polysulfide.xlsx'
        rename_phys = {
            'n_atoms':         'Jumlah Atom (C)',
            'fmax':            'Gaya Maks (eV/Å)',
            'gap_indirect':    'Band Gap Indirect (eV)',
            'gap_direct':      'Band Gap Direct (eV)',
            'E_form_eV_atom':  'Energi Formasi Ef (eV/atom)',
            'K_VRH':           'Bulk Modulus K (GPa)',
            'G_VRH':           'Shear Modulus G (GPa)',
            'E_young':         'Young Modulus E (GPa)',
            'poisson':         'Poisson Ratio',
            'Eads_mean':       'E_ads Rerata (eV)',
            'Eads_best':       'E_ads Terkuat (eV)',
        }
        rename_ads  = {f'Eads_{p}': f'E_ads {p} (eV)' for p in POLYSULFIDES}

        df_main = summary[[c for c in cols_phys + ['Eads_mean', 'Eads_best']
                           if c in summary.columns]].rename(columns=rename_phys)
        df_main.index.name = 'Struktur TPMS'

        ads_only = [c for c in [f'Eads_{p}' for p in POLYSULFIDES] if c in summary.columns]
        df_ads   = summary[ads_only].rename(columns=rename_ads)
        df_ads.index.name = 'Struktur TPMS'

        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            df_main.to_excel(writer, sheet_name='5 Sifat Utama')
            if not df_ads.empty:
                df_ads.to_excel(writer, sheet_name='Adsorpsi 5 Polisulfida')
        print(f"📊 Laporan Excel tersimpan → results/hasil_pure_dft_tpms_polysulfide.xlsx")
    except Exception as exc:
        print(f"   ⚠️  Gagal ekspor Excel: {exc}")

    # ── Grafik Profil Adsorpsi → figures/ ──
    ads_cols_avail = [c for c in [f'Eads_{m}' for m in POLYSULFIDES] if c in summary.columns]
    if ads_cols_avail:
        fig, ax = plt.subplots(figsize=(9, 5))
        lbls    = [c.replace('Eads_', '') for c in ads_cols_avail]
        colors  = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        for (idx, row), color in zip(summary.iterrows(), colors):
            vals = [row[c] for c in ads_cols_avail]
            ax.plot(lbls, vals, marker='o', linewidth=2.0, label=str(idx).upper(), color=color)
        ax.axhline(0, color='gray', linestyle='--', alpha=0.6)
        ax.set_xlabel('Spesies Polisulfida (Siklus Litiasi)', fontweight='bold')
        ax.set_ylabel('Energi Adsorpsi $E_{ads}$ (eV)', fontweight='bold')
        ax.set_title('Profil Afinitas Penjeratan Polisulfida pada Graphene TPMS', fontweight='bold')
        ax.grid(True, linestyle=':', alpha=0.5)
        ax.legend(loc='lower right')
        plt.tight_layout()
        out = FIG_DIR / 'adsorption_profile_5species.png'
        plt.savefig(str(out), dpi=300)
        plt.close()
        print(f"📊 Grafik profil adsorpsi → figures/adsorption_profile_5species.png")

    # ── Grafik DOS Komposit (semua struktur) → figures/ ──
    dos_files = [(nm, RESULTS_DIR / f'{nm}_dos.dat') for nm in RUN
                 if (RESULTS_DIR / f'{nm}_dos.dat').exists()]
    if dos_files:
        n_panel = len(dos_files)
        fig, axes = plt.subplots(n_panel, 1, figsize=(7, 2.5 * n_panel), sharex=True)
        if n_panel == 1:
            axes = [axes]
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        for ax, (nm, dat_path), color in zip(axes, dos_files, colors):
            data = np.loadtxt(str(dat_path))
            e_, d_ = data[:, 0], data[:, 1]
            ax.fill_between(e_, d_, alpha=0.2, color=color)
            ax.plot(e_, d_, lw=1.5, color=color, label=f'Graphene {nm.capitalize()}')
            ax.axvline(0, color='red', linestyle='--', lw=1.0, alpha=0.8)
            ax.set_xlim(-5, 5)
            ax.set_ylim(bottom=0)
            ax.set_ylabel('DOS (st/eV)', fontsize=9)
            ax.legend(loc='upper right', fontsize=9)
            ax.grid(True, linestyle=':', alpha=0.4)
        axes[-1].set_xlabel(r'$E - E_F$ (eV)', fontweight='bold')
        plt.suptitle('Density of States — 5 Graphene TPMS', fontweight='bold', y=1.01)
        plt.tight_layout()
        out = FIG_DIR / 'dos_all.png'
        plt.savefig(str(out), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 Grafik DOS komposit → figures/dos_all.png")

# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    print("=" * 80)
    print("🚀 MASTER DFT GPAW RUNNER [SINGLE-POINT — TANPA RELAKSASI]")
    print(f"   Struktur  : {RUN}")
    print(f"   Mode DFT  : {DFT_MODE.upper()}")
    print(f"   Results   : {RESULTS_DIR}")
    print(f"   Figures   : {FIG_DIR}")
    print("=" * 80)

    for nm in RUN:
        print(f"\n{'>'*5} MEMPROSES: {nm.upper()} {'<'*5}")
        scf_ground_state(nm)
        gap_dos(nm)
        calc_formation_energy(nm)
        adsorb_all_polysulfides(nm)
        calc_elastic(nm)

    summarize_and_plot()
    print("\n🏁 SEMUA KALKULASI DFT 5 TPMS & 5 POLISULFIDA SELESAI!")
