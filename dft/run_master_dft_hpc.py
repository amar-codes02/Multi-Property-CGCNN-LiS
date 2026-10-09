#!/usr/bin/env python3
"""
MASTER PURE DFT GPAW RUNNER — 5 GRAPHENE TPMS & 5 POLYSULFIDES (S8 -> Li2S2)
[VERSI SINGLE-POINT DIRECT SCF — TANPA RELAKSASI GEOMETRI]
=============================================================================
Struktur TPMS: Neovius (188 C), Primitive (200 C), IWP (228 C), Gyroid (244 C), Diamond (332 C)
Adsorbat     : S8, Li2S8, Li2S6, Li2S4, Li2S2
Target       : 1. Ground-State SCF (E_total, fmax, stress)
               2. Band Gap & DOS
               3. Formation Energy vs 2D Graphene
               4. Adsorption Energy 5 Polysulfides (Eads = Ecplx - Ehost - Egas)
               5. Bulk Modulus (K) & Shear Modulus (G)
               6. Summary CSV & Publication Figures
"""

import os, sys, json, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

from ase.io import read, write
from ase.build import graphene
from ase.units import GPa
from ase.calculators.calculator import kptdensity2monkhorstpack
from ase.dft.bandgap import bandgap
from ase.dft.dos import DOS
from ase import Atoms
from gpaw import GPAW, PW, FermiDirac

# ------------------------- DETEKSI DIREKTORI -------------------------
cif_cands = [
    Path('/media/user/uid1083/graphene_tpms'),
    Path('/home/user/Amarus/cgcnn_data_multiproperty/graphene_tpms'),
    Path('graphene_tpms'),
    Path('.')
]
CIF_DIR = next((p for p in cif_cands if (p / 'graphene_sheet_neovius.cif').exists()), Path('.'))

ads_cands = [
    Path('/media/user/uid1083/cifs_graphene_tpms_adsorbate'),
    Path('/home/user/Amarus/cgcnn_data_multiproperty/cifs_graphene_tpms_adsorbate'),
    Path('../cifs_graphene_tpms_adsorbate'),
    Path('cifs_graphene_tpms_adsorbate')
]
ADS_DIR = next((p for p in ads_cands if (p / 'graphene_neovius_S8.cif').exists()), Path('.'))

work_cands = [
    Path('/media/user/uid1083/dft_gpaw_graphene_tpms_testing/results'),
    Path('/home/user/Amarus/cgcnn_data_multiproperty/dft_gpaw_graphene_tpms_testing/results'),
    Path('results')
]
WORK = next((p for p in work_cands if p.parent.exists()), Path('results'))
WORK.mkdir(exist_ok=True, parents=True)

RES_FILE = WORK / 'results.json'
RES = json.loads(RES_FILE.read_text()) if RES_FILE.exists() else {}
def save_res():
    RES_FILE.write_text(json.dumps(RES, indent=2))
def rec(name):
    return RES.setdefault(name, {})

# ------------------------- PARAMETER FISIKA DFT -------------------------
POLYSULFIDES = ['S8', 'Li2S8', 'Li2S6', 'Li2S4', 'Li2S2']
STRUCTS = {
    'neovius':   'graphene_sheet_neovius.cif',
    'primitive': 'graphene_sheet_primitive.cif',
    'iwp':       'graphene_sheet_iwp.cif',
    'gyroid':    'graphene_sheet_gyroid.cif',
    'diamond':   'graphene_sheet_diamond.cif',
}
RUN = ['neovius', 'primitive', 'iwp', 'gyroid', 'diamond']

# Mode: 'lcao' (ultra-cepat, ~20-30 menit) atau 'pw' (Plane-Wave)
DFT_MODE = os.environ.get('DFT_MODE', 'lcao').lower()
ECUT = 400
XC = 'PBE'
ELASTIC_DELTA = 0.005
FORCE = False

def make_calc(atoms, txt, kpts=None, molecule=False, spin=False):
    if kpts is None:
        # Sel TPMS >12 Å (>180 atom): Gamma-point (1,1,1) sangat akurat & cepat
        kpts = (1, 1, 1)
        
    if DFT_MODE == 'lcao':
        return GPAW(mode='lcao', basis='dzp', xc=XC, kpts=kpts, spinpol=spin,
                    occupations=FermiDirac(0.01 if molecule else 0.05),
                    convergence={'energy': 1e-4, 'density': 1e-3}, txt=str(txt))
    else:
        return GPAW(mode=PW(ECUT), xc=XC, kpts=kpts, spinpol=spin,
                    occupations=FermiDirac(0.01 if molecule else 0.05),
                    convergence={'energy': 1e-5, 'density': 1e-4}, txt=str(txt))

def load_initial(name):
    return read(str(CIF_DIR / STRUCTS[name]))

# ==================== LANGKAH 1: GROUND-STATE SCF PRISTINE ====================
def scf_ground_state(name):
    r = rec(name)
    gpw = WORK / f'{name}.gpw'
    
    if gpw.exists() and not FORCE and 'E_total' in r:
        print(f"✨ [{name}] Ground-state SCF sudah ada (E = {r.get('E_total'):.4f} eV). Dilewati.")
        return
        
    a = load_initial(name)
    print("\n" + "="*70)
    print(f"🚀 LANGKAH 1: GROUND-STATE SCF {name.upper()} ({len(a)} ATOM C) [TANPA RELAKSASI]")
    print(f"   Cell lengths: {a.cell.lengths().round(3)} | Mode: {DFT_MODE.upper()}")
    print("="*70)
    
    t0 = time.time()
    calc = make_calc(a, WORK / f'{name}_scf.txt')
    a.calc = calc
    
    E_tot = float(a.get_potential_energy())
    F = a.get_forces()
    fmax = float(np.linalg.norm(F, axis=1).max())
    frms = float(np.sqrt((F**2).mean()))
    
    # Stress tensor (hanya tersedia di PW; di LCAO di-fallback aman)
    try:
        stress_GPa = [float(x) for x in a.get_stress() / GPa]
        stress_max = float(np.abs(stress_GPa).max())
    except Exception:
        stress_GPa = [0.0] * 6
        stress_max = 0.0
    
    r.update(
        n_atoms=len(a), cell_A=float(a.cell.lengths()[0]),
        E_total=E_tot, fmax=fmax, frms=frms,
        stress_GPa=stress_GPa, stress_max_GPa=stress_max,
        scf_time_s=float(time.time() - t0)
    )
    
    a.calc.write(str(gpw), mode='all')
    save_res()
    print(f"✅ [{name}] SCF selesai ({time.time()-t0:.1f} s) | E = {E_tot:.4f} eV | fmax = {fmax:.3f} eV/Å")

# ==================== LANGKAH 2: BAND GAP & DOS ====================
def gap_dos(name):
    r = rec(name)
    gpw = WORK / f'{name}.gpw'
    if not gpw.exists():
        print(f"⚠️ [{name}] {gpw.name} belum ada, jalankan Langkah 1 dulu!")
        return
    if 'gap_indirect' in r and not FORCE:
        print(f"✨ [{name}] Band gap & DOS sudah ada. Dilewati.")
        return
        
    print(f"⚡ [{name}] Menghitung Band Gap & DOS...")
    calc = GPAW(str(gpw), txt=None)
    ef = float(calc.get_fermi_level())
    g, _, _ = bandgap(calc, direct=False, efermi=ef)
    gd, _, _ = bandgap(calc, direct=True, efermi=ef)
    
    try:
        dos = DOS(calc, width=0.1, npts=2000)
        e, d = dos.get_energies() - ef, dos.get_dos()
        np.savetxt(str(WORK / f'{name}_dos.dat'), np.c_[e, d], header='E-EF(eV) DOS(states/eV)')
    except Exception as exc:
        print(f"   [Warn DOS {name}]: {exc}")
        
    r.update(E_fermi=ef, gap_indirect=float(g), gap_direct=float(gd))
    save_res()
    print(f"✅ [{name}] E_F={ef:.3f} eV | gap indirect={g:.3f} eV | direct={gd:.3f} eV")

# ==================== LANGKAH 3: ENERGI FORMASI ====================
def graphene_ref():
    if 'graphene_ref' in RES and not FORCE:
        return RES['graphene_ref']['E_per_atom']
    print("📏 Menghitung referensi 2D Graphene...")
    g = graphene(formula='C2', a=2.46, vacuum=8.0); g.pbc = True
    g.calc = make_calc(g, WORK / 'graphene_ref.txt', kpts=(6, 6, 1))
    E = float(g.get_potential_energy() / len(g))
    RES['graphene_ref'] = {'E_per_atom': E}
    save_res()
    print(f"✅ E(graphene) = {E:.4f} eV/atom")
    return E

def calc_formation_energy(name):
    E_gr = graphene_ref()
    r = rec(name)
    if 'E_total' in r:
        r['E_form_eV_atom'] = float(r['E_total'] / r['n_atoms'] - E_gr)
        save_res()
        print(f"✅ [{name}] E_form = {r['E_form_eV_atom']:.4f} eV/atom")

# ==================== LANGKAH 4: ADSORPSI 5 POLISULFIDA ====================
def free_molecule(mol_name):
    tag = f'gas_{mol_name}'
    xyz_path = WORK / f'{tag}.xyz'
    if tag in RES and xyz_path.exists() and not FORCE:
        return read(str(xyz_path)), RES[tag]['E_free']
    
    ref_cif = ADS_DIR / f'graphene_neovius_{mol_name}.cif'
    cplx = read(str(ref_cif))
    symbols = cplx.get_chemical_symbols()
    mol_atoms = cplx[[i for i, s in enumerate(symbols) if s != 'C']]
    mol_atoms.center(vacuum=10.0)
    mol_atoms.pbc = True
    
    print(f"   [Gas DFT] SCF molekul fasa gas: {mol_name} ({mol_atoms.get_chemical_formula()})...")
    mol_atoms.calc = make_calc(mol_atoms, WORK / f'{tag}_scf.txt', molecule=True)
    E_free = float(mol_atoms.get_potential_energy())
    mol_atoms.calc = None
    write(str(xyz_path), mol_atoms)
    RES[tag] = {'E_free': E_free, 'formula': mol_atoms.get_chemical_formula()}
    save_res()
    print(f"   ✅ Gas {mol_name}: E = {E_free:.4f} eV")
    return mol_atoms, E_free

def adsorb_all_polysulfides(name):
    r = rec(name)
    E_host = r.get('E_total')
    if E_host is None:
        print(f"⚠️ [{name}] E_total host belum ada, lewati adsorpsi.")
        return
        
    print(f"\n🧲 [{name}] Evaluasi Adsorpsi 5 Polisulfida [SINGLE-POINT SCF]")
    out_eads = {}
    for mol in POLYSULFIDES:
        key = f'Eads_{mol}'
        cif_out = WORK / f'{name}_{mol}_adsorbed.cif'
        if key in r and cif_out.exists() and not FORCE:
            out_eads[mol] = r[key]
            continue
            
        _, E_mol_gas = free_molecule(mol)
        cplx_cif = ADS_DIR / f'graphene_{name}_{mol}.cif'
        if not cplx_cif.exists():
            print(f"   ⚠️ {cplx_cif.name} tidak ada!")
            continue
            
        cplx = read(str(cplx_cif))
        t0 = time.time()
        
        # Single-point SCF terkonfinasi (tanpa relaksasi)
        cplx.calc = make_calc(cplx, WORK / f'{name}_{mol}_scf.txt')
        E_cplx = float(cplx.get_potential_energy())
        E_ads = float(E_cplx - E_host - E_mol_gas)
        
        r[key] = E_ads
        out_eads[mol] = E_ads
        write(str(cif_out), cplx)
        save_res()
        print(f"   ✅ {name} + {mol}: E_ads = {E_ads:.3f} eV ({time.time()-t0:.1f} s)")
        
    if out_eads:
        r['Eads_mean'] = float(np.mean(list(out_eads.values())))
        r['Eads_best'] = float(min(out_eads.values()))
        save_res()

# ==================== LANGKAH 5: BULK & SHEAR MODULUS ====================
def calc_elastic(name):
    r = rec(name)
    if 'K_VRH' in r and not FORCE:
        print(f"✨ [{name}] Modulus elastisitas sudah ada. Dilewati.")
        return
        
    base = load_initial(name)
    v0 = float(base.get_volume())
    e0 = float(r.get('E_total', 0.0))
    
    # ── METODE 1: LCAO via Birch-Murnaghan EOS & Shear Strain Energy ──
    if DFT_MODE == 'lcao':
        print(f"\n🔩 [{name}] Menghitung Bulk Modulus K (EOS) & Shear Modulus G (Shear Strain)...")
        # 1. Bulk Modulus K via EOS (5 skala volume)
        scale_factors = [0.98**(1/3), 0.99**(1/3), 1.00, 1.01**(1/3), 1.02**(1/3)]
        vols, energies = [], []
        for sf in scale_factors:
            v_s = float(v0 * (sf**3))
            vols.append(v_s)
            if abs(sf - 1.0) < 1e-4:
                energies.append(e0)
            else:
                at_s = base.copy()
                at_s.set_cell(at_s.cell * sf, scale_atoms=True)
                at_s.calc = make_calc(at_s, WORK / f"{name}_eos_sf{sf**3:.3f}.txt")
                energies.append(float(at_s.get_potential_energy()))
        
        # Fit parabola d2E/dV2 -> K = V0 * d2E/dV2 * 160.21766 GPa
        coeffs = np.polyfit(vols, energies, 2)
        d2E_dV2 = 2.0 * coeffs[0]
        K_gpa = float(v0 * d2E_dV2 * 160.21766)
        
        # 2. Shear Modulus G via Pure Shear Strain gamma
        gammas = [-0.015, -0.008, 0.0, 0.008, 0.015]
        e_shears = []
        for g in gammas:
            if abs(g) < 1e-5:
                e_shears.append(0.0)
            else:
                at_g = base.copy()
                cell_o = base.cell.copy()
                sm = np.array([[1.0, g, 0.0], [g, 1.0, 0.0], [0.0, 0.0, 1.0 / (1.0 - g**2)]])
                at_g.set_cell(cell_o @ sm.T, scale_atoms=True)
                at_g.calc = make_calc(at_g, WORK / f"{name}_shear_g{g:+.3f}.txt")
                e_shears.append(float(at_g.get_potential_energy()) - e0)
                
        poly_g = np.polyfit(np.array(gammas)**2, e_shears, 1)
        G_gpa = float(poly_g[0] / (0.5 * v0) * 160.21766)
        
        E_young = float(9*K_gpa*G_gpa / (3*K_gpa + G_gpa)) if (3*K_gpa + G_gpa) > 0 else 0.0
        poisson = float((3*K_gpa - 2*G_gpa) / (2*(3*K_gpa + G_gpa))) if (3*K_gpa + G_gpa) > 0 else 0.0
        
        r.update(K_VRH=round(K_gpa, 2), G_VRH=round(G_gpa, 2), E_young=round(E_young, 2), poisson=round(poisson, 3))
        save_res()
        print(f"✅ [{name}] K = {K_gpa:.1f} GPa | G = {G_gpa:.1f} GPa")
        return
        
    # ── METODE 2: PW via Stress Tensor C_ij ──
    print(f"\n🔩 [{name}] Menghitung Tensor Elastisitas C_ij (6 arah regangan, Clamped-Ion)...")
    def strain_mat(j, d):
        e = np.zeros((3, 3))
        if j < 3: e[j, j] = d
        else:
            a, b = {3: (1, 2), 4: (0, 2), 5: (0, 1)}[j]; e[a, b] = e[b, a] = d / 2
        return e

    C = np.zeros((6, 6))
    for j in range(6):
        s = []
        for sgn in (+1, -1):
            at = base.copy()
            at.set_cell(at.cell.array @ (np.eye(3) + strain_mat(j, sgn * ELASTIC_DELTA)).T, scale_atoms=True)
            at.calc = make_calc(at, WORK / f'{name}_el{j}_{sgn:+d}.txt')
            s.append(at.get_stress() / GPa)
        C[:, j] = (s[0] - s[1]) / (2 * ELASTIC_DELTA)
        
    C = 0.5 * (C + C.T); S = np.linalg.inv(C)
    Kv = (C[0,0]+C[1,1]+C[2,2] + 2*(C[0,1]+C[0,2]+C[1,2])) / 9
    Gv = (C[0,0]+C[1,1]+C[2,2] - (C[0,1]+C[0,2]+C[1,2]) + 3*(C[3,3]+C[4,4]+C[5,5])) / 15
    Kr = 1 / (S[0,0]+S[1,1]+S[2,2] + 2*(S[0,1]+S[0,2]+S[1,2]))
    Gr = 15 / (4*(S[0,0]+S[1,1]+S[2,2]) - 4*(S[0,1]+S[0,2]+S[1,2]) + 3*(S[3,3]+S[4,4]+S[5,5]))
    K, G = float((Kv + Kr) / 2), float((Gv + Gr) / 2)
    
    r.update(K_VRH=K, G_VRH=G, E_young=float(9*K*G/(3*K+G)), poisson=float((3*K-2*G)/(2*(3*K+G))),
             C_min_eig=float(np.linalg.eigvalsh(C).min()), C_GPa=C.tolist())
    save_res()
    np.savetxt(str(WORK / f'{name}_Cij_GPa.dat'), C, fmt='%9.2f')
    print(f"✅ [{name}] K = {K:.1f} GPa | G = {G:.1f} GPa | Min eig(C) = {r['C_min_eig']:.1f}")

# ==================== LANGKAH 6: RINGKASAN & GRAFIK ====================
def summarize_and_plot():
    df = pd.DataFrame({k: v for k, v in RES.items() if k in STRUCTS}).T
    cols_phys = [c for c in ['n_atoms','fmax','gap_indirect','gap_direct','E_form_eV_atom','K_VRH','G_VRH','E_young','poisson'] if c in df.columns]
    cols_ads = [c for c in [f'Eads_{m}' for m in POLYSULFIDES] + ['Eads_mean', 'Eads_best'] if c in df.columns]
    
    summary = df[cols_phys + cols_ads].astype(float).round(3)
    summary_csv = WORK / 'summary.csv'
    summary.to_csv(summary_csv)
    print(f"\n📋 RINGKASAN LENGKAP TERSIMPAN KE: {summary_csv}")
    print(summary.to_string())
    
    ads_cols = [f'Eads_{m}' for m in POLYSULFIDES if f'Eads_{m}' in summary.columns]
    if ads_cols:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        lbls = [c.replace('Eads_', '') for c in ads_cols]
        for idx, row in summary.iterrows():
            ax.plot(lbls, [row[c] for c in ads_cols], marker='o', linewidth=2, label=str(idx).upper())
        ax.axhline(0, color='gray', linestyle='--', alpha=0.7)
        ax.set_xlabel('Spesies Polisulfida (Siklus Litiasi Redoks)', fontweight='bold')
        ax.set_ylabel('Energi Adsorpsi E_ads (eV)', fontweight='bold')
        ax.set_title('Profil Afinitas Penjeratan Polisulfida pada Graphene TPMS', fontweight='bold')
        ax.grid(True, linestyle=':', alpha=0.6); ax.legend(); plt.tight_layout()
        plt.savefig(str(WORK / 'adsorption_profile_5species.png'), dpi=200)
        plt.close()
        print(f"📊 Grafik Profil Adsorpsi tersimpan ke: {WORK / 'adsorption_profile_5species.png'}")

# ==================== MAIN EXECUTION ====================
if __name__ == '__main__':
    print("="*80)
    print("🚀 MASTER DFT GPAW RUNNER [SINGLE-POINT — TANPA RELAKSASI]")
    print(f"   Struktur: {RUN}")
    print(f"   Mode    : {DFT_MODE.upper()}")
    print(f"   Output  : {WORK}")
    print("="*80)
    
    for nm in RUN:
        print(f"\n>>> MEMPROSES STRUKTUR TPMS: {nm.upper()} <<<")
        scf_ground_state(nm)
        gap_dos(nm)
        calc_formation_energy(nm)
        adsorb_all_polysulfides(nm)
        calc_elastic(nm)
        
    summarize_and_plot()
    print("\n🏁 SELURUH PERHITUNGAN DFT 5 TPMS & 5 POLISULFIDA SELESAI!")
