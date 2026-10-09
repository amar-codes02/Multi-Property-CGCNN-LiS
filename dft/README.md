# Density Functional Theory (DFT) — Graphene TPMS & Polysulfide Adsorption

Folder ini berisi seluruh **kode simulasi mekanika kuantum DFT GPAW**, skrip eksekusi terminal HPC, dan **hasil kalkulasi ground truth** untuk 5 arsitektur Graphene TPMS (*Diamond, Gyroid, IWP, Neovius, Primitive*) serta 5 spesies litium polisulfida ($S_8 \to Li_2S_2$).

---

## 📁 Struktur Direktori `dft/`

```text
dft/
├── run_master_dft_hpc.py      # Master script Python yang berjalan di terminal HPC (Single-Point LCAO DZP)
├── tpms_gpaw.ipynb            # Notebook Jupyter interaktif DFT
├── results.json               # Data hasil kuantum lengkap dalam format JSON
├── export_to_excel.py         # Skrip otomatis ekspor hasil ke Excel multi-sheet (.xlsx)
├── sync_hpc_results.sh        # Skrip otomatis 1-klik untuk menyinkronkan hasil dari HPC
├── figures/                   # Grafik publikasi kurva DOS (Density of States) per TPMS (300 DPI)
│   ├── dos_diamond.png
│   ├── dos_gyroid.png
│   ├── dos_iwp.png
│   ├── dos_neovius.png
│   └── dos_primitive.png
└── results/                   # Berkas data mentah hasil perhitungan DFT
    ├── results.json
    ├── hasil_pure_dft_tpms_polysulfide.xlsx
    ├── *_dos.dat              # Data numerik DOS per energi
    ├── *_scf.txt              # Log konvergensi SCF
    ├── *_adsorbed.cif         # Struktur koordinat kompleks adsorpsi teroptimasi
    └── *_eos_*.txt            # Data regangan volume (Bulk Modulus) & geser (Shear Modulus)
```

---

## ⚡ Cara Menjalankan & Memonitor di Server HPC

### 1. Menjalankan di Background Terminal HPC
```bash
ssh hpc-ub "cd /media/user/uid1083/dft_gpaw_graphene_tpms_testing && nohup python3 run_master_dft_hpc.py > master_dft.log 2>&1 &"
```

### 2. Memonitor Progres Kalkulasi
```bash
# Cek apakah proses masih aktif (CPU %):
ssh hpc-ub "ps aux | grep run_master"

# Cek log kalkulasi real-time:
ssh hpc-ub "tail -f /media/user/uid1083/dft_gpaw_graphene_tpms_testing/master_dft.log"

# Cek ringkasan hasil yang sudah selesai:
ssh hpc-ub "cat /media/user/uid1083/dft_gpaw_graphene_tpms_testing/results/results.json"
```

### 3. Menyinkronkan Hasil Terbaru ke Laptop
Cukup jalankan skrip berikut di terminal lokal:
```bash
cd dft && ./sync_hpc_results.sh
```

---

## 📊 Parameter Fisika Komputasi DFT
* **Perangkat Lunak**: GPAW (Grid-based Projector Augmented Wave method)
* **Mode & Basis**: LCAO (*Linear Combination of Atomic Orbitals*), Double-Zeta Polarized (`dzp`)
* **Exchange-Correlation (XC)**: PBE (Perdew-Burke-Ernzerhof, GGA)
* **Sampling k-point**: Densitas Monkhorst-Pack $3.0 \text{ \AA}^{-1}$
* **Smearing Elektronik**: Fermi-Dirac ($\sigma = 0.05 \text{ eV}$)
* **Kriteria Konvergensi SCF**: Energi $< 10^{-5} \text{ eV}$, Densitas $< 10^{-4}$
