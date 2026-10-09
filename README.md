# Multi-Property-CGCNN-LiS

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0%2B-ee4c2c.svg)](https://pytorch.org/)
[![GPAW](https://img.shields.io/badge/DFT-GPAW-orange.svg)](https://wiki.fysik.dtu.dk/gpaw/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Simultaneous Multi-Target Material Screening for Lithium-Sulfur (Li-S) Battery Cathode Scaffolds using Multi-Task Crystal Graph Convolutional Neural Networks (CGCNN) and First-Principles DFT (GPAW) Quantum Calculations.**

---

## 📌 Project Overview

Developing durable, high-performance cathode hosts for **Lithium-Sulfur (Li-S) batteries** requires simultaneous optimization across multiple physical dimensions:
1. **Thermodynamic stability** against chemical decomposition.
2. **Electronic conductivity** for rapid charge-transfer kinetics.
3. **Mechanical resilience** against substantial volume changes (~80% expansion) during lithiation.
4. **Interfacial polysulfide immobilization** to prevent the shuttle effect ($S_8 \to Li_2S_2$).

This repository provides an end-to-end computational screening and machine learning framework that evaluates **five physical properties simultaneously**:

$$\mathbf{y} = \begin{bmatrix} E_f & E_g & K_{VRH} & G_{VRH} & E_{ads} \end{bmatrix}^T$$

* **$E_f$ (Formation Energy / Curvature Strain, eV/atom)**: Thermodynamic crystal phase stability and strain energy.
* **$E_g$ (Band Gap, eV)**: Electronic conductivity character (semi-metallic Dirac-like transport).
* **$K_{VRH}$ (Bulk Modulus, GPa)**: Isotropic compression resilience against volumetric expansion.
* **$G_{VRH}$ (Shear Modulus, GPa)**: Resistance against shear stresses in porous channels.
* **$E_{ads}$ (Adsorption Energy, eV)**: Chemical binding affinity for lithium polysulfides.

---

## 📁 Repository Structure

```
Multi-Property-CGCNN-LiS/
├── notebook/                                # [CGCNN RESEARCH NOTEBOOK]
│   └── Thesis_Principle_Discovery.ipynb     # Main primary thesis research notebook
├── data/                                    # Datasets, splits, and atomic embeddings
│   ├── atom_init.json                       # 92-dim CGCNN elemental feature vectors
│   ├── catalyst_principle_discovery_5targets.csv # 1,366 samples training dataset
│   ├── evaluasi_metrik_pd_train_val_test.csv# Metric summary table across all partitions
│   └── *.csv                                # Prediction results and rankings
├── models/                                  # Trained PyTorch model checkpoints
│   └── cgcnn_finetuned_dft_tpms.pt          # Fine-tuned weights aligned to DFT GPAW
├── figures/                                 # High-resolution publication figures (300 DPI)
│   ├── tpms_pd_pristine_characterization.png# Multi-panel intrinsic property charts
│   └── tpms_pd_5d_radar_chart.png           # 5D composite radar utility chart
├── graphene_tpms/                           # 5 Pristine 3D Graphene TPMS unit cells (.cif)
│   ├── graphene_sheet_diamond.cif
│   ├── graphene_sheet_gyroid.cif
│   ├── graphene_sheet_iwp.cif
│   ├── graphene_sheet_neovius.cif
│   └── graphene_sheet_primitive.cif
├── cifs_graphene_tpms_adsorbate/            # Confined polysulfide complexes in TPMS channels
├── cif_files_pd/                            # 1,366 CIF files for the catalyst dataset
├── scripts/                                 # Utility scripts
│   └── restore_cif_files_pd.py              # Script to verify and restore dataset CIFs
│
└── dft/                                     # [FIRST-PRINCIPLES DFT QUANTUM CALCULATIONS]
    ├── run_master_dft_hpc.py                # Kode Python yang dijalankan di HPC
    ├── figures/                             # Folder figure dari DFT (DOS plots)
    │   ├── dos_all.png                      # Visualisasi DOS perbandingan seluruh struktur
    │   └── dos_*.png                        # Individual DOS plots (Diamond, Gyroid, IWP, Neovius, Primitive)
    └── results/                             # Hasil kalkulasi DFT (energi, struktur, data DOS)
        ├── hasil_pure_dft_tpms_polysulfide.xlsx # Rekapitulasi Excel hasil DFT
        ├── results.json                     # Ground-truth DFT results (JSON)
        ├── *_dos.dat                        # DOS numerical data
        └── *_adsorbed.cif                   # Adsorbed geometry complexes
```

---

## 🔬 Computational Methodology

### 1. Multi-Head CGCNN (Machine Learning)
* **Dataset**: 1,366 crystal structures curated from Terence Tao et al. (*Nature Communications*) augmented with **JARVIS-DFT** and **Materials Project**.
* **Features**: 92-dimensional periodic atomic embeddings, 41 Gaussian radial basis edge expansions, and 19 domain-specific physicochemical descriptors.
* **Architecture**: 3-seed ensemble (**42, 101, 777**) processing 3D periodic crystal graphs with gated non-linearities and dedicated regression heads for bulk properties.

### 2. First-Principles DFT GPAW (Quantum Ground-Truth)
* **Method**: Grid-based Projector Augmented Wave (GPAW) in LCAO mode (`dzp` basis set, PBE exchange-correlation).
* **Reference Ground State**: Planar 2D Graphene sheet ($E_{\text{gr}} = -9.141\text{ eV/atom}$).
* **Properties Computed**: Single-point SCF, Birch-Murnaghan EOS bulk modulus ($K$), shear modulus ($G$), electronic band gap & DOS, and polysulfide adsorption energy ($S_8 \to Li_2S_2$).

---

## 🚀 Quick Start

### 1. Installation
Clone the repository and set up a virtual environment:
```bash
git clone https://github.com/amar-codes02/Multi-Property-CGCNN-LiS.git
cd Multi-Property-CGCNN-LiS
conda create -n cgcnn-lis python=3.9 -y
conda activate cgcnn-lis
pip install torch torchvision torchaudio
pip install pymatgen ase scikit-learn pandas matplotlib seaborn py3Dmol openpyxl
```

### 2. Running the Primary Research Notebook
Launch Jupyter Lab or Jupyter Notebook:
```bash
jupyter notebook notebook/Thesis_Principle_Discovery.ipynb
```
Run the notebook sequentially to reproduce all feature extraction, model training, evaluation metrics, and TPMS cathode rankings.

### 3. Monitoring or Syncing DFT Results from HPC
```bash
cd dft
./sync_hpc_results.sh
```

---

## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
