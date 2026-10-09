# Multi-Property-CGCNN-LiS

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0%2B-ee4c2c.svg)](https://pytorch.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Simultaneous Multi-Target Material Screening for Lithium-Sulfur (Li-S) Battery Cathode Scaffolds using Multi-Task Crystal Graph Convolutional Neural Networks (CGCNN) and a Hybrid Physicochemical Ensemble.**

---

## 📌 Project Overview

Developing durable, high-performance cathode hosts for **Lithium-Sulfur (Li-S) batteries** requires simultaneous optimization across multiple physical dimensions:
1. **Thermodynamic stability** against chemical decomposition.
2. **Electronic conductivity** for rapid charge-transfer kinetics.
3. **Mechanical resilience** against substantial volume changes (~80% expansion) during lithiation.
4. **Interfacial polysulfide immobilization** to prevent the shuttle effect ($S_8 \to Li_2S_2$).

This repository provides an end-to-end computational screening and machine learning framework that evaluates **five physical properties simultaneously**:

$$\mathbf{y} = \begin{bmatrix} E_f & E_g & K_{VRH} & G_{VRH} & E_{ads} \end{bmatrix}^T$$

* **$E_f$ (Formation Energy, eV/atom)**: Thermodynamic crystal phase stability.
* **$E_g$ (Band Gap, eV)**: Electronic conductivity character.
* **$K_{VRH}$ (Bulk Modulus, GPa)**: Isotropic compression resilience.
* **$G_{VRH}$ (Shear Modulus, GPa)**: Resistance against shear stresses in porous channels.
* **$E_{ads}$ (Adsorption Energy, eV)**: Chemical binding affinity for lithium polysulfides.

---

## 🔬 Dataset & Methodology

### 1. Training & Benchmark Dataset
* **Principle Discovery Electrocatalysts**: 1,366 crystal structures curated from Terence Tao et al. (*Nature Communications* / UC Berkeley).
* **Augmented Multi-Source DFT**: High-precision bulk elasticity and thermodynamic data integrated from **JARVIS-DFT** and **Materials Project (MEGNet)**.
* **Features**: 92-dimensional periodic atomic embeddings (`atom_init.json`), 41 Gaussian radial basis edge expansions, and 23 domain-specific physicochemical descriptors.

### 2. Machine Learning Architecture
* **Multi-Head CGCNN**: 3-seed ensemble (**42, 101, 777**) processing 3D periodic crystal graphs with gated non-linearities and dedicated regression heads for bulk properties.
* **Hybrid Adsorption Ensemble**: Gradient-boosted and tree-based ensemble (*ExtraTrees*, *Random Forest*, *GBM*) for polysulfide binding predictions.
* **Zero Data Leakage Protocol**: Stratified multi-target splitting into Train (80%, 1,092 samples), Validation (10%, 137 samples), and Held-Out Test (10%, 137 samples).

### 3. Application: Triply Periodic Minimal Surface (TPMS) Graphene
The pipeline characterizes and ranks 5 pristine 3D Graphene TPMS porous scaffolds:
* **Diamond** (332 C atoms)
* **Gyroid** (244 C atoms)
* **IWP** (228 C atoms)
* **Neovius** (188 C atoms)
* **Primitive** (200 C atoms)

---

## 📁 Repository Structure

```
Multi-Property-CGCNN-LiS/
├── data/                                    # Processed datasets, splits, and atomic embeddings
│   ├── atom_init.json                       # 92-dim CGCNN elemental feature vectors
│   ├── catalyst_principle_discovery_5targets.csv # Full dataset with 5 target properties
│   ├── evaluasi_metrik_pd_train_val_test.csv# Metric summary table across all partitions
│   └── *.csv                                # Prediction results and rankings
├── cif_files_pd/                            # 1,366 CIF files for the catalyst dataset
├── graphene_tpms/                           # Pristine 3D Graphene TPMS unit cells (.cif)
├── cifs_graphene_tpms_adsorbate/            # Confined polysulfide complexes in TPMS channels
├── figures/                                 # Visualizations, parity plots, and 5D radar charts
├── notebook/
│   └── Thesis_Principle_Discovery.ipynb     # Main executable research notebook
├── .gitignore
└── README.md
```

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
pip install pymatgen ase scikit-learn pandas matplotlib seaborn py3Dmol
```

### 2. Running the Research Notebook
Launch Jupyter Lab or Jupyter Notebook:
```bash
jupyter notebook notebook/Thesis_Principle_Discovery.ipynb
```
Run the notebook sequentially to reproduce all feature extraction, model training, evaluation metrics, and TPMS cathode rankings.

---

## 📊 Summary of Results

* **High Generalization Accuracy**: Multi-task parity plots demonstrate tight $R^2$ correlation across all 5 target physical properties on the strictly isolated test partition.
* **Porous Carbon Scaffolding**: Pristine 3D Graphene TPMS structures exhibit metallic conductivity ($E_g = 0\text{ eV}$), high structural compliance, and favorable polysulfide confinement.

---

## 📜 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
