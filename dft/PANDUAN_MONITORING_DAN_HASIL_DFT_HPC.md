# 📘 PANDUAN LENGKAP: MONITORING & PENGAMBILAN HASIL DFT HPC

Dokumen ini berisi rangkuman panduan untuk memantau proses komputasi **Pure DFT GPAW** yang sedang berjalan di background server HPC (`nebula`), serta cara melihat dan mengecek seluruh hasil perhitungan yang disinkronkan ke komputer lokal Anda.

---

## 📌 Ringkasan Pekerjaan yang Sedang Dijalankan di HPC

| Kategori | Rincian |
| :--- | :--- |
| **Server Host** | `nebula` (HPC UB via `ssh hpc-ub` port 2222) |
| **Proses Background** | `python3 run_master_dft_hpc.py` |
| **5 Struktur TPMS** | `neovius` (188 C), `primitive` (200 C), `iwp` (228 C), `gyroid` (244 C), `diamond` (332 C) |
| **5 Produk Polisulfida** | $S_8$, $Li_2S_8$, $Li_2S_6$, $Li_2S_4$, $Li_2S_2$ (25 sistem kompleks) |
| **Target Fisika** | 1. Relaksasi Geometri ($f_{max} \le 0.05$ eV/Å)<br>2. Band Gap & DOS<br>3. Energi Formasi ($E_f$ vs Graphene)<br>4. Energi Adsorpsi ($E_{ads}$)<br>5. Modulus Bulk ($K$) & Shear ($G$) |
| **Folder Output HPC** | `/media/user/uid1083/dft_gpaw_graphene_tpms_testing/results/` |
| **Folder Output Lokal** | `/home/user/Amarus/cgcnn_data_multiproperty/dft_gpaw_graphene_tpms_testing/results/` |

---

## 🔍 1. Cara Cek Status Proses di HPC

Buka terminal di komputer Anda, lalu jalankan perintah berikut:

### A. Cek Apakah Program Masih Berjalan
```bash
ssh hpc-ub "ps aux | grep run_master_dft_hpc | grep -v grep"
```
* **Jika Muncul Baris Proses**: Program masih aktif berjalan (memakai CPU multi-core HPC).
* **Jika Kosong**: Seluruh rangkaian 5 TPMS dan 5 polisulfida sudah selesai dihitung.

### B. Pantau Log SCF Real-time (*Live Monitoring*)
Untuk melihat iterasi konvergensi SCF yang sedang berlangsung detik demi detik:
```bash
# Pantau iterasi yang sedang aktif:
ssh hpc-ub "tail -f /media/user/uid1083/dft_gpaw_graphene_tpms_testing/results/neovius_relax_scf.txt"
```
*(Tekan `Ctrl + C` pada keyboard untuk keluar dari tampilan live).*

### C. Cek Log Utama Master Runner
```bash
ssh hpc-ub "tail -n 30 /media/user/uid1083/dft_gpaw_graphene_tpms_testing/master_dft.log"
```

---

## 🔄 2. Cara Mengambil Hasil (Sudah Otomatis)

### A. Auto-Sync Background (Sudah Aktif)
Daemon `sync_hpc_results.sh` sudah berjalan otomatis di background komputer lokal Anda. Setiap 15 detik, daemon ini memeriksa HPC dan langsung mendownload file baru ke:
📁 `/home/user/Amarus/cgcnn_data_multiproperty/dft_gpaw_graphene_tpms_testing/results/`

### B. Cek File yang Sudah Tersinkron di Komputer Lokal
```bash
ls -lth /home/user/Amarus/cgcnn_data_multiproperty/dft_gpaw_graphene_tpms_testing/results/
```

### C. Sinkronisasi Manual Kapan Saja
Jika Anda baru menyalakan laptop dan ingin segera menarik seluruh data terbaru dari HPC sekaligus:
```bash
rsync -avz hpc-ub:/media/user/uid1083/dft_gpaw_graphene_tpms_testing/results/ /home/user/Amarus/cgcnn_data_multiproperty/dft_gpaw_graphene_tpms_testing/results/
```

---

## 📊 3. Cara Membaca & Menganalisis Hasil Akhir

Setelah perhitungan selesai, hasil akan tersimpan dalam beberapa format file di folder `results/`:

| Nama File | Fungsi & Isi Data |
| :--- | :--- |
| **`results.json`** | Database nilai mentah seluruh sifat fisik ($E_f, E_g, K, G, E_{ads}$ untuk seluruh 5 TPMS & 5 polisulfida). |
| **`summary.csv`** | Tabel ringkasan yang siap dibuka di Excel / LibreOffice / Pandas. |
| **`<nama>_relaxed.cif`** | File koordinat 3D struktur TPMS yang sudah setimbang dan rileks ($f_{max} \le 0.05$ eV/Å). |
| **`<nama>_<mol>_adsorbed.cif`** | Struktur 3D kompleks adsorpsi polisulfida di dalam kanal pori TPMS. |
| **`<nama>_dos.dat` & `dos_all.png`** | Data densitas keadaan elektronik (DOS) dan grafik kurva DOS. |
| **`adsorption_profile_5species.png`** | Grafik profil afinitas penjeratan $S_8 \to Li_2S_8 \to Li_2S_6 \to Li_2S_4 \to Li_2S_2$. |
| **`summary_bars.png`** | Grafik perbandingan Band Gap, Energi Formasi, dan Bulk Modulus antar kelima TPMS. |

---

## 📓 4. Cara Menampilkan Hasil di Jupyter Notebook

1. Buka file `tpms_gpaw.ipynb`.
2. Langsung gulir ke sel paling bawah (**Langkah 6: Ringkasan 5 Sifat Fisik & Visualisasi Akhir**).
3. Jalankan sel tersebut. Notebook akan otomatis membaca data dari `results.json` dan menampilkan:
   * Tabel komparasi 5 TPMS secara interaktif.
   * Grafik batang sifat elektronik & mekanik.
   * Grafik kurva penjeratan 5 polisulfida untuk bab hasil & pembahasan tesis Anda.
