# MnIn2Se4 lattice thermal conductivity calculation data (force-constant fitting + phono3py RTA)

System: **MnIn₂Se₄** (rhombohedral R-3m #166, 7-atom primitive cell = Mn₁In₂Se₄; 3×3×3 supercell with 189 atoms, 45 random-displacement configurations + 1 equilibrium reference frame).

Workflow: determine shells and cutoff → pheasy force-constant fitting (LASSO / ALASSO + HARM_DENSE, OLS as control, symfc as control) → κ from phono3py RTA.

> This directory is partly parallel to `Mn2In2Se5`; both belong to the lattice thermal conductivity work on the MnIn₂Se system.

---

## Directory structure

| directory | contents |
|---|---|
| `00_report/` | **main report** `MnIn2Se4_report_20260928.md`, RMSE and κ side-by-side table, shell and cutoff table, dynamic stability check, mesh convergence |
| `01_structures/` | `POSCAR_primitive` (7-atom primitive cell), `SPOSCAR_3x3x3` (189-atom supercell) |
| `02_dataset/` | `dataset_disps.npy` / `dataset_forces.npy` (46 frames), `disp_matrix.pkl` / `force_matrix.pkl` (45 frames), `dataset_alignment.json`, `ndata_total.txt` |
| `03_fc_fits/` | **39 fit directories**: OLS / LASSO / ALASSO × 12 cutoffs (2.85–7.00 Å) + symfc × 3. Each contains `POSCAR`, `fc2.hdf5`, `fc3.hdf5`, `run.log` |
| `04_kappa_summary/` | κ summary for each cutoff, NU regeneration logs, mesh convergence |
| `05_kappa_hdf5/` | **39 sets of raw κ output**: one `kappa-m111111.hdf5` per directory (11³ q-mesh, 300 K, RTA, with isotopes) + `POSCAR` + `kappa_run.log` |
| `06_scripts/` | fit driver, κ solver, symfc script, summary script, official `pheasy_fit.sh`(v6c) and `prepare_dataset.py` |
| `07_holdout_rmse/` | whole-configuration 5-fold held-out force-prediction RMSE logs (`holdout_c3_*.log`) |

---

## Key parameters

- **Fitting**: `pheasy-gpu-gh` HEAD `b91da0f` (3090, conda env `wc`, RTX 3090)
  `FIT_ORDER=3  HARM_DENSE=true  C2_CUTOFF=None  C3_CUTOFF=<cutoff>`
  (HARM_DENSE = FC2 fully dense with no L1 penalty, sparsity applied only to FC3)
- **Cutoffs**: 2.85 / 3.70 / 3.95 / 4.30 / 4.70 / 5.00 / 5.20 / 5.50 / 5.80 / 6.10 / 6.50 / 7.00 Å
  (all of them fall in the gaps between atomic-distance shells of the supercell, see `00_report/shells_and_cutoff_20260928.md`)
- **Independent control**: `symfc` 1.7.3 (2.85 / 4.30 / 7.00 Å)
- **κ**: `phono3py` 4.4.0 (hanhai venv mace_cpu), 11³ q-mesh, 300 K, RTA,
  `is_isotope=True, is_N_U=True, write_kappa=True`, supercell = 3x3x3, primitive = POSCAR

---

## Results at a glance

> The κ values below are from the 2026-09-28 3×3×3 fit and are superseded — see `../REPORT_MnIn2Se4_20261115.md`.

### κ versus the third-order cutoff (300 K, 11×11×11, RTA + isotopes, units W/m/K)

| cutoff (Å) | 2.85 | 3.70 | 3.95 | 4.30 | 4.70 | 5.00 | 5.20 | 5.50 | 5.80 | 6.10 | 6.50 | 7.00 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| κ_xx | 5.2219 | 5.1257 | 3.8368 | 2.7065 | 3.8401 | 3.9725 | 3.9126 | 3.8370 | 3.7306 | 3.3048 | 3.3703 | 3.1752 |
| κ_zz | 4.0167 | 3.5122 | 1.0338 | 0.5038 | 1.4173 | 1.9378 | 2.1302 | 1.8752 | 1.2729 | 1.1477 | 0.8922 | 0.6293 |

(ALASSO + HARM_DENSE, report §4.1; κ is sensitive to the cutoff, especially κ_zz)

### Mesh convergence (ALASSO, 300 K)

| q-mesh | 11×11×11 | 14×14×14 | 18×18×18 |
|---|---|---|---|
| length (Å) | 40.0 | 50.0 | 62.5 |
| κ_xx = κ_yy | 3.8828 | 3.9221 | 3.8912 |
| κ_zz | 1.9471 | 1.7513 | **1.5324** |
| relative change | — | 10.1% | 12.5% |

**Conclusion**: the in-plane κ_xx/κ_yy is converged (0.8–1.0%); **κ_zz is not yet converged** (it drops by 10–12.5% per step, and the scan is truncated by the `MESH_CONV_MAX_LENGTH=62.5 Å` upper bound). Interlayer heat transport in layered compounds requires a larger q-mesh.

---

## Data provenance (important)

The four input files provided by the user are identical to the official test material in `pheasy-gpu-gh/test/MnIn2Se4/`:

| file | md5 | note |
|---|---|---|
| `POSCAR` | `cef6af221344890c582da4cfb638f438` | identical to the official file byte by byte |
| `dataset_disps.npy` | `f784694ae79ce1b170cb35987ce76c95` | identical to the official file byte by byte (it stores **fractional coordinates**, not displacements) |
| `dataset_forces.npy` | `22906b11c42f2f32c0964b6b267c7424` | identical to the official file byte by byte |
| `SPOSCAR` | `c702c54dabe657730899e838656be9ba` | same structure as the official file, but with a different atom order |

The `disp_matrix.pkl` generated with the official `tools/prepare_dataset.py SPOSCAR dataset_disps.npy dataset_forces.npy --frac --ref -1`
differs from the official file of the same name by a maximum numerical difference of **4e-17**, and the md5 of `force_matrix.pkl` is **completely identical** → preprocessing matches the official pipeline.

Structure: rhombohedral |a|=|b|=|c|=13.4029 Å, α=17.463°; the equivalent hexagonal description is **a = 4.069 Å, c = 39.586 Å** (layered, with interlayer van der Waals bonding).

---

## Reproduction

```bash
# 1) dataset (requires ase/scipy)
python3 06_scripts/prepare_dataset.py SPOSCAR dataset_disps.npy dataset_forces.npy --frac --align-reference

# 2) fitting (one directory per cutoff)
bash 06_scripts/pheasy_fit.sh FIT_METHOD=LASSO  HARM_DENSE=true C2_CUTOFF=None C3_CUTOFF=5.20 FIT_ORDER=3 NCPU=4
bash 06_scripts/pheasy_fit.sh FIT_METHOD=ALASSO HARM_DENSE=true C2_CUTOFF=None C3_CUTOFF=5.20 FIT_ORDER=3 NCPU=4

# 3) kappa
python3 06_scripts/kappa_full_solver.py 11
```

---

Generated: 2026-09-28 (archived) / 2026-11-07 (uploaded and organized).
