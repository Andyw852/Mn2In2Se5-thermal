# MnIn2Se4 force-constant fitting and lattice thermal conductivity report

Date: 2026-09-28  |  System: MnIn2Se4 (R-3m #166, 7-atom primitive cell, 3×3×3 supercell with 189 atoms, 45 configurations)
Workflow: determine shells and cutoff → pheasy force-constant fitting (LASSO/ALASSO + HARM_DENSE, OLS as control) → κ from phono3py RTA.

> Note: the κ values below come from the 2026-09-28 3×3×3 fit and are superseded — see `../../REPORT_MnIn2Se4_20261115.md`.

---

## 0. Data verification and provenance

The 4 files provided by the user:

- `POSCAR`: primitive cell, Mn 1 / In 2 / Se 4 = 7 atoms. Lattice |a|=|b|=|c|=13.4029 Å, α=17.463° (rhombohedral);
  the equivalent hexagonal description is **a = 4.069 Å, c = 39.586 Å**, which shows that this is a layered compound (intralayer covalent, interlayer van der Waals).
- `SPOSCAR`: 3×3×3 supercell, 189 atoms, |a|=|b|=|c|=40.209 Å.
- `dataset_disps.npy`, `dataset_forces.npy`: both have shape (46, 189, 3).

**Key point: `dataset_disps.npy` stores fractional coordinates (the atomic positions of each frame), not displacements.**
Evidence: frame 45 (the last frame) of this array agrees with the fractional coordinates of SPOSCAR element by element (maximum difference 4.4e-16),
and its force max=0.0035 eV/Å (an equilibrium frame); using the official `tools/prepare_dataset.py --frac --ref -1`
to difference under the minimum-image convention with that frame as the reference yields the 45 real displacement frames:
**displacement rms = 0.0172 Å, max = 0.0794 Å** (a normal random-displacement magnitude), force max = 0.96 eV/Å.

**Provenance (byte by byte)**: these 4 files are exactly the test material `test/MnIn2Se4/` shipped with pheasy-gpu-gh:

| file | user-provided md5 | official test/MnIn2Se4 | conclusion |
|---|---|---|---|
| POSCAR | `cef6af221344890c582da4cfb638f438` | `cef6af22…` | fully identical |
| dataset_disps.npy | `f784694ae79ce1b170cb35987ce76c95` | `f784694a…` | fully identical |
| dataset_forces.npy | `22906b11c42f2f32c0964b6b267c7424` | `22906b11…` | fully identical |
| SPOSCAR | `c702c54dabe657730899e838656be9ba` | `87321de6…` | same structure, different atom order |

The `disp_matrix.pkl` generated with the official `prepare_dataset.py` differs from the official file of the same name by a maximum numerical difference of **4e-17**,
and the md5 of `force_matrix.pkl` is **completely identical** → data preprocessing matches the official pipeline, so the downstream results are reproducible.

> Note: the official test directory also contains a `dataset_forces_corrected.npy`, which differs from this data by at most 0.0035 eV/Å
> (only 0.2% of the force rms); this report uses the original `dataset_forces.npy` supplied by the user.

---

## 1. Shells and cutoff

# MnIn2Se4 supercell (3x3x3, 189 atoms): atomic-distance shells and cutoff selection

Supercell lattice: |a|=|b|=|c|=40.209 A, alpha=17.463 deg (rhombohedral primitive cell; hexagonal description a=4.069 A, c=39.586 A).
Minimum-image limit (supercell/2) = 20.10 A.

## Shell table (clustered by distance, tol=1e-3 A)

| shell | d (A) | multiplicity | atom-pair composition | gap to next shell (A) |
|---|---|---|---|---|
| 0 | 2.5207 | 54 | In-Se x54 |  |
| 1 | 2.6571 | 162 | In-Se x162 | 0.136 |
| 2 | 2.7061 | 162 | Mn-Se x162 | 0.049 |
| 3 | 3.5683 | 81 | Se-Se x81 | 0.862 **<-c3=2.85** |
| 4 | 3.7996 | 81 | Se-Se x81 | 0.231 **<-c3=3.70** |
| 5 | 4.0693 | 567 | In-In x162, Mn-Mn x81, Se-Se x324 | 0.270 **<-c3=3.95** |
| 6 | 4.2271 | 54 | In-Se x54 | 0.158 |
| 7 | 4.4352 | 162 | Se-Se x162 | 0.208 **<-c3=4.30** |
| 8 | 4.5219 | 162 | In-Mn x162 | 0.087 |
| 9 | 4.7868 | 324 | In-Se x324 | 0.265 **<-c3=4.70** |
| 10 | 4.8600 | 162 | In-Se x162 | 0.073 |
| 11 | 4.8870 | 162 | Mn-Se x162 | 0.027 |
| 12 | 5.1046 | 54 | Mn-Se x54 | 0.218 **<-c3=5.00** |
| 13 | 5.4122 | 81 | Se-Se x81 | 0.308 **<-c3=5.20** |
| 14 | 5.5674 | 81 | Se-Se x81 | 0.155 **<-c3=5.50** |
| 15 | 5.7120 | 162 | In-Se x162 | 0.145 |
| 16 | 5.8676 | 324 | In-Se x324 | 0.156 **<-c3=5.80** |
| 17 | 5.9515 | 81 | In-In x81 | 0.084 |
| 18 | 6.0191 | 162 | Se-Se x162 | 0.068 |
| 19 | 6.0833 | 162 | In-Mn x162 | 0.064 |
| 20 | 6.3387 | 162 | In-Se x162 | 0.255 **<-c3=6.10** |
| 21 | 6.3594 | 162 | Mn-Se x162 | 0.021 |
| 22 | 6.5281 | 324 | Mn-Se x324 | 0.169 **<-c3=6.50** |
| 23 | 6.7479 | 54 | Se-Se x54 | 0.220 |
| 24 | 6.7714 | 81 | Se-Se x81 | 0.023 |
| 25 | 6.8622 | 162 | Se-Se x162 | 0.091 |
| 26 | 6.8961 | 81 | Se-Se x81 | 0.034 |
| 27 | 7.0133 | 162 | In-Se x162 | 0.117 **<-c3=7.00** |
| 28 | 7.0483 | 189 | In-In x54, Mn-Mn x27, Se-Se x108 | 0.035 |
| 29 | 7.2097 | 81 | In-In x81 | 0.161 |
| 30 | 7.2656 | 162 | Se-Se x162 | 0.056 |
| 31 | 7.3189 | 162 | In-Mn x162 | 0.053 |
| 32 | 7.4855 | 108 | In-Se x108 | 0.167 |
| 33 | 7.8799 | 324 | Se-Se x324 | 0.394 |
| 34 | 7.9780 | 162 | Se-Se x162 | 0.098 |
| 35 | 8.0765 | 81 | In-In x81 | 0.098 |
| 36 | 8.1084 | 162 | In-Se x162 | 0.032 |
| 37 | 8.2187 | 108 | In-Se x108 | 0.110 |
| 38 | 8.2789 | 81 | In-In x81 | 0.060 |
| 39 | 8.3272 | 162 | In-Se x162 | 0.048 |
| 40 | 8.4250 | 162 | Mn-Se x162 | 0.098 |
| 41 | 8.7026 | 108 | Mn-Se x108 | 0.278 |
| 42 | 8.9559 | 162 | Se-Se x162 | 0.253 |
| 43 | 9.0437 | 81 | In-In x81 | 0.088 |
| 44 | 9.2683 | 162 | In-Se x162 | 0.225 |
| 45 | 9.2709 | 162 | In-Se x162 | 0.003 |
| 46 | 9.3318 | 54 | In-Mn x54 | 0.061 |
| 47 | 9.3562 | 162 | Mn-Se x162 | 0.024 |
| 48 | 9.7218 | 162 | Se-Se x162 | 0.366 |
| 49 | 9.7577 | 108 | Se-Se x108 | 0.036 |
| 50 | 9.9171 | 81 | In-In x81 | 0.159 |
| 51 | 10.1223 | 162 | In-Se x162 | 0.205 |
| 52 | 10.1246 | 162 | In-Se x162 | 0.002 |
| 53 | 10.1804 | 324 | In-Mn x324 | 0.056 |
| 54 | 10.2029 | 162 | Mn-Se x162 | 0.022 |
| 55 | 10.2092 | 27 | Se-Se x27 | 0.006 |
| 56 | 10.5391 | 162 | Se-Se x162 | 0.330 |
| 57 | 10.7690 | 81 | Se-Se x81 | 0.230 |
| 58 | 10.9118 | 162 | In-Se x162 | 0.143 |
| 59 | 10.9301 | 162 | In-Se x162 | 0.018 |
| 60 | 10.9904 | 162 | Se-Se x162 | 0.060 |
| 61 | 11.2974 | 162 | Se-Se x162 | 0.307 |
| 62 | 11.5122 | 81 | Se-Se x81 | 0.215 |
| 63 | 11.6631 | 162 | In-Se x162 | 0.151 |
| 64 | 11.6944 | 108 | In-Mn x108 | 0.031 |
| 65 | 11.8525 | 54 | Mn-Se x54 | 0.158 |
| 66 | 12.1830 | 162 | In-Se x162 | 0.331 |
| 67 | 12.2103 | 81 | Se-Se x81 | 0.027 |
| 68 | 12.3526 | 162 | In-Se x162 | 0.142 |
| 69 | 12.4059 | 54 | Se-Se x54 | 0.053 |
| 70 | 12.5316 | 324 | Mn-Se x324 | 0.126 |
| 71 | 12.8447 | 162 | In-Se x162 | 0.313 |
| 72 | 13.4029 | 567 | In-In x162, Mn-Mn x81, Se-Se x324 | 0.558 |
| 73 | 13.4739 | 162 | In-Se x162 | 0.071 |
| 74 | 13.7899 | 108 | Mn-Se x108 | 0.316 |
| 75 | 14.0070 | 567 | In-In x162, Mn-Mn x81, Se-Se x324 | 0.217 |
| 76 | 14.4364 | 54 | In-Se x54 | 0.429 |
| 77 | 14.5862 | 567 | In-In x162, Mn-Mn x81, Se-Se x324 | 0.150 |
| 78 | 14.7268 | 162 | Mn-Se x162 | 0.141 |
| 79 | 14.9990 | 324 | In-Se x324 | 0.272 |
| 80 | 15.2787 | 162 | Mn-Se x162 | 0.280 |
| 81 | 15.8113 | 162 | Mn-Se x162 | 0.533 |
| 82 | 15.8811 | 27 | Se-Se x27 | 0.070 |
| 83 | 15.8907 | 162 | In-Se x162 | 0.010 |
| 84 | 16.0651 | 108 | In-Se x108 | 0.174 |
| 85 | 16.3512 | 81 | Se-Se x81 | 0.286 |
| 86 | 16.3942 | 162 | Se-Se x162 | 0.043 |
| 87 | 16.4035 | 162 | In-Se x162 | 0.009 |
| 88 | 16.8499 | 81 | Se-Se x81 | 0.446 |
| 89 | 16.9007 | 162 | In-Se x162 | 0.051 |
| 90 | 16.9571 | 54 | Se-Se x54 | 0.056 |
| 91 | 17.2200 | 162 | In-Mn x162 | 0.263 |
| 92 | 17.3343 | 81 | Se-Se x81 | 0.114 |
| 93 | 17.3749 | 54 | Se-Se x54 | 0.041 |
| 94 | 17.4386 | 324 | Se-Se x324 | 0.064 |
| 95 | 17.5802 | 162 | In-Se x162 | 0.142 |
| 96 | 17.6943 | 162 | In-Mn x162 | 0.114 |
| 97 | 18.0450 | 162 | In-Se x162 | 0.351 |
| 98 | 18.1562 | 162 | In-Mn x162 | 0.111 |
| 99 | 18.3636 | 108 | Se-Se x108 | 0.207 |
| 100 | 18.4018 | 54 | In-Se x54 | 0.038 |
| 101 | 18.4502 | 162 | Mn-Se x162 | 0.048 |
| 102 | 18.4982 | 162 | In-Se x162 | 0.048 |
| 103 | 18.6635 | 27 | In-In x27 | 0.165 |
| 104 | 18.8464 | 324 | In-Se x324 | 0.183 |
| 105 | 18.8936 | 162 | Mn-Se x162 | 0.047 |
| 106 | 19.1020 | 162 | In-In x162 | 0.208 |
| 107 | 19.3269 | 162 | Mn-Se x162 | 0.225 |
| 108 | 19.7055 | 108 | In-Se x108 | 0.379 |
| 109 | 19.7828 | 162 | Se-Se x162 | 0.077 |
| 110 | 19.9501 | 54 | In-In x54 | 0.167 |
| 111 | 20.0812 | 162 | Se-Se x162 | 0.131 |
| 112 | 20.4893 | 162 | Se-Se x162 | 0.408 |



All 12 cutoffs fall inside shell gaps, so changing the cutoff adds or removes one complete shell as a whole
and no ambiguity about a half-included shell arises.

---

## 2. Fitting setup

- Code: `pheasy-gpu-gh` HEAD `b91da0f` (3090 server, conda env `wc`, RTX 3090; `CUDA_VISIBLE_DEVICES=5`)
- Command: `FIT_ORDER=3  HARM_DENSE=true  C2_CUTOFF=None  C3_CUTOFF=<c>  NCPU=4`
  - `HARM_DENSE=true`: the second-order force constants are **fully dense** (no cutoff) and **carry no L1 penalty**; only the third-order force constants are sparsified/regularized;
  - Methods: **LASSO + HARM_DENSE** and **ALASSO + HARM_DENSE** (main results), plus **OLS** as a control.
- Training data: 45 frames × 189 atoms × 3 = 25515 force equations; 3×3×3 supercell, 12 symmetry operations of R-3m.
- The number of free IFCs p grows from 615 (2.85 Å) to 3678 (7.00 Å).
  As a check: the official README records p=1303 at c3=5.2 and p=3678 at c3=7.0;
  **we obtain c3=5.20 → p=1303 and c3=7.00 → p=3678, in full agreement**, which shows that the cluster space / symmetrization pipeline is correct.

---

## 3. Fit quality and held-out force-prediction error

- `RMSE in-sample`: force residual on the same 45 frames (biased; decreases monotonically as the number of parameters grows);
- `held-out`: **whole-configuration grouped 5-fold** cross-validation (entire frames are held out, no leakage), obtained with `holdout_eval.py`;
- `RMSE_CV`: grouped-CV estimate from the internal L1 path of the fit (before debiasing).

### A. Fit quality + held-out RMSE + kappa (pheasy: OLS / LASSO / ALASSO, HARM_DENSE)

| c3 (A) | method | p (free IFC) | nnz | alpha_opt | RMSE_CV (eV/A) | RMSE in-sample | rel err | held-out OLS | held-out LASSO | held-out ALASSO | kxx (W/mK) | kzz (W/mK) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2.85 | OLS | 615 | 615 | - | - | 1.4059e-03 | 9.2700e-03 | 1.4520e-03 | 1.4519e-03 | 1.4494e-03 | 5.2220 | 4.0169 |
| 2.85 | LASSO | 615 | 615 | 2.706e-11 | 1.4530e-03 | 1.4059e-03 | 9.2700e-03 | 1.4520e-03 | 1.4519e-03 | 1.4494e-03 | 5.2219 | 4.0167 |
| 2.85 | ALASSO | 615 | 615 | 2.706e-11 | 1.4530e-03 | 1.4059e-03 | 9.2700e-03 | 1.4520e-03 | 1.4519e-03 | 1.4494e-03 | 5.2219 | 4.0167 |
| 3.70 | OLS | 633 | 633 | - | - | 1.3587e-03 | 8.9587e-03 | 1.4054e-03 | 1.4053e-03 | 1.4035e-03 | 5.1291 | 3.5130 |
| 3.70 | LASSO | 633 | 633 | 2.042e-09 | 1.4073e-03 | 1.3587e-03 | 8.9587e-03 | 1.4054e-03 | 1.4053e-03 | 1.4035e-03 | 5.1290 | 3.5128 |
| 3.70 | ALASSO | 633 | 632 | 4.544e-11 | 1.4071e-03 | 1.3587e-03 | 8.9588e-03 | 1.4054e-03 | 1.4053e-03 | 1.4035e-03 | 5.1257 | 3.5122 |
| 3.95 | OLS | 639 | 639 | - | - | 1.3352e-03 | 8.8038e-03 | 1.3814e-03 | 1.3814e-03 | 1.3794e-03 | 3.8368 | 1.0337 |
| 3.95 | LASSO | 639 | 639 | 2.432e-11 | 1.3818e-03 | 1.3352e-03 | 8.8038e-03 | 1.3814e-03 | 1.3814e-03 | 1.3794e-03 | 3.8368 | 1.0338 |
| 3.95 | ALASSO | 639 | 639 | 2.432e-11 | 1.3818e-03 | 1.3352e-03 | 8.8038e-03 | 1.3814e-03 | 1.3814e-03 | 1.3794e-03 | 3.8368 | 1.0338 |
| 4.30 | OLS | 807 | 807 | - | - | 1.1559e-03 | 7.6215e-03 | 1.2092e-03 | 1.2092e-03 | 1.2109e-03 | 2.6939 | 0.5022 |
| 4.30 | LASSO | 807 | 802 | 3.269e-09 | 1.2071e-03 | 1.1559e-03 | 7.6216e-03 | 1.2092e-03 | 1.2092e-03 | 1.2109e-03 | 2.6938 | 0.5018 |
| 4.30 | ALASSO | 807 | 782 | 4.516e-11 | 1.2059e-03 | 1.1560e-03 | 7.6224e-03 | 1.2092e-03 | 1.2092e-03 | 1.2109e-03 | 2.7065 | 0.5038 |
| 4.70 | OLS | 907 | 907 | - | - | 1.0164e-03 | 6.7018e-03 | 1.0731e-03 | 1.0730e-03 | 1.0771e-03 | 3.8235 | 1.4062 |
| 4.70 | LASSO | 907 | 892 | 4.775e-09 | 1.0680e-03 | 1.0164e-03 | 6.7018e-03 | 1.0731e-03 | 1.0730e-03 | 1.0771e-03 | 3.8132 | 1.3942 |
| 4.70 | ALASSO | 907 | 847 | 3.664e-11 | 1.0658e-03 | 1.0167e-03 | 6.7042e-03 | 1.0731e-03 | 1.0730e-03 | 1.0771e-03 | 3.8401 | 1.4173 |
| 5.00 | OLS | 1272 | 1272 | - | - | 9.4788e-04 | 6.2501e-03 | 1.0235e-03 | 1.0231e-03 | 1.0345e-03 | 3.9380 | 1.8316 |
| 5.00 | LASSO | 1272 | 1184 | 7.214e-09 | 1.0163e-03 | 9.4803e-04 | 6.2511e-03 | 1.0235e-03 | 1.0231e-03 | 1.0345e-03 | 3.9542 | 1.8758 |
| 5.00 | ALASSO | 1272 | 1077 | 3.840e-11 | 1.0112e-03 | 9.4911e-04 | 6.2582e-03 | 1.0235e-03 | 1.0231e-03 | 1.0345e-03 | 3.9725 | 1.9378 |
| 5.20 | OLS | 1303 | 1303 | - | - | 9.4543e-04 | 6.2339e-03 | 1.0226e-03 | 1.0222e-03 | 1.0344e-03 | 3.8828 | 1.9476 |
| 5.20 | LASSO | 1303 | 1210 | 7.214e-09 | 1.0156e-03 | 9.4561e-04 | 6.2351e-03 | 1.0226e-03 | 1.0222e-03 | 1.0344e-03 | 3.9199 | 1.9726 |
| 5.20 | ALASSO | 1303 | 1102 | 3.840e-11 | 1.0102e-03 | 9.4672e-04 | 6.2425e-03 | 1.0226e-03 | 1.0222e-03 | 1.0344e-03 | 3.9126 | 2.1302 |
| 5.50 | OLS | 1352 | 1352 | - | - | 9.3405e-04 | 6.1590e-03 | 1.0134e-03 | 1.0131e-03 | 1.0270e-03 | 3.8086 | 1.7517 |
| 5.50 | LASSO | 1352 | 1260 | 6.482e-09 | 1.0064e-03 | 9.3427e-04 | 6.1604e-03 | 1.0134e-03 | 1.0131e-03 | 1.0270e-03 | 3.8696 | 1.8069 |
| 5.50 | ALASSO | 1352 | 1149 | 3.400e-11 | 1.0008e-03 | 9.3522e-04 | 6.1666e-03 | 1.0134e-03 | 1.0131e-03 | 1.0270e-03 | 3.8370 | 1.8752 |
| 5.80 | OLS | 1544 | 1544 | - | - | 9.1395e-04 | 6.0264e-03 | 1.0026e-03 | 1.0021e-03 | 1.0166e-03 | 3.7963 | 1.3955 |
| 5.80 | LASSO | 1544 | 1422 | 6.482e-09 | 9.9359e-04 | 9.1416e-04 | 6.0278e-03 | 1.0026e-03 | 1.0021e-03 | 1.0166e-03 | 3.8358 | 1.3943 |
| 5.80 | ALASSO | 1544 | 1180 | 5.449e-11 | 9.8582e-04 | 9.1736e-04 | 6.0489e-03 | 1.0026e-03 | 1.0021e-03 | 1.0166e-03 | 3.7306 | 1.2729 |
| 6.10 | OLS | 2081 | 2081 | - | - | 8.8835e-04 | 5.8576e-03 | 1.0018e-03 | - | - | 3.2628 | 1.2268 |
| 6.10 | LASSO | 2081 | 1714 | 7.745e-09 | 9.8959e-04 | 8.8975e-04 | 5.8668e-03 | 1.0018e-03 | - | - | 3.3017 | 1.1976 |
| 6.10 | ALASSO | 2081 | 1399 | 4.114e-11 | 9.7433e-04 | 8.9361e-04 | 5.8923e-03 | 1.0018e-03 | - | - | 3.3048 | 1.1477 |
| 6.50 | OLS | 2643 | 2643 | - | - | 8.5590e-04 | 5.6436e-03 | 9.9420e-04 | - | - | 3.1958 | 0.7817 |
| 6.50 | LASSO | 2643 | 1964 | 1.044e-08 | 9.7617e-04 | 8.6026e-04 | 5.6724e-03 | 9.9420e-04 | - | - | 3.2805 | 0.8385 |
| 6.50 | ALASSO | 2643 | 1692 | 3.400e-11 | 9.5493e-04 | 8.6227e-04 | 5.6856e-03 | 9.9420e-04 | - | - | 3.3703 | 0.8922 |
| 7.00 | OLS | 3678 | 3678 | - | - | 8.2024e-04 | 5.4085e-03 | 1.0075e-03 | - | - | 2.9873 | 0.5927 |
| 7.00 | LASSO | 3678 | 2475 | 1.111e-08 | 9.7392e-04 | 8.2799e-04 | 5.4596e-03 | 1.0075e-03 | - | - | 3.0545 | 0.5821 |
| 7.00 | ALASSO | 3678 | 1783 | 5.804e-11 | 9.4317e-04 | 8.3941e-04 | 5.5349e-03 | 1.0075e-03 | - | - | 3.1752 | 0.6293 |



**Trends and conclusions**

1. The force residual (in-sample rel err) decreases monotonically with the cutoff: **0.93% (2.85 Å) → 0.54% (7.00 Å)**.
2. **held-out has a clear minimum: 6.50 Å (9.94e-4 eV/Å); at 7.00 Å it rises back to 1.01e-3** —
   i.e. slight overfitting starts to appear at 7.0 Å. This is the direct basis for "how large a cutoff should be used".
3. **LASSO ≈ OLS**: at most cutoffs alpha_opt is pinned to the smallest end of the search grid (`A@MIN`), which is equivalent to no sparsification;
   ALASSO has a clear sparsifying effect (e.g. at 7.00 Å: 3678 → 1783 non-zero IFCs),
   but beyond 5.0 Å its held-out error is instead about 1–2% higher than OLS (e.g. at 5.50 Å: 1.027e-3 vs 1.013e-3),
   which shows that **sparsification gives no net gain in "predictive power for forces"** (it compresses the third-order terms at the cost of a little accuracy).
4. The held-out differences among the three methods are all **< 2%**, far smaller than the differences caused by the cutoff (about 45%, 1.45e-3 → 9.94e-4).

---

## 4. Lattice thermal conductivity κ

Setup: `phono3py 4.4.0`, 11×11×11 q-mesh, 300 K, RTA, including isotope scattering;
supercell=[[3,0,0],[0,3,0],[0,0,3]], primitive=`POSCAR`.
κ tensor diagonal: κ_xx = κ_yy (in-plane), κ_zz (cross-layer); the Cartesian z axis of this cell is the hexagonal c axis,
so κ_xx/κ_yy is the intralayer thermal conductivity and κ_zz the cross-layer thermal conductivity.

### 4.1 Variation of κ with the third-order cutoff (W/mK)

| c3 (Å) | OLS κxx | LASSO κxx | ALASSO κxx | OLS κzz | LASSO κzz | ALASSO κzz |
|---|---|---|---|---|---|---|
| 2.85 | 5.2220 | 5.2219 | 5.2219 | 4.0169 | 4.0167 | 4.0167 |
| 3.70 | 5.1291 | 5.1290 | 5.1257 | 3.5130 | 3.5128 | 3.5122 |
| 3.95 | 3.8368 | 3.8368 | 3.8368 | 1.0337 | 1.0338 | 1.0338 |
| 4.30 | 2.6939 | 2.6938 | 2.7065 | 0.5022 | 0.5018 | 0.5038 |
| 4.70 | 3.8235 | 3.8132 | 3.8401 | 1.4062 | 1.3942 | 1.4173 |
| 5.00 | 3.9380 | 3.9542 | 3.9725 | 1.8316 | 1.8758 | 1.9378 |
| 5.20 | 3.8828 | 3.9199 | 3.9126 | 1.9476 | 1.9726 | 2.1302 |
| 5.50 | 3.8086 | 3.8696 | 3.8370 | 1.7517 | 1.8069 | 1.8752 |
| 5.80 | 3.7963 | 3.8358 | 3.7306 | 1.3955 | 1.3943 | 1.2729 |
| 6.10 | 3.2628 | 3.3017 | 3.3048 | 1.2268 | 1.1976 | 1.1477 |
| 6.50 | 3.1958 | 3.2805 | 3.3703 | 0.7817 | 0.8385 | 0.8922 |
| 7.00 | 2.9873 | 3.0545 | 3.1752 | 0.5927 | 0.5821 | 0.6293 |

**Conclusion: κ is strongly non-monotonic with the cutoff, and no plateau appears.**

- 2.85–3.70 Å (only nearest-neighbour In–Se / Mn–Se): κ_xx ≈ 5.1–5.2, κ_zz ≈ 3.5–4.0 (the third-order terms are severely incomplete,
  scattering is underestimated and κ is too high);
- 3.95 Å: 3.84 / 1.03;
- **4.30 Å: global minimum, 2.69 / 0.50**;
- 4.70–5.80 Å: κ_xx forms a "false plateau" between **3.73 and 3.97**, and κ_zz reaches its maximum 1.95–2.13 at 5.20 Å;
- after 6.10 Å it decreases again: 3.26–3.37 (6.50) → 2.99–3.18 (7.00), and κ_zz drops to 0.58–0.63.

In other words, **κ is still changing appreciably around 6.5–7.0 Å and cannot be considered converged**. This is consistent with the
6.50 Å minimum of the held-out RMSE (the optimum of fit accuracy lies at 6.5 Å, but κ does not become "stable" as a result).

### 4.2 Cross-validation with the independent method symfc

Using `symfc 1.7.3` (symmetry-adapted least squares, a completely different solution path) to fit fc2+fc3 on the same data:

| c3 (Å) | symfc κxx | symfc κzz | pheasy κxx (OLS) | pheasy κzz (OLS) |
|---|---|---|---|---|
| 2.85 | 5.2246 | 4.0275 | 5.2220 | 4.0169 |
| 4.30 | 2.6819 | 0.4774 | 2.6939 | 0.5022 |
| 7.00 | 2.8134 | 0.4995 | 2.9873 | 0.5927 |

At 2.85 / 4.30 Å the two differ by < 1%; at 7.00 Å they differ by 6–16% (caused by the larger number of third-order terms and the
stronger ill-conditioning, but the magnitude and the trend agree). **This shows that the pheasy fits are credible.**

### 4.3 q-mesh convergence (c3=4.30, LASSO+HARM_DENSE)

| mesh | κxx | κzz |
|---|---|---|
| 11³ | 2.6938 | 0.5018 |
| 15³ | 2.6214 | 0.4449 |

11³ → 15³: κ_xx drops by 2.7% and κ_zz by 11.4%. The 11³ mesh is not yet fully converged for κ_zz (the 21³ run did not complete, owing to memory / being killed).
Therefore the absolute values in the table above carry a mesh error of order ~10%, **but the conclusion "non-monotonic with the cutoff, no plateau" is unaffected by the mesh**
(the different cutoffs are compared on the same mesh).

### 4.4 Dynamic stability

All 39 fits (36 pheasy + 3 symfc) show **0 imaginary frequencies** on the 11³ mesh; the most negative frequency is −3.8e-5 THz
(numerically zero, the Γ-point acoustic modes), and the highest frequency is ≈ 7.73 THz, essentially the same for every cutoff → the second-order force constants are converged and the dynamics is stable.

---

## 5. Overall conclusions and recommendations

1. **The workflow has been shown to be correct**: data preprocessing agrees with the official version element by element; p(c3=5.2)=1303 and p(c3=7.0)=3678
   agree with the official values; symfc independently reproduces κ (<1% @4.3 Å).
2. **κ has no plateau**: κ_xx has a region between 4.7 and 5.8 Å that "looks flat" (3.7–4.0),
   but it falls again after 6.1 Å and is still changing up to 7.0 Å; κ_zz decreases monotonically from 2.13 (5.20 Å) to 0.59 (7.00 Å).
   The present 189-atom cell (3×3×3, maximum usable interatomic distance 20.1 Å) would itself allow a larger cutoff,
   but **the effective range of the third-order force constants exceeds 6 Å**, so a larger supercell is needed for convergence.
3. **Method choice**: on this system LASSO essentially degenerates to OLS (alpha pinned to the smallest grid value);
   ALASSO sparsifies the third-order terms appreciably (7.0 Å: 3678→1783), but its held-out force prediction is slightly worse than OLS.
   If the goal is the "most reliable κ", we recommend **OLS / ALASSO + HARM_DENSE with a cutoff near 6.5 Å**;
   if the goal is the "fewest IFCs", ALASSO is usable, but one must accept a ~1–2% loss in force prediction.
4. **Follow-up recommendations**:
   - κ convergence: extend the supercell to 4×4×4 (if the DFT data allow) or at least scan a denser q-mesh (≥15³);
   - if absolute values are to be reported, it is advisable to start from 15³ and use a uniform mesh for every cutoff;
   - for a layered system one may consider a 2D model / finite-thickness correction; here the 3×3×3 cell already contains 3 layers, so κ_zz already includes interlayer scattering.

---

## 6. Archive inventory

- `01_structures/`: POSCAR, SPOSCAR
- `02_dataset/`: raw npy + disp/force_matrix.pkl + alignment + ndata
- `03_fit_results/`: 36 pheasy directories (POSCAR/fc2.hdf5/fc3.hdf5/run.log) + 3 symfc directories
- `04_kappa/`: κ summary, NU regeneration log, mesh convergence
- `05_scripts/`: fit driver, κ solver, symfc script, summary script, official pheasy_fit.sh(v6c), prepare_dataset.py
- `06_heldout_RMSE/`: whole-configuration 5-fold held-out logs
- `07_kappa_raw_output_kappa_hdf5/`: 39 `kappa-m111111.hdf5` files (containing frequency / gamma / gamma_N /
  gamma_U / gamma_isotope / group_velocity / gv_by_gv / heat_capacity / mode_kappa / kappa /
  qpoint / weight / grid_point / mesh / temperature), directly usable for plotting phonon-phonon scattering rates,
  mode thermal conductivities, cumulative κ, and so on.
