# Two-method force-constant fit of the 4x4x1 dataset (2026-11-10; conclusion: imaginary frequencies, discarded)

## Configuration

- Skill: fit-fc-thermal (AutoZT), hpc = 3090, the two methods run serially within a single job
- `FIT_METHODS = pheasy:RFE-OLS,ALASSO` (RFE-OLS is the primary method, ALASSO the method control)
- `PHEASY_BIN = pheasy-gpu`, `PHEASY_C3_CUTOFF = 6` (gen scans 4.30 / 4.71 / 5.00 / 5.60 / 6.23 / 6.45 A step by step)
- `SUPERCELL = 4 4 1` (112 atoms), `ENABLE_FC = 3`, `CUT3_SCAN = auto`, `CUT3_BOOTSTRAP = 0`
- Dataset: `dft/s4_disp` (46 frames × 112 atoms); the 0.0005 eV/A residual force of the equilibrium frame has been subtracted
- Supercell symmetry-preservation check: 4 4 1 preserves all 12 point-group operations of R-3m (passed)
- The full configuration is in `fit_config.json`, the method list in `methods.json`

## Results (both methods fit well, but both show imaginary frequencies)

| method | force RMSE | relative error | imaginary-frequency gate |
|---|---|---|---|
| RFE-OLS (primary) | 0.88 meV/A | 0.58 % | **FAIL** min = -0.1388 THz (no-NAC -0.151 THz) |
| ALASSO | 0.84 meV/A | 0.55 % | **FAIL** min = -0.1388 THz |

- Criterion: imag_policy `near_gamma_large` (near-Gamma acoustic branches with abs(q) < 0.05x1.2, |nu| > 0.10 THz)
- Minimum: branch 0 (lowest TA), q = (-0.0357, 0, 0.0357), |q| ~ 0.0505, 145 q-points with negative frequencies
- fc2 symmetry spread = 5.3e-9 THz / 12 operations (the symmetry itself is fine)
- The harmonic parts of the two methods agree (HARM is not penalized, the same set of forces is fitted) → the imaginary frequency is unrelated to the fitting method

## Dispersion diagnosis (computed directly from fc2.hdf5)

- The three acoustic branches at the Gamma point are ~ 0 (ASR satisfied, so translational invariance is not broken)
- Lowest TA branch (soft direction (-1,0,1)): t=0.01 → -0.024, 0.05 → -0.109, 0.10 → -0.149, 0.15 → +0.139 THz
- The same along the in-plane direction: 0.01 → -0.012, 0.10 → -0.096, 0.20 → +0.094 THz
→ this is a soft TA branch over the whole range |q| ~ 0.03-0.12 (slightly negative effective elastic constant), not single-point numerical noise

## Conclusions and handling (corrected 2026-11-10)

**This is not a supercell-size problem**: the perpendicular widths of 4x4x1 in all three directions are all >= 13.2 A (safe cutoff 6.50 A), wider than the old 3x3x3
rhombohedral supercell at about 10.5 A (safe cutoff 5.17 A); Mn2In2Se5 with the same 4x4x1 setup is stable. There are two real causes:

1. **Rotational invariance / the Huang equilibrium condition was not imposed** (`PHEASY_RASR=auto` resolves to none for bulk).
   The three acoustic branches at the Gamma point being 0 only shows that the translational sum rule is satisfied; a negative TA branch at |q| ~ 0.03-0.12 is exactly how a missing constraint of this kind
   typically manifests itself.
2. **about -1.23 kBar of isotropic residual stress in the equilibrium frame (disp-00000)** (diagonal -1.246 / -1.246 / -1.207 kBar,
   deviatoric part < 0.03 kBar): the Pulay stress difference between the S1 relaxation with ENCUT=540 and the S4 force evaluation with 410 softens the long-wavelength acoustic branches
   (it cannot explain -0.14 THz on its own, but it contributes).

**Handling**: this batch of 4x4x1 fit results is discarded, but the data are kept (the only force dataset for this material with the same DFT settings as Mn2In2Se5).
The fit is redone on the same data, changing only `PHEASY_RASR = H` (Huang constraint), with the method and all other parameters unchanged.
If near-Gamma negative frequencies below -0.05 THz remain after RASR=H, the next step is to recompute the equilibrium frame (a single frame) with ENCUT=540
before deciding whether to change the S4 settings.

Only fit_config.json / methods.json are kept in this directory: the fc2.hdf5 / fc3.hdf5 / cutoff_scan/ (25 MB) /
phonon_summary.json on the cluster were overwritten by the subsequent refit in the same directory and were not archived.
(Note: the 3x3x3 / 189-atom dataset `datasets/MnIn2Se4_raw` is a test dataset shipped with pheasy whose DFT settings are unknown;
it cannot be used for results on this material and serves only as a reference.)
