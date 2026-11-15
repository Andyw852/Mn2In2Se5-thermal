# MnIn₂Se₄ — reproducible calculation index

Every calculation of this work, with its **inputs and outputs committed to this
repository**. Paths are relative to this directory (`test/MnIn2Se4/`).
All energies in eV, lengths in Å, κ in W m⁻¹ K⁻¹, temperature 300 K unless stated.

---

## 1. FC2 in a 4×4×2 supercell (the core of this rework)

| | |
|---|---|
| **Why** | The old 4×4×1 FC2 gives an interlayer transverse sound velocity of 322.9 m/s vs 1338.3 m/s from the DFT elastic constants (4.1× too soft) → κ_zz badly underestimated |
| **Method** | phonopy symmetry-reduced finite displacement, 7 frames, 0.02 Å, 224 atoms |
| **Inputs** | `fc2_442_phonon/frames/disp-0000{1..7}/` — INCAR, KPOINTS, POSCAR, CONTCAR for each frame |
| | `fc2_442_phonon/disp_plan.json` — displacement plan (fd_distance 0.02, n_disp 7) |
| | `fc2_442_phonon/POSCAR` — primitive cell |
| | `fc2_442_phonon/kl_params.txt` — parameter snapshot |
| **Outputs** | `fc2_442_phonon/fc2.hdf5` — second-order force constants (224×224×3×3) |
| | `fc2_442_phonon/FORCE_CONSTANTS` — same, phonopy text format |
| | `fc2_442_phonon/band-dft-cpu.yaml` — phonon dispersion (499 q-points × 21 branches) |
| | `fc2_442_phonon/phonon_summary.json` — stability verdict (`stable=true`, min_freq +0.000719 THz) |
| | `fc2_442_phonon/phonon_fit.log` — fit log |
| | `fc2_442_phonon/phonon_band_plot/phonon_band_{full,lowfreq}.png` — figures |
| **Key number** | min_freq **+0.000719 THz**, 0 imaginary q-points (4×4×1 gave −0.006 THz) |
| **Run** | `python3 test/MnIn2Se4/kappa_driver_442.py` (see §3) |

## 2. Interlayer sound velocity and C₄₄ (the diagnostic that found the problem)

| | |
|---|---|
| **Script** | `fc2_442_phonon/sound_velocity.py` — its header documents the two pitfalls |
| **Inputs** | `fc2_442_phonon/fc2.hdf5`, `fc2_442_phonon/POSCAR` |
| **Output** | printed table; values stored in `../REPORT_MnIn2Se4_20261115.md` §4 |
| **Result** | TA₁ = TA₂ = **1348.3 m/s** (DFT 1338.3, +0.75 %), LA = 3118.6 (−4.14 %), C₄₄ = **9.543 GPa** (DFT 9.402, +1.50 %) |
| **Pitfall 1** | phono3py uses the *selected primitive* reciprocal basis: the layer normal is `(t,t,t)`, not `(0,0,t)`. Using the POSCAR b₃ gives a spurious 4179/6256 m/s (error factor exactly 3.760) |
| **Pitfall 2** | Γ-point residual f₀ = 0.000719 THz → fit `f² = f₀² + (v q / 2π / 100)²`, not ω/q |

## 3. κ versus FC2 supercell (4×4×1 → 4×4×2)

| | |
|---|---|
| **Script** | `kappa_driver_442.py` |
| **Inputs** | `fit_cut3_6p45/fc2.hdf5` (4×4×1, old) or `fc2_442_phonon/fc2.hdf5` (4×4×2, new) |
| | `fit_cut3_6p45/fc3.hdf5` (FC3, 4×4×1, cutoff 6.45 Å) |
| | `fit/POSCAR` |
| **Outputs** | `kappa_results_442.json`, `kappa_results_mn441.json` |
| **Result (LBTE, mesh 11)** | in-plane 3.9054 → **4.4075** (+12.9 %); κ_zz 0.3629 → **1.2519 (+245 %)** |
| **Method note** | phono3py accepts `phonon_supercell_matrix=4×4×2` with `supercell_matrix=4×4×1`; **ShengBTE cannot** (it requires both on the same supercell) |

## 4. 4×4×1 versus 4×4×2 under one convention

| | |
|---|---|
| **Script** | `compare_441_vs_442.py` |
| **Inputs** | `fit/POSCAR`, `fit_cut3_6p45/fc2.hdf5`, `fc2_442_phonon/fc2.hdf5` |
| **Output** | `compare_441_vs_442_results.json` |
| **Result** | 4×4×1: TA 322.9 (−75.87 %); 4×4×2: TA 1348.3 (+0.75 %) |

## 5. FC3 cutoff scan (5.60 / 6.23 / 6.45 Å)

| | |
|---|---|
| **Ceiling** | a 4×4×1 cell caps the FC3 cutoff at **6.507 Å**; largest contained shell 6.3557 Å, next shell 6.5311 Å is out of reach — **6.45 Å is the physical maximum** |
| **Inputs** | `fc3_cutoff_variants/fc3_cut5p60.hdf5`, `fc3_cut6p23.hdf5`, `fc3_cut6p45.hdf5` |
| | `fc2_442_phonon/fc2.hdf5` (FC2 held fixed) |
| **Scripts** | `kappa_cutscan.py` (RTA), `kappa_lbte_scan.py` (LBTE) |
| **Outputs** | `kappa_cutoff_scan_results.json`, `kappa_cutoff_convergence.png`, `KAPPA_CUTOFF_SCAN.md` |
| **Result (LBTE, mesh 14)** | 5.60: 5.6473 / 2.7388 · 6.23: 4.7842 / 1.2241 · **6.45: 4.4075 / 1.2519** |
| **Verdict** | **not converged**: in-plane still −7.9 % in the last available step; ceiling reached |

## 6. κ versus q-mesh

| | |
|---|---|
| **Scripts** | `kappa_driver_442.py` (`--mesh`), `kappa_run_18_lbte.py` |
| **Inputs** | `fc2_442_phonon/fc2.hdf5`, `fit_cut3_6p45/fc3.hdf5`, `fit/POSCAR` |
| **Outputs** | `kappa_mesh_scan_results.json`, `kappa_mesh_convergence.png`, `KAPPA_MESH_SCAN.md` |
| **Result (RTA)** | 11³ 3.5375 / 1.2495 · 14³ 3.6474 / 1.0934 · 18³ 3.6945 / **0.9435** |
| **Verdict** | in-plane converged (+3.1 % → +1.3 %); **intrinsic κ_zz does not converge** (−13 % per refinement) |

## 7. Finite-grain κ (Casimir boundary scattering)

| | |
|---|---|
| **Input** | same force constants as §6 |
| **Output** | `kappa_finite_grain.json` |
| **Result (mesh 18, RTA)** | ∞: 3.6945 / 0.9435 · 10 μm: 3.6810 / 0.8849 · **1 μm: 3.5729 / 0.6877** |
| **Meaning** | in-plane is boundary-insensitive (−3.3 %); κ_zz is set by the grain size, which is why an "intrinsic" κ_zz cannot be quoted |
| **Unit** | `boundary_mfp` is in **micrometres** in phono3py |

## 8. Wigner coherence term

| | |
|---|---|
| **Script** | `kappa_driver_442.py` (`conductivity_type="wigner"`), results in `kappa_wigner_mn441.json` |
| **Result (Wigner-RTA, mesh 14, 6.45 Å)** | total **3.7492 / 1.1595**; population 3.6584 / 1.0934; coherence 0.0908 / 0.0661 (**+2.5 % / +6.0 %**) |
| **Comparison, Mn₂In₂Se₅** | coherence contributes +4.5 % / **+12.4 %** — the lower κ_zz, the larger the coherence fraction |

## 9. FC3 → ShengBTE export audit

| | |
|---|---|
| **Writers** | `../../analysis/shengbte_crosscheck/write_fc3_sheng.py` (complete), `write_fc2_sheng.py` |
| **Inputs** | `fit_cut3_6p45/fc3.hdf5`, `phono3py_disp.yaml` |
| **Outputs** | `shengbte_corrected/CONTROL_fixed` (Mn at (1,0,1)), `shengbte_corrected/FORCE_CONSTANTS_3RD_complete` (5586 triplets), `shengbte_corrected/triplet_audit.json` |
| **Findings** | fc3.hdf5 has **5586** non-zero triplets; the old `shengbte_export/FORCE_CONSTANTS_3RD` kept only **3276 (−41.4 %)**; the production `kappa_shengbte_20x20x5/CONTROL` places Mn at **(0,0,0)** instead of **(1,0,1)** (AutoZT [0151]/[0154]) |
| **Consequence** | all ShengBTE κ for MnIn₂Se₄ produced before this fix are invalid |
| **Why no ShengBTE rerun** | with one unit cell along c, the 4×4×1 FC2 has degenerate minimum images (atom 0–2: both images at 8.1479 Å), so the R-vectors cannot be assigned uniquely |

## 10. Legacy ShengBTE and phono3py runs (superseded)

| Path | Contents |
|---|---|
| `fit/` | pheasy FC2+FC3 fit (4×4×1), dataset displacements/forces, fit config and manifest |
| `fit_cut3_6p45/` | 4×4×1 FC2 + FC3 at cutoff 6.45 Å |
| `kappa_shengbte_11x11x11/`, `kappa_shengbte_20x20x5/` | ShengBTE runs — **κ invalid** (see §9) |
| `kappa_phono3py_lbte_11x11x11/` | phono3py LBTE output with the *old* 4×4×1 FC2 |
| `shengbte_export/` | the defective FORCE_CONSTANTS_3RD / CONTROL pair |

> The 3×3×3 test dataset that originally stood here is documented under
> `../../MnIn2Se4/` and is superseded by this 4×4×2 work.

---

## Directory tree

```
test/MnIn2Se4/
├── README.md                        # this index
├── fc2_442_phonon/                  # §1 §2  FC2 4x4x2 + sound velocity
│   ├── frames/disp-0000{1..7}/      #   inputs: INCAR KPOINTS POSCAR CONTCAR
│   ├── fc2.hdf5  FORCE_CONSTANTS    #   outputs
│   ├── band-dft-cpu.yaml            #   phonon dispersion
│   ├── phonon_summary.json  phonon_fit.log
│   ├── phonon_band_plot/*.png
│   ├── disp_plan.json  POSCAR  kl_params.txt
│   └── sound_velocity.py            #   §2
├── fc3_cutoff_variants/             # §5  inputs for the cutoff scan
│   ├── fc3_cut5p60.hdf5
│   ├── fc3_cut6p23.hdf5
│   └── fc3_cut6p45.hdf5
├── shengbte_corrected/              # §9  fixed ShengBTE inputs + audit
│   ├── CONTROL_fixed
│   ├── FORCE_CONSTANTS_3RD_complete
│   └── triplet_audit.json
├── kappa_driver_442.py              # §3 §8  main kappa driver (RTA / LBTE / Wigner)
├── kappa_run_18_lbte.py             # §6     18^3 LBTE
├── kappa_lbte_scan.py               # §5     LBTE cutoff scan
├── kappa_cutscan.py                 # §5     RTA cutoff scan
├── compare_441_vs_442.py            # §4
├── kappa_results_442.json           # §3
├── kappa_results_mn441.json         # aggregate snapshot
├── kappa_cutoff_scan_results.json   # §5
├── kappa_mesh_scan_results.json     # §6
├── kappa_finite_grain.json          # §7
├── kappa_wigner_mn441.json          # §8
├── kappa_cutoff_convergence.png     # §5 figure
├── kappa_mesh_convergence.png       # §6 figure
├── KAPPA_442_README.md              # §3 narrative
├── KAPPA_CUTOFF_SCAN.md             # §5 narrative
├── KAPPA_MESH_SCAN.md               # §6 narrative
├── fit/                             # §10 legacy 4x4x1 fit
├── fit_cut3_6p45/                   # §10 legacy force constants
├── kappa_shengbte_*/                # §10 legacy ShengBTE (invalid kappa)
├── kappa_phono3py_lbte_11x11x11/    # §10 legacy phono3py (old FC2)
└── shengbte_export/                 # §10 defective export
```

---

## Environment

phonopy 2.47.1, phono3py 3.24.0, symfc 1.6.0, Python 3.11.
`sound_velocity.py` and `kappa_*.py` only need numpy + phonopy + phono3py.

## Caveats

- FC3 is still the **4×4×1, 6.45 Å** fit; the 4×4×1 cell caps the cutoff at
  6.507 Å, so the FC3 cutoff **cannot be converged** without a new 4×4×2 FC3
  dataset (~2 days of DFT).
- The **intrinsic** κ_zz is not a well-defined number for this material; quote a
  finite-grain value instead (§7).
- ShengBTE cannot be used with a 4×4×2 FC2 and a 4×4×1 FC3, nor can the
  4×4×1 FC2 be exported cleanly (§9).
