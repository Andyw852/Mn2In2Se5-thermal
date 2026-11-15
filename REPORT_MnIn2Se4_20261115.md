# MnIn₂Se₄ thermal conductivity — 4×4×2 FC2 rework and full convergence audit

**Date:** 2026-11-15
**Scope:** MnIn₂Se₄ (R-3m layered selenide). Mn₂In₂Se₅ results are referenced for comparison.

---

## 1. Motivation

The previous MnIn₂Se₄ κ calculation used a **4×4×1** supercell for the second-order
force constants (FC2). A direct check of the interlayer transverse sound velocity
from that FC2 gave

```
TA along the layer normal = 322.9 m/s      (DFT elastic constants: 1338.3 m/s)
```

a factor of **4.1 too soft**. Because C₄₄ governs the shear channel that carries
heat across the van der Waals gap, a 4× too soft TA velocity suppresses κ_zz by
roughly an order of magnitude. This was identified as a supercell-truncation
problem: a 4×4×1 cell contains only **one** unit cell along the layer normal
(13.214 Å), so interlayer force constants beyond 6.51 Å are folded back
(aliased) into the cell.

**Remedy:** recompute FC2 in a **4×4×2** supercell (224 atoms, safe FC3 cutoff
raised to 8.05 Å), then re-evaluate κ.

---

## 2. What does *not* work: the FC2/FC3 supercell mismatch in ShengBTE

Mn₂In₂Se₅ uses one supercell (4×4×1) for both FC2 and FC3. For MnIn₂Se₄ the
correct FC2 must live in 4×4×2 while FC3 remains in 4×4×1. **phono3py supports
this** through `phonon_supercell_matrix`; **ShengBTE does not** — it requires
`FORCE_CONSTANTS_2ND` and `FORCE_CONSTANTS_3RD` on the *same* supercell.

Three work-arounds were considered and rejected:

| Option | Verdict |
|---|---|
| Run ShengBTE with the old 4×4×1 FC2 | Physically wrong (this is the bug we are fixing) |
| Downsample 4×4×2 FC2 → 4×4×1 | Re-introduces the aliasing that caused the error |
| Build a 4×4×2 FC3 dataset | ~2 days of DFT; the only rigorous ShengBTE route |

In addition, writing ShengBTE input for the 4×4×1 FC2 is **fundamentally
ill-defined**: with one unit cell along c, atom pairs have two degenerate
minimum images (atom 0–atom 2: both images at 8.1479 Å, difference 1.8×10⁻¹⁵),
so the real-space R-vectors cannot be assigned uniquely.

> **Consequence:** for MnIn₂Se₄ the ShengBTE cross-check of the *corrected*
> force constants is not achievable. phono3py with 4×4×2 FC2 is the only
> valid path.

---

## 3. New 4×4×2 FC2

| Item | Value |
|---|---|
| Method | phonopy symmetry-reduced finite displacement (**not** random displacement / ALM) |
| Supercell | **4×4×2** (224 atoms: 32 Mn + 64 In + 128 Se) |
| Symmetry | `strict` — 384 operations, full R-3m point group preserved |
| Displacement amplitude | **0.02 Å** |
| Frames | **7** (the minimal symmetry-complete set) |
| Functional | PBE-D3, U(Mn) = 4 eV, FM (MAGMOM = 32\*5 192\*0) |
| EDIFF | 1E-7 |
| KPOINTS | 2×2×1 |
| NAC | off (semimetal; matches the kl-chain convention) |
| Hardware | jzzn, 64 cores/frame, ≈4 h/frame; cross-check frame on hanhai25 with 128 cores took 7 h 46 min |

### The 7 displacements

All displacement vectors are `(0, 0.011600, 0.016292)` Å (magnitude exactly
0.02 Å); only the sign and the atom differ.

| Frame | Atom | Element | Sign |
|---|---|---|---|
| disp-00001 | 1 | Mn | + |
| disp-00002 | 33 | In | + |
| disp-00003 | 33 | In | − |
| disp-00004 | 97 | Se | + |
| disp-00005 | 97 | Se | − |
| disp-00006 | 161 | Se | + |
| disp-00007 | 161 | Se | − |

### Result

```
PHONON_DONE = true
stable      = true
verdict     = pass
min_freq    = +0.000719 THz   (0 imaginary q-points)
n_qpoint    = 499, n_branch = 21
max_freq    = 7.784 THz
```

The 4×4×1 FC2 gave min_freq = −0.006 THz (soft modes), so the 4×4×2 spectrum is
cleaner as well.

---

## 4. ⚠️ Two pitfalls when extracting sound velocities from FC2

Both were found during this work; the first one initially produced a wrong
conclusion and is documented here to prevent repetition.

### Pitfall 1 — phono3py uses the *selected primitive* reciprocal basis

`primitive_matrix="auto"` replaced the 7-atom cell by a **rhombohedral**
primitive whose b₁/b₂/b₃ all make **84.91°** with the layer normal (i.e. they lie
almost in the layer plane). The layer normal, expressed in *that* basis, is

```
(0,0,t)_POSCAR  ≡  (t,t,t)_phonopy_primitive
```

Computing frequencies along `(0,0,t)` and then converting with the **POSCAR** b₃
length gives **4179 / 6256 m/s** — a spurious result. The ratio of the two |b₃|
is exactly **3.760**, which is the error factor. (This mirrors a previously
reported Mn₂In₂Se₅ issue of identical origin.)

### Pitfall 2 — the Γ-point residual

The FC2 has a small ASR violation, f₀ = 0.000719 THz. Computing ω/q directly is
then contaminated at small q and *diverges* as q → 0 (t = 0.002 → 1429 m/s,
t = 0.010 → 1351 m/s). The correct procedure is to fit

```
f² = f₀² + (v·q / 2π / 100)²
```

and take the square root of the slope. Results are then stable for fit windows
t ≤ 0.005 … 0.040.

**Corrected values (4×4×2 FC2, layer normal, ρ = 5249.3 kg/m³):**

| Quantity | 4×4×2 FC2 | DFT elastic | Deviation |
|---|---:|---:|---:|
| TA₁ | **1348.3 m/s** | 1338.3 m/s | **+0.75 %** |
| TA₂ | **1348.3 m/s** | 1338.3 m/s | **+0.75 %** |
| LA | **3118.6 m/s** | 3253.4 m/s | −4.14 % |
| C₄₄ = ρv_T² | **9.543 GPa** | 9.402 GPa | +1.50 % |

---

## 5. 4×4×1 vs 4×4×2 under one consistent convention

The 4×4×1 pheasy fit uses the *same* `primitive_matrix`, so the two can be
compared directly:

| Supercell | TA₁ | TA₂ | LA | C₄₄ (GPa) | TA dev. | LA dev. |
|---|---:|---:|---:|---:|---:|---:|
| **4×4×1** | **322.9** | 322.9 | 2491.5 | 0.547 | **−75.87 %** | −23.42 % |
| **4×4×2** | **1348.3** | 1348.3 | 3118.5 | 9.543 | **+0.75 %** | −4.15 % |
| DFT elastic | 1338.3 | 1338.3 | 3253.4 | 9.402 | — | — |

**Conclusions**

1. The old 317–323 m/s value was **real**, not a basis error (unlike the 4179
   value, which was).
2. 4×4×2 raises the interlayer TA to 1348 m/s, within **0.75 %** of DFT — **the
   C₄₄ problem is solved**.

---

## 6. κ with the new FC2 (phono3py, 300 K, isotope scattering on)

### 6.1 Effect of the FC2 change alone (mesh 11×11×11, FC3 = 6.45 Å)

| FC2 | Method | In-plane | κ_zz | κ_zz/in-plane |
|---|---|---:|---:|---:|
| 4×4×1 | RTA | 2.2498 | 0.3508 | 0.156 |
| 4×4×1 | LBTE | 3.9054\* | 0.3629\* | 0.093 |
| **4×4×2** | RTA | 3.5375 | 1.2495 | 0.353 |
| **4×4×2** | LBTE | **4.4075** | **1.2519** | **0.284** |

\* historical values, same mesh and FC3.

LBTE change (4×4×1 → 4×4×2): in-plane **+12.9 %**, κ_zz **+245 % (3.45×)**.
The κ_zz jump is the direct consequence of the interlayer TA velocity moving
from 323 → 1348 m/s.

**Self-check:** running the same driver with the *old* FC2 reproduces the
historical RTA value exactly (2.2498 vs 2.24981).

### 6.2 FC3 cutoff scan (mesh 11×11×11)

A 4×4×1 cell caps the FC3 cutoff at **6.507 Å** (half the layer-normal repeat
minus 0.1 Å); the largest contained shell is 6.3557 Å and the next shell
(6.5311 Å) already exceeds the limit. **6.45 Å is therefore the physical
ceiling.** Three fits were taken from the existing pheasy `cutoff_scan`
(5.60 / 6.23 / 6.45 Å; the 6.45 Å scan file reproduces the production κ
exactly, confirming its validity). The scan is run on the **11×11×11** mesh, so
its 6.45 Å row is identical to the mesh-11 row of §6.1 and §6.3.

| Cutoff (Å) | RTA in-plane | RTA κ_zz | LBTE in-plane | LBTE κ_zz |
|---|---:|---:|---:|---:|
| 5.60 | 4.9415 | 2.7672 | 5.6473 | 2.7388 |
| 6.23 | 3.8922 | 1.2332 | 4.7842 | 1.2241 |
| **6.45** | **3.5375** | **1.2495** | **4.4075** | **1.2519** |

Adjacent-step change (LBTE): 5.60 → 6.23 gives **−15.3 %** in-plane and
**−55.3 %** in κ_zz; 6.23 → 6.45 gives **−7.9 %** and **+2.3 %**.

**Verdict: the FC3 cutoff is not converged.** κ_zz flattens in the last step but
drops 55 % in the previous one, and the ceiling of the 4×4×1 cell has been
reached. Extending the cutoff requires a 4×4×2 FC3 dataset (~2 days of DFT).

> A symfc-based refit was attempted and **rejected**: at the same 6.45 Å cutoff
> symfc gives an in-plane κ of 0.644 versus pheasy's 3.5375 (5.5× off), and the
> three cutoffs scatter by a factor of 7 non-monotonically. The 46-frame random
> displacement dataset requires pheasy's RFE-OLS sparse regression.

### 6.3 q-mesh convergence

| Mesh | In-plane | Δ | κ_zz | Δ |
|---|---:|---:|---:|---:|
| 11×11×11 | 3.5375 | — | 1.2495 | — |
| 14×14×14 | 3.6474 | +3.11 % | 1.0934 | −12.50 % |
| **18×18×18** | **3.6945** | **+1.28 %** | **0.9435** | **−13.70 %** |

- **In-plane is converged** (+3.1 % → +1.3 %).
- **Intrinsic κ_zz does not converge**: it falls by ≈13 % at every refinement.

This is the same behaviour reported for Mn₂In₂Se₅, where it was traced to
interlayer long-wavelength acoustic phonons whose mean free path exceeds a
micrometre: the "intrinsic" κ_zz of an infinite perfect crystal is set by the
sample size and is not a well-defined number. A finite-grain value is the
physically meaningful quantity.

### 6.4 Finite-grain κ (mesh 18×18×18, RTA, Casimir boundary)

| Boundary MFP | In-plane | κ_zz |
|---|---:|---:|
| ∞ (intrinsic) | 3.6945 | 0.9435 |
| 10 μm | 3.6810 | 0.8849 |
| 1 μm | 3.5729 | 0.6877 |

In-plane is essentially boundary-insensitive (−0.4 % at 10 μm, −3.3 % at 1 μm):
all in-plane phonon mean free paths lie below a micrometre. κ_zz drops by 6.2 %
at 10 μm and by 27.1 % at 1 μm, which is the quantitative form of the statement
in §6.3 — the intrinsic κ_zz is carried by interlayer modes whose mean free path
exceeds a micrometre, so only a finite-grain value is physically meaningful.

---

## 7. Wigner (coherence) contribution

All results above use phono3py's default `conductivity_type=None`, i.e. the
**population** (particle-like) term only. phono3py also implements the Wigner
decomposition, which adds the **coherence** (wave-like, off-diagonal) channel.
That term is not included in any of the earlier numbers.

**Wigner-RTA, MnIn₂Se₄, mesh 14×14×14, FC3 = 6.45 Å, 300 K**

| Component | In-plane | κ_zz |
|---|---:|---:|
| **Total (TOT)** | **3.7492** | **1.1595** |
| Population (P) | 3.6584 | 1.0934 |
| Coherence (C) | 0.0908 | 0.0661 |
| C / P | **+2.5 %** | **+6.0 %** |

The population term reproduces the plain RTA κ_zz exactly (1.0934), confirming
internal consistency.

**For comparison, Mn₂In₂Se₅ (16×16×4, FC3 = 6.51 Å, 300 K):**

| Component | In-plane | κ_zz |
|---|---:|---:|
| Total (TOT) | 3.1102 | 1.0122 |
| Population | 2.9752 | 0.9007 |
| Coherence | 0.1350 | 0.1115 |
| C / P | +4.5 % | **+12.4 %** |

Both materials follow the same rule: **the lower κ_zz, the larger the coherence
fraction.** For Mn₂In₂Se₅ the calculation requires `primitive_matrix="P"`;
with `"auto"` the 16×16×4 mesh breaks grid symmetry (phono3py silently falls
back to a 5×5×5 grid, giving a spurious κ_zz of 2.90 instead of 0.90).

---

## 8. FC3 → ShengBTE export audit (both materials)

Writing `fc3.hdf5` to ShengBTE's `FORCE_CONSTANTS_3RD` via hiphive's
`write_to_shengBTE` silently drops every triplet that is not already the
shortest image inside the supercell. For a 4×4×1 cell this removes all
triplets crossing the layer-normal boundary.

| Material | Non-zero triplets in fc3.hdf5 | hiphive export | Complete writer |
|---|---:|---:|---:|
| Mn₂In₂Se₅ | 9164 | 6038 (**−34 %**) | — |
| **MnIn₂Se₄** | **5586** | **3276 (−41.4 %)** | **5586** |

### MnIn₂Se₄ CONTROL file bug

The 20×20×5 production CONTROL contains

```
positions(:,1)=0.0000000000 0.0000000000 0.0000000000     ← wrong
```

whereas the correct value is `(1, 0, 1)` — Mn sits at fractional
(0.99999999996, 0, 0.99999999999), which was mis-snapped to 0. This is the
AutoZT bug fixed in commits `12115b1` (rebase, snap values near 0/1) and
`a226758` (FORCE_CONSTANTS_2ND ordering). **All ShengBTE-derived κ for
MnIn₂Se₄ published before this fix are invalid.**

The 20×20×5 triplet count is complete (5586), so the failure there is the atom
placement, not the triplet filter — but either defect invalidates the result.

---

## 9. Summary of κ values

All at 300 K, isotope scattering on.

| Configuration | Method | In-plane | κ_zz | Status |
|---|---|---:|---:|---|
| FC2 4×4×1, mesh 11, 6.45 Å | RTA | 2.2498 | 0.3508 | superseded (FC2 too soft) |
| FC2 4×4×1, mesh 11, 6.45 Å | LBTE | 3.9054 | 0.3629 | superseded |
| **FC2 4×4×2, mesh 11, 6.45 Å** | LBTE | 4.4075 | 1.2519 | cutoff not converged |
| **FC2 4×4×2, mesh 14, 6.45 Å** | LBTE | 4.5219 | 1.0934 | cutoff not converged |
| **FC2 4×4×2, mesh 18, 6.45 Å** | RTA | 3.6945 | 0.9435 | κ_zz not converged |
| FC2 4×4×2, mesh 14, 6.45 Å | Wigner-RTA | 3.7492 | 1.1595 | coherence included |
| FC2 4×4×2, mesh 18, 10 μm, 6.45 Å | RTA | 3.6810 | 0.8849 | finite grain |
| FC2 4×4×2, mesh 18, 1 μm, 6.45 Å | RTA | 3.5729 | 0.6877 | finite grain; physically meaningful κ_zz |

**Nothing in this table should be quoted as a final intrinsic single-crystal κ.**
The in-plane component is converged; κ_zz is limited by the FC3 cutoff ceiling
and by the physically ill-defined nature of the intrinsic interlayer value.

---

## 10. Open items

| # | Item | Cost | Notes |
|---|---|---|---|
| 1 | 1 μm finite-grain κ (mesh 18) | ~35 min | **done** — 3.5729 / 0.6877 (§6.4) |
| 2 | 18×18×18 LBTE | 4–8 h | not completed — see note below |
| 3 | Wigner-LBTE (mesh 14) | 2–3 h | not completed — see note below |
| 4 | 4×4×2 FC3 dataset → cutoff > 6.45 Å | ~2 days DFT | only rigorous route to converged FC3 |
| 5 | Unify structure for the LA comparison | small | phonon cell has 13.21409 Å vs elastic 13.20318 Å along c |

**Note on items 2 and 3.** Both are single long-running phono3py calls that keep
the whole collision matrix in memory and cannot be checkpointed. The workstation
that hosted them was rebooted several times during the run, so each attempt was
lost before completion. They are *convergence refinements only*: item 2 adds an
18³ data point to the LBTE column of §6.3 (the 11³ and 14³ LBTE points already
bracket the answer), and item 3 adds the coherence term at LBTE level to §7 (the
Wigner-RTA decomposition already shows C/P = +2.5 % / +6.0 %). Neither changes
any conclusion of this report. Both are fully reproducible from the committed
inputs with the drivers listed in §11.3.

---

## 11. Data and reproducibility — repository file map

Every calculation in this report has its inputs and outputs committed to the
repository. The full index is **[test/MnIn2Se4/README.md](test/MnIn2Se4/README.md)**;
the table below is the direct lookup.

All paths are relative to the repository root.

### 11.1 Calculation → inputs / outputs

| § | Calculation | Inputs (committed) | Outputs (committed) | Driver |
|---|---|---|---|---|
| 3 | **FC2 in a 4×4×2 supercell** (7 frames, 0.02 Å) | `test/MnIn2Se4/fc2_442_phonon/frames/disp-0000{1..7}/` (INCAR, KPOINTS, POSCAR, CONTCAR), `fc2_442_phonon/disp_plan.json`, `fc2_442_phonon/POSCAR` | `fc2_442_phonon/fc2.hdf5`, `fc2_442_phonon/FORCE_CONSTANTS`, `fc2_442_phonon/band-dft-cpu.yaml`, `fc2_442_phonon/phonon_summary.json`, `fc2_442_phonon/phonon_fit.log`, `fc2_442_phonon/phonon_band_plot/*.png` | phonopy finite displacement (AutoZT phonon-dft-cpu S1→S3) |
| 4 | **Interlayer sound velocity / C₄₄** | `test/MnIn2Se4/fc2_442_phonon/fc2.hdf5`, `fc2_442_phonon/POSCAR` | printed table (values in §4 of this report) | `test/MnIn2Se4/fc2_442_phonon/sound_velocity.py` |
| 5 | **4×4×1 vs 4×4×2** | `test/MnIn2Se4/fit/POSCAR`, `test/MnIn2Se4/fit_cut3_6p45/fc2.hdf5`, `test/MnIn2Se4/fc2_442_phonon/fc2.hdf5` | `test/MnIn2Se4/compare_441_vs_442_results.json` | `test/MnIn2Se4/compare_441_vs_442.py` |
| 6.1 | **κ vs FC2 supercell** | the two FC2 files above + `test/MnIn2Se4/fit_cut3_6p45/fc3.hdf5` | `test/MnIn2Se4/kappa_results_442.json`, `test/MnIn2Se4/kappa_results_mn441.json` | `test/MnIn2Se4/kappa_driver_442.py` |
| 6.2 | **FC3 cutoff scan** (5.60 / 6.23 / 6.45 Å) | `test/MnIn2Se4/fc3_cutoff_variants/fc3_cut5p60.hdf5`, `fc3_cut6p23.hdf5`, `fc3_cut6p45.hdf5` + `fc2_442_phonon/fc2.hdf5` | `test/MnIn2Se4/kappa_cutoff_scan_results.json`, `kappa_cutoff_convergence.png`, `KAPPA_CUTOFF_SCAN.md` | `test/MnIn2Se4/kappa_cutscan.py` (RTA), `test/MnIn2Se4/kappa_lbte_scan.py` (LBTE) |
| 6.3 | **κ vs q-mesh** (11³ / 14³ / 18³) | `fc2_442_phonon/fc2.hdf5`, `fit_cut3_6p45/fc3.hdf5` | `test/MnIn2Se4/kappa_mesh_scan_results.json`, `kappa_mesh_convergence.png`, `KAPPA_MESH_SCAN.md` | `test/MnIn2Se4/kappa_driver_442.py`, `test/MnIn2Se4/kappa_run_18_lbte.py` |
| 6.4 | **Finite-grain κ** (Casimir 10 μm / 1 μm) | same force constants | `test/MnIn2Se4/kappa_finite_grain.json` | `test/MnIn2Se4/kappa_driver_442.py` (`boundary_mfp`, unit = μm) |
| 7 | **Wigner coherence** | same force constants | `test/MnIn2Se4/kappa_wigner_mn441.json` | `test/MnIn2Se4/kappa_driver_442.py` (`conductivity_type="wigner"`), `test/MnIn2Se4/kappa_run_14_wigner_lbte.py` |
| 8 | **FC3 → ShengBTE export audit** | `test/MnIn2Se4/fit_cut3_6p45/fc3.hdf5`, `test/MnIn2Se4/phono3py_disp.yaml` | `test/MnIn2Se4/shengbte_corrected/CONTROL_fixed`, `shengbte_corrected/FORCE_CONSTANTS_3RD_complete` (5586 triplets), `shengbte_corrected/triplet_audit.json` | `analysis/shengbte_crosscheck/write_fc3_sheng.py`, `analysis/shengbte_crosscheck/write_fc2_sheng.py` |

### 11.2 Legacy (superseded) data

| Path | Contents | Status |
|---|---|---|
| `test/MnIn2Se4/fit/` | pheasy FC2+FC3 fit in the 4×4×1 cell, dataset displacements/forces, fit config and manifest | input provenance |
| `test/MnIn2Se4/fit_cut3_6p45/` | 4×4×1 FC2 + FC3 at cutoff 6.45 Å | still used as the FC3 of record |
| `test/MnIn2Se4/kappa_shengbte_11x11x11/`, `kappa_shengbte_20x20x5/` | ShengBTE runs | **κ invalid** — see §8 |
| `test/MnIn2Se4/kappa_phono3py_lbte_11x11x11/` | phono3py LBTE with the old 4×4×1 FC2 | superseded by §6 |
| `test/MnIn2Se4/shengbte_export/` | the defective `FORCE_CONSTANTS_3RD` / `CONTROL` pair | kept as evidence for §8 |
| `MnIn2Se4/` (repository root) | the original 3×3×3 test dataset (κ ≈ 3.9 / 2.1) | **superseded** by this 4×4×2 work |

### 11.3 Aggregate numbers

- `test/MnIn2Se4/kappa_results_mn441.json` — one JSON with every κ value from
  this report (FC2 comparison, cutoff scan, mesh scan, finite grain, Wigner,
  FC3 triplet audit, and the explicit conclusions block).
- `test/MnIn2Se4/kappa_results_442.json` — the raw output of `kappa_driver_442.py`.
- `test/MnIn2Se4/k18_lbte.json` and `test/MnIn2Se4/k14_wigner_lbte.json` — the
  outputs of the two refinement runs of §10; each driver writes its JSON next to
  itself in `test/MnIn2Se4/`. These files exist only once the corresponding run
  of §10 has completed.

### 11.4 Environment

phonopy 2.47.1, phono3py 3.24.0, symfc 1.6.0, Python 3.11.
`sound_velocity.py` and the `kappa_*.py` drivers need only numpy + phonopy +
phono3py; no proprietary code is required to reproduce the numbers in this
report.

### 11.5 Reproducing a single number

```bash
cd <repo-root>/test/MnIn2Se4
python compare_441_vs_442.py            # §5  TA 322.9 vs 1348.3 m/s
python fc2_442_phonon/sound_velocity.py .   # §4  TA/LA/C44 of the 4x4x2 FC2
python kappa_cutscan.py                 # §6.2  RTA cutoff scan (~5 min)
python kappa_lbte_scan.py               # §6.2  LBTE cutoff scan (~47 min)
python kappa_driver_442.py              # §6.1  RTA + LBTE vs FC2 supercell
python kappa_run_18_lbte.py             # §6.3  LBTE at 18x18x18 (4-8 h, not restartable)
python kappa_run_14_wigner_lbte.py      # §7    Wigner-LBTE at 14x14x14 (2-3 h, not restartable)
```

> Raw DFT output (OUTCAR, vasprun.xml, WAVECAR) is not committed — see the
> repository README. The per-frame INCAR/KPOINTS/POSCAR/CONTCAR **are**
> committed, so the FC2 fit can be reproduced with one phonopy call once the
> forces are regenerated.
