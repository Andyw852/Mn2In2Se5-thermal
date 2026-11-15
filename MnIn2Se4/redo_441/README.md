# MnIn₂Se₄ recalculated with the same protocol as Mn₂In₂Se₅ (AutoZT kl-dft-cpu, 4×4×1 rebased supercell)

## Why recalculate

The existing force constants in `MnIn2Se4/` come from the official pheasy test dataset: the DFT settings are unknown (whether +U, spin polarization or D3 was used is not known),
and the safe cutoff of the 3×3×3 rhombohedral supercell is only about 5.2 Å (the interlayer third-order force constants are cut off, and κ_zz varies erratically between 0.5 and 4.0 with the cutoff).
It is not the same set of calculations as Mn₂In₂Se₅ (PBE-D3 + U(Mn)=4 + FM, 144 atoms, cutoff 6.51 Å), so the two materials cannot be compared directly,
and the difference cannot be attributed to the materials themselves.

## Instructions for the executing AI

**Copy [`../../protocol/kl-dft-cpu_441_Mn2In2Se5_instructions.md`](../../protocol/kl-dft-cpu_441_Mn2In2Se5_instructions.md) verbatim**,
changing only the items in the table below. Use the current AutoZT master (already containing [0105]–[0107]: spin/U inheritance, symmetry-preserving primitive cell, force-parsing fix, automatic memory increase on OOM).
The iron rules still apply: ask separately about each `rerun`, `clean`, `-f`, `-y`; do not write your own monitoring loop — use `autozt ... monitor -d`.

| item | Mn₂In₂Se₅ (original instructions) | **MnIn₂Se₄ (this run)** |
|---|---|---|
| material name | `Mn2In2Se5_441` | **`MnIn2Se4_441`** |
| material-root POSCAR | 9 atoms | this directory's `POSCAR`: 7 atoms, order Mn In In Se×4, already rebased |
| self-check output | `R-3m (#166), 9 atoms, 12 point-group operations` | `R-3m (#166), 7 atoms, 12 point-group operations`, and `4 4 1 preserves all point-group operations, no rebasing needed` |
| lattice after S1 | |a|=|b|≈4.03, third basis vector ≈18.7 Å | |a|=|b|≈4.07, third basis vector ≈16.2 Å (it will change slightly after the PBE-D3+U relaxation) |
| magnetic order / U | FM, Mn U=4 | **same**: `MAGMOM_OVERRIDE="Mn:5 In:0 Se:0"`, `U_OVERRIDE="Mn:4 In:0 Se:0"` (the primitive cell contains only 1 Mn, so FM is the only collinear state) |
| expected S2 total moment | ≈ 9.98 μB/cell | **≈ 5.0 μB/cell** (denoted M₀; if it deviates appreciably from 5, stop and report) |
| KPOINTS | primitive cell 8×8×2, supercell 2×2×2 | **same** (the two KPOINTS files in this directory; with KSPACING=0.038, |b|/2π = 0.3012/0.0758 → 8 8 2) |
| S4 supercell | 4 4 1 = 144 atoms | 4 4 1 = **112 atoms**, compare the lattice against `SPOSCAR_441_reference` in this directory |
| S4 moment per frame | 16 × M₀ ≈ 159.75 | **16 × M₀ ≈ 80** |
| safe cutoff | 6.61 Å | **6.50 Å** (`CUT3_MAX` stays at 7.0; cutoffs exceeding it are dropped automatically by S4) |
| NAC | S3 failed → S6 `NAC=off` | S3 still runs; **S6 always uses `NAC=off`**, consistent with Mn₂In₂Se₅ so that the two materials can be compared directly |
| S1 force gate | `CELL_FORCE_TOL_EV_A=0.002`, `CELL_STAGE_EDIFFG=1E-6` | **same** |
| S5 | `EQ_FORCE_MAX=0.01`, `FC_IMPORT_DIR=off` | **same** |

**S6 mesh (decided 2026-11-08, identical for both materials)**: `MESH = auto` + `MESH_CONV_MAX_LENGTH = 50`
(the mesh is refined automatically, capped at 16×16×4). **Do not set an explicit `MESH`**: the ShengBTE iterative solve already used 750 GB of memory at 20×20×5,
and 16×16×16 would only be more. κ_xx is already converged at 16×16×4 (20×20×5 differs by only 0.9%); the interlayer mesh refinement for κ_zz is done on the repository side with phono3py,
with the results in `analysis/shengbte_crosscheck/README.md`.

## Reporting and upload

When done, upload to `MnIn2Se4/redo_441/` in this repository following the structure of `Mn2In2Se5/`: `data/kappa`, `data/phonon`, `figures`,
`fc/` (fc2.hdf5, fc3.hdf5, phono3py_disp.yaml for the selected cutoff).
Also write a `dft_summary.txt`: the final maximum force of S1, the S2 band gap and total magnetic moment, the S4 equilibrium-frame residual force, the number of frames, and the per-frame magnetic-moment distribution.

---

## DFT-stage archive (2026-11-09, AutoZT kl-dft-cpu)

S1–S3 are complete and the S4 inputs have been generated (**the calculation was not submitted**); archived in [`dft/`](dft/):

- `dft_summary.txt` — S1 final-state maximum force / S2 magnetic moment and the **semimetal determination** / S3 invalidation note / S4 input checklist and memory configuration
- `s1_std_opt/`, `s2_static/`, `s3_nac/` — inputs (INCAR/POSCAR/KPOINTS/submit.sh), structures (CONTCAR), energies (OSZICAR)
- `s4_disp_inputs/` — SPOSCAR (112 atoms, 4×4×1), sample INCAR/KPOINTS/submit.sh for the displacement frames

Following the repository `.gitignore` conventions, the raw `OUTCAR` files are not committed (they can be re-fetched from the cluster).
**`s2_static/EIGENVAL` is kept as an exception** (via `git add -f`): it is the direct evidence for the "semimetal vs insulator" determination.

**Key results**:
- **S1 ✅** max|F| = 0.00064 eV/Å, R-3m, 7 atoms, a=b=4.0738 Å, c=16.2213 Å.
- **S2 ⚠️ compensated semimetal**: magnetic moment 4.9844 = 5 − 1/64; weighted electron counts up/down = 33.9922 / 29.0078.
  **both spin channels have pockets**: the spin-up valence-band top at 3.9945 is a hole pocket (8 k-points with occupations 0.91–0.99, together short by 1/128),
  compensated by two electron pockets of 1/128 each at the Γ point: spin-up (bands 34=35=3.9371, a degenerate pair) and spin-down (band 30=3.8039),
  giving a spin-crossing overlap of ≈ 0.19 eV. Mn₂In₂Se₅ is completely analogous (up/down 45.9922/36.0078, spin-crossing overlap ≈ 0.13 eV).
  **No band gap is reported**; the uniform statement is "semimetallic under PBE-D3 + U(4 eV) + FM". The real sample is very likely a semiconductor, which is caused by DFT+U underestimating the gap
  combined with the FM assumption; an HSE/AFM calculation would be needed before a conclusion can be drawn (κ_L does not need this).
- **S3 ❌ entirely invalid**: ε∞ is identically 1, the piezoelectric tensor is zero and the Born charges are **negative** (Mn −0.94, In −0.62, whereas cations should be positive) —
  the LEPSILON linear response was never actually computed (every step `E=0, ncg=0`), and the root cause is precisely that the Γ-point occupation rearrangement in S2 violates the DFPT premise.
  **This does not affect the main line** (S6 uses `NAC=off` and S5 does not use the dielectric tensor); if NAC is ever turned on, this step must be fixed first.
- **S4 ✅ complete (2026-11-10, 46/46 frames)**: SPOSCAR 112 atoms, KPOINTS Γ 2×2×2, MAGMOM 16*5 96*0;
  the core count finally used was 48 (`--mem` has no effect on this cluster: `RealMemory=1`, DefMemPerNode=UNLIMITED).
  all 46 frames converged; the equilibrium-frame residual force is **0.000508 eV/Å**; displacement frames have max|F| 0.62–0.83 eV/Å;
  the total magnetic moment of all 46 frames lies within **79.74994–79.74999 μB/cell** (target 16×4.9844 = 79.750, deviation <1e-4),
  i.e. the displacements triggered no occupation/magnetic-moment jump. The raw per-frame data are archived in [`dft/s4_disp/`](dft/s4_disp/):
  46 `disp-*/vasprun.xml` (`git add -f`, including forces) + per-frame energy/force/magnetic moment in `frame_summary.csv`.

### Imaginary-frequency findings and handling for the S4 (4×4×1) fit (2026-11-10)

Two pheasy-gpu methods were fitted to the 46 frames above (`FIT_METHODS = pheasy:RFE-OLS,ALASSO`,
PHEASY_C3_CUTOFF=6, stepwise scan 4.30–6.45 Å, CUT3_BOOTSTRAP=0):

| method | force RMSE | relative error | imaginary-frequency gate |
|---|---|---|---|
| RFE-OLS (primary) | 0.88 meV/Å | 0.58 % | **FAIL** min = **−0.1388 THz (−4.63 cm⁻¹)**, no-NAC −0.151 THz |
| ALASSO | 0.84 meV/Å | 0.55 % | **FAIL** min = −0.1388 THz |

- Minimum: the lowest TA branch (branch 0), q = (−0.0357, 0, 0.0357), |q|≈0.0505; 145 q-points with negative frequencies;
  classified as `near_gamma_large`; the fc2 symmetry spread is 5.3e-9 (the symmetry itself is fine).
- The **harmonic parts of the two methods are identical** (HARM is not penalized) and the results are nearly the same → the imaginary frequency is not a fitting-method problem.
- Dispersion probe (reading fc2.hdf5 directly): the three acoustic branches at Γ are ≈ 0 (ASR is fine), but the lowest TA branch is negative over the whole range |q|≈0.03–0.12
  and only turns positive at 0.15 → a **broad soft TA branch** with a slightly negative effective elastic constant, not single-point noise.
- Dispersion probe direction: the three acoustic branches at Γ are ≈ 0 (which only shows that the **translational** sum rule is satisfied), but the lowest TA branch is negative over the whole range |q|≈0.03–0.12
  and only turns positive at 0.15 → the long-wavelength TA branch is soft, not single-point noise.
- **Root causes (corrected 2026-11-10; it is not the supercell size)**:
  1. **Rotational invariance / the Huang equilibrium condition was not imposed** (`PHEASY_RASR=auto` resolves to none for bulk) —
     a normal Γ point together with negative frequencies at |q|≈0.03–0.12 is exactly the typical signature of a missing constraint of this kind;
  2. **the equilibrium frame carries about −1.23 kBar of isotropic residual stress** (diagonal −1.246 / −1.246 / −1.207 kBar, deviatoric part <0.03 kBar),
     coming from the Pulay stress difference between the S1 relaxation with ENCUT=540 and the S4 force evaluation with 410; this governs the same physics as the Huang condition
     and softens the long-wavelength acoustic branches (it cannot explain −0.14 THz on its own, but it contributes).
  Geometrically, 4×4×1 is not at a disadvantage: the perpendicular widths in all three directions are ≥13.2 Å (safe cutoff 6.50 Å), wider than the old 3×3×3 rhombohedral supercell at
  about 10.5 Å (5.17 Å); Mn₂In₂Se₅ with the same 4×4×1 setup is stable.
- **Handling**: redo the fit on the same 4×4×1 data, changing only `PHEASY_RASR` from none to **H** (Huang constraint),
  with the method/cutoff/all other parameters unchanged. If near-Γ negative frequencies below −0.05 THz remain after RASR=H, recompute the equilibrium frame with ENCUT=540
  (a single frame) to judge whether the residual stress is the dominant cause, and only then consider changing the S4 settings and recalculating.
  (The 3×3×3 / 189-atom dataset `datasets/MnIn2Se4_raw` is a **test dataset shipped with pheasy whose DFT settings are unknown**;
  it can only serve as a reference and cannot be used for results on this material.)

- **Note on NAC (corrected 2026-11-10)**: one cannot claim that "R-3m has an inversion center → the Born effective charges are identically 0" — NaCl and Bi₂Te₃
  both have an inversion center, yet their Born charges are non-zero and they show LO-TO splitting (inversion symmetry only makes the **piezoelectric tensor** zero). The legitimate reasons for switching NAC off in this system are:
  (1) under this DFT setup the system is a semimetal and free carriers screen the LO-TO splitting; (2) DFPT cannot be computed on this ground state
  (S3 is invalid); (3) consistency with Mn₂In₂Se₅, which makes the κ of the two materials directly comparable.
