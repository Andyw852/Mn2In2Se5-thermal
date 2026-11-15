# Mn₂In₂Se₅ lattice thermal conductivity κ_L(T) (first principles, 4×4×1 rebased supercell)

Intrinsic lattice thermal conductivity of the R-3m layered Mn₂In₂Se₅ (9-atom primitive cell): DFT finite displacements → pheasy fits the second- and third-order force constants → ShengBTE iteratively solves the phonon Boltzmann transport equation.
Produced by the AutoZT `kl-dft-cpu` skill chain (S1 optimization → S2 static → S3 NAC → S4 displacements → S5 force constants → S6 κ); the execution protocol is given in
[`../protocol/kl-dft-cpu_441_Mn2In2Se5_instructions.md`](../protocol/kl-dft-cpu_441_Mn2In2Se5_instructions.md).

> This is the intrinsic value for a **single crystal, free of defects, with only three-phonon + natural-isotope scattering**. Comparing with hot-pressed polycrystalline experiments additionally requires texture, antisite disorder and grain boundaries;
> see [`../analysis/extrinsic/`](../analysis/extrinsic/).

> ## ⚠️ The table below is superseded (2026-11-08)
>
> The `FORCE_CONSTANTS_3RD` fed to ShengBTE by the cluster S6 was exported by hiphive and lost 34% of the third-order terms (interlayer triplets crossing the a₃′ boundary),
> which makes κ_xx about 1.7× too low and κ_zz about 4–5× too low. With the same fc2/fc3 run through phono3py, or with the completed 3RD run through ShengBTE, the results agree:
> **300 K RTA κ_xx ≈ 2.96, κ_zz ≈ 0.84 (16×16×4) / 1.11 (16×16×8, not yet converged)**.
> Details in [`../analysis/shengbte_crosscheck/`](../analysis/shengbte_crosscheck/); AutoZT [0127] has been fixed, so S6 must be rerun.
> The table below is kept only as a record.

## Original production results (superseded: ShengBTE iterative solution, q-mesh 16×16×4, third-order cutoff 6.51 Å)

| T (K) | κ_xx = κ_yy | κ_zz |
|---:|---:|---:|
| 30 | 26.862 | 2.1153 |
| 60 | 11.610 | 1.0243 |
| 90 | 7.2043 | 0.6770 |
| 100 | 6.4058 | 0.6099 |
| 120 | 5.2563 | 0.5101 |
| 150 | 4.1569 | 0.4107 |
| 200 | 3.0958 | 0.3106 |
| 300 | **2.0590** | **0.2090** |
| 400 | 1.5452 | 0.1575 |
| 500 | 1.2373 | 0.1264 |
| 600 | 1.0319 | 0.1055 |
| 700 | 0.8852 | 0.0906 |
| 800 | 0.7750 | 0.0793 |

Units W m⁻¹ K⁻¹. At 300 K the in-plane/out-of-plane ratio is about 9.9. The 30–150 K and 100–800 K segments use the same force constants and the same q-mesh.

**Cross-check results (2026-11-08)**:

1. **Solver cross-check** → revealed the export problem above, so the production κ is superseded (see the top of this file).
2. **q-mesh for κ_zz**: automatic convergence only raised the third number from 3 to 4. In phono3py RTA, going from 16×16×4 to 16×16×8 raises κ_zz by 33%
   (κ_xx changes by 0.6%); 12 and 16 continue to rise to 1.304 and 1.541, so the intrinsic κ_zz is not converged (long-mean-free-path interlayer phonons),
   and it converges to ≈ 0.97 for a 1 μm grain. S6 was rerun with `MESH=auto` + `MESH_CONV_MAX_LENGTH=50` (capped at 16×16×4);
   the recommended wording for the report is at the end of the cross-check page. The extrinsic-scattering model used for the comparison with experiment is insensitive to the mesh.

## Methods

**DFT (VASP 6.4.3)**

| item | setting |
|---|---|
| functional | PBE + D3 dispersion (AutoZT automatically selects pbe-d3 for layered multi-component systems) |
| magnetism | **ferromagnetic (FM) assumed**, Mn²⁺ high-spin starting point MAGMOM = 5 μB/Mn; AFM was not compared |
| DFT+U | Mn 3d U_eff = 4 eV (Dudarev); none applied to In or Se |
| electronic state | S2: band gap 0.254 eV, total magnetic moment 9.9844 μB/cell (changing SIGMA from 0.05 to 0.01 leaves the moment unchanged → there is a small genuine partial occupation) |
| k-points | KSPACING = 0.038 (2π/Å): primitive cell Γ 8×8×2, supercell Γ 2×2×2; the two are the same physical k-point set |
| structure optimization | three-stage (ISIF=2 → ISIF=3 ×2); final-state maximum-force gate 0.002 eV/Å, EDIFFG = 1E-6 eV in the variable-cell stages; **actual final-state maximum force 0.000484 eV/Å**; ENCUT 540 (1.5×ENMAX, to suppress Pulay stress) |
| force evaluation | ENCUT 410, PREC=Accurate, EDIFF=1E-7, LREAL=.FALSE., ISYM=0; **equilibrium-frame residual force 0.000933 eV/Å** (the old 3×3×3 gave 0.103) |
| displacement supercell | 4×4×1 = **144 atoms** under the rebased vectors (a, b, a₃′), preserving all 12 point-group operations of R-3m; 1 equilibrium frame + 54 displacement frames |
| NAC | **off**: the S3 DFPT failed (ε∞ came out as the identity matrix, Mn Born charge −0.90) and the result is invalid; `data/phonon/BORN_invalid_not_used` is kept only as a record |

**Force constants (pheasy, GPU 3090)**: ALASSO, fully dense second order without a penalty (HARM_DENSE); after the third-order cutoff scan 6.51 Å was chosen
(supercell safe cutoff 6.61 Å); relative fit error 0.415%, 2440 independent parameters; no imaginary frequencies in the phonon spectrum, and the fc2 symmetry breaking over the 12 operations is 1.5×10⁻¹⁴ THz.

**κ (ShengBTE)**: iterative solution (RTA is also written), isotope scattering on (natural abundance), NAC off;
the cutoff was chosen as 6.51 Å by three criteria — 1-SE, plateau and data volume (`data/kappa/cutoff_scan_S6.json`, `figures/kappa_vs_cutoff.png`).

## Files

```
POSCAR                          input structure for S1 (the 9-atom primitive cell after rebasing, before relaxation)
POSCAR_relaxed                  primitive cell after the S1 relaxation (from fc/phono3py_disp.yaml; the 4×4×1 supercell is generated from it)
data/dft_inputs/                INCAR, KPOINTS and dft_summary.txt for S1/S2/S4
config/step6_kappa.step.conf    S6 parameters (solver, mesh, cutoff criteria)
data/kappa/                     κ(T): CSV, summary JSON, raw ShengBTE BTE.KappaTensorVsT_CONV, cutoff scan
data/phonon/                    phonon_summary.json (stability, fit error), fit_manifest.json, BORN_invalid_not_used
data/provenance/                fingerprints of the S1–S6 input files (containing internal cluster paths; the repository stays private)
figures/                        κ(T), anisotropy, cutoff convergence, phonon spectrum (band_nac uses the invalid BORN; ignore it)
fc/                             force constants (fc2/fc3.hdf5, ShengBTE text; the 3RD exists in incomplete and completed versions, see fc/README.md)
```

## Older versions

The old 3×3×1 (81 atoms, symmetry lost, κ jumping erratically between 0.74 and 0.96 with the mesh) and 3×3×3 (equilibrium-frame residual force 0.103 eV/Å, imaginary frequency −0.44 THz) versions
have been removed from the working tree and remain in the git history (commit `5b5c05b`).
