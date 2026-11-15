# MnIn2Se4 κ (phono3py) — recomputed after switching fc2 to 4×4×2

**2026-11-14** | computed directly with phono3py, bypassing AutoZT (AutoZT's kl-dft-cpu / fit-fc-thermal
do not support using different supercells for fc2 and fc3).

## Why recompute

The old 4×4×1 fc2 gives a transverse sound velocity along the layer normal of only **322.9 m/s**, whereas the DFT elastic constants give **1338.3 m/s**
(see `../fc2_442_phonon/README.md` and `../compare_441_vs_442_results.json`).
The phonon dispersion is severely softened → **κzz is underestimated by about 3.5×**.

The same branch from the 4×4×2 fc2 gives **1348.3 m/s** (+0.75%), and the dispersion is correct.

## Calculation method

| item | value |
|---|---|
| solver | phono3py 3.24.0 |
| fc3 | `fit_cut3_6p45/fc3.hdf5`, **4×4×1 (112 atoms), cutoff 6.45 Å** |
| fc2 (new) | `fc2_442_phonon/fc2.hdf5`, **4×4×2 (224 atoms)** |
| fc2 (old, control) | `fit_cut3_6p45/fc2.hdf5`, 4×4×1 |
| mixed supercells | use phono3py's `phonon_supercell_matrix=4×4×2` (fc3 goes through `supercell_matrix=4×4×1`) |
| primitive-cell matrix | `[[1,1,2],[-1,-2,-1],[1,1,1]]` (consistent with both fits) |
| mesh | **11×11×11** (isotropic; a 4×4×2 mesh cannot be used with a rhombohedral primitive cell) |
| temperature | 300 K |
| isotope | on |
| NAC | off |

## ✅ Results (mesh 11×11×11, 300 K)

| fc2 | method | in-plane (W/m·K) | κzz (W/m·K) | κzz/in-plane | wall time |
|---|---|---:|---:|---:|---:|
| 4×4×1 (old) | RTA | 2.2498 | 0.3508 | 0.156 | 111 s |
| 4×4×1 (old) | LBTE | 3.9054\* | 0.3629\* | 0.093 | — |
| **4×4×2 (new)** | **RTA** | **3.5375** | **1.2495** | **0.353** | 107 s |
| **4×4×2 (new)** | **LBTE** | **4.4075** | **1.2519** | **0.284** | 939 s |

\* Historical values (same mesh, fc2=4×4×1) from `../kappa_phono3py_lbte_11x11x11/kappa_summary.json`.

### Effect of switching fc2 (LBTE, same mesh, same fc3)

| quantity | 4×4×1 | 4×4×2 | change |
|---|---:|---:|---:|
| in-plane | 3.9054 | **4.4075** | **+12.9%** |
| κzz | 0.3629 | **1.2519** | **+245% (3.45×)** |
| κzz/in-plane | 0.093 | **0.284** | 3.05× |

**κzz increases by a factor of 3.45**, exactly the direct consequence of the interlayer transverse sound velocity going from 323 → 1348 m/s.
The in-plane value rises by only 12.9%, because the in-plane sound velocity was already essentially correct.

### Self-check

Running RTA with the old fc2 in my driver gives **2.2498**, identical to the historical **2.24981** — the driver itself is trustworthy.

## ⚠️ These numbers cannot yet be treated as final conclusions

1. **fc3 is unchanged**: still 4×4×1 / cutoff 6.45 Å, previously marked `not_converged`.
   The absolute κ is governed by the fc3 cutoff; this work only **isolates the effect of fc2**.
2. **The mesh is not converged**: 11×11×11. Historical experience: Mn₂In₂Se₅ rose by +7.65% from 12×12×3 → 16×16×4.
3. Producing a final κ requires: (1) a convergence check on the fc3 cutoff (larger cutoff / thicker supercell); (2) a mesh convergence check.

## Reproduction

```bash
cd test/MnIn2Se4
python kappa_driver_442.py            # run both RTA and LBTE (about 20 minutes)
python kappa_driver_442.py --rta      # RTA only (about 4 minutes)
```

Results are written to `kappa_results_442.json`.
