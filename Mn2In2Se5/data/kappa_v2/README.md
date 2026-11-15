# Mn₂In₂Se₅ S6 rerun results (after the 0127 fix, ShengBTE)

> This directory holds the final S6 lattice thermal conductivity results after the [0127] fix to the ShengBTE third-order force-constant export,
> superseding the old kappa (the old results are invalid because the hiphive export of FORCE_CONSTANTS_3RD lost 34% of the interlayer third-order terms).
> Corresponding job: hanhai25 359384 (2026-11-08 12:23→16:25, exit code 0).

## Calculation configuration

| item | value |
|---|---|
| solver | ShengBTE 2.1.0 (hanhai25, OpenMPI 5.0.9-ucx, 256 cores) |
| bte_method | conv (iterative solution, main result); RTA is also written |
| mesh strategy | MESH=auto + MESH_CONV_MAX_LENGTH=50 |
| final mesh | 16×16×4 (mesh_converged=False, stop_reason=MESH_CONV_MAX_LENGTH=50 reached) |
| cutoff | 6.51 Å (16 shells, pick=largest; within the range the data can determine κ does not reach a plateau, so the largest cutoff is taken) |
| temperature | 50–800 K, step 50 (16 points) |
| isotope / NAC | ISOTOPE=true / NAC=off |
| force constants | source_fc = cutoff_scan/cut3_6p51 (see ../fc/) |

## Key results (300 K)

| method | κ_xx=κ_yy | κ_zz | unit |
|---|---|---|---|
| CONV (iterative solution) | 3.1412 | 0.84092 | W/m/K |
| RTA | 2.9535 | 0.83455 | W/m/K |

## Mesh convergence (conv, 300 K, selected cutoff 6.51 Å)

| mesh | κ_xx=κ_yy | κ_zz | relative change (xx/zz) |
|---|---|---|---|
| 12 12 3 (L=40 Å) | 2.8984 | 0.81543 | — |
| 16 16 4 (L=50 Å) | 3.1412 | 0.84092 | +8.38% / +3.13% |

> converged=False: xx/yy differs by 8.38%, which does not meet the 5% threshold, and the MESH_CONV_MAX_LENGTH=50 cap has been reached;
> the next step, 20×20×5, runs out of memory (750 GB of RAM, measured rc=137).
> The mesh refinement for κ_zz is not in this directory; it is supplied by phono3py (the intrinsic κ_zz keeps growing with the mesh and is not converged, see analysis/shengbte_crosscheck/).

## Cutoff scan (starting mesh 12×12×3, 300 K)

> All three cutoffs were run on the **starting mesh 12×12×3** (37 q-points), **not** on the final 16×16×4.
> The cutoff was chosen using the **CONV iterative solution**; after 6.51 Å was selected, that case was converged to 16×16×4
> (see the previous section: at 16×16×4 the same case gives CONV = 3.1412 / 0.84092, 8.4% higher than at 12×12×3).
> **Do not treat the values in the table below as 16×16×4 results.**

| cutoff (Å) | shells | CONV in-plane κ | CONV κ_zz | RTA in-plane κ | RTA κ_zz |
|---|---|---|---|---|---|
| 5.68 | 12 | 3.0229 | 0.96837 | 2.9336 | 0.97466 |
| 6.12 | 14 | 2.9422 | 0.88014 | 2.7954 | 0.88671 |
| **6.51 (selected)** | 16 | **2.8984** | **0.81543** | 2.7257 | 0.82116 |

## 0127 re-export verification

- cut3_6p51: the old export was missing 3126/9164 triplets → after re-export there are 9164 blocks / 9164 triplets, and the fold-back self-check passes.
- FORCE_CONSTANTS_3RD md5 = af86b14cb8013e6ac49426b5a36aefbc (7887357 bytes).
- The old incomplete version has been renamed stale_0127 on the cluster (not uploaded).

## File list

| file | description |
|---|---|
| kappa_summary.json | final κ summary (16×16×4, full 50–800 K table, conv as the main result) |
| mesh_convergence.json | mesh convergence record (12 12 3 → 16 16 4, not converged) |
| mesh_convergence_reference.json | mesh convergence reference |
| cutoff_selection.json | cutoff selection (chosen_cut=6.51 Å) |
| BTE.KappaTensorVsT_CONV | iterative-solution κ(T), 16×16×4 |
| BTE.KappaTensorVsT_RTA | RTA κ(T), 16×16×4 |
| FORCE_CONSTANTS_3RD | 6.51 Å third-order force constants (after the 0127 re-export, 9164 blocks) |
| FORCE_CONSTANTS_2ND | 6.51 Å second-order force constants |
| shengbte.log | ShengBTE run log (including the 0127 re-export record) |

## Notes

- RTA vs CONV: `source=BTE.KappaTensorVsT_CONV` in kappa_summary.json; the main result is conv, and BTE.KappaTensorVsT_RTA is the relaxation-time approximation on the same mesh.
- Cutoff not converged: κ decreases slowly and monotonically with the cutoff; the 6.12–6.51 plateau span is below the 0.5 Å threshold, so pick=largest gives 6.51 Å.
- Uploaded and organized: 2026-11-08.
