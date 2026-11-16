# MnIn₂Se₄ — LBTE and Wigner-LBTE at the largest FC3 cutoff (6.45 Å)

Closes the last gap in the quotable in-plane κ: the report quotes
**LBTE + a coherence term**, and until now the coherence term came from a
Wigner-**RTA** run and the LBTE value from a different (older) run.

| item | setting |
|---|---|
| FC2 | 4×4×2 (`fc2_442_phonon/fc2.hdf5`) |
| FC3 | 4×4×1, **6.45 Å** (`fit_cut3_6p45/fc3.hdf5`) — the ceiling of this cell |
| q-mesh | 14×14×14 |
| temperature | 300 K, isotope scattering **on** |
| solutions | (1) `is_LBTE=True`; (2) `is_LBTE=True, conductivity_type="wigner"` |

## LBTE — done, and it reproduces the committed value

HanHai job **412627**, 1 h 19 min, 3.4 GB peak RSS. The LBTE solve completed and
wrote `kappa-m141414.hdf5`; the numbers read from it are

| | in-plane | κ_zz |
|---|---:|---:|
| **LBTE** (this run) | **4.5219** | **1.0934** |
| RTA, same run | 3.6474 | 1.0934 |

These are **identical to the committed 14³ values** (`kappa_mesh_scan_results.json`),
so the HanHai LBTE path is verified against the original workstation run.

## Wigner-LBTE — had to be re-run (job 414006)

The first script crashed *after* the successful LBTE solve, while writing the
JSON: phono3py returns the kappa arrays as `(1, 6)` or `(1, 1, 6)` depending on
the solver, and `snap()` indexed only one axis. The fix is
`np.asarray(k).reshape(-1, 6)[0]`; the crashing version is kept in `jobs/` with
the fix applied, and `jobs_wigner_only/` holds the trimmed script that runs the
Wigner-LBTE solve alone (HanHai job **414006**), since the LBTE part is already
done and costs 1.3 h.

| file | contents |
|---|---|
| `kappa-m141414.hdf5` | per-mode output of the LBTE run (job 412627) |
| `slurm.out`, `slurm.err` | its log and the traceback of the JSON crash |
| `jobs/` | the original two-solution script (with the indexing fixed) and its submit script |
| `jobs_wigner_only/` | the Wigner-LBTE-only script and submit script actually used to finish |

Purpose: replace the mixed-provenance numbers in report §9.1
(4.5219 LBTE at 14³ + 0.091 coherence from Wigner-RTA) with one LBTE run that
yields both components at the same cutoff and mesh.
