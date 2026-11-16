# MnIn₂Se₄ — Wigner-LBTE at the largest FC3 cutoff (14³, 6.45 Å)

Companion to the LBTE run in `../` (job 412627, which gave in-plane **4.5219** /
κ_zz **1.0934**, identical to the committed 14³ values). This job adds the
coherence term from the **same** solution, so the "in-plane κ = LBTE + coherence"
number in report §9.1 is no longer stitched from two different calculations.

## phono3py 4.4.0 API (two wasted runs to find)

| what | 3.24.0 | 4.4.0 |
|---|---|---|
| select the coherence solver | `conductivity_type="wigner"` | `transport_type="SMM19"` |
| total κ | `kappa_TOT_exact` | `kappa` |
| population part | `kappa_P_exact` | `kappa_intra` |
| coherence part | derived as TOT − P | `kappa_inter` |

Both wrappers accept `(1, 6)` or `(1, 1, 6)` kappa arrays, hence the
`reshape(-1, 6)[0]` in `snap()`.

## Jobs

| job | what | outcome |
|---|---|---|
| 412627 | LBTE then Wigner-LBTE | LBTE **done** (data in `../kappa-m141414.hdf5`); then crashed writing the JSON because of the array shape |
| 414024 | Wigner-LBTE with `transport_type="SMM19"` but the 3.24.0 attribute names | ran 45 min, then would have failed on `kappa_TOT_exact` — superseded |
| **414719** | same, with `kappa` / `kappa_intra` / `kappa_inter` | the run of record → `result_wigner_lbte.json` |

| file | contents |
|---|---|
| `result_wigner_lbte.json` | TOT / population / coherence at 14³, 6.45 Å, isotopes on |
| `kappa-m141414.hdf5` | LBTE per-mode output (job 412627) |
| `slurm.out`, `slurm.err` | logs; `slurm.err` holds the JSON-crash traceback |
| `jobs/` | the original two-solution script (indexing fixed) |
| `jobs_wigner_only/` | **what actually runs**: Wigner-LBTE only, 4.4.0 API |

Purpose: replace the mixed-provenance numbers in report §9.1 — 4.5219 LBTE at 14³
plus 0.091 coherence from a *separate* Wigner-**RTA** run — with one LBTE run that
yields both components at the same cutoff and mesh.
