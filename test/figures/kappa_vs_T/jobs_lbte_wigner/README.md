# Plan: LBTE and Wigner coherence for the ordered and the VCA crystal, κ(T)

**Status (2026-11-19): ready to run, not run yet.** Driver tested on 6³ / 4³ meshes
in a phono3py 3.24.0 environment (results in `test_m6/`).

## Why

Every κ(T) curve in the report — ordered and VCA — is the **RTA population term**.
The 300 K values to quote (report §1, §2.4) add two corrections, but only for the
**ordered** crystal and only at **300 K**:

- LBTE instead of RTA (in-plane +24 % for MnIn₂Se₄, +5.9 % for Mn₂In₂Se₅);
- the Wigner coherence term (+2–5 % in-plane, +10–12 % κ_zz).

These factors cannot be reused for the VCA crystal. Strong elastic disorder
scattering is resistive, so the LBTE correction should shrink, and the larger
linewidths make the coherence term grow. The 6³ test already shows both
(MnIn₂Se₄, 300 K, coarse mesh, qualitative only):

| 6³, 300 K | LBTE / RTA, in-plane | coherence / RTA, in-plane | coherence / RTA, κ_zz |
|---|---:|---:|---:|
| ordered | 1.31 | 0.04 | 0.05 |
| VCA (i = 0.44) | **1.10** | **0.16** | **0.13** |

The jobs below give LBTE/RTA and coherence/RTA as a function of T for both
cases. They are then applied to the converged RTA curves, so the figures show
the best treatment at every temperature and the stars of report §4 are no longer
needed.

## What the driver does

`run_lbte_wigner.py --material M --cases C` runs, for each case, on one mesh:

| step | phono3py call | output |
|---|---|---|
| `rta` | `run_thermal_conductivity()` | RTA population (reference on the same mesh) |
| `lbte` | `run_thermal_conductivity(is_LBTE=True)` | direct LBTE solution |
| `coherence` | 3.x: `conductivity_type="wigner"`; 4.x: `transport_type="IBDB19"` | Wigner-RTA population + coherence |

- Isotopes: ordered = natural isotopes; VCA = averaged site masses + mass variances.
- VCA site model (same as the existing 20³ runs and the L3 model; x = fraction of all Mn on In sites):
  MnIn₂Se₄ x = 0.44 (Range et al. 1991) → 81.29 / 101.65 amu, g = 0.1337 / 0.0595;
  Mn₂In₂Se₅ x = **0.19** (Williams et al. 2024, ESI Table S7) → 66.32 / 103.44 amu, g = 0.1255 / 0.0516.
- **Version rule for the coherence term** (REPORT §7): 4.4.0 `SMM19` is wrong for
  MnIn₂Se₄, so on 4.x the driver uses `IBDB19`; on 3.x it uses `wigner` (SMM19 of 3.24.0).
  Both are accepted. Do not change this to SMM19 on 4.x.
- Writes `result_<MAT>_<case>_m<N>.json`: κ per T for each step and the ratios
  `lbte_over_rta` and `coherence_over_rta`. κ_zz in this file is the plain mesh
  total (it still contains the Γ–Z artifact); only the ratios and the absolute
  coherence term are used later.

## Jobs (HanHai, `submit.sh`)

LBTE holds the collision matrix of **all requested temperatures at once**:
about 2.6 GB per temperature for MnIn₂Se₄ at 14³ and 9.1 GB for Mn₂In₂Se₅ at 16³,
plus a similar amount for the solver. Hence the temperature split for Mn₂In₂Se₅.

| # | job | command | memory | time (64 threads, estimate) |
|---|---|---|---:|---:|
| 0 | smoke test (any machine, first) | `python3 run_lbte_wigner.py --material MnIn2Se4 --mesh 6 --temps 100 300` | < 1 GB | 2 min (4 cores) |
| 1 | MnIn₂Se₄ ordered, 14³, 9 T | `sbatch --export=ALL,MAT=MnIn2Se4,CASE=ordered submit.sh` | ≈ 50 GB | 1–2 h |
| 2 | MnIn₂Se₄ VCA, 14³, 9 T | `sbatch --export=ALL,MAT=MnIn2Se4,CASE=vca submit.sh` | ≈ 50 GB | 1–2 h |
| 3a–c | Mn₂In₂Se₅ ordered, 16³, 3 × 3 T | `sbatch --export=ALL,MAT=Mn2In2Se5,CASE=ordered,TEMPS="30 60 100" submit.sh`, then `TEMPS="150 200 300"`, `TEMPS="400 500 700"` | ≈ 60 GB each | 3–5 h each |
| 4a–c | Mn₂In₂Se₅ VCA (x = 0.19), 16³, 3 × 3 T | same as 3 with `CASE=vca` | ≈ 60 GB each | 3–5 h each |
| 5 | Mn₂In₂Se₅ VCA (x = 0.19), RTA only, 20³, 12 T, per-mode hdf5 | `sbatch --export=ALL,MAT=Mn2In2Se5,CASE=vca,MESH=20,STEPS=rta,WRITE_KAPPA=1,TEMPS="30 60 90 100 120 150 200 300 400 500 600 700" submit.sh` | ≈ 20 GB | 1 h |

Order: 0 → 1, 2 (MnIn₂Se₄; these alone already settle the question for one
material) → 3, 4 → 5. Jobs 1–4 are independent and can run in parallel.

Job 5 replaces the existing Mn₂In₂Se₅ VCA curve, which used x = 0.13 from the
fit, by the literature x = 0.19. It writes the per-mode hdf5 (`WRITE_KAPPA=1`)
for the Γ–Z decomposition with `test/MnIn2Se4/kappa_zz_decompose.py`.

If a job is killed for memory, split its temperatures further (`TEMPS="30 60"`
etc.); results of different temperature sets are combined afterwards.

## Checks before trusting the results

| check | expected | source |
|---|---|---|
| job 0 RTA, MnIn₂Se₄ 6³ 300 K, ordered | in-plane 2.6628, κ_zz 1.2775 | `test_m6/result_MnIn2Se4_ordered_m6.json` |
| job 0 RTA, VCA | in-plane 0.8895, κ_zz 1.1571 | `test_m6/result_MnIn2Se4_vca_m6.json` |
| job 0 LBTE/RTA in-plane, ordered / VCA | 1.309 / 1.102 (identical population in 3.24.0 and 4.4.0) | same |
| job 0 coherence, 4.4.0 IBDB19 vs 3.24.0 | within ≈ 10 % (11³ reference: 0.114 vs 0.106) | REPORT §7 |
| job 1 RTA 300 K, MnIn₂Se₄ 14³ | 3.6474 / 1.0934 | REPORT §6.3 |
| job 1 LBTE 300 K | in-plane 4.5219 (± 0.5 %); κ_zz = RTA to < 0.1 % | REPORT §9.1, job 412627 |
| job 1 coherence 300 K | ≈ 0.09 / 0.066 | REPORT §7 |
| job 3b RTA 300 K, Mn₂In₂Se₅ 16³ | in-plane 2.9599, κ_zz total 1.5432 | `test/Mn2In2Se5/iso_rerun/` |
| job 3b coherence 300 K | ≈ 0.147 / 0.114 (that was SMM19; IBDB19 may differ by a few %) | job 419208 |

If job 1 does not reproduce 4.5219, stop and report; nothing else is meaningful
before that.

## What to send back

Commit the `run_*/result_*.json` files and `run.log` (and, for job 5, the
`kappa_vca_rta_m20.hdf5`, ≈ 40 MB) to the repository, or copy them to the E
drive, and say which jobs finished. Do not edit the report.

## What happens with the results (repository side)

1. For each T and case: in-plane = converged RTA × (LBTE/RTA) + coherence;
   κ_zz = converged RTA term × (LBTE/RTA)_zz + coherence_zz.
2. Replace the RTA curves in report §4–§5 (ordered and VCA) with these; drop the
   stars; update §1 (VCA rows become LBTE + coherence) and §2.4.
3. Mn₂In₂Se₅ VCA at x = 0.19 from job 5 replaces the x = 0.13 curve.
4. Optional, later: the L3 model is built on RTA bases; rebuilding it on LBTE would
   change the fitted L by roughly the LBTE correction. Not part of this plan.
