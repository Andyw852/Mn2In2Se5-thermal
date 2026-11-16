# MnIn₂Se₄ — Wigner (coherence) at the largest FC3 cutoff (14³, 6.45 Å)

Companion to the LBTE run in `../` (job 412627: in-plane **4.5219** / κ_zz
**1.0934**, identical to the committed 14³ values), so the "LBTE + coherence"
composition in report §9.1 has both halves on the same mesh and cutoff.

## Wigner-RTA — done (job 414857, 11.6 s)

| component | in-plane | κ_zz |
|---|---:|---:|
| population | 3.6474 | 1.0934 |
| **coherence** | **0.6539** | 0.5982 |
| TOT | 4.3013 | 1.6916 |

The population part reproduces the standard RTA value (3.6474 / 1.0934), which is
the cross-check that the setup is right.

> ⚠️ This coherence term is **7× larger** than the 0.0908 in
> `../kappa_wigner_mn441.json` (phono3py **3.24.0**,
> `conductivity_type="wigner"`). Same physical setting, population agreeing to
> 0.3 %, so the difference sits in the effective heat-capacity matrix — 4.4.0
> uses `mode_cv_matrix_smm19`, labelled experimental in the source. Local
> compute is forbidden here, so 3.24.0 cannot be re-run for a controlled
> comparison. Report §7 states both numbers and tells the reader to quote
> in-plane κ including coherence as **4.6–5.2** until this is settled.

## phono3py 4.4.0 API map (each wrong guess cost a run)

| what | 3.24.0 | 4.4.0 |
|---|---|---|
| select the coherence solver | `conductivity_type="wigner"` | `transport_type="SMM19"` |
| RTA: total / population / coherence | `kappa_TOT_exact` / `kappa_P_exact` / — | `kappa` / `kappa_intra` / `kappa_inter` |
| LBTE: total / population-exact / population-RTA / coherence | `kappa_TOT_exact` / `kappa_P_exact` / — / — | `kappa` / `kappa_intra_exact` / `kappa_intra_RTA` / `kappa_inter` |

Every kappa array can be `(1, 6)` or `(1, 1, 6)`; hence `reshape(-1, 6)[0]`.

## Job history

| job | what | outcome |
|---|---|---|
| 412627 | LBTE then Wigner-LBTE | LBTE **done** (`kappa-m141414.hdf5`); then crashed writing JSON (array shape) |
| 414024 | Wigner-LBTE, `conductivity_type`→`transport_type` but 3.24.0 attribute names | ran 45 min, failed on `kappa_TOT_exact` |
| 414719 | same, 4.4.0 RTA attribute names | died after 5 min 53 s with exit 1 and an **empty stderr** — node problem, not code (`anode116`) |
| 414857 | Wigner-RTA (saved) + Wigner-LBTE | RTA **done**; LBTE **solved in 1 h 06 min** but then failed on `kappa_intra` (LBTE uses `kappa_intra_exact`) — traceback in `lbte_error_414857.txt` |
| **415323** | Wigner-LBTE with the correct names, TOT written to JSON *before* the part attributes | the run of record → `result_wigner_lbte.json` |

## Files

| file | contents |
|---|---|
| `result_wigner_rta.json` | Wigner-RTA TOT / population / coherence at 14³ |
| `result_wigner_lbte.json` | Wigner-LBTE TOT / population / coherence (job 415323) |
| `kappa-m141414.hdf5` | LBTE per-mode output (job 412627) |
| `lbte_error_414857.txt` | why the first successful Wigner-LBTE solve still produced nothing |
| `jobs/`, `jobs_wigner_only/` | earlier scripts |
| `jobs_wigner_lbte/` | **what actually runs**: Wigner-LBTE with guaranteed TOT capture |
