# Mn₂In₂Se₅ — coherence term and iterative solution at 16³ / 6.51 Å

The in-plane κ quoted in the reports is **RTA × 1.06 (iterative correction
borrowed from ShengBTE) + 0.135 coherence (from a Wigner-RTA run whose data was
never committed)**. This directory replaces both borrowed numbers with phono3py
runs on the same force constants as the κ_zz mesh scans.

| item | setting |
|---|---|
| FC2 / FC3 | `Mn2In2Se5/fc/` (4×4×1), FC3 cutoff **6.51 Å** — the ceiling of this cell |
| q-mesh | 16×16×16 |
| temperature | 300 K, isotopes **off** (same as the κ_zz mesh scans) |
| solvers | Wigner-RTA (`transport_type="SMM19"`), then plain LBTE |

## Why three attempts

1. **412628** (LBTE then Wigner-LBTE in one go) — the script crashed after the
   successful LBTE solve: phono3py returns `(1, 6)` or `(1, 1, 6)` kappa arrays
   and `snap()` indexed one axis. Fixed to `reshape(-1, 6)[0]`.
2. **414467** used `conductivity_type="wigner"`, which **does not exist in
   phono3py 4.4.0** — the coherence solver is selected with
   `transport_type="SMM19"` (variants SMM19 / NJC23 / IBDB19).
3. **412628's rerun** was SIGKILLed by the kernel after 2 h 10 min at only
   11 GB RSS: the node was out of memory (251 GB total, 239 GB used by other
   jobs), not our cgroup limit. Re-ran on a different node.

So the current job (**414467**, `jobs_lbte_rta/`) runs in two steps and writes
its JSON after *each* solve, so a kill can no longer lose everything:

1. **Wigner-RTA** → `result_wigner_rta.json`: TOT / population / coherence.
2. **plain LBTE** → `result_lbte.json`: iterative κ plus the RTA value from the
   same run.

## Files

| file | contents |
|---|---|
| `result_wigner_rta.json` | coherence term: TOT − population, at 16³ |
| `result_lbte.json` | LBTE κ and RTA κ from the same run |
| `kappa-m161616.hdf5` | per-mode output |
| `jobs/` | the original two-solution script (indexing fixed) |
| `jobs_wigner_only/` | Wigner-only variant, superseded by `jobs_lbte_rta/` |
| `jobs_lbte_rta/` | **what actually runs**: incremental Wigner-RTA + LBTE |
