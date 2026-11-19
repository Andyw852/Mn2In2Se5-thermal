# For the local AI: tasks for 2026-11-19

One task: run LBTE and the Wigner coherence term for the ordered and the VCA
crystal as a function of temperature. Report after each step before moving on;
the AGENTS.md iron rules still apply.

## 0. Background (in one sentence)

All κ(T) curves so far are RTA population only; the LBTE correction and the
coherence term exist only for the ordered crystal at 300 K, and the 6³ test shows
that both change strongly with disorder, so the VCA curves need their own runs.

## 1. The plan

Everything is in
[`test/figures/kappa_vs_T/jobs_lbte_wigner/README.md`](../test/figures/kappa_vs_T/jobs_lbte_wigner/README.md):
driver `run_lbte_wigner.py`, HanHai script `submit.sh`, the job list with memory
and time estimates, and the checks.

Steps, in order:

1. Pull the branch with this file onto HanHai (`git pull` in `$HOME/work/Mn2In2Se5-thermal`).
2. Job 0 (smoke test, 2 min): compare with `test_m6/`. Report the numbers.
3. Jobs 1 and 2 (MnIn₂Se₄ ordered and VCA, 14³). **Check job 1 at 300 K: LBTE
   in-plane 4.5219 ± 0.5 %.** If not, stop and report.
4. Jobs 3a–c and 4a–c (Mn₂In₂Se₅, 16³, temperatures split for memory).
5. Job 5 (Mn₂In₂Se₅ VCA at x = 0.19, RTA 20³, with per-mode hdf5).

## 2. What to send back

The `run_*/result_*.json` and `run.log` files (job 5 also `kappa_vca_rta_m20.hdf5`),
committed to the repository or copied to the E drive, and which jobs finished.
Do not change the report or the figures; that is done on the repository side.
