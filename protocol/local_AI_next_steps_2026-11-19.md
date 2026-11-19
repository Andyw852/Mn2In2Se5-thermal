# For the local AI: tasks for 2026-11-19

Two tasks, independent of each other (task 1 on HanHai, task 2 on the 3090):
1. LBTE and the Wigner coherence term for the ordered and the VCA crystal, κ(T).
2. FC3 cutoff beyond 6.45 / 6.51 Å by refitting the existing frames (no new DFT).

Report after each step before moving on; the AGENTS.md iron rules still apply.

## 0. Background (in one sentence)

All κ(T) curves so far are RTA population only; the LBTE correction and the
coherence term exist only for the ordered crystal at 300 K, and the 6³ test shows
that both change strongly with disorder, so the VCA curves need their own runs.

## Task 1 — LBTE and Wigner, κ(T)

### 1.1 The plan

Everything is in
[`test/figures/kappa_vs_T/jobs_lbte_wigner/README.md`](../test/figures/kappa_vs_T/jobs_lbte_wigner/README.md):
driver `run_lbte_wigner.py`, HanHai script `submit.sh`, the job list with memory
and time estimates, and the checks.

Steps, in order:

1. Pull the branch with this file onto HanHai (`git pull` in `$HOME/work/Mn2In2Se5-thermal`).
2. Job 0 (smoke test, 2 min): compare with `test_m6/`. Report the numbers.
3. Jobs 1 and 2 (MnIn₂Se₄ ordered and VCA, 14³). **Check job 1 at 300 K: LBTE
   in-plane 4.5219 ± 0.5 %.** If not, stop and report.
4. **Revised 2026-11-19:** jobs 3L, 4L (Mn₂In₂Se₅ LBTE on 12³, 9 T) and 3C, 4C
   (coherence on 16³, 9 T) replace the former 3a–c / 4a–c. LBTE cost grows
   linearly with the number of temperatures (≈ 45 min per T for job 1), and the
   LBTE/RTA ratio hardly depends on the mesh, so 12³ is enough for it. See the
   job table in the plan.
5. Job 5 (Mn₂In₂Se₅ VCA at x = 0.19, RTA 20³, with per-mode hdf5).

Job 1 taking 7–8 h is expected. Do not cancel running jobs: the driver writes its
JSON only at the end of each case.

### 1.2 What to send back

The `run_*/result_*.json` and `run.log` files (job 5 also `kappa_vca_rta_m20.hdf5`),
committed to the repository or copied to the E drive, and which jobs finished.
Do not change the report or the figures; that is done on the repository side.

## Task 2 — FC3 cutoff extension, no new DFT (decision of 2026-11-18 revised)

Why now: in the experiment Mn₂In₂Se₅ has the lower κ in both directions. The
calculation reproduces this in-plane and in the direction average, but its
Mn₂In₂Se₅ κ_zz (≈ 1.09) is above MnIn₂Se₄ (≈ 0.71), and that κ_zz is only an
upper bound: its FC3-cutoff steps (−11.9 %, −9.9 %) do not shrink. Report §1
("Which material is lower") has the numbers.

The procedure is in [`dft/fc3_cutoff_extend/README.md`](../dft/fc3_cutoff_extend/README.md), steps 0 and 1:

1. **Step 0: done.** pheasy `structure/atoms.py::wigner_seitz()` takes the shortest
   Cartesian image over the 27 neighbouring cells (code check, local side). The
   ±1-cell search is exact for both supercells (repository side, 2026-11-19): of
   11 248 Mn₂In₂Se₅ pairs below 7.96 Å and 8 432 MnIn₂Se₄ pairs below 8.01 Å, none
   differs from a ±2-cell brute-force search. No numerical cross-check through the
   pheasy API is needed.
2. **Step 1, Mn₂In₂Se₅ first:** pheasy refit of the existing 55 frames with the
   same method and settings as the 6.51 Å fit (**ALASSO**, `test/Mn2In2Se5/fit/fit_config.json`;
   the cutoff README says RFE-OLS for both, which is right only for MnIn₂Se₄), at 6.67 / 6.90 / 7.11 / 7.43 Å. Record the
   data/parameter ratio and the train/CV error for each.
3. κ for each new fc3: phono3py RTA, 300 K, isotopes on, 16³, `primitive_matrix="P"`;
   split κ_zz with `test/MnIn2Se4/kappa_zz_decompose.py`. The 6.51 Å reference is
   in-plane 2.8925, κ_zz rest term 0.9850.
4. Then MnIn₂Se₄: existing 46 frames, FC3 only, RFE-OLS as before (`test/MnIn2Se4/fit/fit_config.json`), shells up to ≈ 7.9 Å, 14³ with the
   4×4×2 FC2 (reference at 6.45 Å: in-plane 3.6474, κ_zz rest 0.6444).

Settings (answers of 2026-11-19):
- **Do not change the production project** (`Mn2In2Se5_441/kl-dft-cpu`,
  `cut3_candidates`). Run the refits as a separate scan, e.g. in
  `step5_fc/cutoff_scan_ext/`, from the same dataset and the same pheasy settings.
  If AutoZT can only do this by editing that project's configuration, ask the
  user first (iron rule 3) and then only append the new cutoffs, after a backup
  of the original file.
- **Fits: the 4 new cutoffs plus a refit at 6.51 Å** as an anchor, not all 10.
  The 6.51 Å refit must agree with the existing 6.51 Å fits within the known
  fit-to-fit scatter of ≈ 2 % (κ_zz rest 0.985–1.009, in-plane 2.89–3.00 at 16³).
- **FC2: keep the production 4×4×1 FC2 fixed** (`Mn2In2Se5/fc/fc2.hdf5`, md5
  `b611f139…`) for every cutoff and change only FC3, as the existing scan did
  (job 412029). That isolates the FC3 effect.

Decision rule (from that README):
- data/parameter ratio ≥ 3 and κ_zz steps below ≈ 5 %: converged, done, **no DFT**.
- ratio < 3 above ≈ 7 Å (expected for Mn₂In₂Se₅ at 7.4–7.9 Å from the earlier
  extrapolation, ≈ 2–2.5): **stop and report**. Only then are new DFT frames needed
  (step 2 there: ≈ 30–35 frames in the same 4×4×1 cell, ≈ 4–6 k core-h, ≈ 1 day),
  and that needs the user's go-ahead.

Send back: the fit logs (ratio, errors), the fc3.hdf5 of each cutoff (or their
location on the E drive) and the κ table. Do not change the report.
