# For the local AI: Mn₂In₂Se₅ FC3 cutoff, overfitting or physics? (2026-11-20)

Context: `test/Mn2In2Se5/cutoff_ext/README.md`. κ_zz rest falls by −10.8 % at
7.11 Å and −24.6 % at 7.43 Å, the in-plane value by −5.4 % and −13.0 %. Step 0
(minimum image) passed, so there are two explanations left: real long-range
interlayer terms, or the 7.11/7.43 Å fits absorbing noise (data/parameter ratio
3.55 / 2.69). The training error (0.408 → 0.351 %) cannot tell them apart,
because it always falls as parameters are added.

This round separates the two **with the existing 54 frames only**.

## Ground rules

- **No new DFT.** Step 2 (`make_frames.py`) is still not approved.
- Do not modify the production `step5_fc/` or `step5_fc/cutoff_scan_ext/`.
  Work in a new sibling directory `step5_fc/cutoff_cv/` on 3090, set up the
  same way as `cutoff_scan_ext` (same driver, same config except the dataset
  and the cutoff list, same `submit_fit_ext_gpu.sh` environment).
- Same fit method for every fit: pheasy-gpu, ALASSO, `PHEASY_RASR=BHH`, and
  the same hyperparameter selection as the scan. If ALASSO picks α
  internally, let it do so in every fit. Do not fix α by hand for some fits
  and not others.
- The GPU traps in `dft/fc3_cutoff_extend/README.md` still apply:
  `PHEASY_GPU_MODE=required`, pin the freest card, no
  `PHEASY_GPU_RESIDENT_FALLBACK`.
- Do not cancel 7709 (MnIn₂Se₄ step 1) or 7702. 7709 has priority on the GPU
  queue; these fits go in after it or alongside it when a card is free.
- Fix every random seed and write the frame indices of every subset into the
  output JSON, so each fit can be reproduced.
- Report after step A before submitting the rest.

## A. Set-up and one sanity fit (do first, then report)

1. Find out how `fc_fit_driver.py` reads the frames (the `dataset_dir`, the
   normalised dataset, `disp_matrix.pkl`). Find the least invasive way to fit on
   a subset of frames: a subset dataset directory of symlinks, or an index list
   the driver already accepts. Do not edit pheasy itself. Write down which way
   you chose.
2. Count the frames: 54 displaced frames plus the equilibrium reference
   (the plan says 55 in total). The reference frame stays in **every**
   training set and is never held out.
3. Write `predict_forces.py`. Given the fitted force constants of one fit
   (every order the fit contains, not only FC2+FC3) and a set of frames, it
   predicts
   `F_i = −Σ_j Φ²_ij u_j − ½ Σ_jk Φ³_ijk u_j u_k (− …higher orders if fitted)`
   on the full 144-atom supercell. Use the full (not compact) FC3, or expand it.
4. **Sanity check, mandatory:** run `predict_forces.py` on the **training**
   frames of the existing full 54-frame 7.43 Å fit. It must reproduce pheasy's
   reported train relative error 0.351 % to about 0.001 %. If it does not, the
   prediction or the error definition differs from pheasy's: stop and report.
   Do not go on with a metric that disagrees with pheasy's.
5. Run one hold-out fit (7.43 Å, fold 0, see B) to get the time per fit.
   Full-data fits took ≈ 50–70 min each at 6.67–7.43 Å.

## B. Test 1: hold-out validation

- Cutoffs: **6.51, 6.90, 7.11, 7.43 Å**.
- Split the 54 displaced frames into fixed random folds, with a fixed seed.
  Use **3 folds of 10 held-out frames**: each fit trains on 44 frames plus the
  reference and is tested on its 10. The three held-out sets are disjoint.
  Use the same three splits for every cutoff.
- That is 4 × 3 = **12 fits**.
- For each fit, record on its held-out frames:
  - the total-force RMSE (eV/Å) and the relative error as pheasy defines it;
  - **the anharmonic relative error**:
    `ε_anh = ‖F_DFT − F_pred‖ / ‖F_DFT − F_harm,prod‖`, where `F_harm,prod`
    uses the production FC2 (`b611f139…`) for every fit, so the denominator
    is the same everywhere. The total-force error is dominated by FC2 and is
    not sensitive to FC3. This one is the decision quantity;
  - the same two numbers on its training frames, for comparison.

**Reading:** compare the mean ε_anh(held-out) over the 3 folds between cutoffs,
with the fold-to-fold spread as the error bar.

- ε_anh(held-out) at 7.11/7.43 **higher** than at 6.90 by more than the
  spread, while the training error keeps falling → **overfitting**.
- ε_anh(held-out) flat or still falling → the extra parameters carry real
  information, not noise.

## C. Test 2: learning curve in κ

- The 12 fits of B already give κ at N = 44 frames (3 subsets per cutoff), and
  the scan gives N = 54.
- Add **N = 36** (36 displaced frames + reference, 2 random subsets, fixed
  seeds) at **6.51, 7.11 and 7.43 Å**: **6 more fits**. The data/parameter
  ratio at 7.43 Å is then ≈ 1.8, which is intended: it shows what an
  under-determined fit does.
- κ for every fit of B and C (18 in total): phono3py 4.4.0 RTA, 16³, 300 K,
  isotopes **on**, `primitive_matrix="P"`, production FC2 fixed and only FC3
  varying. Exactly the settings of `test/Mn2In2Se5/cutoff_ext/`. Split κ_zz
  with `test/MnIn2Se4/kappa_zz_decompose.py` and use **κ_zz rest**, never the
  total.
- Run this κ check first: the existing 54-frame κ at 6.51 Å must come out at
  2.9348 / 0.9801 (in-plane / κ_zz rest) again.

**Reading:** at each cutoff, put κ_zz rest and in-plane κ against N = 36 / 44
/ 54, with the subset-to-subset spread.

- At 7.43 Å, κ stable in N (the change from 44 to 54 and the subset spread
  both **< 5 %**) → the −24.6 % is determined by the data.
- κ still drifting with N, or a subset spread **≳ 10 %**, comparable to the
  drop itself → not determined by the data. Then more frames (step 2) are
  needed.
- 6.51 Å is the control: it should be stable in N. If it is not, report that
  before interpreting anything else.

## D. Verdict

| Test 1 (hold-out) | Test 2 (κ vs N) | verdict | next |
|---|---|---|---|
| error rises at 7.11/7.43 | κ drifts with N | **overfitting** | quote 6.90 Å as the last trusted cutoff, with the 6.51–6.90 spread as its uncertainty; no new DFT |
| flat or falling | κ stable in N | **real long-range terms**, cutoff not converged inside the 8 Å reach of 4×4×1 | stop; the user decides on 5×5×1 |
| anything else | | **undecided** | report; step 2 (more frames) is the candidate, the user decides |

Write the verdict as the table row it matches, with the numbers. Do not soften
or strengthen it.

## Output

- Scripts: `dft/fc3_cutoff_cv/` (subset maker, `predict_forces.py`, submit
  scripts, README with the subset method and the sanity-check result).
- Results: `test/Mn2In2Se5/cutoff_cv/`
  - `holdout.json`: per cutoff and fold, the frame indices and the train and
    held-out errors (total and ε_anh);
  - `learning_curve.json`: per cutoff and N, the frame indices, in-plane κ,
    κ_zz total, κ_zz rest, and the first Γ–Z point;
  - `README.md`: both tables, the verdict row, and anything that did not go to
    plan.
- The per-mode κ hdf5 and the fitted fc3.hdf5 files go to the E drive, not to
  git: `E:\backup\joint_research_project\xia\raw_hdf5\cutoff_cv\`.
- Commit, push to master, mirror to E, and report the commit hash.
