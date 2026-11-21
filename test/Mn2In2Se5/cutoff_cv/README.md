# FC3 cutoff: overfitting or real long-range terms? (Mn₂In₂Se₅)

Protocol: [`protocol/local_AI_next_steps_2026-11-20.md`](../../../protocol/local_AI_next_steps_2026-11-20.md).
Scripts: `dft/fc3_cutoff_cv/`. Work dirs: 3090 `step5_fc/cutoff_cv/`, A800
`~/work/mn_cv/` (A800 abandoned, see below). No new DFT was run.

## What was done

* Subsets are **array slices** of the normalised dataset (the `arrays` signature the
  driver accepts). Frames: 0–53 displaced, 54 = equilibrium (zero displacement, zero
  force), equilibrium in every training set and never held out.
* Test 1: 3 folds of 10 held-out frames, seed 20261120, same splits for every cutoff
  → 4 cutoffs × 3 folds, 44 training frames each.
* Test 2: N = 36 subsets, seeds 20261122/20261123, shared by 6.51 / 7.11 / 7.43 Å.
* Every fit: pheasy-gpu ALASSO, `PHEASY_RASR=BHH`, `PHEASY_GPU_MODE=required`,
  freest card pinned. The cutoff is set through **`pheasy_c3_cutoff`** (the driver
  only uses `cut3_candidates` in scan mode; a first attempt that set only the
  candidate list silently fitted every run at 6 Å and was discarded).
* κ: 16³, 300 K, isotopes on, `primitive_matrix="P"`, production FC2
  `b611f139…` fixed, only FC3 varying. κ_zz split with
  `test/MnIn2Se4/kappa_zz_decompose.py`, **rest term only**.

## Test 1 — hold-out error (the decision quantity)

`ε_anh = ‖F_pred − F_DFT‖ / ‖F_DFT − F_harm,prod‖`, harmonic denominator from the
production FC2 for every fit.

| cutoff | ε_anh held-out (mean of folds) | fold spread | training re |
|---|---:|---:|---:|
| 6.51 | 0.07528 | 0.00345 | 0.00403 |
| 6.90 | 0.07611 | 0.00343 | 0.00391 |
| 7.11 | 0.07806 | 0.00356 | 0.00365 |
| 7.43 | 0.08099 (2 folds) | 0.00496 | 0.00334 |

The training error **falls** monotonically with the cutoff while the held-out error
**rises**. The N=36 subsets confirm the direction: 7.43 Å gives 0.0871 / 0.0831 at
36 frames against 0.0810 at 44.

**Paired by fold** (added 2026-11-21). The three folds are the same frame splits at
every cutoff, so the right comparison is the step within each fold, not the mean
against the fold-to-fold spread (that spread is mostly which frames are held out,
and it is common to all cutoffs):

| held-out ε_anh step | fold 0 | fold 1 | fold 2 |
|---|---:|---:|---:|
| 6.51 → 6.90 | +0.00060 | +0.00097 | +0.00095 |
| 6.90 → 7.11 | +0.00216 | +0.00177 | +0.00190 |
| 7.11 → 7.43 | +0.00466 | +0.00261 | — |

All eight steps have the same sign. The held-out error already rises from 6.51 to
6.90 Å, in every fold.

## Test 2 — κ versus number of frames

κ_zz rest term (in-plane in the same cell). N=54 is the full-frame reference from
`test/Mn2In2Se5/cutoff_ext`.

| cutoff | N = 36 | N = 44 | N = 54 | change 44→54 | subset spread at N=44 |
|---|---:|---:|---:|---:|---:|
| **6.51 (control)** | 0.888 (2.739) | 0.976 (2.909) | **0.980 (2.935)** | **+0.4 %** | 8.4 % |
| 6.90 | — | 0.831 (2.772) | 0.948 (2.873) | +14.1 % | 7.2 % |
| 7.11 | 0.701 (2.507) | 0.802 (2.670) | 0.875 (2.777) | +9.1 % | 3.8 % |
| 7.43 | 0.586 (2.204) | 0.620 (2.355) | 0.739 (2.555) | **+19.2 %** | **18.6 %** |

**Anchor reproduced**: the 6.51 Å N=44 mean (0.976 / 2.909) matches the N=54 values
(0.980 / 2.935) to 0.4 % / 0.9 %, inside the subset spread. The control is stable in
N (44→54 = +0.4 %), as the protocol requires before anything else is interpreted.

At 7.43 Å the κ moves +19 % between 44 and 54 frames and the subset spread is 18.6 %
— both far outside the 5 % stability band. The −24.6 % drop reported in
`cutoff_ext` is therefore **not determined by the data**.

6.90 Å fails the same test: +14.1 % from 44 to 54 frames, also outside the 5 %
band. Paired by fold at N = 44, 6.51 → 6.90 Å moves κ_zz rest by −19.9 / −13.7 /
−10.7 % and in-plane by −4.7 / −5.0 / −4.5 %; with all 54 frames the same step is
only −3.3 % / −2.1 %. The step shrinks as frames are added, so at 6.90 Å it is
still set by the amount of data.

## Verdict (revised 2026-11-21)

Protocol section D, row 1:

| Test 1 (hold-out) | Test 2 (κ vs N) | verdict |
|---|---|---|
| error rises from 6.90 Å on, in every fold | κ drifts with N at 6.90, 7.11 and 7.43 Å | **overfitting** |

**Mn₂In₂Se₅: every FC3 cutoff above 6.51 Å is overfitted with the 54 frames.**
6.51 Å is the only cutoff that passes both tests, so it stays the production cutoff
and the quoted values do not change. The first version of this section named
6.90 Å as the last trusted cutoff; that followed the wording pre-written in the
protocol, not the data, and is withdrawn.

**Why every drop is negative.** Fitting noise in Φ₃ enters the three-phonon rates
through |V₃|², so it raises the scattering on average and biases κ **down**. That is
what the data show: at every cutoff, fewer frames give a lower κ (6.51 Å: 0.888 at
N = 36 → 0.976 at 44 → 0.980 at 54). The 54-frame value at 6.90 Å (0.948) is
therefore a **lower bound** on what a converged fit at 6.90 Å would give.

**Cutoff uncertainty to quote** (one-sided):

| quantity (16³, 300 K, isotopes on) | value at 6.51 Å | change at 6.90 Å |
|---|---:|---:|
| κ_zz rest of zone | **0.980** | no more negative than **−3.3 %** |
| in-plane κ | **2.935** | no more negative than **−2.1 %** |

Above 6.90 Å the existing frames cannot resolve anything. The −10.8 % (7.11 Å) and
−24.6 % (7.43 Å) values in `cutoff_ext` must not be quoted. No new DFT is justified
by this test.

The cutoff steps that the data do resolve shrink: κ_zz rest −11.8 % (5.68 → 6.12 Å),
−9.9 % (6.12 → 6.51 Å, data/parameter ratio 9.1 → 6.1), then at most −3.3 % at
6.90 Å.

Optional check, no DFT: 3 hold-out folds at 6.12 Å (about 3 GPU hours) would show
whether 6.51 Å is at the minimum of the held-out error or already slightly past it.

## Caveats and what did not go to plan

1. **7.43 Å has 2 folds, not 3.** The third fold fit was rejected by pheasy's
   FISTA convergence gate, and the retry hit `GPU mode is required but CUDA is
   unavailable` because an AutoZT job took the card. Its two folds already straddle
   far more than the stability band, so the verdict does not depend on it.
2. **The cutoff is set by `pheasy_c3_cutoff`.** The first four "completed" fits were
   silently all at 6 Å (identical `fc3.hdf5` md5) because only `cut3_candidates` was
   set; they were deleted and everything was refitted.
3. **A800 was set up and abandoned** at the user's request once 3090's GPU freed;
   5 of the missing κ runs were finished on HanHai instead of waiting for the 3090
   `cpu192` node held by a long AutoZT job.
4. **MnIn₂Se₄ is not covered by this protocol.** Its cutoff scan (6.45/6.79/7.19/7.74,
   46 frames) gave κ_zz rest +8.0 % / +1.5 % / −27.7 % against the 6.45 Å anchor, with
   a data/parameter ratio of 2.39 at 7.74 Å. The first version of this caveat called
   7.19 Å "the last trusted MnIn₂Se₄ cutoff". That does not follow: the Mn₂In₂Se₅
   test shows overfitting from the first step above the production cutoff, MnIn₂Se₄
   has fewer frames, and a non-monotonic +8.0 % / +1.5 % sequence looks like noise.
   MnIn₂Se₄ stays at 6.45 Å with its earlier ≈ −10 % estimate. Two things are still
   open: (a) the extension results are not in the repository yet; (b) the production
   6.45 Å FC3 is an RFE-OLS fit while the extension used ALASSO. If the anchor
   3.6474 / 0.6444 is the production fit and not an ALASSO refit, the extension points
   are compared against a different fit method.
5. **Frame numbering.** "0–53" above are array indices; the JSON files use the frame
   ids 1–54 (`disp-00001` … `disp-00054`, id = index + 1), and the equilibrium frame
   is `disp-00000`. Both conventions describe the same frames.

## Files

| file | contents |
|---|---|
| `holdout.json` | per fit: cutoff, split, frame ids, `re_total` and `ε_anh` on train and held-out frames |
| `learning_curve.json` | per fit: cutoff, N, split, in-plane κ, κ_zz total, κ_zz rest, first Γ–Z point; plus the N=54 reference |
| per-mode κ hdf5 | backup drive `E:\backup\joint_research_project\xia\raw_hdf5\cutoff_cv\` |
| fitted fc3.hdf5 (17) | same backup directory |

Unicode note: the Python source of `kappa_one_fc3.py` deliberately avoids non-ASCII.
