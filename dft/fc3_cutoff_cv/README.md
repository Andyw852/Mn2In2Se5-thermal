# FC3 cutoff overfitting test — hold-out validation and learning curve

Protocol: [`protocol/local_AI_next_steps_2026-11-20.md`](../../protocol/local_AI_next_steps_2026-11-20.md).
Work directory on 3090: `step5_fc/cutoff_cv/` (production `step5_fc/` and
`step5_fc/cutoff_scan_ext/` untouched).

## Subset method (A.1)

`fc_fit_driver.py` normalises whatever sits in `dataset_dir` (symlinks are fine:
it only opens `disp-*/vasprun.xml`). So each subset is a directory of **symlinks**
to the original frame folders, plus the equilibrium frame `disp-00000`, which is in
every training set and never held out. `make_subsets.py` builds them:

| subset | frames | used by |
|---|---|---|
| `fold0/1/2` | 44 training (+ equilibrium); 10 held out | 4 cutoffs × 3 folds = 12 fits |
| `n36_0/1` | 36 training (+ equilibrium) | 3 cutoffs × 2 subsets = 6 fits |

No pheasy internals were changed. Frame indices, seeds and the exact splits are in
`subsets_manifest.json` next to the runs on the cluster.

## Sanity check (A.4) — passed

`predict_forces.py` predicts `F_i = -Σ_j Φ²_ij u_j - ½ Σ_jk Φ³_ijk u_j u_k` on the
144-atom supercell from a fit's own FC2 and FC3.

On the **training frames of the existing full 54-frame 7.43 Å fit**:

| quantity | value |
|---|---:|
| prediction time | 12.5 s for 54 frames |
| this script, `re = ‖F_pred − F_dft‖/‖F_dft‖` | **0.003511 = 0.3511 %** |
| pheasy's reported training relative error | **0.3512 %** |

Agreement to 0.0001 %, i.e. the prediction and the error definition match pheasy.
(Subtracting the equilibrium frame's forces gives 70 %, which confirms that the
driver's normalised dataset already has that residual removed — do not subtract it
again.)

## A800 fallback (2026-11-20 evening)

The 3090 `gpu` partition shares its single node with `cpu192`, and a long AutoZT
job holds it, so the whole CV study was blocked. The same 18 runs were therefore
rebuilt on A800:

- `subsets/` are **array slices** of the normalised dataset (`dataset_disps.npy` /
  `dataset_forces.npy`), which the driver accepts as the `arrays` signature — no
  vasprun frames needed, and the same splits as on 3090 (from
  `subsets_manifest_3090.json`, the authoritative one).
- **Dataset convention (found the hard way):** in `dataset_*.npy` the **last**
  frame is the equilibrium (zero displacement, zero force); indices 0–53 are the
  54 displaced frames. The 3090 manifest names frames by `disp-XXXXX` id (1–54),
  so frame id *f* maps to array index *f−1*, and the equilibrium is appended last.
  With that fixed, `prep` on A800 reproduces the production normalisation exactly
  (44 training frames, RMS displacement 0.0172 Å, equilibrium residual 0).
- `submit_cv_fit_a800.sh` is the same job on the `a800` partition (env `mace`,
  `--gres=gpu:a800:1`), with the same GPU rules.
- `rebuild_a800.py` (on A800) rebuilds the subsets/runs from the 3090 manifest.

Status: A800's queue is congested by other users (many pending jobs), so the first
fit is queued rather than running.

## Runs

`submit_cv_fit.sh` runs `prep → fit → post` for one run directory on the GPU
partition, with the same GPU rules as `dft/fc3_cutoff_extend/`:
`PHEASY_GPU_MODE=required`, the freest card pinned, no
`PHEASY_GPU_RESIDENT_FALLBACK`, `PHEASY_FIT_METHOD=ALASSO`, `PHEASY_RASR=BHH`.

```bash
sbatch --export=ALL,RUN=c0743_f0 submit_cv_fit.sh
```

Results land in `runs/<RUN>/` and the fitted `fc2.hdf5`/`fc3.hdf5` are copied to
`fits/<RUN>/`.
