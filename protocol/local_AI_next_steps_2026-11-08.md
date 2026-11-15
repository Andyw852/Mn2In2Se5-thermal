# For the local AI: tasks for 2026-11-08

Steps 1–3 of the 11-07 list are done (repository merged, force constants uploaded, the two old repositories deleted). The checking turned up an export bug that affects every ShengBTE result,
so the order is adjusted as follows. Report after each step before moving on; the AGENTS.md iron rules still apply.

## 0. Background (in one sentence)

When S6 used hiphive to convert `fc3.hdf5` into ShengBTE's `FORCE_CONSTANTS_3RD`, it lost 34% of the third-order terms (the interlayer terms crossing the a₃′ boundary in the 4×4×1 supercell),
so the production κ of Mn₂In₂Se₅ (2.059 / 0.209) is superseded. With the same force constants the correct values are about 300 K RTA κ_xx ≈ 2.96 and κ_zz ≈ 0.84–1.1 (the κ_zz mesh is still converging).
The evidence is in `analysis/shengbte_crosscheck/README.md` in this repository.

## 1. Merge the AutoZT [0127] fix

Branch `claude/jolly-curie-jsvosv` (Andyw852/AutoZT, commit `25a67fd`). It changes only 8 files relative to master
(`fc_common.py`, `kappa_driver.py`, `fc_fit_driver.py`, two mlff files, `lattice_kappa.py`, CHANGELOG, a new test).

```bash
cd ~/software/AutoZT
git fetch origin
git checkout master
git merge --no-ff origin/claude/jolly-curie-jsvosv
python3 -m pytest -q tests/test_0127_shengbte_fc3_export.py
git push origin master
```

- The branch's original OOM auto-memory patch (old id [0107]) has been reverted: master later classified signal 9 as `mpi_killed` and routed it to human handling,
  and the two logics conflict, so it is no longer merged. The branch history contains one revert commit and one merge commit; the net change is only [0127].
- After merging, the AutoZT on the cluster must also be updated to the same version (the S5/S6 scripts are pushed by gen and come from the local repository files).

## 2. Stop the superseded S6 job (requires user approval)

The S6 job 357276 (jzzn) of `Mn2In2Se5_441` still uses the incomplete 3RD, so its result is superseded as well. First run `autozt -tt kl-dft-cpu -p Mn2In2Se5_441 brief`
to see what it is computing, then request a `stop` per iron rule 2.

## 3. Rerun S6 for Mn₂In₂Se₅ (revised 2026-11-08: no explicit MESH)

The originally written `MESH="16 16 16"` has been withdrawn: the ShengBTE iterative solve ran out of memory at 20×20×5 (750 GB), and 16×16×16 would be larger. The current settings
(identical for both materials): `MESH = auto`, `MESH_CONV_MAX_LENGTH = 50` (capped at 16×16×4), `T_MIN/T_MAX/T_STEP = 50/800/50`.
κ_xx is converged at 16×16×4; the κ_zz mesh refinement is done on the repository side with phono3py, see `analysis/shengbte_crosscheck/README.md`.

After each `conf --set`, re-check `autozt conf` and the file itself; after `retry` and before `start`, confirm that the generated `kappa_config.json`
has `mesh`, `mesh_conv_max_length` and the temperature list matching the above; if not, stop.

## 4. Cluster: recalculate MnIn₂Se₄ with the same protocol (can start now)

Follow `MnIn2Se4/redo_441/README.md` (material name `MnIn2Se4_441`). Two additions:

- **Write the memory request before submitting S4**, do not wait for it to blow up: use the `[submit]` memory settings from the Mn₂In₂Se₅ S4 run that finally succeeded (a whole node, or the value used at the time),
  and put them into the step4_disp step.conf of `MnIn2Se4_441` before `init` generates S4; after generation, verify that `disp-00000/submit.sh` contains `--mem`.
- S5/S6 must run only after [0127] is merged in step 1 (otherwise ShengBTE's 3RD will be incomplete again).

## 5. Reporting

- Step 1: the merge commit hash, the test results, and the AutoZT version on the cluster.
- Step 2: what job 357276 is doing and whether it has been stopped.
- Step 4: the S1 final force, the S2 band gap and total magnetic moment (expected ≈ 5 μB/cell).
