# For the local AI: tasks after 2026-11-07

Work through these in order. Report the result after each step before moving on.

## 1. Merge the tidied-up repository

The tidying-up lives on the branch `claude/jolly-curie-jsvosv` of `Andyw852/Mn2In2Se5-thermal`; fast-forward it into master first:

```bash
git clone https://github.com/Andyw852/Mn2In2Se5-thermal.git   # if you already have a local copy, cd into it
cd Mn2In2Se5-thermal
git fetch origin
git checkout master
git merge --ff-only origin/claude/jolly-curie-jsvosv
git push origin master
```

If `--ff-only` fails, something new has been pushed to master: stop and report, do not force-push.

After merging, check: the root contains only `Mn2In2Se5/ MnIn2Se4/ experiment/ analysis/ protocol/ README.md .gitignore`;
and `Mn2In2Se5/data/kappa/kappa_main_T100-800K.csv` has the same md5 as the file of the same name in the original `Mn2In2Se5_441-thermal`.

## 2. Upload the Mn₂In₂Se₅ force constants and DFT inputs (read-only copy; change nothing on the cluster)

The list is in `Mn2In2Se5/fc/README.md`. S5 ran on 3090 and S6 on hanhai25, so the files may be split across the two machines; find them with `ssh … ls`.

```bash
cd Mn2In2Se5-thermal/Mn2In2Se5/fc
# example: <S5 dir> = Mn2In2Se5_441/kl-dft-cpu/step5_fc on 3090
scp 3090:<S5 dir>/cutoff_scan/cut3_6p51/fc2.hdf5 .
scp 3090:<S5 dir>/cutoff_scan/cut3_6p51/fc3.hdf5 .
scp 3090:<S5 dir>/phono3py_disp.yaml .          # if absent, take it from step4_disp
# the ShengBTE run: <S6 dir>/shengbte_runs/cut3_6p51_m16x16x4/
scp hanhai25:<S6 dir>/shengbte_runs/cut3_6p51_m16x16x4/FORCE_CONSTANTS_2ND .
scp hanhai25:<S6 dir>/shengbte_runs/cut3_6p51_m16x16x4/FORCE_CONSTANTS_3RD .
scp hanhai25:<S6 dir>/shengbte_runs/cut3_6p51_m16x16x4/BTE.KappaTensorVsT_RTA .
scp jzzn:<S4 dir>/SPOSCAR .
ls -la                                         # if any single file is > 95 MB, stop and report; do not push
md5sum fc2.hdf5 fc3.hdf5 > md5.txt
```

Put the DFT inputs in `Mn2In2Se5/data/dft_inputs/` (all small text files):
`step1_std_opt/{INCAR,KPOINTS,workflow_method.txt}`, `step2_static/{INCAR,KPOINTS}`, `step4_disp/disp-00000/{INCAR,KPOINTS}`,
Also write a `dft_summary.txt`: the final maximum force of S1, the S2 band gap and total magnetic moment, the S4 equilibrium-frame residual force, the number of S4 frames, and the per-frame magnetic-moment distribution (the output of `sort | uniq -c`).
POTCAR is not uploaded (copyrighted).

```bash
git add Mn2In2Se5/fc Mn2In2Se5/data/dft_inputs
git commit -m "Mn2In2Se5: force constants (cut3 6.51 A) + DFT inputs"
git push origin master
```

## 3. Delete the two redundant repositories (their contents are now merged into this one)

| repository | destination |
|---|---|
| `Andyw852/Mn2In2Se5_441-thermal` | everything merged into `Mn2In2Se5/` (`config/` kept only the step6 step.conf; cluster configuration such as hpc.yaml was not moved) |
| `Andyw852/MnIn2Se4-thermal` | everything merged into `MnIn2Se4/archive_2026-10-28_fit-fc-thermal/` |

Delete only after step 1 has been checked; deletion is irreversible:

```bash
gh auth refresh -h github.com -s delete_repo
gh repo delete Andyw852/Mn2In2Se5_441-thermal --yes
gh repo delete Andyw852/MnIn2Se4-thermal --yes
```

Other repositories (taskflow, taskflow-v2.0, etc.) are out of scope; do not touch them.

## 4. Cluster: recalculate MnIn₂Se₄ with the same protocol

Follow `MnIn2Se4/redo_441/README.md` (material name `MnIn2Se4_441`; the protocol differs from Mn₂In₂Se₅ only in the items in that table).
Human-check the S1→S2→S3 and S4 generation stages; before submitting S4, verify SPOSCAR (112 atoms), KPOINTS (2 2 2), and ISPIN/MAGMOM in `disp-00000/INCAR`.

## 5. What not to do

- Do not rerun S6 for Mn₂In₂Se₅: the κ_zz mesh refinement and the phono3py cross-check are done on the repository side with the uploaded force constants;
  if the cross-check turns up problems, then decide how to follow up on the cluster.
