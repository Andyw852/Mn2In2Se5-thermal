# For the local AI: Mn₂In₂Se₅ scattering-rate figure at x = 0.19 (2026-11-22)

Context: `test/Mn2In2Se5/vca_m20_x019/README.md`. Since 2026-11-21 the Mn₂In₂Se₅ VCA
curves use the literature x = 0.19. Their numbers come from HanHai job 5 (JSON) plus
one grid point recomputed in the cloud, so no per-mode file was needed. One item is
still at x = 0.13: the Mn₂In₂Se₅ scattering-rate figure
(`test/figures/scattering/mn2in2se5/`, report §6, bottom row). It needs the per-mode
hdf5 that job 5 wrote on HanHai (`WRITE_KAPPA=1`), which was never uploaded.

This round copies that file in and regenerates the figure and the report. It involves
**no DFT and no new phonon calculation**, unless step 1 finds the file missing.

## Ground rules

- Do not edit any README, report text or script. The HTML caption reads x from the
  data (`vca_x.json`), so no text needs to change. If something does not fit, stop
  and report.
- Do not `git add -A`: the HanHai clone holds large untracked run directories. Add
  exactly the files listed in step 6.
- If any `assert` fires, stop. Do not commit, and send the full output.
- Python: the one job 5 used (`$HOME/software/pybin/python3`, phono3py 4.4.0), or any
  environment with h5py, phonopy and matplotlib. Step 5 also needs `seekpath`.

## Step 0 — repository state

The user merges `claude/gifted-clarke-dbnv8f` into `master` first. Then, in the clone
you work in (on HanHai: `$HOME/work/Mn2In2Se5-thermal`, the `REPO` of `submit.sh`):

```bash
git pull origin master
grep -c vca_x.json test/Mn2In2Se5/vca_m20_x019/postprocess.py      # must be >= 1
cat test/figures/scattering/mn2in2se5/vca_x.json                    # "x": 0.13 before this round
```

## Step 1 — find and check the job-5 file (HanHai)

`submit.sh` put job 5 in
`test/figures/kappa_vs_T/jobs_lbte_wigner/run_Mn2In2Se5_vca_m20_rta_T30-60-90-100-120-150-200-300-400-500-600-700/`
of the HanHai clone (the repository copy of that directory was renamed to
`run_Mn2In2Se5_vca_m20_rta`). The driver renamed phono3py's `kappa-m202020.hdf5` to
`kappa_vca_rta_m20.hdf5`.

```bash
cd $HOME/work/Mn2In2Se5-thermal/test/figures/kappa_vs_T/jobs_lbte_wigner
ls -ld run_Mn2In2Se5_vca_m20_rta*/ && ls -l run_Mn2In2Se5_vca_m20_rta*/*.hdf5
```

Check that it is the right file. Expect 12 temperatures, mesh 20³, 781 irreducible
points, and at 300 K the job-5 JSON values in-plane 1.46109 and κ_zz total 1.20477:

```bash
$HOME/software/pybin/python3 - <<'EOF'
import glob, h5py, numpy as np
fn, = glob.glob("run_Mn2In2Se5_vca_m20_rta*/kappa_vca_rta_m20.hdf5")
with h5py.File(fn) as f:
    T, k = list(f["temperature"][()]), f["kappa"][()]
    print(fn, "mesh", list(f["mesh"][()]), "mode_kappa", f["mode_kappa"].shape)
    print("T", T)
    i = T.index(300.0)
    print("300 K in-plane %.5f (expect 1.46109)  kzz %.5f (expect 1.20477)" % (0.5 * (k[i, 0] + k[i, 1]), k[i, 2]))
EOF
```

Expected: `mode_kappa (12, 781, 27, 6)`, T = 30 60 90 100 120 150 200 300 400 500 600
700, and both 300 K values within ±1×10⁻⁵ of the expected numbers.

**If the file is missing, or only `kappa-m202020.hdf5` is there**, report that. A
`kappa-m202020.hdf5` that passes the check above can be used under the new name. If
there is no file at all, resubmit job 5 unchanged. It is the row "5" in
`test/figures/kappa_vs_T/jobs_lbte_wigner/README.md`, about 5 min on 64 threads:

```bash
cd $HOME/work/Mn2In2Se5-thermal/test/figures/kappa_vs_T/jobs_lbte_wigner
sbatch --export=ALL,MAT=Mn2In2Se5,CASE=vca,MESH=20,STEPS=rta,WRITE_KAPPA=1,TEMPS="30 60 90 100 120 150 200 300 400 500 600 700" submit.sh
```

Its new `result_Mn2In2Se5_vca_m20.json` must equal the committed one
(`run_Mn2In2Se5_vca_m20_rta/result_Mn2In2Se5_vca_m20.json`) in `rta.population`. Then
run the check above again.

## Step 2 — copy it into the repository

```bash
cd $HOME/work/Mn2In2Se5-thermal
cp test/figures/kappa_vs_T/jobs_lbte_wigner/run_Mn2In2Se5_vca_m20_rta*/kappa_vca_rta_m20.hdf5 \
   test/Mn2In2Se5/vca_m20_x019/
md5sum test/figures/kappa_vs_T/jobs_lbte_wigner/run_Mn2In2Se5_vca_m20_rta*/kappa_vca_rta_m20.hdf5 \
       test/Mn2In2Se5/vca_m20_x019/kappa_vca_rta_m20.hdf5
```

The two md5 sums must agree. Size about 16 MB (the x = 0.13 12-T file in
`test/figures/kappa_vs_T/raw/` is 15.8 MB).

## Step 3 — post-process

```bash
cd test/Mn2In2Se5/vca_m20_x019 && $HOME/software/pybin/python3 postprocess.py
```

Expected output, in this order:

```
wrote test/figures/kappa_vs_T/data_mn2in2se5_vca.json
300 K: in-plane 1.4611 (baseline 2.998, 2.05x), kzz rest 0.7534 (baseline 0.9837, 1.31x)
full hdf5 vs single-point route: max abs deviation ...e-06 W/mK (JSON totals have 5 decimals)
wrote kappa_vca_300K.hdf5 (300 K slice), copied it to test/figures/scattering/mn2in2se5/kappa_vca.hdf5 and set vca_x.json there to 0.19
```

The deviation line is the real test. It decomposes the full file at all 12
temperatures and compares it with the single-point route. It must be ≤ 2×10⁻⁵; the
script asserts this. In the dry run with the x = 0.13 file it was 8×10⁻⁶.

Afterwards `git status` should show:

- **modified**: `test/figures/scattering/mn2in2se5/kappa_vca.hdf5` and `test/figures/scattering/mn2in2se5/vca_x.json` (now `"x": 0.19`);
- **new**: `test/Mn2In2Se5/vca_m20_x019/kappa_vca_rta_m20.hdf5` and `kappa_vca_300K.hdf5`;
- **unchanged**: `test/figures/kappa_vs_T/data_mn2in2se5_vca.json` and `test/Mn2In2Se5/vca_m20_x019/decomposition.json`. If either shows a diff, send it and stop.

## Step 4 — scattering figure

```bash
cd $HOME/work/Mn2In2Se5-thermal/test/figures/scattering/mn2in2se5
$HOME/software/pybin/python3 plot_scattering_figure.py
```

It prints "Temperature point: baseline 300 K, vca 300 K" and rewrites
`scattering_compare.png / .pdf`. Expected look: the left panel (3-phonon) is nearly
unchanged. In the right panel the red VCA disorder points sit slightly higher than
before, because g rises from 0.103 / 0.035 to 0.125 / 0.052.

## Step 5 — HTML report

Run this where the report has been built before. `export_plot_data.py` needs
`seekpath` (`pip install seekpath`), and the Chinese version needs the CJK font used
so far.

```bash
cd $HOME/work/Mn2In2Se5-thermal
python3 report/make_report_html.py && python3 report/make_report_html_zh.py
cat report/plot_data/scattering_vca_x.json          # Mn2In2Se5 must be 0.19
grep -o "bottom row Mn₂In₂Se₅ (x = [^)]*)" report/calc_vs_experiment_report.html
```

The grep must print `bottom row Mn₂In₂Se₅ (x = 0.19)`, without the "not yet in the
repository" note. The Chinese page reads `下排 Mn₂In₂Se₅（x = 0.19）`.

If you build the report on a different machine from HanHai, first move the step 2–4
results there: the files in step 6, without the `report/` ones.

## Step 6 — commit and push

```bash
git add test/Mn2In2Se5/vca_m20_x019/kappa_vca_rta_m20.hdf5 \
        test/Mn2In2Se5/vca_m20_x019/kappa_vca_300K.hdf5 \
        test/figures/scattering/mn2in2se5/kappa_vca.hdf5 \
        test/figures/scattering/mn2in2se5/vca_x.json \
        test/figures/scattering/mn2in2se5/scattering_compare.png \
        test/figures/scattering/mn2in2se5/scattering_compare.pdf \
        report/
git status --short      # nothing else staged; send this output
git commit -m "data: Mn2In2Se5 VCA x = 0.19 per-mode file (job 5); scattering figure and report at x = 0.19"
git push origin master
```

`report/` contains only generated files (HTML, `figures/`, `plot_data/`). PDFs whose
only change is the timestamp may be committed with the rest.

## What to send back

1. The output of the step 1 check (file path, shapes, T list, the two 300 K values).
   If the file was missing, say so, and give the job ID if job 5 was resubmitted.
2. The two md5 sums from step 2.
3. The full output of `postprocess.py`, in particular the deviation line.
4. The `git status --short` before the commit, and the commit hash after the push.
5. Anything that did not match the expectations above, verbatim.
