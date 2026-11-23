# For the local AI: Mn₂In₂Se₅ paramagnetic (spin-disordered) FC2 and its κ (2026-11-23)

Context, method, literature and decision rules: `dft/spin_disorder_fc2/README.md`.

In short: the Γ screen (`dft/spin_phonon_gamma/results/`) found strong spin–lattice
coupling. We now compute the harmonic force constants of the paramagnet in the
production 4×4×1 supercell:
- 8 spin-disordered configurations, whose pair correlations are exactly zero in the
  first five Mn–Mn shells;
- an FM control and an AFM bound.

Afterwards, the cloud side reruns the intrinsic κ with each FC2 and the production fc3.

This round needs **new DFT**: 43 VASP single points of 144 atoms. At the 2–3 h × 48
cores per frame measured for the same cell (`dft/fc3_cutoff_extend/README.md`), that is
about **4–6 k core-hours**, possibly up to ~8 k if the spin-disordered SCF is slower.
The user has approved plan A.

## Ground rules

- Use the production settings and POTCAR. Do not change any INCAR tag that `prepare.py`
  writes, and do not change `spin_configs.json` or `frames.npz`. If something fails,
  report it; do not tune.
- Work directory **outside** the repository, e.g. `$HOME/work/spin_disorder_fc2_run/`
  on jzzn. WAVECAR and OUTCAR are not committed. Afterwards, back the OUTCARs up to the
  E: drive as for the Γ test.
- POTCAR: the same file as the Mn₂In₂Se₅ 4×4×1 S4 run (species order Mn, In, Se).
- Python: an environment with phonopy, phono3py, **symfc** and h5py (e.g. `$HOME/software/pybin/python3`) for `prepare.py` and
  `analyse.py`. If `import symfc` fails, install it into that
  environment (`pip install symfc`) and report the versions.
- Stop at every failed check and report. Do not continue to the next step.
- Do not `git add -A`. Add exactly the files listed in step 6. Do not edit any README,
  report text or script.

## Step 0 — repository and self-tests

```bash
cd <clone of Mn2In2Se5-thermal> && git pull origin master
python -c "import phonopy, symfc; print(phonopy.__version__, symfc.__version__)"
python dft/spin_disorder_fc2/analyse.py --repo . --selftest
python dft/spin_disorder_fc2/prepare.py --repo . --list-configs
```

The self-test must end with `selftest passed`. Its gate table shows
`G1 ... FAIL`, followed by the line saying this is expected on synthetic data. The other
two gates must pass. Expected lines (the last digits may differ):

```
selftest FM_ctl: fitted (read back from hdf5) vs exact synthetic FC2: rms |dw/w| 0.05 %, below 3 THz 0.08 %
selftest AFM   : fitted (read back from hdf5) vs exact synthetic FC2: rms |dw/w| 0.08 %, below 3 THz 0.12 %
selftest PM    : fitted (read back from hdf5) vs exact synthetic FC2: rms |dw/w| 0.20 %, below 3 THz 0.32 %
selftest: synthetic PM vs FM shift (rms < 3 THz) 1.96 % (must be resolved)
selftest passed
```

`--list-configs` must print `192 space-group permutations`. Then, for PM01–PM08, it must
print `<s s> = +0.000` in all five shells and `E_J - E_FM = -5.74 meV/Mn`, with AFM at
`-1.000 +1.000 -1.000 -1.000 +1.000` and `-10.44`. It takes about a minute.

## Step 1 — prepare the runs

```bash
python dft/spin_disorder_fc2/prepare.py --repo . --out $HOME/work/spin_disorder_fc2_run \
       --potcar <path to the production Mn2In2Se5 S4 POTCAR>
```

Expected: `wrote ...: 43 VASP runs of 144 atoms in 10 sets` and the same correlation
table. The script stops by itself in three cases:
- an INCAR tag differs from the production S4 INCAR;
- KPOINTS is not Γ 2×2×2;
- the POTCAR species are not Mn, In, Se in that order.

## Step 2 — the ten reference cells first

```bash
cd $HOME/work/spin_disorder_fc2_run
cp <clone>/dft/spin_disorder_fc2/run_set.sh .
for s in FM_ctl AFM PM01 PM02 PM03 PM04 PM05 PM06 PM07 PM08; do sbatch run_set.sh $s ref; done
```

When all ten are done:

```bash
python <clone>/dft/spin_disorder_fc2/analyse.py --repo <clone> --work $HOME/work/spin_disorder_fc2_run --refs-only \
       | tee refs.log
```

It must end with `ref checks passed`:
- every run finished, below NELM;
- all 32 Mn moments have |m| > 4 μB and the intended sign;
- the total magnetization is ≈ 159.7 for FM_ctl and ≈ 0 for the others.

Expected, not gates (report the values):
- AFM ≈ −10.4 meV/Mn (the 9-atom value);
- PM sets near the Heisenberg column (−5.74);
- FM_ctl static force < 0.005 eV/Å;
- AFM and PM static forces of a few 0.01 eV/Å.

If a PM or AFM ref has moments with the wrong sign, do not change its INCAR. Stop and
report which set and which ions. If the job hits the 2-day limit and the queue rejects
`--time=2-00:00:00`, set 1 day and resubmit: finished runs are skipped.

## Step 3 — the displaced cells

```bash
for s in FM_ctl AFM PM01 PM02 PM03 PM04 PM05 PM06 PM07 PM08; do sbatch run_set.sh $s disp; done
```

FM_ctl has 9 frames, AFM 8, each PM 2. Every frame starts from its set's `ref/WAVECAR`.
A stopped job can be resubmitted with the same command.

## Step 4 — fit and analysis

```bash
python <clone>/dft/spin_disorder_fc2/analyse.py --repo <clone> --work $HOME/work/spin_disorder_fc2_run \
       | tee $HOME/work/spin_disorder_fc2_run/analyse.log
```

If it prints `PROBLEMS`, fix nothing; report the list. It must end with
`ALL GATES PASSED`. If it prints `GATE FAILED`, still do steps 5 and 6 (the
fc2 files are needed to find the cause), but say so in the commit message and send the
whole log. It writes:
- `results/{results.json, fc2_FM_ctl.hdf5, fc2_AFM.hdf5, fc2_PM.hdf5}`;
- `jackknife/fc2_PM_without_PM0k.hdf5` (8 files).

## Step 5 — back up

Copy every OUTCAR, the `manifest.json`, the `displacements.json` files and the logs to
`E:\backup\joint_research_project\xia\raw_hdf5\spin_disorder_fc2\`, as for the Γ test.
Spot-check the md5 of three files. WAVECARs are not needed.

## Step 6 — commit the results

The κ runs (`kappa_extrinsic.py bases` with each FC2 and the production fc3, then
`kappa_summary.py`) are done on the cloud side from these files. Do not run them.

```bash
W=$HOME/work/spin_disorder_fc2_run
cd <clone>
mkdir -p dft/spin_disorder_fc2/results/jackknife
cp $W/results/{results.json,fc2_FM_ctl.hdf5,fc2_AFM.hdf5,fc2_PM.hdf5} dft/spin_disorder_fc2/results/
cp $W/jackknife/fc2_PM_without_PM0?.hdf5 dft/spin_disorder_fc2/results/jackknife/
cp $W/{analyse.log,refs.log,manifest.json} dft/spin_disorder_fc2/results/
git add dft/spin_disorder_fc2/results/
git commit -m "dft: Mn2In2Se5 paramagnetic spin-disordered FC2" && git push origin master
```

That is 11 fc2 files of about 1.5 MB each.

## What to send back

1. The step 0 output: versions, the self-test tail and the `--list-configs` table.
2. `refs.log`.
3. The full `analyse.log`, including the gates line.
4. The commit hash on master.
5. Core-hours used, wall time, the number of SCF steps of a PM ref compared with the
   FM_ctl ref, and anything that did not match the expectations above.
