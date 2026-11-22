# For the local AI: Mn₂In₂Se₅ FM vs AFM Γ phonons, a spin–lattice coupling screen (2026-11-22)

Context and decision rule: `dft/spin_phonon_gamma/README.md`. In short, we compute the
Γ phonons of the 9-atom primitive cell with ferromagnetic (FM, production) and
antiferromagnetic (AFM: the two Mn layers antiparallel) order, at the same geometry and
with the production DFT settings. We also compute three total energies in a 2×1×1 cell
with in-plane stripe order. The analysis script compares the two phonon sets and gives
a verdict by rules fixed in advance.

This round needs **new DFT**: 113 small VASP runs (110 with 9 atoms, 3 with 18 atoms),
roughly 100–300 core-hours in total. The user has approved this test.

## Ground rules

- Use the production settings and POTCAR. Do not change any INCAR tag that `prepare.py`
  writes. If something fails, report it; do not tune.
- Work directory **outside** the repository, e.g. `$HOME/work/spin_phonon_gamma_run/`
  on jzzn. WAVECARs and OUTCARs are not committed.
- POTCAR: the same file as the Mn₂In₂Se₅ 4×4×1 S4 run (species order Mn, In, Se). This
  is the POTCAR that `dft/fc3_cutoff_extend/submit_frames.sh` also asks for.
- Python: any environment with phonopy and h5py, e.g. the one used for job 5
  (`$HOME/software/pybin/python3`). Use the same one for `prepare.py` and `analyse.py`.
- Stop at every failed check below and report. Do not continue to the next step.
- `prepare.py` and `analyse.py` can run on any machine with a clone; only VASP needs jzzn.
  If the clone is elsewhere, copy the work directory to jzzn after step 1. For step 4,
  copy it back; only the OUTCARs and the JSON files are read.

## Step 0 — repository and self-test

```bash
cd <clone of Mn2In2Se5-thermal> && git pull origin master
python dft/spin_phonon_gamma/analyse.py --repo . --selftest
```

Expected output (the numbers may differ in the last digits):

```
selftest FM: max |w - w_production| = 1.3e-08 THz
selftest AFM-like: max |w - direct diagonalisation| = 1.5e-07 THz
selftest AFM-like perturbation: rms delta 0.92 %, max 2.09 %, noise 9.1e-13 %
selftest passed
```

If it does not end with `selftest passed`, the phonopy in that environment behaves
differently: report the full output.

## Step 1 — prepare the runs

```bash
python dft/spin_phonon_gamma/prepare.py --repo . --out $HOME/work/spin_phonon_gamma_run \
       --potcar <path to the production Mn2In2Se5 S4 POTCAR>
```

Expected output: "113 VASP runs" and these MAGMOM lines:

```
  FM           MAGMOM = 5 5 7*0
  AFM          MAGMOM = 5 -5 7*0
  2x1_FM       MAGMOM = 5 5 5 5 14*0
  2x1_stripeA  MAGMOM = 5 -5 5 -5 14*0
  2x1_stripeB  MAGMOM = 5 -5 -5 5 14*0
```

The script stops by itself if the INCAR tags differ from the production S2 static INCAR,
or if the POTCAR species are not Mn, In, Se in that order.

## Step 2 — the two reference cells first

```bash
cd $HOME/work/spin_phonon_gamma_run
cp <clone>/dft/spin_phonon_gamma/run_config.sh .
sbatch run_config.sh FM ref
sbatch run_config.sh AFM ref
sbatch run_config.sh energy_2x1 all
```

When the two ref runs are done, check them:

```bash
for c in FM AFM; do echo "== $c"; grep "energy(sigma->0)" $c/ref/OUTCAR | tail -1
  grep "number of electron" $c/ref/OUTCAR | tail -1
  grep -A 5 "magnetization (x)" $c/ref/OUTCAR | tail -2
  awk '/TOTAL-FORCE/{f=1;n=0;m=0;next} f&&/---/{if(n>0)f=0;next} f{n++;for(i=4;i<=6;i++){v=$i<0?-$i:$i;if(v>m)m=v}} END{print "max |F| =",m," eV/A"}' $c/ref/OUTCAR
done
```

| check | FM | AFM |
|---|---|---|
| total magnetization | 9.98 ± 0.05 (production S1: 9.9844) | \|M\| < 0.05 |
| Mn moments (ions 1, 2) | both ≈ +4.5 … +4.7 | ion 1 > +4, ion 2 < −4 |
| max residual force (`TOTAL-FORCE` block) | < 0.005 eV/Å | report the value (it is part of the result) |

If FM fails, stop: the setup differs from production. If the AFM run has fallen back to
FM, set `NUPDOWN = 0` in `AFM/ref/INCAR` only, rerun it and report that you did.

## Step 3 — the displaced cells

```bash
sbatch run_config.sh FM disp
sbatch run_config.sh AFM disp
```

Each job runs 54 cells one after another. Each cell starts from `../ref/WAVECAR`, so
it stays in the intended spin state and converges faster. A job that stops can be
resubmitted with the same command; finished runs are skipped.

## Step 4 — analysis

```bash
python <clone>/dft/spin_phonon_gamma/analyse.py --repo <clone> --work $HOME/work/spin_phonon_gamma_run \
       | tee $HOME/work/spin_phonon_gamma_run/analyse.log
```

If it prints `PROBLEMS`, it lists the affected runs: unfinished, NELM reached, or a Mn
moment with the wrong sign or below 4 μB. Fix those runs or report them before going on.

## Step 5 — commit the results

Copy into the repository and commit only these files:

```bash
mkdir -p <clone>/dft/spin_phonon_gamma/results
cp $HOME/work/spin_phonon_gamma_run/{results.json,analyse.log,manifest.json} <clone>/dft/spin_phonon_gamma/results/
cd <clone> && git add dft/spin_phonon_gamma/results/ && git commit -m "dft: Mn2In2Se5 FM vs AFM Gamma phonons (spin-lattice screen)" && git push origin master
```

Do not edit any README or report text.

## What to send back

1. Step 0 output.
2. The step 2 checks (energies, total magnetization, Mn moments, max residual force) for FM and AFM.
3. The full `analyse.log`, including the mode table, the noise, the energies, the VERDICT
   and any NOTE or PROBLEMS lines.
4. Total core-hours used, and anything that did not match the expectations above.
