#!/bin/bash
# hanhai_symfc_run.sh -- independent symfc (fc2+fc3) fit for MnIn2Se4 at selected cutoffs.
BASE=$HOME/work/mninse4_189_kappa
S=$BASE/symfc_src
for c3 in 2.85 4.30 7.00; do
  D=$BASE/symfc_c3_$c3
  mkdir -p "$D"
  cp -f "$S/POSCAR" "$S/SPOSCAR" "$S/dataset_disps.npy" "$S/dataset_forces.npy" "$D/"
  sed "s/^FC3_CUTOFF = .*/FC3_CUTOFF = $c3                # scan value/" "$S/symfc_fit_test.py" > "$D/fit_test.py"
  grep -m1 "^FC3_CUTOFF" "$D/fit_test.py"
  ( cd "$D" && OMP_NUM_THREADS=8 timeout 21600 ~/venvs/mace_cpu/bin/python -u fit_test.py > symfc_fit.log 2>&1; echo "rc=$?" >> symfc_fit.log )
  echo "[symfc] c3=$c3 done $(date)"
  ls -la "$D/fc2.hdf5" "$D/fc3.hdf5" 2>/dev/null
done
echo SYMFC_ALL_DONE $(date)
