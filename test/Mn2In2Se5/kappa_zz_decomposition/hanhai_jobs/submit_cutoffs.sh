#!/bin/bash
#SBATCH -J m215cut
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=01:00:00
#SBATCH --mem=120G
#SBATCH -o log/slurm-%x-%j.out
#SBATCH -e log/slurm-%x-%j.err
#
# Redo the Mn2In2Se5 kappa_zz FC3-cutoff comparison on an isotropic 16^3 mesh,
# so the decomposition can be applied to the "rest of the zone" term instead of
# to the contaminated totals (same idea as MnIn2Se4 section 6.6).
#
# Three cases, each in its own output directory so nothing is overwritten:
#   c568 -> in/fc3_cut5p68.hdf5
#   c612 -> in/fc3_cut6p12.hdf5
#   c651 -> in/fc3_cut6p51_scan.hdf5   (the scan copy of 6.51 A; the production
#                                        copy is a different file, md5 ac7bf6b2...,
#                                        and is compared against this run)
set -u
BASE=$HOME/work/mn2in2se5_cutoff
PY=$HOME/software/pybin/python3
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
cd "$BASE" || exit 1
mkdir -p log out
echo "host=$(hostname) threads=$OMP_NUM_THREADS"
for spec in "c568 5p68" "c612 6p12" "c651 6p51_scan"; do
  set -- $spec
  tag=$1; f=$2
  d="$BASE/out/${tag}_m161616"
  h5="$d/kappa-m161616.hdf5"
  mkdir -p "$d"
  if [ -f "$h5" ]; then echo "[skip] $tag"; continue; fi
  echo "[run ] $tag (fc3_cut$f)"
  ( cd "$d" && "$PY" "$BASE/run_kzz_cut.py" --mesh 16 16 16 --fc3 "fc3_cut$f.hdf5" \
      --out "$d/result.json" >> "$d/run.log" 2>&1 )
  if [ -f "$h5" ]; then echo "[done] $tag"; else echo "[FAIL] $tag"; tail -3 "$d/run.log"; fi
done
echo "ALL DONE"