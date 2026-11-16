#!/bin/bash
#SBATCH -J m215fc-kzz
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=02:00:00
#SBATCH --mem=120G
#SBATCH -o log/slurm-%x-%j.out
#SBATCH -e log/slurm-%x-%j.err
#
# Decomposition-support meshes for Mn2In2Se5 kappa_zz. One job loops over the
# meshes and skips the ones that are already done, so it can be re-run.
set -u
BASE=$HOME/work/mn2in2se5_kzz_fcset
PY=$HOME/software/pybin/python3
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
cd "$BASE" || exit 1
mkdir -p out log
echo "host=$(hostname) threads=$OMP_NUM_THREADS"
for m in "16 16 4" "16 16 8" "16 16 12" "16 16 16" "24 24 24"; do
  tag=$(echo "$m" | tr -d ' ')
  d="out/m$tag"
  h5="$d/kappa-m$tag.hdf5"
  mkdir -p "$d"
  if [ -f "$h5" ]; then echo "[skip] $m (already done)"; continue; fi
  echo "[run ] mesh $m"
  ( cd "$d" && "$PY" "$BASE/run_kzz.py" --mesh $m --out "$BASE/out/m$tag/result.json" \
      >> "$BASE/out/m$tag/run.log" 2>&1 )
  [ -f "$h5" ] && echo "[done] $m" || echo "[FAIL] $m"
done
echo "ALL DONE"
