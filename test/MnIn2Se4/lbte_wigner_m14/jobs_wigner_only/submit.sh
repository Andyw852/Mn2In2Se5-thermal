#!/bin/bash
#SBATCH -J mn441-wig
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=06:00:00
#SBATCH --mem=200G
#SBATCH -o log/slurm-%x-%j.out
#SBATCH -e log/slurm-%x-%j.err
#SBATCH --exclude=anode132
set -u
BASE=$HOME/work/mn441_wigner_only
PY=$HOME/software/pybin/python3
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
cd "$BASE" || exit 1
mkdir -p log
if [ -f result.json ]; then echo "[skip] already done"; exit 0; fi
echo "host=$(hostname) threads=$OMP_NUM_THREADS job=$SLURM_JOB_NAME"
"$PY" "$BASE/run.py"
echo DONE