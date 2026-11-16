#!/bin/bash
#SBATCH -J mn215-lb
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=08:00:00
#SBATCH --mem=240G
#SBATCH --exclusive
#SBATCH -o log/slurm-%x-%j.out
#SBATCH -e log/slurm-%x-%j.err
#SBATCH --exclude=anode107,anode132
set -u
BASE=$HOME/work/mn215_lbte
PY=$HOME/software/pybin/python3
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
cd "$BASE" || exit 1
mkdir -p log
echo "host=$(hostname) threads=$OMP_NUM_THREADS job=$SLURM_JOB_NAME"
"$PY" "$BASE/run.py"
echo DONE