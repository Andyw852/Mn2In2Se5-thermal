#!/bin/bash
#SBATCH -J mn215-vca
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=04:00:00
#SBATCH --mem=200G
#SBATCH --exclusive
#SBATCH -o log/slurm-%x-%j.out
#SBATCH -e log/slurm-%x-%j.err
#SBATCH --exclude=anode107,anode116,anode118,anode132
set -u
BASE=$HOME/work/mn215_vca
PY=$HOME/software/pybin/python3
export OMP_NUM_THREADS=64
export OPENBLAS_NUM_THREADS=64
export MKL_NUM_THREADS=64
cd "$BASE" || exit 1
mkdir -p log
echo "host=$(hostname)"
"$PY" -u run.py 2>&1 | tee -a run.log
echo DONE rc=$?