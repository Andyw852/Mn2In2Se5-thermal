#!/bin/bash
#SBATCH -J mn215i-wig
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --mem=100G
#SBATCH -o log/slurm-%x-%j.out
#SBATCH -e log/slurm-%x-%j.err
#SBATCH --exclude=anode107,anode116,anode118,anode132
set -u
BASE=$HOME/work/mn215_iso
PY=$HOME/software/pybin/python3
export OMP_NUM_THREADS=32
export OPENBLAS_NUM_THREADS=32
export MKL_NUM_THREADS=32
cd "$BASE" || exit 1
mkdir -p log
echo "host=$(hostname)"
"$PY" -u run_wigner.py 2>&1 | tee -a run_wigner.log
echo DONE rc=$?