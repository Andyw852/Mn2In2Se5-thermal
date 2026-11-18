#!/bin/bash
#SBATCH -J cut215
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --mem=100G
#SBATCH -o log/cut-%j.out
#SBATCH -e log/cut-%j.err
#SBATCH --exclude=anode107,anode116,anode118,anode132
set -u
BASE=$HOME/work/coh_chk
PY=$HOME/software/pybin/python3
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 MKL_NUM_THREADS=32
cd "$BASE" || exit 1
mkdir -p log
"$PY" -u run_cutoff.py 16 0.3 2>&1 | grep -v Warning
echo DONE