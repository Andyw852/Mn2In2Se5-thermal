#!/bin/bash
#SBATCH -J mn215i-32
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=08:00:00
#SBATCH --mem=200G
#SBATCH --exclusive
#SBATCH -o log/slurm-%x-%j.out
#SBATCH -e log/slurm-%x-%j.err
#SBATCH --exclude=anode107,anode116,anode118,anode132
set -u
BASE=$HOME/work/mn215_iso
PY=$HOME/software/pybin/python3
export OMP_NUM_THREADS=64
export OPENBLAS_NUM_THREADS=64
export MKL_NUM_THREADS=64
cd "$BASE" || exit 1
mkdir -p log
echo "host=$(hostname)"
"$PY" -u run.py --mesh 32 32 32 --fc3 in/fc3.hdf5 --temps 300 --tag m32 | tee out_m32.log
echo DONE rc=$?