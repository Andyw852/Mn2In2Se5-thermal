#!/bin/bash
#SBATCH -J mn215i-T
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
mkdir -p out_T100_700 out_Tlow
( cd out_T100_700 && "$PY" -u ../run.py --mesh 16 16 16 --fc3 in/fc3.hdf5 --temps 100 200 300 400 500 600 700 --tag T100_700 | tee ../log/T100_700.log )
( cd out_Tlow && "$PY" -u ../run.py --mesh 16 16 16 --fc3 in/fc3.hdf5 --temps 30 60 90 120 150 --tag Tlow | tee ../log/Tlow.log )
echo DONE rc=$?