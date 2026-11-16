#!/bin/bash
#SBATCH -J mn215i-cut
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
for c in 5p68 6p12; do
  d=out_cut$c; mkdir -p $d
  ( cd $d && "$PY" -u ../run.py --mesh 16 16 16 --fc3 in/fc3_cut$c.hdf5 --temps 300 --tag cut$c | tee ../log/cut$c.log )
done
d=out_cut6p51; mkdir -p $d
( cd $d && "$PY" -u ../run.py --mesh 16 16 16 --fc3 in/fc3_cut6p51_scan.hdf5 --temps 300 --tag cut6p51 | tee ../log/cut6p51.log )
echo DONE rc=$?