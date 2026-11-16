#!/bin/bash
#SBATCH -J kzz
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=04:00:00
#SBATCH --mem=120G
#SBATCH -o log/slurm-%x-%j.out
#SBATCH -e log/slurm-%x-%j.err
#
# usage: sbatch submit_case.sh <case> <mesh-N> <fc3-file>
#   e.g. sbatch submit_case.sh c560_m14 14 fc3_cut5p60.hdf5
set -eu
CASE=$1
MESH=$2
FC3=$3
BASE=$HOME/work/mn441_kappa_zz
mkdir -p "$BASE/out/$CASE" "$BASE/log"
cd "$BASE/out/$CASE"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
echo "host=$(hostname) case=$CASE mesh=$MESH fc3=$FC3 threads=$OMP_NUM_THREADS"
$HOME/software/pybin/python3 "$BASE/run_kzz.py" \
    --mesh "$MESH" "$MESH" "$MESH" --fc3 "$FC3" --out "$BASE/out/$CASE/result.json"
ls -la "$BASE/out/$CASE"
