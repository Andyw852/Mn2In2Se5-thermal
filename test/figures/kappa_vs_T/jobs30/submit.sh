#!/bin/bash
# HanHai: one job per material, ordered + VCA at 30^3, 12 temperatures.
#   sbatch --export=ALL,MAT=MnIn2Se4  submit.sh
#   sbatch --export=ALL,MAT=Mn2In2Se5 submit.sh
#SBATCH -J kvsT30
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=12:00:00
#SBATCH --mem=200G
#SBATCH --exclusive
#SBATCH -o slurm-%x-%j.out
#SBATCH --exclude=anode107,anode116,anode118,anode132
set -u
: "${MAT:?set MAT=MnIn2Se4 or MAT=Mn2In2Se5}"
REPO=${REPO:-$HOME/work/Mn2In2Se5-thermal}      # clone of the repository on HanHai
PY=${PY:-$HOME/software/pybin/python3}
export OMP_NUM_THREADS=64 OPENBLAS_NUM_THREADS=64 MKL_NUM_THREADS=64
WD=$REPO/test/figures/kappa_vs_T/jobs30/run_$MAT
mkdir -p "$WD" && cd "$WD" || exit 1
echo "host=$(hostname) MAT=$MAT"
"$PY" -u ../run_kappa_vs_T.py --material "$MAT" --mesh 30 --cases ordered vca 2>&1 | tee run.log
echo "DONE rc=${PIPESTATUS[0]}"
