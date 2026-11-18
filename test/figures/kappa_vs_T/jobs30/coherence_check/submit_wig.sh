#!/bin/bash
#SBATCH -J wig-chk
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --mem=100G
#SBATCH -o log/wig-%j.out
#SBATCH -e log/wig-%j.err
#SBATCH --exclude=anode107,anode116,anode118,anode132
set -u
BASE=$HOME/work/coh_chk
PY=$HOME/software/pybin/python3
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 MKL_NUM_THREADS=32
cd "$BASE" || exit 1
mkdir -p log
W=test/figures/kappa_vs_T/jobs30/wig.py
"$PY" -c 'import phono3py;print("phono3py", phono3py.__version__)'
for V in SMM19 NJC23 IBDB19; do "$PY" -u $W 11 $V 2>&1 | tail -1; done
echo DONE