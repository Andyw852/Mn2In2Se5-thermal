#!/bin/bash
#SBATCH -J m215kzz
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=08:00:00
#SBATCH --mem=160G
#SBATCH -o log/slurm-%x-%j.out
#SBATCH -e log/slurm-%x-%j.err
#
# One isotropic mesh of Mn2In2Se5 kappa_zz per job.
#   sbatch submit_one_mesh.sh 28
#   sbatch submit_one_mesh.sh 32
#
# Uses $HOME/work/mn2in2se5_kzz_fcset/, i.e. the force constants of
# Mn2In2Se5/fc/ (the set the 16^3 and 24^3 runs used and the set whose C44 was
# verified in test/Mn2In2Se5/fc_sound_velocity/). primitive_matrix="P",
# isotopes off -- the same physics as the 16^3/24^3 runs, so the rest-term
# series stays on one footing.
set -u
N=$1
BASE=$HOME/work/mn2in2se5_kzz_fcset
PY=$HOME/software/pybin/python3
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
cd "$BASE" || exit 1
d="$BASE/out/m$N$N$N"
h5="$d/kappa-m$N$N$N.hdf5"
mkdir -p "$d" log
if [ -f "$h5" ]; then echo "[skip] mesh $N already done"; exit 0; fi
echo "host=$(hostname) mesh=$N threads=$OMP_NUM_THREADS isotope=False fc=Mn2In2Se5/fc"
# phono3py writes kappa-*.hdf5 into the CWD, so run from the case directory
cd "$d" || exit 1
"$PY" "$BASE/run_kzz.py" --mesh "$N" "$N" "$N" --out "$d/result.json" \
    >> "$d/run.log" 2>&1
if [ -f "$h5" ]; then echo "[done] mesh $N"; ls -la "$h5"; else echo "[FAIL] mesh $N"; tail -5 "$d/run.log"; fi
