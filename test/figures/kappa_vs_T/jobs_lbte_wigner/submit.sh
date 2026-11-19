#!/bin/bash
# HanHai: RTA + LBTE + Wigner coherence for one material, one case, one set of temperatures.
#   sbatch --export=ALL,MAT=MnIn2Se4,CASE=ordered            submit.sh
#   sbatch --export=ALL,MAT=Mn2In2Se5,CASE=vca,TEMPS="30 60 100" submit.sh
# Optional: MESH (default 14 / 16 from the driver), STEPS (default "rta lbte coherence"), X (VCA antisite fraction),
#           WRITE_KAPPA=1 (RTA step also writes the per-mode hdf5).
# LBTE keeps the collision matrix of ALL requested temperatures in memory at once:
#   MnIn2Se4 14^3 ~2.6 GB per T, Mn2In2Se5 16^3 ~9.1 GB per T (x2 for the solver) -> split TEMPS for Mn2In2Se5.
#SBATCH -J lbtewig
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --mem=230G
#SBATCH --exclusive
#SBATCH -o slurm-%x-%j.out
#SBATCH --exclude=anode107,anode116,anode118,anode132
set -u
: "${MAT:?set MAT=MnIn2Se4 or MAT=Mn2In2Se5}"
: "${CASE:?set CASE=ordered or CASE=vca}"
REPO=${REPO:-$HOME/work/Mn2In2Se5-thermal}      # clone of the repository on HanHai
PY=${PY:-$HOME/software/pybin/python3}
export OMP_NUM_THREADS=64 OPENBLAS_NUM_THREADS=64 MKL_NUM_THREADS=64
TAG=${TEMPS:+_T$(echo "$TEMPS" | tr ' ' '-')}
WD=$REPO/test/figures/kappa_vs_T/jobs_lbte_wigner/run_${MAT}_${CASE}${TAG}
mkdir -p "$WD" && cd "$WD" || exit 1
ARGS=(--material "$MAT" --cases "$CASE")
[ -n "${MESH:-}" ]  && ARGS+=(--mesh "$MESH")
[ -n "${X:-}" ]     && ARGS+=(--x "$X")
[ -n "${TEMPS:-}" ] && ARGS+=(--temps $TEMPS)
[ -n "${STEPS:-}" ] && ARGS+=(--steps $STEPS)
[ -n "${WRITE_KAPPA:-}" ] && ARGS+=(--write-kappa)
echo "host=$(hostname) MAT=$MAT CASE=$CASE ${ARGS[*]}"
"$PY" -c "import phono3py; print('phono3py', phono3py.__version__)"
/usr/bin/time -v "$PY" -u ../run_lbte_wigner.py "${ARGS[@]}" 2>&1 | tee run.log
echo "DONE rc=${PIPESTATUS[0]}"
