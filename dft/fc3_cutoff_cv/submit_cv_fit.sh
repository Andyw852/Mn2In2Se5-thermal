#!/bin/bash
# One hold-out / learning-curve fit of the FC3-cutoff CV study.
#   sbatch --export=ALL,RUN=c0743_f0 submit_cv_fit.sh
# Runs in $BASE/runs/$RUN, driver and modules in $BASE, subset data in $BASE/subsets.
# GPU handling is the same as dft/fc3_cutoff_extend/submit_fit_ext_gpu.sh:
# PHEASY_GPU_MODE=required, pin the freest card, no resident fallback.
#SBATCH -J cvfit
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --time=08:00:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
set -u
: "${RUN:?set RUN=<run name>}"
BASE=/data/wangchaoyue852/taskflow/work/Mn2In2Se5_441/kl-dft-cpu/step5_fc/cutoff_cv
NEED_MIB=${NEED_MIB:-6000}
cd "$BASE/runs/$RUN" || exit 1
source /home/wangchaoyue852/miniconda3/etc/profile.d/conda.sh
conda activate mace-gpu
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
export PHEASY_USE_GPU=1
export PHEASY_GPU_MODE=required
export PHEASY_GPU_SM=1
export PHEASY_GPU_SM_NGPU=1
export PHEASY_GPU_SM_DEVICES=0
export PHEASY_GPU_LASSO_RESIDENT=1
export PHEASY_GPU_MEM_FRACTION=0.9
export PHEASY_TWOLEVEL_CACHE_T=0
export PHEASY_FIT_METHOD=ALASSO
export PHEASY_RASR=BHH
BEST=""; BESTF=0
for i in 0 1 2 3 5 4 6 7; do
  F=$(nvidia-smi -i $i --query-gpu=memory.free --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
  case "$F" in ''|*[!0-9]*) continue;; esac
  if [ "$F" -gt "$BESTF" ]; then BESTF=$F; BEST=$i; fi
done
if [ -z "$BEST" ] || [ "$BESTF" -lt "$NEED_MIB" ]; then echo "[gpu][FATAL] no card with >=${NEED_MIB}MiB free" >&2; exit 3; fi
export CUDA_VISIBLE_DEVICES=$BEST
echo "[cv] host=$(hostname) run=$RUN gpu=$BEST free=${BESTF}MiB"
T0=$(date +%s)
python -u "$BASE/fc_fit_driver.py" prep fit_config.json > prep.log 2>&1 || { echo '[FAIL] prep'; tail -5 prep.log; exit 1; }
python -u "$BASE/fc_fit_driver.py" fit  fit_config.json > fit.log  2>&1 || { echo '[FAIL] fit';  tail -20 fit.log; exit 1; }
python -u "$BASE/fc_fit_driver.py" post fit_config.json > post.log 2>&1 || echo '[warn] post returned non-zero (fc2/fc3 still written)'
T1=$(date +%s)
echo "[cv] $RUN done in $((T1-T0)) s"
ls -la fc2.hdf5 fc3.hdf5 2>/dev/null
cp fc2.hdf5 fc3.hdf5 "$BASE/fits/$RUN/" 2>/dev/null || { mkdir -p "$BASE/fits/$RUN"; cp fc2.hdf5 fc3.hdf5 "$BASE/fits/$RUN"; }
echo DONE