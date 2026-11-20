#!/bin/bash
# One CV fit on A800.  sbatch --export=ALL,RUN=c0743_f0 submit_cv_fit_a800.sh
#SBATCH -J cvfit
#SBATCH -p a800
#SBATCH --gres=gpu:a800:1
#SBATCH -N1
#SBATCH -n1
#SBATCH -c 8
#SBATCH --time=08:00:00
#SBATCH -o log/cv-%j.out
#SBATCH -e log/cv-%j.err
set -u
: "${RUN:?set RUN=<run name>}"
B=$HOME/work/mn_cv
cd "$B/runs/$RUN" || exit 1
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mace
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
echo "[cv] host=$(hostname) run=$RUN gpu=${CUDA_VISIBLE_DEVICES:-?}"
python -c "import torch;print('[gpu] cuda=%s n=%d'%(torch.cuda.is_available(),torch.cuda.device_count()))"
T0=$(date +%s)
python -u "$B/fc_fit_driver.py" prep fit_config.json > prep.log 2>&1 || { echo '[FAIL] prep'; tail -12 prep.log; exit 1; }
python -u "$B/fc_fit_driver.py" fit  fit_config.json > fit.log  2>&1 || { echo '[FAIL] fit'; tail -20 fit.log; exit 1; }
python -u "$B/fc_fit_driver.py" post fit_config.json > post.log 2>&1 || echo '[warn] post non-zero'
T1=$(date +%s)
echo "[cv] $RUN done in $((T1-T0)) s"
mkdir -p "$B/fits/$RUN" && cp fc2.hdf5 fc3.hdf5 "$B/fits/$RUN/" && echo '[cv] copied fits'
echo DONE