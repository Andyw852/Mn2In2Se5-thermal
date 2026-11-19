#!/bin/bash
# FC3 cutoff-extension refit on the 3090 cluster: pheasy-gpu + ALASSO.
#
# This is the single kept copy. Earlier drafts of this script were wrong and are
# deleted -- do not reintroduce them:
#   * a CPU version (the task requires the GPU fit),
#   * PHEASY_GPU_MODE=prefer  -> this pheasy-gpu build only accepts
#     auto | cpu | required, and a wrong value raises ValueError deep inside
#     optimizer.fit (see the note in README.md).
#
# GPU handling on this node (5x RTX 3090, one physically dead at 84:00.0) is
# what makes this script non-trivial:
#   * SLURM's --gres=gpu:1 does not protect the card from other jobs of the same
#     account that grab "the first visible GPU" without a GRES request, so a
#     freshly started job can find <1 GiB free and the dense solve in
#     pheasy_gpu/core/optimizer.py aborts with "GPU dense solve exceeds the
#     configured VRAM budget";
#   * therefore probe every card and pin CUDA_VISIBLE_DEVICES to the one with the
#     most free memory, and require a hard minimum (NEED_MIB, default 6000).
#
# Usage:  sbatch --export=ALL submit_fit_ext_gpu.sh
#         (or let gpu_watch.sh submit it when a card frees up)
#SBATCH -J fc3ext-gpu
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --time=08:00:00
#SBATCH -o log/fc3ext-%j.out
#SBATCH -e log/fc3ext-%j.err
set -u
X=${X:-/data/wangchaoyue852/taskflow/work/Mn2In2Se5_441/kl-dft-cpu/step5_fc/cutoff_scan_ext}
NEED_MIB=${NEED_MIB:-6000}
mkdir -p "$X/log" && cd "$X" || exit 1

source /home/wangchaoyue852/miniconda3/etc/profile.d/conda.sh
conda activate mace-gpu
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8

export PHEASY_USE_GPU=1
export PHEASY_GPU_MODE=required
export PHEASY_GPU_SM=1
export PHEASY_GPU_SM_NGPU=1
export PHEASY_GPU_SM_DEVICES=0
export PHEASY_GPU_LASSO_RESIDENT=1
# Do NOT set PHEASY_GPU_RESIDENT_FALLBACK=1 here: that opts into a substitute
# backend when the VRAM pre-flight fails, which is exactly the silent downgrade
# the required-GPU rule is meant to prevent. With the flag unset the run raises
# instead, and gpu_watch.sh retries when a card is actually free.
export PHEASY_GPU_MEM_FRACTION=0.9
export PHEASY_TWOLEVEL_CACHE_T=0
export PHEASY_FIT_METHOD=ALASSO
export PHEASY_RASR=BHH

echo "[env] host=$(hostname) SLURM_CVD=${CUDA_VISIBLE_DEVICES:-unset}"
BEST=""; BESTF=0
for i in 0 1 2 3 5 4 6 7; do
    F=$(nvidia-smi -i $i --query-gpu=memory.free --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
    case "$F" in ''|*[!0-9]*) continue;; esac
    echo "[probe] gpu $i free=${F}MiB"
    if [ "$F" -gt "$BESTF" ]; then BESTF=$F; BEST=$i; fi
done
if [ -z "$BEST" ] || [ "$BESTF" -lt "$NEED_MIB" ]; then
    echo "[gpu][FATAL] no card with >=${NEED_MIB}MiB free (best=${BEST:-none} ${BESTF}MiB); refusing to run" >&2
    exit 3
fi
export CUDA_VISIBLE_DEVICES=$BEST
echo "[gpu] pinned to physical GPU $BEST (${BESTF}MiB free)"
python - <<'PYG'
import torch
print('[gpu] cuda=%s count=%d' % (torch.cuda.is_available(), torch.cuda.device_count()))
if torch.cuda.is_available():
    f, t = torch.cuda.mem_get_info()
    print('[gpu] free %.1f GiB / total %.1f GiB' % (f / 2**30, t / 2**30))
PYG

echo '=== prep'; python -u fc_fit_driver.py prep fit_config.json 2>&1 | tail -8
echo '=== fit';  python -u fc_fit_driver.py fit  fit_config.json 2>&1 | tail -50
echo '=== post'; python -u fc_fit_driver.py post fit_config.json 2>&1 | tail -10
echo DONE
