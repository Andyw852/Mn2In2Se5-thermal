#!/bin/bash
# mni4_holdout_driver.sh -- 5-fold (whole-configuration) held-out force RMSE for
# MnIn2Se4 (45 configs x 567 rows) for OLS / LASSO / ALASSO at each C3 cutoff.
# PHEASY_HARM_DENSE=1 so LASSO/ALASSO match the HARM_DENSE production fits.
# Runs up to 3 cutoffs concurrently (small SM, plenty of cores).
source /home/wangchaoyue852/miniconda3/etc/profile.d/conda.sh
conda activate wc
export CUDA_VISIBLE_DEVICES=5
export PHEASY_HARM_DENSE=1
R=/data/wangchaoyue852/pheasy_mninse4_189
G=/home/wangchaoyue852/software/pheasy-gpu-gh
OUT=$R/holdout_results
mkdir -p "$OUT"; : > "$OUT/holdout.log"
run_one() {
  local c3=$1
  local d="$R/ghsc_LASSO_c3_${c3}"
  [ -f "$d/sm_prime.npz" ] || { echo "[skip-nodata] c3=$c3" >> "$OUT/holdout.log"; return; }
  if [ -s "$OUT/holdout_c3_${c3}.log" ]; then return; fi
  echo "[run ] c3=$c3 $(date)" >> "$OUT/holdout.log"
  ( cd "$d" && OMP_NUM_THREADS=4 timeout 28800 python3 "$G/holdout_eval.py" . \
      --methods "OLS LASSO ALASSO" --n-configs 45 --rows-per-config 567 --n-splits 5 \
      --paired-refs "OLS" > "$OUT/holdout_c3_${c3}.log" 2>&1 )
  echo "[done] c3=$c3 rc=$? $(date)" >> "$OUT/holdout.log"
}
N=3
for c3 in 2.85 3.70 3.95 4.30 4.70 5.00 5.20 5.50 5.80 6.10 6.50 7.00; do
  run_one "$c3" &
  while [ "$(jobs -r | wc -l)" -ge "$N" ]; do wait -n; done
done
wait
echo "MNI4_HOLDOUT_DONE $(date)" > "$R/mni4_holdout_done.txt"
echo "[holdout driver] done $(date)"
