#!/bin/bash
# mni4_after_scan.sh -- wait for the 3090 cutoff scan, then launch the held-out
# RMSE driver on 3090 and the kappa sync+run on hanhai.
T=/home/wangchao/software/taskflow-v2.0/tmp/mn2inse4
R=/data/wangchaoyue852/pheasy_mninse4_189
for i in $(seq 1 300); do
  if ssh wangchao_3090 "test -f $R/mni4_ghsc_done.txt" 2>/dev/null; then
    echo "[after] scan DONE detected $(date)"; break
  fi
  sleep 60
done
echo "[after] launching holdout driver on 3090 $(date)"
scp -q "$T/mni4_holdout_driver.sh" wangchao_3090:$R/
ssh wangchao_3090 "cd $R && nohup bash mni4_holdout_driver.sh > mni4_holdout_driver.out 2>&1 & echo HOLDOUT_PID=\$!"
echo "[after] running kappa sync $(date)"
bash "$T/mni4_kappa_sync.sh"
echo "[after] ALL DONE $(date)"
