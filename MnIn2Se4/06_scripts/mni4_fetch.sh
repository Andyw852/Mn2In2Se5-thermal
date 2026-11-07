#!/bin/bash
# mni4_fetch.sh -- pull fit logs (3090) and kappa/holdout lines (hanhai) into results/
T=/home/wangchao/software/taskflow-v2.0/tmp/mn2inse4
RES=$T/results
mkdir -p "$RES"
scp -q wangchao_3090:/data/wangchaoyue852/pheasy_mninse4_189/mni4_scan_results.txt wangchao_3090:/data/wangchaoyue852/pheasy_mninse4_189/mni4_scan_results_b.txt wangchao_3090:/data/wangchaoyue852/pheasy_mninse4_189/mni4_finish_results.txt "$RES/" 2>/dev/null
for d in $(ssh wangchao_3090 'ls -d /data/wangchaoyue852/pheasy_mninse4_189/ghsc_*/ 2>/dev/null | xargs -n1 basename'); do
  m=$(echo "$d" | sed -E 's/^ghsc_([A-Z]+)_c3_(.*)$/\1/')
  c=$(echo "$d" | sed -E 's/^ghsc_([A-Z]+)_c3_(.*)$/\2/')
  scp -q wangchao_3090:/data/wangchaoyue852/pheasy_mninse4_189/$d/run.log "$RES/run_${m}_c3_${c}.log" 2>/dev/null
done
SOCK="$HOME/.ssh/tf-wangchao217@114.214.255.25:22.sock"
/home/wangchao/bin/hanhai25-connect >/dev/null 2>&1 || true
ssh -o ControlPath="$SOCK" -o BatchMode=yes hanhai25 'cd ~/work/mninse4_189_kappa && for d in ghsc_*/ symfc_c3_*/; do d=${d%/}; [ -f "$d/kappa_run.log" ] && echo "$d $(grep -a "^KAPPA" $d/kappa_run.log | tail -1)"; done' > "$RES/kappa_all.txt"
# holdout RMSE logs live on 3090 (ran there)
scp -q wangchao_3090:'/data/wangchaoyue852/pheasy_mninse4_189/holdout_results/holdout_c3_*.log' "$RES/" 2>/dev/null
scp -q wangchao_3090:/data/wangchaoyue852/pheasy_mninse4_189/holdout_results/holdout.log "$RES/holdout_driver.log" 2>/dev/null
echo "run logs: $(ls $RES/run_*.log 2>/dev/null | wc -l)  kappa: $(wc -l < $RES/kappa_all.txt)  holdout: $(ls $RES/holdout_c3_*.log 2>/dev/null | wc -l)"
