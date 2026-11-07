#!/bin/bash
# mni4_finish.sh -- after driver B finishes, (re)run every missing (method,c3) fit
# cleanly, then signal driver-A style completion for the after-scan watcher.
source /home/wangchaoyue852/miniconda3/etc/profile.d/conda.sh
conda activate wc
export CUDA_VISIBLE_DEVICES=5
R=/data/wangchaoyue852/pheasy_mninse4_189
V6=$R/pheasy_fit.sh
cd "$R" || exit 1
echo "[finish] waiting for driver B ... $(date)"
for i in $(seq 1 180); do [ -f "$R/mni4_ghsc_b_done.txt" ] && break; sleep 30; done
echo "[finish] B done; filling gaps $(date)"
SUM=$R/mni4_finish_results.txt; : > "$SUM"
runone() {
  local method=$1 c3=$2 extra=$3
  local d="$R/ghsc_${method}_c3_${c3}"
  if [ -s "$d/fc3.hdf5" ]; then echo "[ok] $d"; return; fi
  echo "[redo] $method c3=$c3 $(date)"
  rm -rf "$d"; mkdir -p "$d"
  cp "$R/POSCAR" "$R/SPOSCAR" "$R/disp_matrix.pkl" "$R/force_matrix.pkl" "$V6" "$d"/
  ( cd "$d" && timeout 28800 bash pheasy_fit.sh FIT_METHOD=$method $extra C2_CUTOFF=None C3_CUTOFF=$c3 FIT_ORDER=3 NCPU=4 > run.log 2>&1 )
  local rc=$?
  printf "%-8s c3=%-5s rc=%-3s p=%-6s rel=%s\n" "$method" "$c3" "$rc" \
      "$(grep -m1 'Free IFC terms' "$d/run.log" | awk '{print $NF}')" \
      "$(grep -m1 'Relative error' "$d/run.log" | awk '{print $NF}')" | tee -a "$SUM"
}
for c3 in 2.85 3.70 3.95 4.30 4.70 5.00 5.20 5.50 5.80 6.10 6.50 7.00; do
  runone OLS "$c3" ""
  runone LASSO "$c3" "HARM_DENSE=true"
  runone ALASSO "$c3" "HARM_DENSE=true"
done
echo "MNI4_ALL_FITS_DONE $(date)" > "$R/mni4_ghsc_done.txt"
echo "[finish] done $(date)"
