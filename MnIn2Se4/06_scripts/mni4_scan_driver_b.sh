#!/bin/bash
# mni4_scan_driver_b.sh -- parallel tail group (large cutoffs) of the MnIn2Se4 scan.
source /home/wangchaoyue852/miniconda3/etc/profile.d/conda.sh
conda activate wc
export CUDA_VISIBLE_DEVICES=5
R=/data/wangchaoyue852/pheasy_mninse4_189
V6=$R/pheasy_fit.sh
cd "$R" || exit 1
SUM=$R/mni4_scan_results_b.txt
: > "$SUM"
runone() {
  local method=$1 c3=$2 extra=$3
  local d="$R/ghsc_${method}_c3_${c3}"
  if [ -s "$d/fc3.hdf5" ]; then echo "[skip] $d"; return; fi
  rm -rf "$d"; mkdir -p "$d"
  cp "$R/POSCAR" "$R/SPOSCAR" "$R/disp_matrix.pkl" "$R/force_matrix.pkl" "$d"/
  cp "$V6" "$d/pheasy_fit.sh"
  echo "[run ] $method c3=$c3 $extra $(date)"
  ( cd "$d" && timeout 28800 bash pheasy_fit.sh FIT_METHOD=$method $extra C2_CUTOFF=None C3_CUTOFF=$c3 FIT_ORDER=3 NCPU=4 > run.log 2>&1 )
  local rc=$?
  local fc er cv nz f3 aopt
  fc=$(grep -m1 "Free IFC terms" "$d/run.log" | awk '{print $NF}')
  er=$(grep -m1 "Relative error" "$d/run.log" | awk '{print $NF}')
  cv=$(grep -m1 "RMSE_CV" "$d/run.log" | awk '{print $NF}')
  nz=$(grep -m1 "Non-zero IFC terms" "$d/run.log" | awk '{print $NF}')
  f3=$(grep -m1 "fc3       max" "$d/run.log" | awk '{print $NF}')
  aopt=$(grep -m1 "alpha_opt" "$d/run.log" | awk '{print $NF}')
  printf "%-8s c3=%-5s rc=%-3s p=%-6s rmse_cv=%-12s rel=%-12s nz=%-6s alpha=%-14s fc3max=%s\n" \
      "$method" "$c3" "$rc" "$fc" "$cv" "$er" "$nz" "$aopt" "$f3" | tee -a "$SUM"
  echo "[done] $method c3=$c3 rc=$rc $(date)"
}
for c3 in 5.50 5.80 6.10 6.50 7.00; do
  runone OLS "$c3" ""
  runone LASSO "$c3" "HARM_DENSE=true"
  runone ALASSO "$c3" "HARM_DENSE=true"
done
echo "MNI4_SCAN_B_DONE $(date)" > "$R/mni4_ghsc_b_done.txt"
echo "[driver B] done $(date)"
