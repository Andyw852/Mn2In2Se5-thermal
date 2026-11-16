#!/bin/bash
# 30x30x30 RTA -- the overnight refinement. Runs only after the 24^3 case of
# kappa_zz_scan.sh has finished, so the two never compete for the cores.
# Like kappa_zz_scan.sh it skips work that is already done, so it can be
# re-invoked (e.g. from cron after a reboot) without repeating a run.
set -u
# ── 2026-11-16, explicit user instruction: no phono3py on this workstation ──
# All 12 cores get saturated and the Windows host plus the dsh web UI freeze.
# To re-enable local runs: rm ~/.no-local-compute
if [ -f "$HOME/.no-local-compute" ]; then
  echo "[BLOCKED] local computation is disabled ($HOME/.no-local-compute exists); run phono3py on a cluster"
  exit 0
fi
# ─────────────────────────────────────────────────────────────────────────────
PY=/home/wangchao/miniconda3/envs/atomate2_p_a/bin/python3
command -v "$PY" >/dev/null || { echo "[FAIL] $PY not found"; exit 1; }
# ── 2026-11-16, fixes from the dsh diagnosis ────────────────────────────────
# 1) This workstation has only 12 cores; when phono3py/OpenMP/BLAS take all of
#    them, Windows and the dsh web UI become unresponsive (the 12:33 freeze
#    happened during the all-core 18^3 run). Default 6 threads, override with
#    KAPPA_THREADS=N.
# 2) Run the calculation with nice 10 so the foreground and Windows win.
# 3) Cap virtual memory (KAPPA_MEM_MAX_GB, default 32): beyond it Python raises
#    MemoryError instead of pushing the whole machine into swap. 0 disables.
THREADS=${KAPPA_THREADS:-6}
export OMP_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS \
       OPENBLAS_NUM_THREADS=$THREADS NUMEXPR_NUM_THREADS=$THREADS
MEM_MAX_GB=${KAPPA_MEM_MAX_GB:-32}
case "$MEM_MAX_GB" in
  ''|*[!0-9]*) : ;;                      # not a number: do not set
  0) : ;;                                # 0: disabled
  *) ulimit -v $((MEM_MAX_GB * 1024 * 1024)) 2>/dev/null || true ;;
esac
# ─────────────────────────────────────────────────────────────────────────────
cd "$(dirname "$(readlink -f "$0")")" || exit 1
SCAN=kappa_zz_scan
[ -f "$SCAN/m24/kappa-m242424.hdf5" ] || { echo "[wait] 24^3 not finished yet"; exit 0; }
d="$SCAN/m30"
h5="$d/kappa-m303030.hdf5"
mkdir -p "$d"
exec 9>"$SCAN/.lock30"
flock -n 9 || { echo "[skip] 30^3 already running"; exit 0; }
[ -f "$h5" ] && { echo "[skip] 30^3 done"; exit 0; }
echo "[run ] 30^3"
( cd "$d" && nice -n 10 "$PY" ../../kappa_driver_442.py --mesh 30 30 30 --rta --fc2 4x4x2 \
    --write-kappa --out "$(pwd)/result.json" >> run.log 2>&1 )
[ -f "$h5" ] && echo "[done] 30^3" || echo "[FAIL] 30^3 (see $d/run.log)"
