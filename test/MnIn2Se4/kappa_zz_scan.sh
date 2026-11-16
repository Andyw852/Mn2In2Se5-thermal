#!/bin/bash
# Converge kappa_zz on the q-mesh -- CPU only, no new DFT.
#
# NOTE (2026-11-16): this workstation is no longer allowed to run these
# calculations. The script is kept as the definition of the scan and as the
# reference driver; run it on a cluster (the 24^3 / 30^3 / cutoff cases of this
# study were run on HanHai -- see ../hanhai_jobs/ and section 11.6 of the
# report). The cases are skipped when their kappa-*.hdf5 already exists, so it
# is safe to re-run.
#
# Each case runs phono3py RTA with write_kappa=True in its own directory, so the
# kappa-*.hdf5 files that kappa_zz_decompose.py needs do not collide. A case is
# skipped when its kappa hdf5 already exists, so the script can be re-run (or
# resumed after an interruption) without losing finished work.
#
# Cases
#   m11, m14, m18             mesh scan, fc3 = 6.45 A, fc2 = 4x4x2 (workstation)
#   m11_b1um                  same as m11 plus a 1 um Casimir boundary,
#                             used ONLY to validate the a-posteriori boundary
#                             reconstruction in kappa_zz_decompose.py
#   c560_m14, c623_m14        fc3 cutoff 5.60 / 6.23 A at 14^3, to redo the
#                             cutoff comparison on the converged "rest" term
#                             instead of on the 11^3 total
set -u
# ── 2026-11-16, explicit user instruction: no phono3py on this workstation ──
# All 12 cores get saturated and the Windows host plus the dsh web UI freeze.
# To re-enable local runs: rm ~/.no-local-compute
if [ -f "$HOME/.no-local-compute" ]; then
  echo "[BLOCKED] local computation is disabled ($HOME/.no-local-compute exists); run phono3py on a cluster"
  exit 0
fi
# ─────────────────────────────────────────────────────────────────────────────
# cron does not inherit the conda environment, so pin the interpreter that has
# phono3py instead of relying on PATH
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
mkdir -p "$SCAN"
exec 9>"$SCAN/.lock"
flock -n 9 || { echo "[skip] another scan is already running"; exit 0; }

case_run() {          # name  mesh  extra-args...
  local name=$1 mesh=$2; shift 2
  local d="$SCAN/$name"
  local d_abs h5
  d_abs="$(pwd)/$d"     # --out must be absolute: we cd into $d below, so a relative path would become $d/$d/result.json
  h5="$d/kappa-m$mesh$mesh$mesh.hdf5"
  mkdir -p "$d"
  if [ -f "$h5" ]; then echo "[skip] $name (hdf5 exists)"; return 0; fi
  echo "[run ] $name  mesh=$mesh  $*"
  ( cd "$d" && nice -n 10 "$PY" ../../kappa_driver_442.py --mesh "$mesh" "$mesh" "$mesh" \
      --rta --fc2 4x4x2 --write-kappa --out "$d_abs/result.json" "$@" \
      >> run.log 2>&1 )
  if [ -f "$h5" ]; then echo "[done] $name"; else echo "[FAIL] $name (see $d/run.log)"; fi
}

case_run m11 11
case_run m11_b1um 11 --boundary-mfp 1.0
case_run m14 14
case_run m18 18
case_run c560_m14 14 --fc3 "$(pwd)/fc3_cutoff_variants/fc3_cut5p60.hdf5"
case_run c623_m14 14 --fc3 "$(pwd)/fc3_cutoff_variants/fc3_cut6p23.hdf5"
# 24^3 last: it is the long one (hours) and the one a reboot is most likely to kill
case_run m24 24
echo "ALL DONE"
