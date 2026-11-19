#!/bin/bash
# Submit submit_fit_ext_gpu.sh as soon as one GPU on this node has NEED_MIB free.
# Bookkeeping only (no computation); bounded to 24 h so it cannot become an orphan.
# Use it because the 3090 cards are shared with other jobs of the same account.
set -u
X=${X:-/data/wangchaoyue852/taskflow/work/Mn2In2Se5_441/kl-dft-cpu/step5_fc/cutoff_scan_ext}
NEED_MIB=${NEED_MIB:-6000}
START=$(date +%s); MAX=$((24*3600))
cd "$X" || exit 1
log() { echo "$(date '+%F %T') $*" >> "$X/gpu_watch.log"; }
log "watcher start, need ${NEED_MIB} MiB free, max 24 h"
while :; do
    NOW=$(date +%s)
    [ $((NOW-START)) -gt $MAX ] && { log "timeout, giving up"; exit 0; }
    BEST=""; BESTF=0
    for i in 0 1 2 3 5 4 6 7; do
        F=$(nvidia-smi -i $i --query-gpu=memory.free --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
        case "$F" in ''|*[!0-9]*) continue;; esac
        if [ "$F" -gt "$BESTF" ]; then BESTF=$F; BEST=$i; fi
    done
    if [ -n "$BEST" ] && [ "$BESTF" -ge "$NEED_MIB" ]; then
        rm -f "$X/.fcfit.lock"
        JOB=$(sbatch --parsable submit_fit_ext_gpu.sh 2>>"$X/gpu_watch.log")
        log "submitted $JOB (gpu=$BEST free=${BESTF}MiB)"
        exit 0
    fi
    log "waiting: best gpu=${BEST:-none} free=${BESTF}MiB"
    sleep 180
done
