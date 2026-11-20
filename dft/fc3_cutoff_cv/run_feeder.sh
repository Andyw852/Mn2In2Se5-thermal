#!/bin/bash
# Feed the 18 FC3-cutoff CV fits, keeping at most MAXINFLIGHT of them queued or
# running. Bookkeeping only; bounded to 48 h so it cannot become an orphan.
#   RUNQ=<queue file> SENT=<sent file> SUB=<submit script> sbatch_cmd=<...>
# On A800 the partition is congested by other users, so this only tops up as
# slots free; the caller may run the same script on several clusters.
set -u
B=${B:-$HOME/work/mn_cv}
Q=${RUNQ:-$B/run_queue.txt}
S=${SENT:-$B/runs_sent.txt}
L=${LOG:-$B/run_feeder.log}
SUB=${SUB:-$B/submit_cv_fit_a800.sh}
MAXINFLIGHT=${MAXINFLIGHT:-2}
JOBNAME=${JOBNAME:-cvfit}
START=$(date +%s); MAX=$((48*3600))
log() { echo "$(date '+%F %T') $*" >> "$L"; }
inflight() { squeue -u "$USER" -n "$JOBNAME" -h -o %i 2>/dev/null | wc -l; }
log "feeder start on $(hostname): maxinflight=$MAXINFLIGHT queue=$Q"
while read -r RUN; do
  [ -z "${RUN:-}" ] && continue
  case "$RUN" in \#*) continue;; esac
  grep -qxF "$RUN" "$S" 2>/dev/null && continue
  while :; do
    NOW=$(date +%s)
    [ $((NOW-START)) -gt $MAX ] && { log 'timeout, giving up'; exit 0; }
    [ "$(inflight)" -lt "$MAXINFLIGHT" ] && break
    sleep 180
  done
  JOB=$(sbatch --parsable --export=ALL,RUN="$RUN" "$SUB" 2>>"$L")
  echo "$RUN" >> "$S"
  log "submitted $JOB  RUN=$RUN"
done < "$Q"
log 'queue exhausted'