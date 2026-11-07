#!/bin/bash
# mni4_kappa_sync.sh -- pull fitted POSCAR/fc2/fc3 for every scan dir on 3090,
# stage them on hanhai, and run phono3py RTA (11^3, 300 K, isotope).
set -u
T=/home/wangchao/software/taskflow-v2.0/tmp/mn2inse4
ST=$T/kappa_stage
mkdir -p "$ST"
SOCK="$HOME/.ssh/tf-wangchao217@114.214.255.25:22.sock"
/home/wangchao/bin/hanhai25-connect >/dev/null 2>&1 || true
RS=/data/wangchaoyue852/pheasy_mninse4_189
DIRS=$(ssh wangchao_3090 "ls -d $RS/ghsc_* 2>/dev/null | xargs -n1 basename" 2>/dev/null)
echo "[sync] dirs: $DIRS"
for n in $DIRS; do
  if ssh -o ControlPath="$SOCK" -o BatchMode=yes hanhai25 "[ -f ~/work/mninse4_189_kappa/$n/kappa-m111111.hdf5 ] && echo yes" 2>/dev/null | grep -q yes; then
    echo "[skip-kappa] $n"; continue
  fi
  ssh wangchao_3090 "[ -s $RS/$n/fc3.hdf5 ]" 2>/dev/null || { echo "[no-fc3] $n"; continue; }
  mkdir -p "$ST/$n"
  scp -q wangchao_3090:"$RS/$n/POSCAR" wangchao_3090:"$RS/$n/fc2.hdf5" wangchao_3090:"$RS/$n/fc3.hdf5" "$ST/$n/" || { echo "[pull-fail] $n"; continue; }
  ssh -o ControlPath="$SOCK" -o BatchMode=yes hanhai25 "mkdir -p ~/work/mninse4_189_kappa/$n"
  scp -q -o ControlPath="$SOCK" -o BatchMode=yes "$ST/$n/POSCAR" "$ST/$n/fc2.hdf5" "$ST/$n/fc3.hdf5" "hanhai25:~/work/mninse4_189_kappa/$n/"
  echo "[push] $n"
done
ssh -o ControlPath="$SOCK" -o BatchMode=yes hanhai25 'cd ~/work/mninse4_189_kappa && for d in ghsc_*/; do d=${d%/}; [ -f "$d/fc3.hdf5" ] || continue; [ -f "$d/kappa-m111111.hdf5" ] && { echo "[skip] $d"; continue; }; echo "[kappa] $d $(date +%H:%M:%S)"; ( cd "$d" && OMP_NUM_THREADS=16 timeout 7200 ~/venvs/mace_cpu/bin/python ../kappa_full_solver.py 11 > kappa_run.log 2>&1 ); echo "  rc=$? $(grep -a KAPPA $d/kappa_run.log | tail -1)"; done; echo KAPPA_ALL_DONE $(date)'
echo "[kappa sync] done $(date)"
