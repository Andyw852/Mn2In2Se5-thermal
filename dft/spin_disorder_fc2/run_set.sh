#!/bin/bash
# One job per set directory written by prepare.py (FM_ctl, AFM, PM01 ... PM08); its VASP runs one after another.
#   sbatch run_set.sh <work>/PM01 ref      # undisplaced cell only (check its moments before the rest)
#   sbatch run_set.sh <work>/PM01 disp     # the displaced cells, each started from ../ref/WAVECAR
# Runs that already have a finished OUTCAR ("General timing") are skipped, so a job can simply be resubmitted.
# Partition, modules, MPI settings and VASP binary are those of dft/fc3_cutoff_extend/submit_frames.sh
# (jzzn, cpu192, 48 cores), which ran the same 144-atom cell in about 2-3 h per frame.
#SBATCH --partition=cpu192
#SBATCH --job-name=spinfc2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --output=spinfc2-%j.out
#SBATCH --qos=regular
#SBATCH --time=2-00:00:00
set -u
D=$(cd "${1:?set directory}" && pwd); STAGE=${2:?ref or disp}
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited
export OMPI_MCA_btl="self,vader,tcp"
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl_vader_single_copy_mechanism=none
export OMPI_MCA_osc=pt2pt
VASP="mpirun -np $SLURM_NTASKS /public/home/wangchao/software/vasp.6.4.3/bin/vasp_std"

run() {   # run VASP in directory $1 unless it already finished
  if grep -q "General timing" "$1/OUTCAR" 2>/dev/null; then echo "done   $1"; return 0; fi
  [ -f "$1/POTCAR" ] || { echo "no POTCAR in $1"; exit 1; }
  echo "start  $1  $(date)"
  (cd "$1" && $VASP > vasp.log 2>&1)
  grep -q "General timing" "$1/OUTCAR" || { echo "FAILED $1 (see $1/vasp.log)"; exit 1; }
}

case $STAGE in
  ref)  run "$D/ref" ;;
  disp)
    [ -f "$D/ref/WAVECAR" ] || { echo "run the ref stage first"; exit 1; }
    for d in "$D"/disp-*; do
      grep -q "General timing" "$d/OUTCAR" 2>/dev/null || cp "$D/ref/WAVECAR" "$d/"
      run "$d"
      rm -f "$d/WAVECAR"
    done ;;
  *)    echo "stage must be ref or disp"; exit 1 ;;
esac
echo "finished $D $STAGE $(date)"
