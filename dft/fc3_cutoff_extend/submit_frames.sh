#!/bin/bash
# Submit VASP force runs for frames written by make_frames.py, on jzzn, one job per frame.
#   bash submit_frames.sh new_frames 55 84      # Mn2In2Se5 frames disp-00055 .. disp-00084
# Skips frames that already have a complete vasprun.xml or a queue.out.
# Needs POTCAR in the frame directory root (the same set as the existing 4x4x1 S4 run).
set -u
ROOT=$(cd "${1:?frame directory}" && pwd); first=${2:?first frame}; last=${3:?last frame}
[ -f "$ROOT/POTCAR" ] || { echo "copy the 4x4x1 S4 POTCAR to $ROOT/POTCAR first"; exit 1; }
for k in $(seq "$first" "$last"); do
  d=$ROOT/$(printf "disp-%05d" "$k")
  [ -d "$d" ] || { echo "skip $d (not generated)"; continue; }
  if grep -q "</modeling>" "$d/vasprun.xml" 2>/dev/null; then echo "done $d"; continue; fi
  [ -f "$d/queue.out" ] && { echo "queued/started $d"; continue; }
  cp "$ROOT/POTCAR" "$d/"
  cat > "$d/job.sh" <<JOB
#!/bin/bash
#SBATCH --partition=cpu192
#SBATCH --job-name=fc3x-$(printf "%05d" "$k")
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --output=queue.out
#SBATCH --error=queue-%j.err
#SBATCH --qos=regular
#SBATCH --time=1-00:00:00
cd \$SLURM_SUBMIT_DIR
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited
export OMPI_MCA_btl="self,vader,tcp"
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl_vader_single_copy_mechanism=none
export OMPI_MCA_osc=pt2pt
mpirun -np \$SLURM_NTASKS /public/home/wangchao/software/vasp.6.4.3/bin/vasp_std
JOB
  (cd "$d" && sbatch job.sh)
done
