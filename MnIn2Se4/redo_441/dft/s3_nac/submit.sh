#!/bin/bash
# =====================================================================
# submit_std_3d.tpl -- vasp_std submission template (3D, step1/step2; stock VASP)
# Cluster parameters are hard-coded in this template; to change machine/queue/VASP edit only this file.
# The script only fills the JOBNAME placeholder (built from system and step; override with --jobname).
# =====================================================================
#SBATCH --partition=cpu192   # [goal R256] revert the mem-64G of R29: the large kl supercells GaP(250 atoms ~78GB)/CaO(105GB)/h-BeO(95GB) always OOM in 64GB
#SBATCH --job-name=MnIn2Se4_441-kl-dft-cpu-S3_nac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32   # 2026-10-29 used for the kl S4 force frames (128-atom small supercell); 32 cores lets several jobs share a 192-core node
#SBATCH --output=queue.out
#SBATCH --error=queue-%j.err
#SBATCH --qos=regular
#SBATCH --time=1-00:00:00     # the regular QOS wall-time limit is 1 day, so it must be set

cd $SLURM_SUBMIT_DIR

# load the same AOCC environment used at compile time (do not use Intel)
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited

# [2026-11-02] on some jzzn nodes InfiniBand (OpenFabrics mlx5) init fails -> MPI deadlock, VASP hangs at 0 CPU.
#   Force TCP to avoid openib; same origin as the xpmem workaround on hanhai25.
export OMPI_MCA_btl="self,vader,tcp"
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl_vader_single_copy_mechanism=none
export OMPI_MCA_osc=pt2pt

# standard collinear build: step1/step2 use vasp_std
export AUTOZT_CELL_CONSTRAINT=none
mpirun -np $SLURM_NTASKS /public/home/wangchao/software/vasp.6.4.3/bin/vasp_std
