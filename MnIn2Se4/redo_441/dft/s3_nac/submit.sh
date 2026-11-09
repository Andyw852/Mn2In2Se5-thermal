#!/bin/bash
# =====================================================================
# submit_std_3d.tpl —— vasp_std 提交模板（3D，step1/step2；官方 VASP）
# 集群参数在此模板中直接写死；换机器/换队列/换 VASP 时只改本文件。
# 脚本只填充 JOBNAME 占位符（按体系与步骤自动生成，可用 --jobname 覆盖）。
# =====================================================================
#SBATCH --partition=cpu192   # [goal R256] 撤销 R29 的 mem-64G：kl 大超胞 GaP(250原子~78GB)/CaO(105GB)/h-BeO(95GB) 在 64GB 必 OOM
#SBATCH --job-name=MnIn2Se4_441-kl-dft-cpu-S3_nac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32   # 2026-10-29 用来跑 kl S4 取力帧（128 原子小超胞），32 核可在 192 核节点上并排多作业
#SBATCH --output=queue.out
#SBATCH --error=queue-%j.err
#SBATCH --qos=regular
#SBATCH --time=1-00:00:00     # regular QOS 墙时上限 1 天，必须写

cd $SLURM_SUBMIT_DIR

# 加载与编译时一致的 AOCC 环境（不要用 Intel）
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited

# [2026-11-02] jzzn 部分节点 InfiniBand(OpenFabrics mlx5) 初始化失败 → MPI 死锁、VASP 0 CPU 挂死。
#   强制走 TCP 规避 openib，与 hanhai25 的 xpmem 规避同源。
export OMPI_MCA_btl="self,vader,tcp"
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl_vader_single_copy_mechanism=none
export OMPI_MCA_osc=pt2pt

# 标准共线版本：step1/step2 使用 vasp_std
export AUTOZT_CELL_CONSTRAINT=none
mpirun -np $SLURM_NTASKS /public/home/wangchao/software/vasp.6.4.3/bin/vasp_std
