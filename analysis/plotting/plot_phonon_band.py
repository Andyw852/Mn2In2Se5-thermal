#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_phonon_band.py — 用拟合得到的 fc2.hdf5 计算并绘制声子谱 (+虚频体检)

用法: 在拟合完成的目录 (含 POSCAR / SPOSCAR / fc2.hdf5) 下直接运行
      python3 plot_phonon_band.py
输出: band.png  band.yaml  (终端打印虚频检查结果)
超胞矩阵自动从 SPOSCAR 反推 (含 pymatgen→phonopy 约定转换与手性修正, 与 v4 拟合脚本一致)
"""
import os, sys
import numpy as np
import h5py

PRIM     = "POSCAR"
SPOSCAR  = "SPOSCAR"
FC2_H5   = "fc2.hdf5"
BORN_F   = "BORN"
NPOINTS  = 101          # 每段路径的采样点数
MESH     = [11, 11, 11] # 虚频体检用的网格
IMAG_TOL = -0.05        # THz, 低于此值视为虚频 (容忍数值噪声)

from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

for f in (PRIM, SPOSCAR, FC2_H5):
    if not os.path.isfile(f):
        sys.exit(f"❌ 缺少 {f}")

# ── 超胞矩阵自动检测 (同 v4) ──
uc = read_vasp(PRIM)
sc_vasp = read_vasp(SPOSCAR)
R = sc_vasp.cell @ np.linalg.inv(uc.cell)
M_pmg = np.round(R).astype(int)
if not np.allclose(R, M_pmg, atol=1e-4):
    sys.exit("❌ SPOSCAR 与 POSCAR 不配套")
M = M_pmg.T
if np.linalg.det(M) < 0:
    M = M[:, [1, 0, 2]]
print(f"超胞矩阵 (phonopy约定) = {M.tolist()}")

# ── 建 Phonopy 对象, 读 fc2 ──
ph = Phonopy(uc, supercell_matrix=M, primitive_matrix=np.eye(3))
with h5py.File(FC2_H5) as f:
    fc2 = f["force_constants"][:]
natom_sc = len(ph.supercell.numbers)
if fc2.shape[0] == fc2.shape[1] == natom_sc:
    pass                                    # full
elif fc2.shape[0] == len(ph.primitive.numbers) and fc2.shape[1] == natom_sc:
    pass                                    # compact, phonopy 直接接受
else:
    sys.exit(f"❌ fc2 形状 {fc2.shape} 与超胞({natom_sc})不符")
ph.force_constants = fc2

if os.path.isfile(BORN_F):
    from phonopy.file_IO import parse_BORN
    ph.nac_params = parse_BORN(ph.primitive, filename=BORN_F)
    print(f"✅ NAC: 已读入 {BORN_F}")
else:
    print(f"ℹ️  无 {BORN_F}, 不含 NAC (极性材料 LO-TO 劈裂缺失, 建议补 DFPT)")

# ── 高对称路径: 优先 seekpath 自动, 失败则退回菱方通用路径 ──
try:
    from phonopy.phonon.band_structure import get_band_qpoints_by_seekpath
    bands, labels, connections = get_band_qpoints_by_seekpath(ph.primitive, NPOINTS)
    print(f"路径: seekpath 自动 ({' '.join(labels)})")
except Exception as e:
    print(f"ℹ️  seekpath 不可用 ({type(e).__name__}), 用菱方通用路径 Γ-L-F-Γ-Z")
    paths = [[[0,0,0],[0.5,0,0]], [[0.5,0,0],[0.5,0.5,0]],
             [[0.5,0.5,0],[0,0,0]], [[0,0,0],[0.5,0.5,0.5]]]
    bands, connections = get_band_qpoints_and_path_connections(paths, npoints=NPOINTS)
    labels = ["$\\Gamma$", "L", "F", "$\\Gamma$", "Z"]

ph.run_band_structure(bands, path_connections=connections, labels=labels,
                      with_eigenvectors=False)
ph.write_yaml_band_structure(filename="band.yaml")

import matplotlib
matplotlib.use("Agg")
plt = ph.plot_band_structure()
plt.savefig("band.png", dpi=200, bbox_inches="tight")
print("✅ band.png / band.yaml 已写出")

# ── 虚频体检 ──
bs = ph.get_band_structure_dict()
fmin_band = min(fr.min() for fr in bs["frequencies"])
ph.run_mesh(MESH)
fmin_mesh = ph.get_mesh_dict()["frequencies"].min()
# Γ 点声学支
ph.run_qpoints([[0, 0, 0]])
f_gamma = np.sort(ph.get_qpoints_dict()["frequencies"][0])[:3]
print(f"\n[体检] Γ 点声学支: {np.round(f_gamma, 4)} THz (应≈0, ASR 检查)")
print(f"[体检] 路径最低频率 = {fmin_band:.4f} THz | {MESH} 网格最低频率 = {fmin_mesh:.4f} THz")
fmin = min(fmin_band, fmin_mesh)
if fmin < IMAG_TOL:
    print(f"⚠️  存在虚频 ({fmin:.3f} THz)! 先排查: 结构未充分弛豫 / 数据不足 / "
          f"位移振幅不当, 再考虑算热导率")
else:
    print("✅ 无显著虚频, 可以进入热导率计算")
