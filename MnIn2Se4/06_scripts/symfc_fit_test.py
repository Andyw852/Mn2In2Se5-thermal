#!/usr/bin/env python3
"""
独立的 phono3py + symfc 拟合脚本 (fc2 + fc3)
仅做拟合，不提交任何作业，输出 fc2.hdf5, fc3.hdf5, phono3py_params.yaml
"""

import os
import sys
import numpy as np
from phonopy.interface.vasp import read_vasp
from phonopy.structure.cells import get_supercell
from phono3py import Phono3py
from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5

# ---------- 参数配置 ----------
PRIM_FILE = "POSCAR"            # 原胞
SPOSCAR_FILE = "SPOSCAR"        # 超胞 (VASP 用的胞)
DISPS_FILE = "dataset_disps.npy"
FORCES_FILE = "dataset_forces.npy"
FC3_CUTOFF = 5.5                # 截断半径 (Å)，设为 None 表示不截断
BORN_FILE = "BORN"              # 可选，如存在则读入 NAC

# ----- 抽帧控制（可减少用于拟合的结构数）-----
NDATA = "all"                   # 可设为整数，比如 20，只取前 N 帧或均匀抽样
MODE = "even"                   # "even" 或 "head"

# ---------------------------------------------

def main():
    # 1. 读入
    uc = read_vasp(PRIM_FILE)
    sc_vasp = read_vasp(SPOSCAR_FILE)
    disps_in = np.load(DISPS_FILE)      # (nframe+1, natom, 3) 分数坐标
    forces_in = np.load(FORCES_FILE)    # (nframe+1, natom, 3) eV/Å
    natom = len(sc_vasp.numbers)
    assert disps_in.shape[1] == natom == forces_in.shape[1], \
        "原子数不匹配"

    # 2. 自动反推超胞矩阵 (pymatgen → phonopy 约定)
    R = sc_vasp.cell @ np.linalg.inv(uc.cell)
    M_pmg = np.round(R).astype(int)
    if not np.allclose(R, M_pmg, atol=1e-4):
        sys.exit(f"SPOSCAR 晶格不是 POSCAR 的整数倍组合\nR=\n{R}")
    M = M_pmg.T   # 转置为 phonopy 列约定
    if np.linalg.det(M) < 0:
        M = M[:, [1, 0, 2]]   # 修手性
    print(f"超胞矩阵 M (phonopy约定):\n{M}\ndet = {int(round(np.linalg.det(M)))}")

    # 3. 位移和力：扣除平衡帧
    eq_frac = disps_in[-1]
    d_frac = disps_in[:-1] - eq_frac[None]
    d_frac -= np.round(d_frac)          # 最小镜像
    u_vasp = d_frac @ sc_vasp.cell      # 笛卡尔位移
    f_vasp = forces_in[:-1] - forces_in[-1][None]
    ntot = len(u_vasp)
    print(f"总训练帧数: {ntot}, 平衡残余力 max = {np.linalg.norm(forces_in[-1], axis=1).max():.5f} eV/Å")

    # 4. 抽帧（可选）
    if NDATA == "all":
        idx = np.arange(ntot)
    else:
        n = int(NDATA)
        if MODE == "head":
            idx = np.arange(n)
        else:   # even
            idx = np.round(np.linspace(0, ntot-1, n)).astype(int)
            idx = np.unique(idx)
            while len(idx) < n:
                rest = np.setdiff1d(np.arange(ntot), idx)
                idx = np.sort(np.concatenate([idx, rest[:n-len(idx)]]))
    u, f = u_vasp[idx], f_vasp[idx]
    nuse = len(u)
    print(f"参与拟合: {nuse} 帧 → {nuse*natom*3} 个力方程")
    with open("ndata_used.txt", "w") as fo:
        fo.write(str(nuse))

    # 5. 原子映射：保证 VASP 超胞和 phonopy 超胞原子顺序一致
    sc_ph = get_supercell(uc, M)
    assert len(sc_ph.numbers) == natom
    # 验证两胞晶格是否同一超晶格
    T = sc_vasp.cell @ np.linalg.inv(sc_ph.cell)
    if not (np.allclose(T, np.round(T), atol=1e-5) and abs(abs(np.linalg.det(T))-1.0)<1e-6):
        sys.exit("两个超胞不描述同一超晶格，检查矩阵")
    # 根据坐标寻找置换
    cart_va = sc_vasp.scaled_positions @ sc_vasp.cell
    fp_va = cart_va @ np.linalg.inv(sc_ph.cell)
    fp_ph = sc_ph.scaled_positions
    perm = np.full(natom, -1, dtype=int)
    used = np.zeros(natom, dtype=bool)
    for i in range(natom):
        d = fp_va - fp_ph[i]
        d -= np.round(d)
        dist = np.linalg.norm(d @ sc_ph.cell, axis=1)
        dist[used] = np.inf
        j = int(dist.argmin())
        if dist[j] > 1e-3 or sc_vasp.numbers[j] != sc_ph.numbers[i]:
            sys.exit(f"原子映射失败 i={i}, min dist={dist[j]:.4f}")
        perm[i] = j; used[j] = True
    u = u[:, perm]
    f = f[:, perm]
    print("原子映射完成")

    # 6. 拟合
    ph3 = Phono3py(uc, supercell_matrix=M, primitive_matrix=np.eye(3))
    ph3.dataset = {"displacements": u, "forces": f}
    # fc2
    print("\n拟合 fc2 ...")
    ph3.produce_fc2(fc_calculator="symfc")
    # fc3
    opts = None if FC3_CUTOFF is None else f"cutoff = {FC3_CUTOFF}"
    print(f"拟合 fc3 (options: {opts}) ...")
    ph3.produce_fc3(fc_calculator="symfc", fc_calculator_options=opts)

    # 7. 写出
    p2s = ph3.primitive.p2s_map
    write_fc2_to_hdf5(ph3.fc2, filename="fc2.hdf5", p2s_map=p2s)
    write_fc3_to_hdf5(ph3.fc3, filename="fc3.hdf5", p2s_map=p2s)

    # 8. NAC（如果有 BORN）
    if os.path.isfile(BORN_FILE):
        from phonopy.file_IO import parse_BORN
        ph3.nac_params = parse_BORN(ph3.primitive, filename=BORN_FILE)
        print(f"已读入 {BORN_FILE}")
    ph3.save("phono3py_params.yaml")
    print("✅ 已输出 fc2.hdf5, fc3.hdf5, phono3py_params.yaml")

    # 9. 验证力误差
    fc2 = ph3.fc2 if ph3.fc2.shape[0] == natom else None
    fc3 = ph3.fc3 if ph3.fc3.shape[0] == natom else None
    if fc2 is not None and fc3 is not None:
        f_pred = -np.einsum("iajb,njb->nia", fc2.transpose(0,2,1,3), u)
        for n in range(len(u)):
            f_pred[n] -= 0.5 * np.einsum("ijkabc,jb,kc->ia", fc3, u[n], u[n])
        err = np.linalg.norm(f_pred - f) / np.linalg.norm(f)
        print(f"\nRelative error: {err:.6f}")
        with open("rel_err.txt", "w") as fo:
            fo.write(f"{err:.6f}")
    else:
        print("ℹ️  compact fc, 跳过力验证")

if __name__ == "__main__":
    main()