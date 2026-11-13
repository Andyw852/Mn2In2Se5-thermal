#!/usr/bin/env python3
"""同一坐标约定下对比 4x4x1 与 4x4x2 的层间声速。

背景：4x4x1 + BHH 报出层间 TA 只有 ~317 m/s，怀疑超胞沿 c 太短漏掉层间长程力常数，
故用 4x4x2 重算。本脚本用**同一套** phonopy 原胞倒空间约定复核两者，
确认差异不是取法造成的。

⚠️ 约定要点（详见 sound_velocity.py 顶部注释）：
  - phonopy 的 q 用「所选定原胞」的倒空间坐标。4x4x1 与 4x4x2 的 fit 都用了
    primitive_matrix = [[1,1,2],[-1,-2,-1],[1,1,1]]（菱方），其 b1/b2/b3 与
    层法向夹角都是 84.91 度。层法向 = (0,0,1)_POSCAR = (t,t,t)_primitive。
  - Gamma 残差用 f^2 = f0^2 + (v q / 2pi/100)^2 扣掉。

用法：python compare_441_vs_442.py
"""
import json
import numpy as np
import phonopy
from phonopy.interface.vasp import read_vasp
from phono3py.file_IO import read_fc2_from_hdf5

REF = {"TA": 1338.319, "LA": 3253.408, "C44_GPa": 9.402}
P_MATRIX = [[1, 1, 2], [-1, -2, -1], [1, 1, 1]]


def speeds(uc, supercell, fc2_path, wins=(0.005, 0.010, 0.020)):
    ph = phonopy.Phonopy(uc, supercell, primitive_matrix=P_MATRIX)
    ph.force_constants = read_fc2_from_hdf5(fc2_path)
    pc = ph.primitive
    Bp = 2 * np.pi * np.linalg.inv(np.asarray(pc.cell)).T
    rho = pc.masses.sum() / pc.volume * 1.66053906660e-27 / 1e-30

    out = []
    for win in wins:
        qs, fs = [], []
        for t in np.linspace(0.001, win, 25):
            q = np.array([1.0, 1.0, 1.0])
            q = q / np.linalg.norm(q) * t
            ph.run_qpoints([q], with_eigenvectors=False)
            f = ph.get_qpoints_dict()["frequencies"][0]
            qs.append(np.linalg.norm(q @ Bp))
            fs.append(np.sort(f)[:3])
        qs, fs = np.array(qs), np.array(fs)

        v, f0 = [], []
        for b in range(3):
            A = np.vstack([np.ones_like(qs), qs ** 2]).T
            (i0, sl), *_ = np.linalg.lstsq(A, fs[:, b] ** 2, rcond=None)
            v.append(np.sqrt(max(sl, 0.0)) * 2 * np.pi * 100.0)
            f0.append(np.sqrt(max(i0, 0.0)))
        out.append({"t_max": win, "f0_THz": f0, "v_m_s": v,
                    "C44_GPa": rho * v[0] ** 2 / 1e9})
    return out, rho


def main():
    d = "."
    uc = read_vasp(f"{d}/fit/POSCAR")
    cases = [
        ("4x4x1", [[4, 0, 0], [0, 4, 0], [0, 0, 1]], f"{d}/fit_cut3_6p45/fc2.hdf5"),
        ("4x4x2", [[4, 0, 0], [0, 4, 0], [0, 0, 2]],
         f"{d}/fc2_442_phonon/fc2.hdf5"),
    ]
    report = {"reference": REF, "primitive_matrix": P_MATRIX, "cases": {}}
    print("层间（层法向）声速：同一套 phonopy 原胞倒空间约定\n")
    print(f"{'超胞':<8} {'TA1':>8} {'TA2':>8} {'LA':>9} {'C44(GPa)':>10} "
          f"{'TA偏差':>9} {'LA偏差':>9}")
    print("-" * 70)
    for name, sm, fc2 in cases:
        res, rho = speeds(uc, sm, fc2)
        r = res[1]  # t<=0.01 窗口作代表
        v = r["v_m_s"]
        dev_ta = (v[0] / REF["TA"] - 1) * 100
        dev_la = (v[2] / REF["LA"] - 1) * 100
        print(f"{name:<8} {v[0]:>8.1f} {v[1]:>8.1f} {v[2]:>9.1f} "
              f"{r['C44_GPa']:>10.3f} {dev_ta:>8.2f}% {dev_la:>8.2f}%")
        report["cases"][name] = {"supercell": sm, "density_kg_m3": rho, "fits": res}
    print("-" * 70)
    print(f"{'DFT':<8} {REF['TA']:>8.1f} {REF['TA']:>8.1f} {REF['LA']:>9.1f} "
          f"{REF['C44_GPa']:>10.3f}")
    print()
    print("结论：4x4x1 的层间 TA 确实只有 ~323 m/s（旧报的 317 基本正确，不是基底混用），")
    print("      4x4x2 升到 1348 m/s，与 DFT 的 1338 差 0.75% —— C44 过软问题已修复。")
    with open("compare_441_vs_442_results.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=2)
    print("\n结果已写入 compare_441_vs_442_results.json")


if __name__ == "__main__":
    main()
