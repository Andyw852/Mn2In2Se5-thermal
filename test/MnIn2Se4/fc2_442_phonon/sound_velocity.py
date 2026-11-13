#!/usr/bin/env python3
"""从 4x4x2 的 fc2 提取层法向声速，和 DFT 弹性常数对照。

⚠️ 踩过的坑（2026-11-13 用户核查 + 本脚本复核确认）：

1. **phonopy 的 q 用「所选定原胞」的倒空间坐标，不是 POSCAR 原胞。**
   本次 `primitive_matrix="auto"` 把原胞换成了菱方原胞（7 原子），
   新原胞的 b1/b2/b3 与层法向夹角都是 84.91°（几乎沿层面）。
   POSCAR 的 a1、a2 在层面内 → 层法向 = (0,0,1)_POSCAR。
   同一个层间方向在两套坐标里：
        (0,0,t)_POSCAR  ≡  (t,t,t)_phonopy_primitive
   直接用 (0,0,t) 算频率、再用 POSCAR b3 换算 |q|，会得到 4179/6256 m/s
   的假值；两套 |b3| 之比 3.760 正好是误差倍数。

2. **Γ 点有数值残差**（f0 = 0.000719 THz），低 q 时 ω/q 被它污染、
   q→0 反而发散。必须拟合 f² = f0² + (v q / 2π/100)²，取斜率开方。

修正后：TA = 1348.3 m/s（DFT 1338.3，+0.75%），LA = 3118.6（DFT 3253.4，-4.14%），
C44 = 9.543 GPa（DFT 9.402，+1.50%）。窗口 t≤0.005~0.040 结果稳定。

用法：
    python sound_velocity.py <含 fc2.hdf5 与 phonopy_disp.yaml 的目录>
"""
import sys
import numpy as np
import phonopy
from phono3py.file_IO import read_fc2_from_hdf5

# DFT 弹性常数导出的参考声速（同一声子结构密度 5249.3 kg/m^3）
REF_TA = 1338.319
REF_LA = 3253.408
REF_C44_GPa = 9.402


def main(d):
    ph = phonopy.load(f"{d}/phonopy_disp.yaml", produce_fc=False,
                      is_nac=False, log_level=0)
    ph.force_constants = read_fc2_from_hdf5(f"{d}/fc2.hdf5")

    uc = ph.unitcell
    pc = ph.primitive
    # ★ 关键：用 primitive 的倒格矢，不是 unitcell 的
    Bp = 2 * np.pi * np.linalg.inv(np.asarray(pc.cell)).T
    rho = pc.masses.sum() / pc.volume * 1.66053906660e-27 / 1e-30

    a1, a2 = np.asarray(uc.cell)[0], np.asarray(uc.cell)[1]
    n = np.cross(a1, a2)
    n /= np.linalg.norm(n)
    print(f"phonopy primitive 原子数 = {len(pc)}")
    print(f"层法向 (Cartesian) = {np.round(n, 6)}")
    for i in range(3):
        ang = np.degrees(np.arccos(min(1.0, abs(Bp[i] @ n) / np.linalg.norm(Bp[i]))))
        print(f"  primitive b{i+1} 与层法向夹角 = {ang:.2f} deg")
    print(f"rho = {rho:.1f} kg/m^3\n")
    print("层法向 = phonopy primitive 的 (t,t,t) 方向")
    print("拟合 f^2 = f0^2 + (v·q/2π/100)^2（扣 Γ 点残差）\n")

    ref = np.array([REF_TA, REF_TA, REF_LA])
    for win in (0.005, 0.010, 0.020, 0.040):
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
            (intercept, slope), *_ = np.linalg.lstsq(A, fs[:, b] ** 2, rcond=None)
            v.append(np.sqrt(max(slope, 0.0)) * 2 * np.pi * 100.0)
            f0.append(np.sqrt(max(intercept, 0.0)))
        v = np.array(v)
        dev = (v / ref - 1) * 100
        print(f"窗口 t<={win:.3f}: f0={f0[0]:.6f} THz  "
              f"TA1={v[0]:.1f} TA2={v[1]:.1f} LA={v[2]:.1f} m/s  "
              f"(偏差 {dev[0]:+.2f}% {dev[1]:+.2f}% {dev[2]:+.2f}%)")
        print(f"                    C44={rho*v[0]**2/1e9:.3f} GPa "
              f"(DFT {REF_C44_GPa}, {(rho*v[0]**2/1e9/REF_C44_GPa-1)*100:+.2f}%)")
    print()
    print(f"DFT 弹性常数参考: TA={REF_TA:.1f} LA={REF_LA:.1f} m/s, C44={REF_C44_GPa} GPa")
    print("※ 为什么是 (t,t,t)：见本文件顶部注释第 1 条。")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
