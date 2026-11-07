#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kappa_extrinsic.py —— 本征 κ（phono3py RTA）+ Mn/In 反位质量无序 + 晶界散射 + 织构，对实验 κ 拟合。

思路：散射率可加
    Γ_λ = Γ_ph-ph,λ(T) + Σ_元素 g_元素(x)·B_元素,λ + Γ_晶界,λ(L)
  * Γ_ph-ph：一次 phono3py RTA（不开同位素）；
  * B_元素：Tamura 同位素/质量无序散射对 g 是线性的 —— 每种元素的位点单独设 g = 1，
    用 read_gamma 只算质量无序那一项（很便宜）；
  * Γ_晶界 = |v| / (4π L)，与 phono3py --bmfp 同一个式子（直接调用 phono3py 的函数）。
  于是任意 (x, L) 的 κ 只是对模式求和：
      κ(x, L) = Σ_λ κ_λ^(0) · Γ0_λ / (Γ0_λ + Γ_iso,λ(x) + Γ_b,λ(L))
  κ_λ^(0) 是只有 ph-ph 散射时的 RTA 模式热导率。

无序参数 x 的定义（两个材料通用）：
  x = 跑到 In 位上的 Mn 原子占全部 Mn 的比例（Mn↔In 成对交换，化学计量不变）。
  Mn 位：(1-x) Mn + x In；In 位：Mn 的比例 = x·N_Mn/N_In。
  完全无序上限 x_rand = N_In/(N_Mn+N_In)：Mn2In2Se5 = 0.5，MnIn2Se4 = 2/3。
  只算了质量差（下限）；键强/局域应变的无序会让真实散射更强，拟合出的 x 是“等效 x”。

织构参数 s = <sin²θ>（θ = 晶粒 c 轴与热压方向的夹角；s = 0 完全织构，s = 2/3 完全随机）：
  κ_Para = s·κ_xx + (1-s)·κ_zz            （平行热压方向 ≈ 面外）
  κ_Perp = (1-s/2)·κ_xx + (s/2)·κ_zz       （垂直热压方向 ≈ 面内）
  (2·Perp + Para)/3 = (2κ_xx + κ_zz)/3 与 s 无关。这是并联（Voigt）平均，是上限。

子命令
  bases  算 Γ0、κ_λ^(0)、B_元素、群速度，存 <out>/bases.npz
  eval   给定 x、L，打印 κ_xx、κ_zz(T)
  fit    对实验 CSV 拟合 (x, L, s)，写 fit_<material>.json

示例
  python3 kappa_extrinsic.py bases --poscar POSCAR --dim 4 4 1 --fc2 fc2.hdf5 --fc3 fc3.hdf5 \\
          --mesh 16 16 4 --out work_Mn2In2Se5
  python3 kappa_extrinsic.py eval --bases work_Mn2In2Se5/bases.npz --x 0.1 --L 1.0
  python3 kappa_extrinsic.py fit --bases work_Mn2In2Se5/bases.npz \\
          --exp ../../experiment/kappa_tot_digitized.csv --material Mn2In2Se5 --tmin 290
"""
import argparse
import csv
import json
import os
import sys
import time

import numpy as np

DEFAULT_T = [10, 20, 30, 40, 50, 60, 80, 100, 125, 150, 200, 250, 300, 325, 350, 375, 400,
             450, 500, 550, 600, 625, 650]
A_EL, B_EL = "Mn", "In"           # 交换的两种阳离子


# ---------------------------------------------------------------------------
# 质量无序 g
# ---------------------------------------------------------------------------
def _natural(el):
    try:
        from phonopy.atomic_data import get_atomic_data
        iso = get_atomic_data().isotope_data[el]
    except Exception:                         # 旧 phonopy
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from phonopy.structure.atoms import isotope_data
        iso = isotope_data[el]
    return [(m, a) for _, m, a in iso]


def _g(mix):
    mb = sum(m * f for m, f in mix)
    return sum(f * (1.0 - m / mb) ** 2 for m, f in mix)


def site_g(symbols, x):
    """每种元素位点的 g（含天然同位素 + Mn↔In 反位），返回 {元素: g}。"""
    n_a = sum(s == A_EL for s in symbols)
    n_b = sum(s == B_EL for s in symbols)
    fa_on_b = x * n_a / n_b if n_b else 0.0
    out = {}
    for el in sorted(set(symbols)):
        if el == A_EL:
            mix = [(m, (1 - x) * a) for m, a in _natural(A_EL)] + [(m, x * a) for m, a in _natural(B_EL)]
        elif el == B_EL:
            mix = [(m, (1 - fa_on_b) * a) for m, a in _natural(B_EL)] + \
                  [(m, fa_on_b * a) for m, a in _natural(A_EL)]
        else:
            mix = _natural(el)
        out[el] = _g(mix)
    return out


def x_random(symbols):
    n_a = sum(s == A_EL for s in symbols)
    n_b = sum(s == B_EL for s in symbols)
    return n_b / float(n_a + n_b)


# ---------------------------------------------------------------------------
# bases
# ---------------------------------------------------------------------------
def _load_cell(a):
    from phonopy.interface.vasp import read_vasp
    if a.disp_yaml:
        import phono3py
        p = phono3py.load(os.path.abspath(a.disp_yaml), produce_fc=False, log_level=0)
        return p.unitcell, np.array(p.supercell_matrix, int)
    if not (a.poscar and a.dim):
        sys.exit("[ERROR] 需要 --disp-yaml，或 --poscar + --dim")
    return read_vasp(os.path.abspath(a.poscar)), np.diag(a.dim)


def cmd_bases(a):
    from phono3py import Phono3py
    from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5
    from phono3py.conductivity.scattering_solvers import compute_bulk_boundary_scattering

    uc, smat = _load_cell(a)
    fc2p, fc3p = os.path.abspath(a.fc2), os.path.abspath(a.fc3)
    T = np.array(a.temps, float)
    os.makedirs(a.out, exist_ok=True)
    os.chdir(a.out)
    t0 = time.time()
    ph3 = Phono3py(uc, supercell_matrix=smat, primitive_matrix="P", log_level=0)
    ph3.fc2 = read_fc2_from_hdf5(fc2p)
    ph3.fc3 = read_fc3_from_hdf5(fc3p)
    ph3.mesh_numbers = a.mesh
    ph3.init_phph_interaction()
    symbols = list(ph3.primitive.symbols)
    print("[..] 原胞 %s，超胞 %d 原子，q 网格 %s，%d 个温度" %
          ("".join("%s%d" % (e, symbols.count(e)) for e in dict.fromkeys(symbols)),
           len(ph3.supercell), list(a.mesh), len(T)), flush=True)

    # 1) 本征：只有 ph-ph（同位素关），写 kappa-m*.hdf5 供下面 read_gamma
    ph3.run_thermal_conductivity(temperatures=T, is_isotope=False, write_kappa=True, log_level=0)
    tc = ph3.thermal_conductivity
    kappa0 = np.asarray(tc.kappa)[0]                 # (nT, 6)
    mk0 = np.asarray(tc.mode_kappa)[0]               # (nT, ngp, nb, 6)
    g0 = np.asarray(tc.gamma)[0]                     # (nT, ngp, nb)
    gv = np.asarray(tc.group_velocities)             # (ngp, nb, 3)
    freq = np.asarray(tc.frequencies)
    w = np.asarray(tc.grid_weights)
    # phono3py：kappa = Σ_{q,ν} mode_kappa / N（N = 网格点总数 = Σ weights）
    sum_err = float(np.abs(mk0.sum(axis=(1, 2)) / w.sum() - kappa0).max())
    print("[OK] 本征 RTA 完成（%.0f s）；Σ mode_kappa/N 与 kappa 最大差 %.2e" % (time.time() - t0, sum_err),
          flush=True)

    # 2) 每种元素位点 g = 1 的质量无序散射率（读回 ph-ph gamma，不重算）
    B = {}
    for el in dict.fromkeys(symbols):
        mv = [1.0 if s == el else 0.0 for s in symbols]
        t1 = time.time()
        ph3.run_thermal_conductivity(temperatures=T, is_isotope=True, mass_variances=mv,
                                     read_gamma=True, log_level=0)
        gi = np.asarray(ph3.thermal_conductivity.gamma_isotope)[0]
        B[el] = gi
        print("[OK] B_%s（%.0f s）" % (el, time.time() - t1), flush=True)

    gb_1um = compute_bulk_boundary_scattering(gv, 1.0)      # L = 1 μm 时的 Γ_b；∝ 1/L
    np.savez_compressed(
        "bases.npz", T=T, mesh=np.array(a.mesh), symbols=np.array(symbols), kappa0=kappa0,
        mk0=mk0[..., [0, 1, 2]], g0=g0, gb_1um=gb_1um, freq=freq, weights=w,
        **{"B_" + el: B[el] for el in B})
    meta = {"poscar": a.poscar, "disp_yaml": a.disp_yaml, "supercell_matrix": smat.tolist(),
            "fc2": fc2p, "fc3": fc3p, "mesh": list(a.mesh), "T": T.tolist(), "symbols": symbols,
            "seconds": round(time.time() - t0), "mode_kappa_sum_err": sum_err}
    json.dump(meta, open("bases_meta.json", "w"), indent=1, ensure_ascii=False)
    i300 = int(np.argmin(abs(T - 300)))
    print("[DONE] bases.npz；本征 κ(%g K) xx/zz = %.4f / %.4f" %
          (T[i300], kappa0[i300, 0], kappa0[i300, 2]), flush=True)


# ---------------------------------------------------------------------------
# 模型
# ---------------------------------------------------------------------------
class Model:
    def __init__(self, path):
        z = np.load(path)
        self.T = z["T"]
        self.symbols = [str(s) for s in z["symbols"]]
        # phono3py 的 mode_kappa 求和后还要除以网格点总数 N 才是 κ
        self.mk0 = z["mk0"] / float(z["weights"].sum())   # (nT, ngp, nb, 3) : xx, yy, zz
        err = np.abs(self.mk0.sum(axis=(1, 2)) - z["kappa0"][:, :3]).max()
        if err > 1e-6 * max(1.0, np.abs(z["kappa0"]).max()):
            sys.exit("[ERROR] Σ mode_kappa/N 与 kappa 不一致（差 %.3g），bases.npz 不对" % err)
        self.g0 = z["g0"]
        self.gb1 = z["gb_1um"]
        self.B = {k[2:]: z[k] for k in z.files if k.startswith("B_")}
        self.kappa0 = z["kappa0"]
        self.xr = x_random(self.symbols)
        self._mk0_inplane = 0.5 * (self.mk0[..., 0] + self.mk0[..., 1])
        self._mk0_zz = self.mk0[..., 2]
        with np.errstate(divide="ignore", invalid="ignore"):
            self._ok = self.g0 > 0

    def gamma_iso(self, x):
        g = site_g(self.symbols, x)
        return sum(g[el] * self.B[el] for el in self.B), g

    def kappa(self, x, L_um):
        """返回 (κ_inplane(T), κ_zz(T))，T = self.T。L_um = None/inf 表示无晶界。"""
        gi, _ = self.gamma_iso(x)
        gt = self.g0 + gi[None]
        if L_um is not None and np.isfinite(L_um):
            gt = gt + self.gb1[None] / L_um
        r = np.where(self._ok, self.g0 / np.where(gt > 0, gt, 1.0), 0.0)
        return (self._mk0_inplane * r).sum(axis=(1, 2)), (self._mk0_zz * r).sum(axis=(1, 2))


def _interp_log(Tg, y, Tq):
    return np.exp(np.interp(Tq, Tg, np.log(np.maximum(y, 1e-12))))


def cmd_eval(a):
    m = Model(a.bases)
    L = None if a.L in (None, 0) else a.L
    xx, zz = m.kappa(a.x, L)
    g = site_g(m.symbols, a.x)
    print("# x = %.3f（x_rand = %.3f），L = %s μm，g = %s" %
          (a.x, m.xr, "∞" if L is None else a.L, {k: round(v, 6) for k, v in g.items()}))
    print("# T(K)   κ_xx     κ_zz     (2xx+zz)/3   本征 xx   本征 zz")
    for i, t in enumerate(m.T):
        print("%6.1f %8.4f %8.4f %10.4f %10.4f %8.4f" %
              (t, xx[i], zz[i], (2 * xx[i] + zz[i]) / 3, m.kappa0[i, 0], m.kappa0[i, 2]))


# ---------------------------------------------------------------------------
# 拟合
# ---------------------------------------------------------------------------
def read_exp(path, material, tmin, tmax):
    rows = []
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(l for l in f if not l.startswith("#")):
            if r["material"] != material:
                continue
            t = float(r["T_K"])
            if t < tmin or t > tmax:
                continue
            ke = float(r.get("kappa_e_WmK") or 0.0)
            rows.append((r["direction"], t, float(r["kappa_tot_WmK"]) - ke))
    if not rows:
        sys.exit("[ERROR] %s 里没有 %s 在 %g–%g K 的数据" % (path, material, tmin, tmax))
    return rows


def cmd_fit(a):
    m = Model(a.bases)
    rows = read_exp(a.exp, a.material, a.tmin, a.tmax)
    d_para = np.array([r[0] == "Para" for r in rows])
    t_exp = np.array([r[1] for r in rows])
    k_exp = np.array([r[2] for r in rows])
    if t_exp.min() < m.T.min() or t_exp.max() > m.T.max():
        sys.exit("[ERROR] 实验温度 %g–%g K 超出 bases 的温度范围 %g–%g K" %
                 (t_exp.min(), t_exp.max(), m.T.min(), m.T.max()))
    xs = np.round(np.arange(0.0, m.xr + 1e-9, a.dx), 4)
    Ls = np.concatenate([[np.inf], np.logspace(np.log10(a.lmin), np.log10(a.lmax), a.nl)])
    ss = np.linspace(0.0, 2.0 / 3.0, 201)
    best = None
    grid = np.full((len(xs), len(Ls)), np.nan)
    sgrid = np.full_like(grid, np.nan)
    for i, x in enumerate(xs):
        for j, L in enumerate(Ls):
            xx, zz = m.kappa(x, None if not np.isfinite(L) else L)
            kx = _interp_log(m.T, xx, t_exp)
            kz = _interp_log(m.T, zz, t_exp)
            # (n_s, n_pts)
            para = ss[:, None] * kx[None] + (1 - ss[:, None]) * kz[None]
            perp = (1 - ss[:, None] / 2) * kx[None] + (ss[:, None] / 2) * kz[None]
            mod = np.where(d_para[None], para, perp)
            rms = np.sqrt(np.mean(np.log(mod / k_exp[None]) ** 2, axis=1))
            k = int(np.argmin(rms))
            grid[i, j], sgrid[i, j] = rms[k], ss[k]
            if best is None or rms[k] < best[0]:
                best = (rms[k], x, L, ss[k], kx, kz)
    rms0, x0, L0, s0, kx0, kz0 = best
    tol = max(rms0 * a.tol_factor, rms0 + a.tol_abs)
    ok = grid <= tol
    xi, lj = np.where(ok)
    region = {"rms_threshold": round(float(tol), 4),
              "x_range": [float(xs[xi].min()), float(xs[xi].max())],
              "L_um_range": [float(Ls[lj].min()), float(Ls[lj].max())],
              "s_range": [float(sgrid[ok].min()), float(sgrid[ok].max())]}
    # 本征（x=0, L=∞）对同一批实验点的最佳织构
    xx_i, zz_i = m.kappa(0.0, None)
    res = {
        "material": a.material, "exp_file": a.exp, "T_range": [a.tmin, a.tmax], "n_points": len(rows),
        "model": "RTA; Γ = Γ_phph + Σ g(x)·B + |v|/(4πL); texture Voigt (s = <sin²θ>)",
        "best": {"x": float(x0), "x_random_limit": round(m.xr, 4),
                 "L_um": None if not np.isfinite(L0) else round(float(L0), 4),
                 "s": round(float(s0), 4), "rms_log_err": round(float(rms0), 4),
                 "g": {k: round(v, 6) for k, v in site_g(m.symbols, x0).items()}},
        "acceptable_region": region,
        "intrinsic_rms_log_err": round(float(grid[0, 0]), 4),
        "profile_x": [{"x": float(x), "rms": round(float(np.nanmin(grid[i])), 4),
                       "L_um": (None if not np.isfinite(Ls[int(np.nanargmin(grid[i]))])
                                else round(float(Ls[int(np.nanargmin(grid[i]))]), 4))}
                      for i, x in enumerate(xs)],
        "profile_L": [{"L_um": (None if not np.isfinite(L) else round(float(L), 4)),
                       "rms": round(float(np.nanmin(grid[:, j])), 4),
                       "x": float(xs[int(np.nanargmin(grid[:, j]))])}
                      for j, L in enumerate(Ls)],
        "points": [],
        "model_curve": [],
    }
    for r, kx, kz in zip(rows, kx0, kz0):
        mod = s0 * kx + (1 - s0) * kz if r[0] == "Para" else (1 - s0 / 2) * kx + (s0 / 2) * kz
        res["points"].append({"direction": r[0], "T": r[1], "exp": round(r[2], 4),
                              "model": round(float(mod), 4), "ratio": round(float(mod / r[2]), 3)})
    xx, zz = m.kappa(x0, None if not np.isfinite(L0) else L0)
    for i, t in enumerate(m.T):
        res["model_curve"].append({"T": float(t), "xx": round(float(xx[i]), 4), "zz": round(float(zz[i]), 4),
                                   "xx_intrinsic": round(float(xx_i[i]), 4),
                                   "zz_intrinsic": round(float(zz_i[i]), 4)})
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(a.bases)), "fit_%s.json" % a.material)
    json.dump(res, open(out, "w"), indent=1, ensure_ascii=False)
    b = res["best"]
    print("[FIT] %s  %d 点 %g–%g K" % (a.material, len(rows), a.tmin, a.tmax))
    print("  最优 x = %.3f（上限 %.3f），L = %s μm，s = %.3f，rms(ln) = %.3f；本征 rms = %.3f" %
          (b["x"], m.xr, b["L_um"] if b["L_um"] is not None else "∞", b["s"], b["rms_log_err"],
           res["intrinsic_rms_log_err"]))
    print("  可接受区（rms ≤ %.3f）：x %s，L %s μm，s %s" %
          (tol, region["x_range"], region["L_um_range"], region["s_range"]))
    for p in res["points"]:
        print("   %-4s %6.1f K  exp %.3f  model %.3f  (%.2f)" %
              (p["direction"], p["T"], p["exp"], p["model"], p["ratio"]))
    print("[OK] 写出 %s" % out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    b = sp.add_parser("bases")
    b.add_argument("--poscar")
    b.add_argument("--dim", type=int, nargs=3)
    b.add_argument("--disp-yaml", help="phono3py_disp.yaml（给了就从它读单胞与超胞矩阵）")
    b.add_argument("--fc2", required=True)
    b.add_argument("--fc3", required=True)
    b.add_argument("--mesh", type=int, nargs=3, required=True)
    b.add_argument("--temps", type=float, nargs="+", default=DEFAULT_T)
    b.add_argument("--out", required=True)
    e = sp.add_parser("eval")
    e.add_argument("--bases", required=True)
    e.add_argument("--x", type=float, default=0.0)
    e.add_argument("--L", type=float, default=None, help="晶界平均自由程 μm；不给 = 无晶界")
    f = sp.add_parser("fit")
    f.add_argument("--bases", required=True)
    f.add_argument("--exp", required=True)
    f.add_argument("--material", required=True)
    f.add_argument("--tmin", type=float, default=290.0)
    f.add_argument("--tmax", type=float, default=700.0)
    f.add_argument("--dx", type=float, default=0.01)
    f.add_argument("--lmin", type=float, default=0.005, help="μm")
    f.add_argument("--lmax", type=float, default=100.0, help="μm")
    f.add_argument("--nl", type=int, default=61)
    f.add_argument("--tol-factor", type=float, default=1.25)
    f.add_argument("--tol-abs", type=float, default=0.02)
    f.add_argument("--out")
    a = ap.parse_args()
    {"bases": cmd_bases, "eval": cmd_eval, "fit": cmd_fit}[a.cmd](a)


if __name__ == "__main__":
    main()
