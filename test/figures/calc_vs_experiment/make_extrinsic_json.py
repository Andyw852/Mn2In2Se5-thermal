#!/usr/bin/env python3
"""L3 curves for the figures with the antisite fraction x fixed at the literature value.

The free (x, L, s) fit in analysis/extrinsic/kappa_extrinsic.py leaves x undetermined
(x and L are degenerate). Here x is fixed from a structure refinement and only the
boundary mean free path L and the texture s are fitted, on the same grid and with the
same metric (rms of ln(model/exp), 290-700 K) as `kappa_extrinsic.py fit`.

  MnIn2Se4   x = 0.44  Range et al., Z. Naturforsch. 46b, 1122 (1991), single crystal
  Mn2In2Se5  x = 0.19  Williams et al., J. Mater. Chem. C 12, 1753 (2024), ESI Table S7

Writes extrinsic_<material>.json (rows T, kxx, kzz, intrinsic_kxx, intrinsic_kzz and
the fitted s), the format read by plot_calc_vs_experiment.py and
../paper/make_paper_figures.py.

  python3 make_extrinsic_json.py        # needs numpy only (bases.npz is in the repository)
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "analysis" / "extrinsic"))
import kappa_extrinsic as ke  # noqa: E402

CASES = {
    "MnIn2Se4": {"bases": "MnIn2Se4_442_m24", "x": 0.44,
                 "x_source": "Range et al., Z. Naturforsch. 46b, 1122 (1991)"},
    "Mn2In2Se5": {"bases": "Mn2In2Se5_m24", "x": 0.19,
                  "x_source": "Williams et al., J. Mater. Chem. C 12, 1753 (2024), ESI Table S7"},
}
LS = np.concatenate([[np.inf], np.logspace(np.log10(0.005), np.log10(100.0), 61)])   # um
SS = np.linspace(0.0, 2.0 / 3.0, 201)


def rms_over_s(m, x, L, dp, te, kexp):
    xx, zz = m.kappa(x, None if not np.isfinite(L) else L)
    kx = ke._interp_log(m.T, xx, te)
    kz = ke._interp_log(m.T, zz, te)
    para = SS[:, None] * kx[None] + (1 - SS[:, None]) * kz[None]
    perp = (1 - SS[:, None] / 2) * kx[None] + (SS[:, None] / 2) * kz[None]
    return np.sqrt(np.mean(np.log(np.where(dp[None], para, perp) / kexp[None]) ** 2, axis=1))


def main():
    for mat, c in CASES.items():
        m = ke.Model(str(REPO / "analysis/extrinsic/results" / c["bases"] / "bases.npz"))
        rows = ke.read_exp(str(REPO / "experiment/kappa_tot_digitized.csv"), mat, 290.0, 700.0)
        dp = np.array([r[0] == "Para" for r in rows])
        te = np.array([r[1] for r in rows])
        kexp = np.array([r[2] for r in rows])
        grid = np.array([rms_over_s(m, c["x"], L, dp, te, kexp) for L in LS])   # (nL, ns)
        j, k = np.unravel_index(np.argmin(grid), grid.shape)
        rms0, L0, s0 = float(grid[j, k]), float(LS[j]), float(SS[k])
        tol = max(1.25 * rms0, rms0 + 0.02)
        okL, okS = np.where(grid <= tol)
        xx, zz = m.kappa(c["x"], L0)
        xi, zi = m.kappa0[:, 0], m.kappa0[:, 2]        # three-phonon only, as in the earlier files
        out = {"material": mat, "s": round(s0, 4), "x": c["x"], "x_source": c["x_source"],
               "L_nm": round(1000 * L0, 1), "rms_log_err": round(rms0, 4),
               "acceptable": {"rms_threshold": round(tol, 4),
                              "L_nm": [round(1000 * LS[okL].min(), 1), round(1000 * LS[okL].max(), 1)],
                              "s": [round(float(SS[okS].min()), 4), round(float(SS[okS].max()), 4)]},
               "note": "L3 with x fixed at the literature value; L and s fitted (make_extrinsic_json.py)",
               "rows": [{"T": float(t), "kxx": round(float(a), 4), "kzz": round(float(b), 4),
                         "intrinsic_kxx": round(float(p), 4), "intrinsic_kzz": round(float(q), 4)}
                        for t, a, b, p, q in zip(m.T, xx, zz, xi, zi)]}
        fn = HERE / f"extrinsic_{mat.lower()}.json"
        fn.write_text(json.dumps(out, indent=1, ensure_ascii=False))
        i = int(np.argmin(abs(m.T - 300)))
        print(f"{mat}: x = {c['x']} fixed, L = {out['L_nm']} nm, s = {out['s']}, rms = {rms0:.4f}; "
              f"acceptable L {out['acceptable']['L_nm']} nm, s {out['acceptable']['s']}; "
              f"300 K xx/zz = {xx[i]:.4f} / {zz[i]:.4f} -> {fn.name}")


if __name__ == "__main__":
    main()
