#!/usr/bin/env python3
"""Compare the intrinsic RTA kappa of the FC2 sets (bases.npz written by analysis/extrinsic/kappa_extrinsic.py).

  python kappa_summary.py --production <repo>/analysis/extrinsic/results/Mn2In2Se5_m16164 \
         --fm-ctl B/FM_ctl --pm B/PM --afm B/AFM --jackknife B/jk_PM01 ... B/jk_PM08 --out kappa_summary.json

All sets must use the same fc3 (production), q mesh and temperatures. Reported per T: kappa_xx and
kappa_zz of each set, R = kappa(PM)/kappa(FM_ctl) and R_AFM, the jackknife error of R, and the method
check kappa(FM_ctl)/kappa(production). Decision rule (fixed before the runs, see README):
  |R - 1| <= 5 % at 300 K in both directions  ->  "small": the mean (spin-averaged) harmonic effect is
                                                  within the other uncertainties; keep the FM-based L1
  otherwise                                   ->  "rebase": L1 and the L3 fit move to the PM FC2
"""
import argparse
import json
from pathlib import Path

import numpy as np

T_REPORT = [100, 200, 300, 400, 500, 600]


def load(d):
    b = np.load(Path(d) / "bases.npz")
    meta = json.loads((Path(d) / "bases_meta.json").read_text())
    return b["T"], b["kappa0"][:, [0, 2]], meta     # xx, zz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--production", required=True)
    ap.add_argument("--fm-ctl", required=True)
    ap.add_argument("--pm", required=True)
    ap.add_argument("--afm", required=True)
    ap.add_argument("--jackknife", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sets = {"production": a.production, "FM_ctl": a.fm_ctl, "PM": a.pm, "AFM": a.afm}
    K, T0, metas = {}, None, {}
    for k, d in sets.items():
        T, kap, meta = load(d)
        if T0 is None:
            T0, mesh0, fc3 = T, meta["mesh"], Path(meta["fc3"]).name
        assert np.allclose(T, T0) and meta["mesh"] == mesh0, f"{k}: T or mesh differ"
        if Path(meta["fc3"]).name != fc3:
            print(f"WARNING {k}: fc3 file {meta['fc3']}")
        K[k], metas[k] = kap, meta
    jk = np.array([load(d)[1] for d in a.jackknife])
    n = len(jk)
    R = K["PM"] / K["FM_ctl"]
    Rjk = jk / K["FM_ctl"]
    errR = np.sqrt((n - 1) / n * np.sum((Rjk - Rjk.mean(axis=0)) ** 2, axis=0))
    rows = []
    print(f"intrinsic RTA kappa (W/mK), mesh {mesh0}; R = PM/FM_ctl with jackknife error ({n} sets)")
    print("   T    prod xx/zz      FM_ctl xx/zz     PM xx/zz        AFM xx/zz      R_xx            R_zz            R_AFM xx/zz   FM_ctl/prod xx/zz")
    for t in T_REPORT:
        i = int(np.argmin(np.abs(T0 - t)))
        r = {"T": float(T0[i])}
        for k in K:
            r[k] = K[k][i].tolist()
        r["R_PM"] = R[i].tolist()
        r["R_PM_err"] = errR[i].tolist()
        r["R_AFM"] = (K["AFM"][i] / K["FM_ctl"][i]).tolist()
        r["FM_ctl_over_production"] = (K["FM_ctl"][i] / K["production"][i]).tolist()
        rows.append(r)
        print(f"  {T0[i]:4.0f}  " + "  ".join(f"{K[k][i][0]:6.3f}/{K[k][i][1]:5.3f}" for k in K)
              + f"  {R[i][0]:.3f}±{errR[i][0]:.3f}  {R[i][1]:.3f}±{errR[i][1]:.3f}"
              + f"  {r['R_AFM'][0]:.3f}/{r['R_AFM'][1]:.3f}  {r['FM_ctl_over_production'][0]:.3f}/{r['FM_ctl_over_production'][1]:.3f}")
    r300 = next(r for r in rows if abs(r["T"] - 300) < 1)
    dec = "small" if max(abs(x - 1) for x in r300["R_PM"]) <= 0.05 else "rebase"
    print(f"\nDECISION (300 K, |R-1| <= 5 % in xx and zz): {dec}")
    Path(a.out).write_text(json.dumps({"rows": rows, "decision_300K": dec, "mesh": mesh0, "fc3": fc3,
                                       "n_jackknife": n, "inputs": {k: str(v) for k, v in sets.items()}}, indent=1))
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
