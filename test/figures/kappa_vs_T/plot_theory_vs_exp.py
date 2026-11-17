#!/usr/bin/env python3
"""Experimental kappa_tot (hot-pressed) with the calculated curves overlaid:
ordered (natural isotopes only) and VCA antisite disorder.

Inputs (all in this repository):
  experiment/kappa_tot_digitized.csv        digitized experimental points
  test/figures/kappa_vs_T/data_*_vca.json   20^3 phono3py RTA, baseline + VCA

Mapping, stated on the figure as well:
  Perp (perpendicular to the pressing direction) <-> calculated in-plane kappa
  Para (parallel to the pressing direction)      <-> calculated kappa_zz (rest of the zone,
                                                    i.e. Gamma-Z quadrature error removed)
This is the fully-textured limit. A real hot-pressed sample mixes the two
components through its texture, and the experiment also contains grain/defect
scattering and kappa_e, so the curves are a reference, not a fit.

Outputs:
  theory_vs_exp.png/.pdf      same two-panel layout and axis ranges as the
                              experimental figure (curves that leave the frame
                              are clipped, which is the point: ordered is far above)
  theory_vs_exp_log.png/.pdf  2 x 2 (material x direction), log scale,
                              30-700 K, everything visible
"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
EXP = ROOT / "experiment" / "kappa_tot_digitized.csv"

# colours follow the experimental figure
COL = {("Mn2In2Se5", "Perp"): "#3a3a3a", ("Mn2In2Se5", "Para"): "#e8282b",
       ("MnIn2Se4", "Perp"): "#9a9a9a", ("MnIn2Se4", "Para"): "#f08080"}
MARK = {"Perp": "h", "Para": "D"}
LABEL = {"Mn2In2Se5": r"In$_2$Mn$_2$Se$_5$", "MnIn2Se4": r"In$_2$MnSe$_4$"}
CALC_KEY = {"Perp": ("baseline_inplane", "vca_inplane"),
            "Para": ("baseline_kzz_rest", "vca_kzz_rest")}
ORDER = [("Mn2In2Se5", "Perp"), ("Mn2In2Se5", "Para"), ("MnIn2Se4", "Perp"), ("MnIn2Se4", "Para")]


def load_exp():
    rows = {}
    with open(EXP, encoding="utf-8") as f:
        lines = [l for l in f if not l.startswith("#")]
    for r in csv.DictReader(lines):
        key = (r["material"], r["direction"], r["panel"])
        rows.setdefault(key, []).append((float(r["T_K"]), float(r["kappa_tot_WmK"])))
    return {k: sorted(v) for k, v in rows.items()}


def load_calc():
    out = {}
    for m, fn in (("MnIn2Se4", "data_mnin2se4_vca.json"), ("Mn2In2Se5", "data_mn2in2se5_vca.json")):
        out[m] = json.loads((HERE / fn).read_text())
    return out


def exp_points(ax, exp, m, d, panel, ms=6, label=True, **kw):
    pts = exp.get((m, d, panel), [])
    if not pts:
        return
    T, k = zip(*pts)
    filled = panel == "highT"
    c = COL[(m, d)]
    ax.plot(T, k, MARK[d], ms=ms, color=c, mfc=c if filled else "none", mew=1.2,
            ls="none", label=f"{LABEL[m]}-{d}" if label else None, **kw)


def calc_lines(ax, calc, m, d, tmin, tmax, label=True):
    dat = calc[m]
    T = dat["T_K"]
    kb, kv = CALC_KEY[d]
    sel = [i for i, t in enumerate(T) if tmin <= t <= tmax]
    Ts = [T[i] for i in sel]
    c = COL[(m, d)]
    ax.plot(Ts, [dat[kb][i] for i in sel], "--", color=c, lw=1.5,
            label="calc., ordered" if label else None)
    ax.plot(Ts, [dat[kv][i] for i in sel], "-", color=c, lw=2.0,
            label="calc., VCA" if label else None)


def style(ax):
    ax.tick_params(direction="in", top=True, right=True, labelsize=11)
    for s in ax.spines.values():
        s.set_linewidth(1.2)


NOTE = ("calc.: phono3py RTA, 20$^3$, natural isotopes; Perp ↔ in-plane, Para ↔ κ$_{zz}$ (Γ–Z term removed); "
        "fully textured limit, no grain boundary, no κ$_e$. VCA: MnIn$_2$Se$_4$ i = 0.44, Mn$_2$In$_2$Se$_5$ x = 0.13")


def fig_same_axes(exp, calc):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.9))
    for m, d in ORDER:
        exp_points(a1, exp, m, d, "highT")
        calc_lines(a1, calc, m, d, 280, 720, label=False)
        exp_points(a2, exp, m, d, "lowT", ms=5, label=False)
        calc_lines(a2, calc, m, d, 0, 160, label=False)
    a1.set(xlim=(280, 650), ylim=(0, 1.2), xlabel="Temperature (K)",
           ylabel=r"$\kappa_{\rm tot}$ (W m$^{-1}$ K$^{-1}$)")
    a2.set(xlim=(0, 150), ylim=(0, 1.5), xlabel="Temperature (K)",
           ylabel=r"$\kappa_{\rm tot}$ (W m$^{-1}$ K$^{-1}$)")
    for a in (a1, a2):
        style(a)
        a.xaxis.label.set_size(13)
        a.yaxis.label.set_size(13)
    from matplotlib.lines import Line2D
    leg1 = a1.legend(fontsize=9, frameon=False, loc="upper right", title="exp.", title_fontsize=9)
    a1.add_artist(leg1)
    a1.legend(handles=[Line2D([], [], color="k", ls="--", lw=1.5, label="calc., ordered"),
                       Line2D([], [], color="k", ls="-", lw=2.0, label="calc., VCA")],
              fontsize=9, frameon=False, loc="lower left")
    a1.text(0.98, 0.02, "in-plane calc. (Perp) runs above the frame",
            transform=a1.transAxes, fontsize=8, color="#555", ha="right")
    a2.text(0.98, 0.02, "calc. only at T ≥ 30 K;\nmost calc. curves are above the frame here",
            transform=a2.transAxes, fontsize=8, color="#555", ha="right")
    fig.text(0.5, 0.005, NOTE, ha="center", fontsize=7.5, color="#444")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"theory_vs_exp.{ext}", dpi=300)
    plt.close(fig)


def fig_log(exp, calc):
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8), sharex=True)
    for (m, d), ax in zip([("MnIn2Se4", "Perp"), ("MnIn2Se4", "Para"),
                           ("Mn2In2Se5", "Perp"), ("Mn2In2Se5", "Para")], axes.ravel()):
        calc_lines(ax, calc, m, d, 0, 1000)
        exp_points(ax, exp, m, d, "highT")
        exp_points(ax, exp, m, d, "lowT", ms=5, label=False)
        ax.set_yscale("log")
        from matplotlib.ticker import ScalarFormatter
        for axis in (ax.yaxis,):
            axis.set_major_formatter(ScalarFormatter())
            axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlim(0, 720)
        comp = "in-plane" if d == "Perp" else r"$\kappa_{zz}$"
        ax.set_title(f"{LABEL[m]} — {d} (calc. {comp})", fontsize=11)
        ax.legend(fontsize=7.5, frameon=False, loc="upper right")
        style(ax)
    for ax in axes[1]:
        ax.set_xlabel("Temperature (K)", fontsize=12)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\kappa$ (W m$^{-1}$ K$^{-1}$)", fontsize=12)
    fig.text(0.5, 0.005, NOTE, ha="center", fontsize=7.5, color="#444")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"theory_vs_exp_log.{ext}", dpi=300)
    plt.close(fig)


def main():
    exp, calc = load_exp(), load_calc()
    fig_same_axes(exp, calc)
    fig_log(exp, calc)
    print("wrote theory_vs_exp.{png,pdf} and theory_vs_exp_log.{png,pdf}")


if __name__ == "__main__":
    main()
