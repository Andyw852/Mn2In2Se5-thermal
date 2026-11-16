#!/usr/bin/env python3
"""Calculated kappa(T) overlaid on the experimental kappa_tot figure.

Style follows the AutoZT kl chain: this script imports
skill/kl-dft-cpu/kappa_plots.py and uses its RC block (Arial, 9 pt, inward ticks,
black thin frame, no grid, 300 dpi) and its direction palette (xx blue / zz
green). The experimental points are the digitised ones already in the repo
(experiment/kappa_tot_digitized.csv); Perp = in-plane (kappa_xx),
Para = out-of-plane (kappa_zz).

Three things are plotted:
  symbols  experiment (hot-pressed polycrystal)
  thin dashed  our intrinsic single-crystal calculation
  thick solid  our extrinsic model (antisite disorder + grain boundary + texture),
               only where a fit exists on the current force constants
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

KL = Path("/home/wangchao/software/AutoZT/skill/kl-dft-cpu")
sys.path.insert(0, str(KL))
import kappa_plots as kp            # noqa: E402

import matplotlib                   # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt      # noqa: E402

plt.rcParams.update(kp.RC)
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]

XX = kp.DIR_STYLE[0][1]             # blue   -> Perp / in-plane
ZZ = kp.DIR_STYLE[2][1]             # green  -> Para / out-of-plane
C_441 = kp.DIR_STYLE[1][1]          # orange -> the second material
INK2 = kp.INK2


def exp_rows():
    rows = []
    with open(REPO / "experiment/kappa_tot_digitized.csv") as f:
        for r in csv.DictReader(l for l in f if not l.startswith("#")):
            rows.append({"material": r["material"], "dir": r["direction"],
                         "T": float(r["T_K"]), "k": float(r["kappa_tot_WmK"]),
                         "panel": r["panel"]})
    return rows


def calc(name):
    return json.loads((REPO / "test/figures/kappa_vs_T" / name).read_text())


def model(fname, s):
    d = json.loads((HERE / fname).read_text())
    T = np.array([r["T"] for r in d["rows"]])
    kxx = np.array([r["kxx"] for r in d["rows"]])
    kzz = np.array([r["kzz"] for r in d["rows"]])
    return T, s * kxx + (1 - s) * kzz, (1 - s / 2) * kxx + (s / 2) * kzz


def panel(ax, exp, d441, d215, xlim, logy, show_model):
    ax.axvspan(*xlim, color="none")
    # --- experiment (colour = material, marker = direction)
    series = [("Mn2In2Se5", "Perp", XX, "o", "In$_2$Mn$_2$Se$_5$: exp, Perp = in-plane"),
              ("Mn2In2Se5", "Para", XX, "D", "In$_2$Mn$_2$Se$_5$: exp, Para = out-of-plane"),
              ("MnIn2Se4", "Perp", C_441, "o", "In$_2$MnSe$_4$: exp, Perp = in-plane"),
              ("MnIn2Se4", "Para", C_441, "D", "In$_2$MnSe$_4$: exp, Para = out-of-plane")]
    lo, hi = xlim
    for mat, dr, col, mk, lab in series:
        sel = [r for r in exp if r["material"] == mat and r["dir"] == dr and lo - 1e-6 <= r["T"] <= hi]
        if not sel:
            continue
        filled = mk == "o"
        ax.plot([r["T"] for r in sel], [r["k"] for r in sel], "none", marker=mk, ms=4.2,
                mfc=col if filled else "white", mec=col, mew=0.8, label=lab, zorder=4)
    # --- our intrinsic calculation (solid = in-plane, dash-dot = kappa_zz)
    for d, col in ((d215, XX), (d441, C_441)):
        T = np.array(d["T_K"])
        keep = (T >= lo) & (T <= hi)
        ax.plot(T[keep], np.array(d["inplane"])[keep], "-", color=col, lw=1.15, zorder=2)
        ax.plot(T[keep], np.array(d["kzz_rest"])[keep], "-.", color=col, lw=1.15, zorder=2)
    # --- extrinsic model
    if show_model:
        T, para, perp = model("extrinsic_mn2in2se5.json", 0.2533)
        keep = (T >= lo) & (T <= hi)
        ax.plot(T[keep], perp[keep], "-", color=XX, lw=2.4, alpha=0.55, zorder=1)
        ax.plot(T[keep], para[keep], "-.", color=XX, lw=2.4, alpha=0.55, zorder=1)
    ax.set_xlim(*xlim)
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$\kappa_{\rm tot}$ (W m$^{-1}$ K$^{-1}$)")


def main():
    exp = exp_rows()
    d441, d215 = calc("data_mnin2se4.json"), calc("data_mnin2se5.json")
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    panel(a, exp, d441, d215, (290, 700), False, True)
    panel(b, exp, d441, d215, (20, 150), True, True)
    a.set_title("(a) 300-700 K", fontsize=9, loc="left")
    b.set_title("(b) 20-150 K (log scale)", fontsize=9, loc="left")
    from matplotlib.lines import Line2D
    handles, labels = a.get_legend_handles_labels()
    handles += [Line2D([], [], color=INK2, lw=1.15, ls="-"),
                Line2D([], [], color=INK2, lw=1.15, ls="-."),
                Line2D([], [], color=INK2, lw=2.4, alpha=0.55, ls="-")]
    labels += ["calc intrinsic, in-plane", "calc intrinsic, $\\kappa_{zz}$",
               "calc extrinsic model, Mn$_2$In$_2$Se$_5$ (both directions)"]
    a.legend(handles, labels, fontsize=6.0, loc="upper left", ncol=1,
             handletextpad=0.3, borderpad=0.2, labelspacing=0.25)
    fig.text(0.01, -0.06,
             "Symbols: experiment (hot-pressed polycrystal, digitised from "
             "experiment/kappa_tot_figure.png). Solid lines: our intrinsic single-crystal RTA, "
             "in-plane; dash-dot: the same for $\\kappa_{zz}$ (MnIn$_2$Se$_4$ 14$^3$/6.45 A, "
             "Mn$_2$In$_2$Se$_5$ 16$^3$/6.51 A; neither FC3 cutoff converged).\n"
             "Thick translucent: extrinsic model for Mn$_2$In$_2$Se$_5$ (x = 0.13 Mn/In antisite "
             "disorder, L = 16 nm grain boundary, texture s = 0.25) -- that is what turns the "
             "intrinsic single-crystal curves into the measured polycrystal values. "
             "The MnIn$_2$Se$_4$ extrinsic fit still rests on the superseded 4x4x1 FC2, so only its "
             "intrinsic curve is shown.",
             fontsize=6, color=INK2, transform=fig.transFigure)
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    for ext in ("png", "pdf"):
        fig.savefig(str(HERE / ("calc_vs_experiment." + ext)), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote calc_vs_experiment")


if __name__ == "__main__":
    main()
