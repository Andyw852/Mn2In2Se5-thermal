#!/usr/bin/env python3
"""Calculated kappa(T) over the experimental kappa_tot figure.

Drawn in the style of the project's own plotting programs
(analysis/plotting/plot_kappa_compare.py): figsize (10, 4.2), two panels,
dotted grid at low alpha, "-o" for in-plane and "-s" with mfc="white" for
out-of-plane, frameon=False legend, log y where the range needs it.

Data
  experiment  experiment/kappa_tot_digitized.csv  (Perp = in-plane, Para = out-of-plane)
  calc        ../kappa_vs_T/data_*.json            (intrinsic single crystal)
  model       extrinsic_mn2in2se5.json             (x = 0.13, L = 16 nm, s = 0.25)
"""
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]

C_215 = "#1f77b4"        # In2Mn2Se5
C_441 = "#d62728"        # In2MnSe4
C_EXP_EDGE = "#333333"


def exp_rows():
    rows = []
    with open(REPO / "experiment/kappa_tot_digitized.csv") as f:
        for r in csv.DictReader(l for l in f if not l.startswith("#")):
            rows.append((r["material"], r["direction"], float(r["T_K"]),
                         float(r["kappa_tot_WmK"])))
    return rows


def calc(name):
    return json.loads((REPO / "test/figures/kappa_vs_T" / name).read_text())


def extrinsic(name):
    """extrinsic model curve: returns T, kappa_Perp, kappa_Para"""
    d = json.loads((HERE / name).read_text())
    T = np.array([r["T"] for r in d["rows"]])
    kxx = np.array([r["kxx"] for r in d["rows"]])
    kzz = np.array([r["kzz"] for r in d["rows"]])
    s = d["s"]
    return T, (1 - s / 2) * kxx + (s / 2) * kzz, s * kxx + (1 - s) * kzz


def draw(ax, lo, hi):
    rows = exp_rows()
    # experiment
    for mat, col in (("Mn2In2Se5", C_215), ("MnIn2Se4", C_441)):
        for dr, mk, filled in (("Perp", "o", True), ("Para", "s", False)):
            sel = [(t, k) for m, d, t, k in rows if m == mat and d == dr and lo <= t <= hi]
            if not sel:
                continue
            ax.plot(*zip(*sel), "none", marker=mk, ms=6.0,
                    mfc=col if filled else "white", mec=col, mew=1.1, zorder=4)
    # intrinsic calculation
    # 30^3 ordered totals: the convention used by the figures from 2026-11-17
    # (kappa_zz = plain mesh total, NOT the rest-of-zone term)
    for d, col in ((calc("data_mn2in2se5_m30.json"), C_215),
                   (calc("data_mnin2se4_m30.json"), C_441)):
        T = np.array(d["T_K"])
        keep = (T >= lo) & (T <= hi)
        ax.plot(T[keep], np.array(d["inplane"])[keep], "-", color=col, lw=1.5,
                marker="o", ms=3.0, mfc=col, zorder=2)
        ax.plot(T[keep], np.array(d["kzz_rest"])[keep], "-", color=col, lw=1.5,
                marker="s", ms=3.0, mfc="white", zorder=2)
    # extrinsic model, both materials
    for name, col in (("extrinsic_mn2in2se5.json", C_215),
                      ("extrinsic_mnin2se4.json", C_441)):
        T, perp, para = extrinsic(name)
        keep = (T >= lo) & (T <= hi)
        ax.plot(T[keep], perp[keep], "--", color=col, lw=1.6, alpha=0.9, zorder=3)
        ax.plot(T[keep], para[keep], "--", color=col, lw=1.6, alpha=0.9, zorder=3)
    ax.set_xlim(lo, hi)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$\kappa_{\rm tot}$ (W m$^{-1}$ K$^{-1}$)")
    ax.grid(True, ls=":", alpha=0.4)


def main():
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.2))
    draw(a, 290, 700)
    draw(b, 0, 150)
    a.set_yscale("log")
    b.set_yscale("log")
    a.set_title("(a) 300-700 K", fontsize=10, loc="left")
    b.set_title("(b) 0-150 K", fontsize=10, loc="left")

    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", ls="none", mfc=C_215, mec=C_215, ms=6,
               label="In$_2$Mn$_2$Se$_5$ exp, Perp (in-plane)"),
        Line2D([], [], marker="s", ls="none", mfc="white", mec=C_215, ms=6,
               label="In$_2$Mn$_2$Se$_5$ exp, Para (out-of-plane)"),
        Line2D([], [], marker="o", ls="none", mfc=C_441, mec=C_441, ms=6,
               label="In$_2$MnSe$_4$ exp, Perp (in-plane)"),
        Line2D([], [], marker="s", ls="none", mfc="white", mec=C_441, ms=6,
               label="In$_2$MnSe$_4$ exp, Para (out-of-plane)"),
        Line2D([], [], color="#555555", lw=1.5, marker="o", ms=3,
               label="calc intrinsic, in-plane"),
        Line2D([], [], color="#555555", lw=1.5, marker="s", ms=3, mfc="white",
               label="calc intrinsic, $\\kappa_{zz}$"),
        Line2D([], [], color=C_215, lw=1.6, ls="--",
               label="calc extrinsic model (both materials)"),
    ]
    fig.legend(handles=handles, fontsize=7, frameon=False, ncol=4,
               loc="lower center", bbox_to_anchor=(0.5, -0.06), handletextpad=0.4,
               columnspacing=1.2)
    fig.text(0.01, -0.15,
             "Symbols: hot-pressed polycrystal experiment (digitised).\n"
             "Solid: ordered single crystal, phono3py RTA on a 30$^3$ mesh, natural isotopes; "
             "$\\kappa_{zz}$ is the plain mesh total (no $\\Gamma$–Z decomposition, so it may still "
             "be ~15 % high). FC2 4x4x2 (In$_2$MnSe$_4$) / 4x4x1 (In$_2$Mn$_2$Se$_5$).\n"
             "Dashed: extrinsic model (Mn/In antisite mass disorder + grain boundary + texture), "
             "fitted to the same CSV: In$_2$Mn$_2$Se$_5$ x = 0.13, L = 16 nm, s = 0.25; "
             "In$_2$MnSe$_4$ x = 0.45, L = 31 nm, s = 0.35.",
             fontsize=6, color="#555555", transform=fig.transFigure, wrap=True)
    fig.tight_layout(rect=(0, 0.14, 1, 0.97))
    for ext in ("png", "pdf"):
        fig.savefig(str(HERE / ("calc_vs_experiment." + ext)), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote calc_vs_experiment")


if __name__ == "__main__":
    main()
