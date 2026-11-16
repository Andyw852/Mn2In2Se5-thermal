#!/usr/bin/env python3
"""kappa(T) for MnIn2Se4 and Mn2In2Se5 from the committed RTA runs.

One figure per material (in-plane and kappa_zz) plus a combined comparison.
kappa_zz plots the REST-OF-THE-ZONE term at every temperature (the raw total is
contaminated by the first Gamma-Z q-point), while the total is kept in the data
files for reference.

The sub-150 K region is shaded: on 14^3 / 16^3 meshes the long-mean-free-path
phonons are not sampled, so those points are lower bounds on the intrinsic kappa
and are NOT mesh-converged. The 300 K column is the converged value.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
INK2 = "#666666"
C_441, C_215 = "#1f4e79", "#c00000"
MESH_LIMITED_BELOW = 150.0


def load(n):
    return json.loads((HERE / n).read_text())


def one(d, color, stem):
    T = d["T_K"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 3.6))
    for ax, vals, title, ylab in (
            (a1, d["inplane"], "(a) in-plane", r"$\kappa_{in}$ (W m$^{-1}$ K$^{-1}$)"),
            (a2, d["kzz_rest"], "(b) cross-plane (rest of the zone)",
             r"$\kappa_{zz}$ (W m$^{-1}$ K$^{-1}$)")):
        ax.axvspan(0, MESH_LIMITED_BELOW, color=INK2, alpha=0.08, lw=0)
        ax.plot(T, vals, "-o", color=color, lw=1.5, ms=4.5, zorder=3)
        i300 = T.index(300)
        ax.plot([300], [vals[i300]], "o", ms=9, mfc="none", mec="#2e7d32", mew=1.6, zorder=4)
        ax.annotate("%.2f" % vals[i300], (300, vals[i300]), xytext=(6, 6),
                    textcoords="offset points", fontsize=7, color="#2e7d32")
        ax.annotate("mesh-limited\n(long MFP)", (75, ax.get_ylim()[1] if False else 0.92),
                    xycoords=("data", "axes fraction"), fontsize=6.5, color=INK2, ha="center")
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=9, loc="left")
        ax.set_xlim(0, 720)
        ax.grid(alpha=0.18, lw=0.5)
        ax.set_yscale("log")
    fig.suptitle("%s  %s, FC3 %s Å, isotopes %s" %
                 (d["material"], d["mesh"], d["fc3_cutoff_A"],
                  "on" if d["isotopes"] else "off"), fontsize=9.5, x=0.01, ha="left")
    fig.text(0.01, -0.04,
             "Shaded: below 150 K the phonon mean free paths exceed the mesh sampling, so those points are "
             "lower bounds, not converged values.\n"
             "Green circle = 300 K, the value quoted in the reports. kappa_zz is the rest-of-the-zone term; "
             "the raw total (contaminated) is in the data JSON.",
             fontsize=6, color=INK2, transform=fig.transFigure)
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    for ext in ("png", "pdf"):
        fig.savefig(str(HERE / (stem + "." + ext)), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", stem)


def both(d1, d2):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    for ax, key, title, ylab in ((a1, "inplane", "(a) in-plane",
                                  r"$\kappa_{in}$ (W m$^{-1}$ K$^{-1}$)"),
                                 (a2, "kzz_rest", "(b) cross-plane (rest of the zone)",
                                  r"$\kappa_{zz}$ (W m$^{-1}$ K$^{-1}$)")):
        ax.axvspan(0, MESH_LIMITED_BELOW, color=INK2, alpha=0.08, lw=0)
        for d, col in ((d1, C_441), (d2, C_215)):
            ax.plot(d["T_K"], d[key], "-o", color=col, lw=1.5, ms=4,
                    label=d["material"])
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=9, loc="left")
        ax.set_xlim(0, 720)
        ax.set_yscale("log")
        ax.grid(alpha=0.18, lw=0.5)
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("MnIn2Se4 vs Mn2In2Se5, 300 K region and its temperature dependence",
                 fontsize=9.5, x=0.01, ha="left")
    fig.text(0.01, -0.05,
             "Shaded below 150 K: mesh-limited. MnIn2Se4 14^3 / 6.45 A / isotopes on; "
             "Mn2In2Se5 16^3 / 6.51 A / isotopes off. Neither FC3 cutoff is converged.",
             fontsize=6, color=INK2, transform=fig.transFigure)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    for ext in ("png", "pdf"):
        fig.savefig(str(HERE / ("both_kappa_vs_T." + ext)), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote both_kappa_vs_T")


if __name__ == "__main__":
    d1, d2 = load("data_mnin2se4.json"), load("data_mnin2se5.json")
    one(d1, C_441, "MnIn2Se4_kappa_vs_T")
    one(d2, C_215, "Mn2In2Se5_kappa_vs_T")
    both(d1, d2)
