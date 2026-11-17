#!/usr/bin/env python3
"""kappa(T) with and without VCA antisite disorder, same mesh for each pair.

Reads data_<material>_vca.json (20^3 baseline + VCA, from the multi-temperature
VCA jobs 424495 / 424496). Solid = ordered baseline (natural isotopes),
dashed = VCA (averaged site masses + mass-variance scattering).

The point of the figure: the disorder does not just scale kappa down uniformly --
it changes the shape of kappa(T), because the mass-variance rate (~g*omega^2) wins
over the anharmonic rate at low temperature.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
C = {"MnIn2Se4": "#d62728", "Mn2In2Se5": "#1f77b4"}
INK2 = "#666666"


def load(name):
    p = HERE / name
    return json.loads(p.read_text()) if p.exists() else None


def main():
    sets = [d for d in (load("data_mnin2se4_vca.json"), load("data_mn2in2se5_vca.json")) if d]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.4))
    for d in sets:
        T = d["T_K"]
        col = C.get(d["material"], "#333333")
        a1.plot(T, d["baseline_inplane"], "-o", color=col, ms=4, mfc=col,
                label="%s, ordered" % d["material"])
        a1.plot(T, d["vca_inplane"], "--s", color=col, ms=4, mfc="white",
                label="%s, VCA" % d["material"])
        a2.plot(T, d["baseline_kzz_rest"], "-o", color=col, ms=4, mfc=col)
        a2.plot(T, d["vca_kzz_rest"], "--s", color=col, ms=4, mfc="white")
    for ax, ylab, title in ((a1, r"$\kappa_{in}$ (W m$^{-1}$ K$^{-1}$)", "(a) in-plane"),
                            (a2, r"$\kappa_{zz}$ rest (W m$^{-1}$ K$^{-1}$)",
                             "(b) cross-plane, rest of the zone")):
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_yscale("log")
        ax.set_xlim(0, 720)
        ax.grid(True, ls=":", alpha=0.4)
    a1.legend(fontsize=7, frameon=False)
    fig.text(0.01, -0.05,
             "Solid symbols: ordered crystal with natural isotopes. Open symbols, dashed: VCA "
             "(site masses averaged over the Mn/In antisite occupancy, disorder returned as "
             "mass-variance scattering).\n"
             "Each pair is computed on the same 20^3 mesh, largest FC3 cutoff (In$_2$MnSe$_4$ "
             "6.45 A, In$_2$Mn$_2$Se$_5$ 6.51 A), so only the disorder differs between the two "
             "curves of a pair. x = 0.44 / 0.13 respectively.",
             fontsize=6, color=INK2, transform=fig.transFigure)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    for ext in ("png", "pdf"):
        fig.savefig(str(HERE / ("vca_vs_T." + ext)), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote vca_vs_T")


if __name__ == "__main__":
    main()
