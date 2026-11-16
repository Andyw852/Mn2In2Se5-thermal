#!/usr/bin/env python3
"""kappa vs FC3 cutoff for MnIn2Se4 and Mn2In2Se5, with the Wigner decomposition.

Data files: data_mnin2se4.json / data_mnin2se5.json (provenance in README.md and
in each file's "note"). One figure per material plus a normalized comparison.

Conventions follow the kl chain's cutoff figure (skill/_common/fcfit/kappa_driver.py
-> kappa_vs_cutoff.png): third-order cutoff on the x axis, kappa in W/m/K, 300 dpi,
png + pdf, step percentages annotated between points.

The kappa_zz column plots the REST-OF-THE-ZONE term, not the total: the total is
contaminated by the single Gamma-Z q-point (see the decomposition READMEs).
The coherence term for kappa_zz was measured on that contaminated total, hence it
is drawn in grey and must not be added to the rest term.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

HERE = Path(__file__).resolve().parent
INK, INK2 = "#222222", "#666666"
C_POP, C_COH, C_TOT = "#1f4e79", "#c00000", "#2e7d32"


def load(name):
    return json.loads((HERE / name).read_text())


def step_pct(v):
    return [100.0 * (v[i + 1] / v[i] - 1.0) for i in range(len(v) - 1)]


def panel(ax, xs, series, ylabel, title, ceils):
    """series: list of (values, label, color, marker, ls, lw, filled, val_dy)"""
    lo, hi = None, None
    for vals, lab, col, mk, ls, lw, filled, dy in series:
        if any(v is None for v in vals):
            continue
        ax.plot(xs, vals, linestyle=ls, color=col, marker=mk, lw=lw, ms=6,
                markerfacecolor=col if filled else "white", markeredgecolor=col,
                label=lab, zorder=3)
        for i in (0, len(xs) - 1):        # only the endpoints carry their value
            x, v = xs[i], vals[i]
            ha = "left" if i == 0 else "right"
            ax.annotate("%.2f" % v, (x, v), textcoords="offset points",
                        xytext=(-8 if i == 0 else 8, dy), ha=ha, fontsize=7, color=col)
        lo = min(vals) if lo is None else min(lo, min(vals))
        hi = max(vals) if hi is None else max(hi, max(vals))
        d = step_pct(vals)
        for i, p in enumerate(d):
            ax.annotate("%+.1f%%" % p, ((xs[i] + xs[i + 1]) / 2,
                        (vals[i] + vals[i + 1]) / 2), textcoords="offset points",
                        xytext=(0, -15), ha="center", fontsize=6.5, color=INK2)
    ax.set_xlabel(r"Third-order cutoff $r_c$ (Å)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9, loc="left")
    ax.set_xticks(xs)
    ax.set_xticklabels(["%.2f" % x for x in xs])
    if lo is not None:
        span = hi - lo if hi > lo else max(hi, 1.0)
        ax.set_ylim(lo - 0.22 * span, hi + 0.18 * span)
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.grid(alpha=0.18, lw=0.5)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")


def one(d, outstem):
    xs = d["cutoff_A"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 3.6))
    panel(a1, xs, [
        (d["inplane_population"], "population (BTE)", C_POP, "s", "-", 1.4, True, -14),
        (d["inplane_coherence"], "coherence (Wigner)", C_COH, "^", (0, (5, 2)), 1.2, False, 10),
        (d["inplane_total"], "total (P + C)", C_TOT, "o", "-", 1.8, True, 10),
    ], r"$\kappa_{in}$ (W m$^{-1}$ K$^{-1}$)", "(a) in-plane", xs)
    panel(a2, xs, [
        (d["kzz_rest"], "population, rest of the zone", C_POP, "s", "-", 1.4, True, 10),
        (d["kzz_coherence"], "coherence (on the contaminated total)", "#999999", "^",
         (0, (5, 2)), 1.0, False, -16),
    ], r"$\kappa_{zz}$ (W m$^{-1}$ K$^{-1}$)", "(b) cross-plane", xs)
    fig.suptitle("%s  %s, %d K, isotopes %s" % (d["material"], d["mesh"], d["T_K"],
                 "on" if d["isotope"] else "off"), fontsize=9.5, x=0.01, ha="left")
    fig.text(0.01, -0.03,
             "Rightmost point = the largest cutoff the 4x4x1 cell allows (6.45 A); the curve "
             "is still falling there, so neither component is cutoff-converged.\n"
             "population/rest: point-by-point RTA; coherence: Wigner-RTA (SMM19). The kappa_zz "
             "coherence is grey because it was measured on the contaminated total and must not "
             "be added to the rest term.",
             fontsize=6, color=INK2, transform=fig.transFigure)
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    for ext in ("png", "pdf"):
        fig.savefig(str(HERE / (outstem + "." + ext)), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", outstem)


def normalized(ds, outstem):
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    styles = {"in-plane": ("o", "-", "#1f4e79"), "kappa_zz (rest)": ("s", (0, (5, 2)), "#c00000")}
    for d in ds:
        xs = d["cutoff_A"]
        for key, lab, vals in (("in-plane", "in-plane", d["inplane_population"]),
                               ("kappa_zz (rest)", "kappa_zz (rest)", d["kzz_rest"])):
            mk, ls, col = styles[key]
            ax.plot(xs, [100.0 * v / vals[-1] for v in vals], linestyle=ls, color=col, marker=mk,
                    lw=1.2, ms=5, markerfacecolor=col if mk == "o" else "white",
                    markeredgecolor=col,
                    label="%s, %s" % (d["material"], lab))
    ax.axhline(100.0, color=INK2, lw=0.7, ls=(0, (4, 2)))
    ax.set_xlabel(r"Third-order cutoff $r_c$ (Å)")
    ax.set_ylabel("κ / κ(largest available cutoff)  (%)")
    ax.set_title("FC3-cutoff sensitivity, 300 K", fontsize=9, loc="left")
    ax.grid(alpha=0.18, lw=0.5)
    ax.legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(str(HERE / (outstem + "." + ext)), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", outstem)


if __name__ == "__main__":
    d1, d2 = load("data_mnin2se4.json"), load("data_mnin2se5.json")
    one(d1, "MnIn2Se4_kappa_vs_cutoff")
    one(d2, "Mn2In2Se5_kappa_vs_cutoff")
    normalized([d1, d2], "both_normalized_vs_cutoff")
