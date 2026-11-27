#!/usr/bin/env python3
"""The three paper figures, from data already in the repository (no recomputation).

  fig1_kappa_vs_T     intrinsic kappa(T), both materials, in-plane and kappa_zz
  fig2_calc_vs_exp    hot-pressed experiment vs calculation, 2 x 2 (material x direction)
  fig3_cutoff         kappa vs FC3 cutoff, both materials (the remaining uncertainty)

One convention throughout, the same as REPORT §9.1:
  * kappa_zz is the mesh-converged rest-of-the-zone term (Gamma-Z quadrature
    error removed). The 30^3 totals are not plotted.
  * curves are RTA population, natural isotopes on.
  * the stars in fig1 are the values quoted in the text at 300 K: in-plane
    LBTE + coherence, kappa_zz converged + coherence (3.24.0 coherence for
    MnIn2Se4).

Data (all under test/figures/ unless noted):
  kappa_vs_T/data_mnin2se4.json   14^3, ordered, rest term per T
  kappa_vs_T/data_mnin2se5.json   16^3, ordered, rest term per T
  kappa_vs_T/data_*_vca.json      20^3, ordered + VCA pair, rest term per T
  calc_vs_experiment/extrinsic_*.json   L3 fit (x, L, s), texture-mixed Perp / Para
  cutoff/data_*.json              cutoff scans, population and rest term
  experiment/kappa_tot_digitized.csv   (repository root)
"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
ROOT = HERE.parents[2]

MAT = ("MnIn2Se4", "Mn2In2Se5")
LABEL = {"MnIn2Se4": r"In$_2$MnSe$_4$", "Mn2In2Se5": r"In$_2$Mn$_2$Se$_5$"}
COL = {"MnIn2Se4": "#d62728", "Mn2In2Se5": "#1f77b4"}       # analysis/plotting convention
QUOTED = {  # 300 K, REPORT §9.1
    "MnIn2Se4": {"inplane": 4.67, "kzz": 0.714},
    # in-plane 3.37 x R_xx(300 K) = 3.37 x 0.8571 = 2.89: spin-disordered FC2
    # (dft/spin_disorder_fc2/results/KAPPA_VERDICT.md); the RTA ratio is applied
    # to the LBTE + coherence value as an approximation
    "Mn2In2Se5": {"inplane": 2.89, "kzz": 1.089},
}
# R_xx(T) = kappa_xx(PM) / kappa_xx(FM_ctl) for Mn2In2Se5, 16x16x4, production fc3.
# Spin disorder lowers the in-plane single-crystal kappa; kappa_zz is NOT scaled
# (not resolved by that dataset). MnIn2Se4 is not scaled (not computed).
SPIN_R = json.loads((ROOT / "dft/spin_disorder_fc2/results/R_xx_vs_T.json").read_text())


def spin_factor(m, T):
    """R_xx(T) for Mn2In2Se5 (linear interpolation, 1 outside the range); 1 for MnIn2Se4."""
    if m != "Mn2In2Se5":
        return np.ones_like(np.asarray(T, float))
    return np.interp(np.asarray(T, float), SPIN_R["T_K"], SPIN_R["R_xx"])
CUTOFF_DOWN = {  # fractional downward FC3-cutoff uncertainty of the quoted values
    "MnIn2Se4": {"inplane": 0.10, "kzz": 0.10},    # estimate, REPORT §9.1
    "Mn2In2Se5": {"inplane": 0.021, "kzz": 0.033},  # next-step bound, cutoff_cv (2026-11-21)
}
# Mn2In2Se5 6.90 A refit (54 frames) relative to 6.51 A, same run (test/Mn2In2Se5/cutoff_ext);
# a lower bound for the converged 6.90 A value (cutoff_cv)
MN215_690_OVER_651 = {"inplane_population": 2.8731 / 2.9348, "kzz_rest": 0.9478 / 0.9801}


def load(rel):
    return json.loads((FIG / rel).read_text())


def style(ax):
    ax.grid(True, ls=":", alpha=0.4)
    ax.tick_params(direction="in", top=True, right=True)
    if ax.get_yscale() == "log":
        from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
        ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())


def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def exp_rows():
    out = {}
    with open(ROOT / "experiment" / "kappa_tot_digitized.csv", encoding="utf-8") as f:
        for r in csv.DictReader(l for l in f if not l.startswith("#")):
            out.setdefault((r["material"], r["direction"], r["panel"]), []).append(
                (float(r["T_K"]), float(r["kappa_tot_WmK"])))
    return {k: np.array(sorted(v)) for k, v in out.items()}


# --------------------------------------------------------------------------- fig 1
def fig1():
    data = {"MnIn2Se4": load("kappa_vs_T/data_mnin2se4.json"),
            "Mn2In2Se5": load("kappa_vs_T/data_mnin2se5.json")}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.2))
    for m in MAT:
        d = data[m]
        T = np.array(d["T_K"], float)
        keep = T >= 30
        # Mn2In2Se5 in-plane carries the spin-disorder factor R_xx(T); kappa_zz not
        a1.plot(T[keep], spin_factor(m, T[keep]) * np.array(d["inplane"])[keep], "-o",
                color=COL[m], ms=4, label=f"{LABEL[m]} ({d['mesh']})")
        a2.plot(T[keep], np.array(d["kzz_rest"])[keep], "-s", color=COL[m], ms=4, mfc="white",
                label=f"{LABEL[m]} ({d['mesh']})")
    # quoted 300 K values with the cutoff uncertainty
    for ax, key in ((a1, "inplane"), (a2, "kzz")):
        for m in MAT:
            v = QUOTED[m][key]
            x = 300 + (12 if m == "MnIn2Se4" else -12)
            ax.plot([x], [v], "*", color=COL[m], ms=13, mec="k", mew=0.6, zorder=5)
            # downward FC3-cutoff uncertainty: MnIn2Se4 ~10 % estimated from shrinking
            # steps; Mn2In2Se5 bound on the next step (6.90 A) from the hold-out test
            # (test/Mn2In2Se5/cutoff_cv/README.md)
            down = CUTOFF_DOWN[m][key]
            ax.errorbar([x], [v], yerr=[[down * v], [0]], color=COL[m], capsize=4, lw=1.6, zorder=4)
    for ax, t in ((a1, "(a) in-plane"), (a2, r"(b) cross-plane, $\kappa_{zz}$")):
        ax.set_yscale("log")
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(r"$\kappa_L$ (W m$^{-1}$ K$^{-1}$)")
        ax.set_title(t, fontsize=11)
        style(ax)
    a1.legend(handles=a1.get_legend_handles_labels()[0] + [
        Line2D([], [], ls="none", marker="*", color="0.4", mec="k", ms=11,
               label="300 K value quoted in the text")], fontsize=8, frameon=False)
    a2.legend(fontsize=8, frameon=False)
    fig.text(0.5, -0.04,
             "Lines: RTA population, natural isotopes, ordered crystal; κ$_{zz}$ = mesh-converged term "
             "(Γ–Z quadrature error removed). Stars: in-plane LBTE + coherence, κ$_{zz}$ + coherence.\n"
             "Bars: FC3-cutoff uncertainty, In$_2$MnSe$_4$ ≈ −10 % (estimate); In$_2$Mn$_2$Se$_5$ −2.1 % / −3.3 % "
             "(bound on the next cutoff step from a hold-out test).\n"
             "In$_2$Mn$_2$Se$_5$ in-plane is scaled by the spin-disorder ratio R$_{xx}$(T) = "
             "κ$_{xx}$(PM)/κ$_{xx}$(FM) = 0.822 (30 K) … 0.857 (300 K) … 0.858 (≥400 K), jackknife $\pm$0.018; "
             "κ$_{zz}$ is not scaled (that ratio is not resolved). In$_2$MnSe$_4$ is not scaled.",
             ha="center", fontsize=7, color="0.3")
    fig.tight_layout()
    save(fig, "fig1_kappa_vs_T")


# --------------------------------------------------------------------------- fig 2
def extrinsic(m):
    d = load(f"calc_vs_experiment/extrinsic_{m.lower()}.json")
    T = np.array([r["T"] for r in d["rows"]])
    kxx = np.array([r["kxx"] for r in d["rows"]])
    kzz = np.array([r["kzz"] for r in d["rows"]])
    s = d["s"]
    return T, {"Perp": (1 - s / 2) * kxx + (s / 2) * kzz, "Para": s * kxx + (1 - s) * kzz}, s


def fig2():
    exp = exp_rows()
    vca = {"MnIn2Se4": load("kappa_vs_T/data_mnin2se4_vca.json"),
           "Mn2In2Se5": load("kappa_vs_T/data_mn2in2se5_vca.json")}
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.6), sharex=True)
    for i, m in enumerate(MAT):
        d = vca[m]
        T = np.array(d["T_K"], float)
        Tl3, l3, s = extrinsic(m)
        for j, (dr, key, comp) in enumerate((("Perp", "inplane", "in-plane"),
                                             ("Para", "kzz_rest", r"$\kappa_{zz}$"))):
            ax = axes[i, j]
            c = COL[m]
            # in-plane carries the spin-disorder factor R_xx(T) for Mn2In2Se5,
            # applied to both the ordered and the VCA curve (assumes the spin and
            # antisite effects multiply); kappa_zz is not scaled
            f = spin_factor(m, T) if key == "inplane" else np.ones_like(T)
            ax.plot(T, f * np.array(d[f"baseline_{key}"]), "--", color=c, lw=1.4,
                    label="calc. ordered (single crystal)")
            ax.plot(T, f * np.array(d[f"vca_{key}"]), "-", color=c, lw=1.8,
                    label="calc. VCA (single crystal)")
            ax.plot(Tl3, l3[dr], ":", color="k", lw=1.6, label=f"calc. L3 fit (texture s = {s:.2f})")
            for panel, filled in (("highT", True), ("lowT", False)):
                p = exp.get((m, dr, panel))
                if p is not None:
                    ax.plot(p[:, 0], p[:, 1], "o" if dr == "Perp" else "s", ms=5, color=c,
                            mfc=c if filled else "white", mec=c, ls="none",
                            label="experiment (hot-pressed)" if filled else None)
            ax.set_yscale("log")
            ax.set_xlim(0, 720)
            ax.set_title(f"{LABEL[m]} — {dr} (calc. {comp})", fontsize=10.5)
            style(ax)
            if j == 0:
                ax.set_ylabel(r"$\kappa$ (W m$^{-1}$ K$^{-1}$)")
            if i == 1:
                ax.set_xlabel("Temperature (K)")
            ax.legend(fontsize=7.5, frameon=False, loc="upper right")
    vp = {m: vca[m]["vca"] for m in MAT}
    fig.text(0.5, -0.075,
             "Ordered / VCA: RTA, 20$^3$ mesh (same mesh for each pair), natural isotopes, κ$_{zz}$ = "
             "mesh-converged term; single-crystal components (Perp ↔ in-plane, Para ↔ κ$_{zz}$), "
             "no boundary or texture.\nAntisite occupancy (VCA and L3): In$_2$MnSe$_4$ i = %s, "
             "In$_2$Mn$_2$Se$_5$ x = %s, literature structure refinements (not of the measured samples).\n"
             "L3: ordered-mass κ + antisite mass disorder + boundary mean free path + texture; x fixed, "
             "L and s fitted to the 300–623 K points. Experiment: κ$_{tot}$, κ$_e$ not subtracted.\n"
             "In$_2$Mn$_2$Se$_5$ in-plane (ordered and VCA) is scaled by the spin-disorder ratio R$_{xx}$(T) = "
             "κ$_{xx}$(PM)/κ$_{xx}$(FM) = 0.857 $\pm$ 0.018 at 300 K, which assumes the spin and antisite "
             "effects multiply; κ$_{zz}$ is not scaled (unresolved). L3 still uses the FM FC2 — switching its "
             "in-plane basis to PM would lower the intrinsic in-plane κ by 14 %% and mainly change L$_{eff}$; "
             "MnIn$_2$Se$_4$ is not scaled (its spin-disorder effect was not computed)."
             % (vp["MnIn2Se4"]["x"], vp["Mn2In2Se5"]["x"]),
             ha="center", fontsize=7.2, color="0.3")
    fig.tight_layout()
    save(fig, "fig2_calc_vs_exp")


# --------------------------------------------------------------------------- fig 3
def fig3():
    cut = {"MnIn2Se4": load("cutoff/data_mnin2se4.json"), "Mn2In2Se5": load("cutoff/data_mnin2se5.json")}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.2))
    for m in MAT:
        d = cut[m]
        x = np.array(d["cutoff_A"])
        for ax, key, mk, mfc in ((a1, "inplane_population", "-o", COL[m]), (a2, "kzz_rest", "-s", "white")):
            y = np.array(d[key])
            ax.plot(x, y / y[-1], mk, color=COL[m], ms=6, mfc=mfc,
                    label=f"{LABEL[m]} ({d['mesh']}, last = {y[-1]:.3f})")
            up = m == "MnIn2Se4"
            for k in range(1, len(x)):
                ax.annotate(f"{100 * (y[k] / y[k - 1] - 1):+.1f} %",
                            ((x[k] + x[k - 1]) / 2, (y[k] + y[k - 1]) / 2 / y[-1]),
                            xytext=(6 if up else -6, 6 if up else -10), textcoords="offset points",
                            color=COL[m], fontsize=8, ha="left" if up else "right")
            if m == "Mn2In2Se5":
                # 6.90 A refit: over-parameterised, its kappa is a lower bound (cutoff_cv)
                r = MN215_690_OVER_651[key]
                ax.plot([6.90], [r], "^", color=COL[m], ms=7, mfc="white", ls="none")
                ax.plot([x[-1], 6.90], [1.0, r], ":", color=COL[m], lw=1)
                ax.annotate(f"≥ {100 * (r - 1):+.1f} %", (6.90, r), xytext=(4, -12),
                            textcoords="offset points", color=COL[m], fontsize=8)
    for ax, t in ((a1, "(a) in-plane, population"), (a2, r"(b) $\kappa_{zz}$, mesh-converged term")):
        ax.set_xlabel("FC3 cutoff (Å)")
        ax.set_ylabel("κ / κ(largest cutoff)")
        ax.set_title(t, fontsize=11)
        style(ax)
        ax.legend(fontsize=8, frameon=False)
    fig.text(0.5, -0.05, "300 K, RTA, natural isotopes. The largest cutoffs (6.45 / 6.51 Å) are the values "
             "used throughout; percentages are the step-to-step changes.\nOpen triangle: In$_2$Mn$_2$Se$_5$ "
             "refit at 6.90 Å on the same frames. A hold-out test shows it is over-parameterised, so it is a "
             "lower bound on the next step.", ha="center", fontsize=7.5, color="0.3")
    fig.tight_layout()
    save(fig, "fig3_cutoff")


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    print("wrote fig1_kappa_vs_T, fig2_calc_vs_exp, fig3_cutoff (.png/.pdf)")
