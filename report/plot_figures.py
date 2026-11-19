#!/usr/bin/env python3
"""Draw the four report figures from the CSV / JSON files in a plot_data folder.

Needs only numpy + matplotlib. Works inside the repository and in the exported
bundle:

  python3 plot_figures.py                      # data: ./plot_data or ../data, output: ./figures
  python3 plot_figures.py --data DIR --out DIR [--lang en|zh|both]

Writes fig1_phonon_dispersion, fig2_kappa_intrinsic, fig3_calc_vs_experiment and
fig4_scattering_300K, each as .png (300 dpi) and .pdf; the Chinese-label versions end in _zh.
"""
import argparse
import csv
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

MAT = ("MnIn2Se4", "Mn2In2Se5")
LABEL = {"MnIn2Se4": r"MnIn$_2$Se$_4$", "Mn2In2Se5": r"Mn$_2$In$_2$Se$_5$"}
COL = {"MnIn2Se4": "#d62728", "Mn2In2Se5": "#1f77b4"}      # analysis/plotting convention
LANG = "en"
ZH = {
    "branches": "支", "max": "最高", "Frequency (THz)": "频率 (THz)", "300 K value to quote": "300 K 推荐值",
    "(a) in-plane": "(a) 面内", "(b) cross-plane": "(b) 层间（沿 c 轴）", "Temperature (K)": "温度 (K)",
    "calc. ordered single crystal": "计算：有序单晶", "calc. VCA single crystal": "计算：VCA 单晶",
    "calc. L3 fit": "计算：L3 拟合", "experiment 300–623 K (fitted)": "实验 300–623 K（参与拟合）",
    "experiment 55–145 K (not fitted)": "实验 55–145 K（未参与拟合）", "calc.": "计算", "in-plane": "面内",
    "baseline 3-phonon": "基线 三声子", "VCA 3-phonon": "VCA 三声子", "3-phonon scattering": "三声子散射",
    "VCA disorder (mass variance)": "VCA 无序（质量方差）", "baseline natural isotope": "基线 天然同位素",
    "disorder elastic scattering": "无序弹性散射", "diffuson κ": "diffuson κ",
}
CJK_FONTS = ["Microsoft YaHei", "SimHei", "PingFang SC", "Heiti SC", "Noto Sans CJK SC", "Source Han Sans SC",
             "WenQuanYi Zen Hei"]


UNI = {"MnIn2Se4": "MnIn₂Se₄", "Mn2In2Se5": "Mn₂In₂Se₅"}


def lab(m):
    """Material label: mathtext in English, Unicode subscripts in Chinese (mixed with CJK text)."""
    return UNI[m] if LANG == "zh" else LABEL[m]


def t(text):
    return ZH.get(text, text) if LANG == "zh" else text


def set_lang(lang):
    """'en' or 'zh'. For 'zh' a CJK font is added as fallback (the first one installed is used)."""
    global LANG
    LANG = lang
    if lang == "zh":
        # DejaVu Sans first keeps digits and the minus sign; Chinese glyphs fall back to the first
        # installed CJK font (per-glyph fallback over the font.family list, matplotlib >= 3.6)
        from matplotlib import font_manager
        have = {f.name for f in font_manager.fontManager.ttflist}
        cjk = [f for f in CJK_FONTS if f in have]
        if not cjk:
            print("warning: no CJK font found; install one of", CJK_FONTS)
        matplotlib.rcParams["font.family"] = ["DejaVu Sans"] + cjk[:1]
        matplotlib.rcParams["axes.unicode_minus"] = False
        # mathtext has no font fallback, so Chinese labels never contain $math$ (see lab() and ZH)


STEMS = {"bands": "fig1_phonon_dispersion", "intrinsic": "fig2_kappa_intrinsic",
         "calc_vs_exp": "fig3_calc_vs_experiment", "scattering": "fig4_scattering_300K"}


def read_csv(path):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.reader(l for l in f if not l.startswith("#")))
    return rows[0], rows[1:]


def columns(path):
    head, rows = read_csv(path)
    out = {}
    for j, h in enumerate(head):
        col = [r[j] for r in rows]
        try:
            out[h] = np.array(col, float)
        except ValueError:
            out[h] = np.array(col)
    return out


def style(ax):
    ax.grid(True, ls=":", alpha=0.4)
    ax.tick_params(direction="in", top=True, right=True)


# --------------------------------------------------------------------------- figures
def fig_bands(D):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, m in zip(axes, MAT):
        b = columns(D / f"bands_{m}.csv")
        tk = columns(D / f"bands_{m}_ticks.csv")
        bands = [k for k in b if k.startswith("band")]
        for s in np.unique(b["segment"]):
            k = b["segment"] == s
            for name in bands:
                ax.plot(b["distance"][k], b[name][k], color=COL[m], lw=1.0)
        for x in tk["distance"]:
            ax.axvline(x, color="0.7", lw=0.7)
        ax.set_xticks(tk["distance"], tk["label"])
        ax.set_xlim(tk["distance"][0], tk["distance"][-1])
        ax.set_ylim(0, 8.2)
        fmax = max(b[n].max() for n in bands)
        ax.set_title(f"{lab(m)}  ({len(bands)} {t('branches')}, {t('max')} {fmax:.2f} THz)", fontsize=11)
        ax.tick_params(direction="in", right=True)
        ax.grid(True, axis="y", ls=":", alpha=0.4)
    axes[0].set_ylabel(t("Frequency (THz)"))
    fig.tight_layout()
    return fig


def fig_intrinsic(D):
    q = {}
    for r in read_csv(D / "kappa_quoted_300K.csv")[1]:
        q[(r[0], r[1])] = float(r[2])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for m in MAT:
        d = columns(D / f"kappa_intrinsic_{m}.csv")
        k = d["T_K"] >= 30
        a1.plot(d["T_K"][k], d["kappa_inplane"][k], "-o", color=COL[m], ms=4, label=lab(m))
        a2.plot(d["T_K"][k], d["kappa_zz"][k], "-s", color=COL[m], ms=4, mfc="white", label=lab(m))
    for ax, key in ((a1, "inplane"), (a2, "kzz")):
        for m in MAT:
            v = q[(m, key)]
            x = 300 + (12 if m == "MnIn2Se4" else -12)
            ax.plot([x], [v], "*", color=COL[m], ms=13, mec="k", mew=0.6, zorder=5)
            if m == "MnIn2Se4":                        # FC3-cutoff estimate, about -10 %
                ax.errorbar([x], [v], yerr=[[0.10 * v], [0]], color=COL[m], capsize=4, lw=1.6, zorder=4)
            else:                                       # upper bound
                ax.annotate("", xy=(x, v * 0.72), xytext=(x, v * 0.95),
                            arrowprops=dict(arrowstyle="-|>", color=COL[m], lw=1.6), zorder=4)
    for ax, ttl in ((a1, t("(a) in-plane")), (a2, t("(b) cross-plane") + ("" if LANG == "zh" else r" $\kappa_{zz}$"))):
        ax.set_yscale("log")
        ax.set_xlabel(t("Temperature (K)"))
        ax.set_ylabel(r"$\kappa_L$ (W m$^{-1}$ K$^{-1}$)")
        ax.set_title(ttl, fontsize=11)
        style(ax)
    star = Line2D([], [], ls="none", marker="*", color="0.4", mec="k", ms=11, label=t("300 K value to quote"))
    a1.legend(handles=a1.get_legend_handles_labels()[0] + [star], fontsize=8.5, frameon=False)
    a2.legend(fontsize=8.5, frameon=False)
    fig.tight_layout()
    return fig


def fig_calc_vs_exp(D):
    exp = {}
    for r in read_csv(D / "experiment_kappa_tot_digitized.csv")[1]:
        exp.setdefault((r[0], r[1], r[5]), []).append((float(r[2]), float(r[3])))
    km = json.loads((D / "independent_checks.json").read_text())["kmin_exp_velocities"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for i, m in enumerate(MAT):
        d = columns(D / f"kappa_ordered_vca_{m}.csv")
        l3 = columns(D / f"kappa_L3_{m}.csv")
        s = float(l3["s"][0])
        for j, (dr, key, comp) in enumerate((("Perp", "inplane", t("in-plane")),
                                             ("Para", "kzz", "层间" if LANG == "zh" else r"$\kappa_{zz}$"))):
            ax, c = axes[i, j], COL[m]
            ax.plot(d["T_K"], d[f"ordered_{key}"], "--", color=c, lw=1.4, label=t("calc. ordered single crystal"))
            ax.plot(d["T_K"], d[f"vca_{key}"], "-", color=c, lw=1.8, label=t("calc. VCA single crystal"))
            ax.plot(l3["T_K"], l3[dr], ":", color="k", lw=1.6, label=f"{t('calc. L3 fit')} (s = {s:.2f})")
            if m == "MnIn2Se4":                 # minimum kappa from the measured sound velocities
                ax.axhline(km["cahill"], color="0.45", lw=1.0, ls="-.", label=f"Cahill κ_min = {km['cahill']:.2f}")
                ax.axhline(km["diffuson"], color="0.65", lw=1.0, ls="-.", label=f"diffuson κ = {km['diffuson']:.2f}")
            for panel, filled in (("highT", True), ("lowT", False)):
                p = exp.get((m, dr, panel))
                if p:
                    p = np.array(sorted(p))
                    ax.plot(p[:, 0], p[:, 1], "o" if dr == "Perp" else "s", ms=5, color=c,
                            mfc=c if filled else "white", mec=c, ls="none",
                            label=t("experiment 300–623 K (fitted)") if filled else t("experiment 55–145 K (not fitted)"))
            ax.set_yscale("log")
            ax.set_xlim(0, 720)
            ax.set_title(f"{lab(m)} — {dr} ({t('calc.')} {comp})", fontsize=10.5)
            style(ax)
            if j == 0:
                ax.set_ylabel(r"$\kappa$ (W m$^{-1}$ K$^{-1}$)")
            if i == 1:
                ax.set_xlabel(t("Temperature (K)"))
            ax.legend(fontsize=7.5, frameon=False, loc="upper right")
    fig.tight_layout()
    return fig


def fig_scattering(D):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.2))
    for i, m in enumerate(MAT):
        f = next(D.glob(f"scattering_{m}_*K.csv"))
        T = f.stem.rsplit("_", 1)[1]
        d = columns(f)
        a1, a2 = axes[i]
        a1.scatter(d["freq_baseline_THz"], d["tau_inv_3ph_baseline_ps"], s=5, alpha=0.4, color="#1f77b4",
                   label=t("baseline 3-phonon"), rasterized=True)
        a1.scatter(d["freq_vca_THz"], d["tau_inv_3ph_vca_ps"], s=5, alpha=0.4, color="#d62728",
                   label=t("VCA 3-phonon"), rasterized=True)
        a1.set_title(f"{lab(m)} — {t('3-phonon scattering')} (~{T})", fontsize=10.5)
        fv, gv = d["freq_vca_THz"], d["tau_inv_disorder_vca_ps"]
        a2.scatter(fv, gv, s=5, alpha=0.5, color="#d62728", label=t("VCA disorder (mass variance)"), rasterized=True)
        a2.scatter(d["freq_baseline_THz"], d["tau_inv_disorder_baseline_ps"], s=5, alpha=0.5, color="#1f77b4",
                   label=t("baseline natural isotope"), rasterized=True)
        k = (fv > 1) & (gv > 0)
        sc = np.median(gv[k] / fv[k] ** 2)
        ff = np.linspace(0.1, fv.max(), 60)
        a2.plot(ff, sc * ff ** 2, "k--", lw=1, alpha=0.6, label=r"$\propto\omega^2$")
        a2.set_title(f"{lab(m)} — {t('disorder elastic scattering')}", fontsize=10.5)
        for ax in (a1, a2):
            ax.set_yscale("log")
            ax.set_ylabel(r"$\tau^{-1}$ (ps$^{-1}$)")
            style(ax)
            ax.legend(fontsize=8, frameon=False, loc="lower right", markerscale=2)
            if i == 1:
                ax.set_xlabel(t("Frequency (THz)"))
        a1.set_ylim(1e-4, 2)
        a2.set_ylim(1e-21, 2)
    fig.tight_layout()
    return fig


MAKERS = {"bands": fig_bands, "intrinsic": fig_intrinsic, "calc_vs_exp": fig_calc_vs_exp,
          "scattering": fig_scattering}


def make_all(data_dir, out_dir=None, png_for_html=False, lang="en"):
    """Draw every figure. Saves .png/.pdf to out_dir (if given; Chinese labels get a _zh suffix); returns
    {key: PNG bytes at 150 dpi} when png_for_html is set (used by the make_report_html scripts)."""
    set_lang(lang)
    suffix = "_zh" if lang == "zh" else ""
    data_dir = Path(data_dir)
    if out_dir:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
    out = {}
    for key, fn in MAKERS.items():
        fig = fn(data_dir)
        if out_dir:
            for ext in ("png", "pdf"):
                fig.savefig(Path(out_dir) / f"{STEMS[key]}{suffix}.{ext}", dpi=300, bbox_inches="tight")
        if png_for_html:
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
            out[key] = buf.getvalue()
        plt.close(fig)
    return out


def main():
    here = Path(__file__).resolve().parent
    default = next((p for p in (here / "plot_data", here.parent / "data") if p.is_dir()), here / "plot_data")
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=str(default))
    ap.add_argument("--lang", choices=("en", "zh", "both"), default="both")
    ap.add_argument("--out", default=str(here / "figures" if (here / "plot_data").is_dir() else here.parent / "figures"))
    a = ap.parse_args()
    for lang in (("en", "zh") if a.lang == "both" else (a.lang,)):
        make_all(a.data, a.out, lang=lang)
    print(f"wrote {len(MAKERS)} figures (.png + .pdf, languages: {a.lang}) to {a.out}")


if __name__ == "__main__":
    main()
