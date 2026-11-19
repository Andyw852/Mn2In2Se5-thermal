#!/usr/bin/env python3
"""Static HTML report: calculation vs experiment for MnIn2Se4 and Mn2In2Se5.

Writes report/calc_vs_experiment_report.html, a single self-contained file
(figures embedded as PNG, no JavaScript, opens offline). Every number comes from
files already in the repository; nothing is recomputed except the Mn2In2Se5
phonon dispersion, which is interpolated from its committed fc2.hdf5.

Sections
  1  300 K values to quote vs experiment (PROGRESS_REPORT_20261118.md §2)
  2  phonon dispersion, both materials
  3  intrinsic kappa(T)            (= paper fig. 1 data)
  4  calculation vs hot-pressed experiment, 2 x 2 (= paper fig. 2 data)
  5  scattering rates at 300 K, both materials (3-phonon and disorder)
  6  independent checks against experiment (lattice, Raman, sound velocity,
     minimum kappa); numbers from report/independent_checks.json
  7  what the L3 fit is
  8  caveats

Needs: numpy, matplotlib, h5py, pyyaml, phonopy (+ seekpath for the band path).

  python3 report/make_report_html.py
"""
import base64
import csv
import io
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "test" / "figures"
OUT = ROOT / "report" / "calc_vs_experiment_report.html"

MAT = ("MnIn2Se4", "Mn2In2Se5")
LABEL = {"MnIn2Se4": r"MnIn$_2$Se$_4$", "Mn2In2Se5": r"Mn$_2$In$_2$Se$_5$"}
HTML_LABEL = {"MnIn2Se4": "MnIn₂Se₄", "Mn2In2Se5": "Mn₂In₂Se₅"}
COL = {"MnIn2Se4": "#d62728", "Mn2In2Se5": "#1f77b4"}      # analysis/plotting convention
QUOTED = {"MnIn2Se4": {"inplane": 4.67, "kzz": 0.714},     # REPORT §9.1, 300 K
          "Mn2In2Se5": {"inplane": 3.37, "kzz": 1.089}}
CHK = json.loads((ROOT / "report" / "independent_checks.json").read_text())
SC_DIR = {"MnIn2Se4": FIG / "scattering", "Mn2In2Se5": FIG / "scattering" / "mn2in2se5"}


def load(rel):
    return json.loads((FIG / rel).read_text())


def style(ax):
    ax.grid(True, ls=":", alpha=0.4)
    ax.tick_params(direction="in", top=True, right=True)


def png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def tex_label(s):
    return (s.replace("$", "").replace("\\mathrm{", "").replace("}", "")
             .replace("\\Gamma", "Γ").replace("_", ""))


# --------------------------------------------------------------------------- bands
def bands_mnin2se4():
    b = yaml.safe_load((ROOT / "test/MnIn2Se4/fc2_442_phonon/band-dft-cpu.yaml").read_text())
    segs, i = [], 0
    for n in b["segment_nqpoint"]:
        blk = b["phonon"][i:i + n]
        i += n
        segs.append((np.array([p["distance"] for p in blk]),
                     np.array([[bb["frequency"] for bb in p["band"]] for p in blk])))
    labels = []
    for k, (a, c) in enumerate(b["labels"]):
        a, c = tex_label(a), tex_label(c)
        if k == 0:
            labels.append(a)
        elif a != labels[-1]:
            labels[-1] += "|" + a
        labels.append(c)
    return segs, labels


def bands_mn2in2se5():
    import phonopy
    d = ROOT / "Mn2In2Se5" / "fc"
    ph = phonopy.load(str(d / "phono3py_disp.yaml"), force_constants_filename=str(d / "fc2.hdf5"),
                      produce_fc=False, is_nac=False)
    ph.auto_band_structure(npoints=101, write_yaml=False)
    bs = ph.band_structure
    segs = [(np.array(x), np.array(f)) for x, f in zip(bs.distances, bs.frequencies)]
    lab, conn = [tex_label(x) for x in bs.labels], bs.path_connections
    labels, j = [lab[0]], 1
    for k, c in enumerate(conn):
        labels.append(lab[j])
        j += 1
        if not c and k < len(conn) - 1:
            labels[-1] += "|" + lab[j]
            j += 1
    return segs, labels


def fig_bands():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, m, (segs, labels) in zip(axes, MAT, (bands_mnin2se4(), bands_mn2in2se5())):
        ticks = [s[0][0] for s in segs] + [segs[-1][0][-1]]
        for x, f in segs:
            ax.plot(x, f, color=COL[m], lw=1.0)
        for t in ticks:
            ax.axvline(t, color="0.7", lw=0.7)
        ax.set_xticks(ticks, labels)
        ax.set_xlim(ticks[0], ticks[-1])
        ax.set_ylim(0, 8.2)
        nb = segs[0][1].shape[1]
        ax.set_title(f"{LABEL[m]}  ({nb} branches, max {max(s[1].max() for s in segs):.2f} THz)", fontsize=11)
        ax.tick_params(direction="in", right=True)
        ax.grid(True, axis="y", ls=":", alpha=0.4)
    axes[0].set_ylabel("Frequency (THz)")
    fig.tight_layout()
    return png(fig)


# --------------------------------------------------------------------------- kappa(T)
def fig_intrinsic():
    data = {"MnIn2Se4": load("kappa_vs_T/data_mnin2se4.json"),
            "Mn2In2Se5": load("kappa_vs_T/data_mnin2se5.json")}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for m in MAT:
        d = data[m]
        T = np.array(d["T_K"], float)
        k = T >= 30
        a1.plot(T[k], np.array(d["inplane"])[k], "-o", color=COL[m], ms=4, label=f"{LABEL[m]} ({d['mesh']})")
        a2.plot(T[k], np.array(d["kzz_rest"])[k], "-s", color=COL[m], ms=4, mfc="white",
                label=f"{LABEL[m]} ({d['mesh']})")
    for ax, key in ((a1, "inplane"), (a2, "kzz")):
        for m in MAT:
            v = QUOTED[m][key]
            x = 300 + (12 if m == "MnIn2Se4" else -12)
            ax.plot([x], [v], "*", color=COL[m], ms=13, mec="k", mew=0.6, zorder=5)
            if m == "MnIn2Se4":
                ax.errorbar([x], [v], yerr=[[0.10 * v], [0]], color=COL[m], capsize=4, lw=1.6, zorder=4)
            else:
                ax.annotate("", xy=(x, v * 0.72), xytext=(x, v * 0.95),
                            arrowprops=dict(arrowstyle="-|>", color=COL[m], lw=1.6), zorder=4)
    for ax, t in ((a1, "(a) in-plane"), (a2, r"(b) cross-plane $\kappa_{zz}$")):
        ax.set_yscale("log")
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(r"$\kappa_L$ (W m$^{-1}$ K$^{-1}$)")
        ax.set_title(t, fontsize=11)
        style(ax)
    star = Line2D([], [], ls="none", marker="*", color="0.4", mec="k", ms=11, label="300 K value to quote")
    a1.legend(handles=a1.get_legend_handles_labels()[0] + [star], fontsize=8.5, frameon=False)
    a2.legend(fontsize=8.5, frameon=False)
    fig.tight_layout()
    return png(fig)


def exp_rows():
    out = {}
    with open(ROOT / "experiment" / "kappa_tot_digitized.csv", encoding="utf-8") as f:
        for r in csv.DictReader(l for l in f if not l.startswith("#")):
            out.setdefault((r["material"], r["direction"], r["panel"]), []).append(
                (float(r["T_K"]), float(r["kappa_tot_WmK"])))
    return {k: np.array(sorted(v)) for k, v in out.items()}


def extrinsic(m):
    d = load(f"calc_vs_experiment/extrinsic_{m.lower()}.json")
    T = np.array([r["T"] for r in d["rows"]])
    kxx = np.array([r["kxx"] for r in d["rows"]])
    kzz = np.array([r["kzz"] for r in d["rows"]])
    s = d["s"]
    return T, {"Perp": (1 - s / 2) * kxx + (s / 2) * kzz, "Para": s * kxx + (1 - s) * kzz}, s


def fig_calc_vs_exp():
    exp = exp_rows()
    vca = {"MnIn2Se4": load("kappa_vs_T/data_mnin2se4_vca.json"),
           "Mn2In2Se5": load("kappa_vs_T/data_mn2in2se5_vca.json")}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for i, m in enumerate(MAT):
        d = vca[m]
        T = np.array(d["T_K"], float)
        Tl3, l3, s = extrinsic(m)
        for j, (dr, key, comp) in enumerate((("Perp", "inplane", "in-plane"),
                                             ("Para", "kzz_rest", r"$\kappa_{zz}$"))):
            ax = axes[i, j]
            c = COL[m]
            ax.plot(T, d[f"baseline_{key}"], "--", color=c, lw=1.4, label="calc. ordered single crystal")
            ax.plot(T, d[f"vca_{key}"], "-", color=c, lw=1.8, label="calc. VCA single crystal")
            ax.plot(Tl3, l3[dr], ":", color="k", lw=1.6, label=f"calc. L3 fit (s = {s:.2f})")
            if m == "MnIn2Se4":                 # minimum kappa from the measured sound velocities
                km = CHK["kmin_exp_velocities"]
                ax.axhline(km["cahill"], color="0.45", lw=1.0, ls="-.", label=f"Cahill κ_min = {km['cahill']:.2f}")
                ax.axhline(km["diffuson"], color="0.65", lw=1.0, ls="-.", label=f"diffuson κ = {km['diffuson']:.2f}")
            for panel, filled in (("highT", True), ("lowT", False)):
                p = exp.get((m, dr, panel))
                if p is not None:
                    ax.plot(p[:, 0], p[:, 1], "o" if dr == "Perp" else "s", ms=5, color=c,
                            mfc=c if filled else "white", mec=c, ls="none",
                            label="experiment 300–623 K (fitted)" if filled else "experiment 55–145 K (not fitted)")
            ax.set_yscale("log")
            ax.set_xlim(0, 720)
            ax.set_title(f"{LABEL[m]} — {dr} (calc. {comp})", fontsize=10.5)
            style(ax)
            if j == 0:
                ax.set_ylabel(r"$\kappa$ (W m$^{-1}$ K$^{-1}$)")
            if i == 1:
                ax.set_xlabel("Temperature (K)")
            ax.legend(fontsize=7.5, frameon=False, loc="upper right")
    fig.tight_layout()
    return png(fig)


# --------------------------------------------------------------------------- scattering
def load_sc(fn):
    with h5py.File(fn, "r") as f:
        fr = f["frequency"][:].ravel()
        g3 = f["gamma"][0].ravel() * 2.0          # tau^-1 = 2 Gamma (ps^-1)
        gi = f["gamma_isotope"][:].ravel() * 2.0
        T = float(f["temperature"][0])
    o = np.argsort(fr)
    return fr[o], g3[o], gi[o], T


def fig_scattering():
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.2))
    for i, m in enumerate(MAT):
        fb, g3b, gib, T = load_sc(SC_DIR[m] / "kappa_baseline.hdf5")
        fv, g3v, giv, _ = load_sc(SC_DIR[m] / "kappa_vca.hdf5")
        a1, a2 = axes[i]
        a1.scatter(fb, g3b, s=5, alpha=0.4, color="#1f77b4", label="baseline 3-phonon", rasterized=True)
        a1.scatter(fv, g3v, s=5, alpha=0.4, color="#d62728", label="VCA 3-phonon", rasterized=True)
        a1.set_title(f"{LABEL[m]} — 3-phonon scattering (~{T:.0f} K)", fontsize=10.5)
        a2.scatter(fv, giv, s=5, alpha=0.5, color="#d62728", label="VCA disorder (mass variance)", rasterized=True)
        a2.scatter(fb, gib, s=5, alpha=0.5, color="#1f77b4", label="baseline natural isotope", rasterized=True)
        mask = (fv > 1) & (giv > 0)
        sc = np.median(giv[mask] / fv[mask] ** 2)
        ff = np.linspace(0.1, fv.max(), 60)
        a2.plot(ff, sc * ff ** 2, "k--", lw=1, alpha=0.6, label=r"$\propto\omega^2$")
        a2.set_title(f"{LABEL[m]} — disorder elastic scattering", fontsize=10.5)
        for ax in (a1, a2):
            ax.set_yscale("log")
            ax.set_ylabel(r"$\tau^{-1}$ (ps$^{-1}$)")
            style(ax)
            ax.legend(fontsize=8, frameon=False, loc="lower right", markerscale=2)
            if i == 1:
                ax.set_xlabel("Frequency (THz)")
        a1.set_ylim(1e-4, 2)
        a2.set_ylim(1e-21, 2)
    fig.tight_layout()
    return png(fig)


# --------------------------------------------------------------------------- HTML
CSS = """
:root{--fg:#1b2129;--muted:#5b6573;--rule:#dde2e8;--bg:#f6f7f9;--panel:#fff;--accent:#2b5d8a;--red:#c8372d;--blue:#1f6fb4}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.6 "IBM Plex Sans","Segoe UI",Helvetica,Arial,sans-serif}
.wrap{max-width:1040px;margin:0 auto;padding:32px 16px 64px}
h1{font-size:28px;line-height:1.25;margin:0 0 8px}
h2{font-size:21px;margin:44px 0 6px}
h2 .n{color:var(--muted);font:500 14px ui-monospace,Menlo,Consolas,monospace;margin-right:10px}
.eyebrow{font:12px ui-monospace,Menlo,Consolas,monospace;letter-spacing:.08em;text-transform:uppercase;color:var(--accent)}
p.lead,.cap{color:var(--muted);max-width:82ch}
.cap{font-size:13.5px}
.cap b{color:var(--fg)}
figure{margin:12px 0 0;background:var(--panel);border:1px solid var(--rule);border-radius:6px;padding:10px}
figure img{width:100%;height:auto;display:block}
.tbl{overflow-x:auto;background:var(--panel);border:1px solid var(--rule);border-radius:6px;margin-top:12px}
table{border-collapse:collapse;width:100%;font-size:14px;font-variant-numeric:tabular-nums}
th,td{padding:7px 12px;text-align:right;border-bottom:1px solid var(--rule)}
th:first-child,td:first-child{text-align:left}
thead th{font:500 12px ui-monospace,Menlo,Consolas,monospace;color:var(--muted)}
tbody tr:last-child td{border-bottom:none}
.r{color:var(--red)}.b{color:var(--blue)}
.x{font:12px ui-monospace,Menlo,Consolas,monospace;color:var(--muted)}
code{font:13px ui-monospace,Menlo,Consolas,monospace}
.eq{background:var(--panel);border:1px solid var(--rule);border-radius:6px;padding:12px 16px;overflow-x:auto;
    font:13.5px/1.7 ui-monospace,Menlo,Consolas,monospace;white-space:pre}
ul{max-width:82ch;padding-left:20px}
li{margin:6px 0}
"""


def checks_html():
    L, c = HTML_LABEL, CHK
    lc, le = c["lattice"]["calc"], c["lattice"]["exp"]
    vc = c["velocity_calc"]
    e1, e2 = c["velocity_exp"]
    pct = lambda a, b: f"{100 * (a / b - 1):+.1f} %"
    raman = ", ".join(f"{m['irrep']} {m['cm-1']:.0f}" for m in c["raman_calc"])
    km, kc = c["kmin_exp_velocities"], c["kmin_calc_velocities"]
    return f"""<h2><span class="n">6</span>Independent checks against experiment ({L['MnIn2Se4']})</h2>
<p class="lead">Experimental data that were not used in any fit. They test the force constants and give a lower
reference for κ. Computed by <code>report/independent_checks.py</code> from the 4×4×2 FC2.</p>
<div class="tbl"><table>
<thead><tr><th>quantity</th><th>experiment</th><th>calculation</th><th>calc. / exp.</th></tr></thead>
<tbody>
<tr><td>lattice a (Å, hexagonal)</td><td>{le['a']:.2f}</td><td>{lc['a']:.3f}</td><td>{pct(lc['a'], le['a'])}</td></tr>
<tr><td>lattice c (Å)</td><td>{le['c']:.2f}</td><td>{lc['c']:.2f}</td><td>{pct(lc['c'], le['c'])}</td></tr>
<tr><td>Raman peaks (cm⁻¹)</td><td>≈ 78, ≈ 177 (read off the slide)</td><td>A1g 75, A1g 177</td><td>within the reading error</td></tr>
<tr><td>v_L (m/s)</td><td>{e1['vL_ms']:.0f} / {e2['vL_ms']:.0f}</td><td>{vc['avg_vL']} (c axis {vc['c_axis'][2]:.0f}, in-plane {vc['in_plane'][2]:.0f})</td><td>{pct(vc['avg_vL'], (e1['vL_ms'] + e2['vL_ms']) / 2)}</td></tr>
<tr><td>v_T (m/s)</td><td>{e1['vT_ms']:.0f} / {e2['vT_ms']:.0f}</td><td>{vc['avg_vT']} (range {min(vc['in_plane'][0], vc['c_axis'][0]):.0f}–{max(vc['in_plane'][1], vc['c_axis'][1]):.0f})</td><td>{pct(vc['avg_vT'], (e1['vT_ms'] + e2['vT_ms']) / 2)}</td></tr>
<tr><td>Debye mean v_m (m/s)</td><td>{e1['v_mean_ms']:.0f} / {e2['v_mean_ms']:.0f}</td><td>{vc['debye_mean']}</td><td>{pct(vc['debye_mean'], (e1['v_mean_ms'] + e2['v_mean_ms']) / 2)}</td></tr>
<tr><td>Cahill κ_min (W m⁻¹ K⁻¹)</td><td>{km['cahill']:.2f} (measured v)</td><td>{kc['cahill']:.2f} (calculated v)</td><td></td></tr>
<tr><td>diffuson κ (W m⁻¹ K⁻¹)</td><td>{km['diffuson']:.2f} (measured v)</td><td>{kc['diffuson']:.2f} (calculated v)</td><td></td></tr>
</tbody></table></div>
<ul>
<li><b>Lattice and Raman</b>: the relaxed structure is within 0.7 % of the single-crystal XRD values, and the two
strong Raman peaks match the calculated A1g modes (all Raman-active modes: {raman} cm⁻¹). Both support the force constants.</li>
<li><b>Sound velocity</b>: ultrasonic, hot-pressed pellets, two samples (agree within 2 %). The measured v_L is below
the calculated v_L in <i>every</i> direction (c axis {vc['c_axis'][2]:.0f} m/s is the lowest), so texture cannot explain
it; porosity and defect softening of the pellet are the likely reasons. The pellet density is not known yet
(theoretical {c['theoretical_density_g_cm3']:.2f} g/cm³), so the two cannot be separated. The propagation direction
relative to the pressing axis is also not known.</li>
<li><b>Minimum κ</b> (Cahill, Watson &amp; Pohl, PRB 46, 6131 (1992); diffuson model, Agne et al., Energy Environ. Sci.
11, 609 (2018)), n = {c['number_density_A-3']:.4f} Å⁻³: the measured κ_tot (Perp 0.93, Para 0.49 at 303 K) lies between
the diffuson value {km['diffuson']:.2f} and the Cahill value {km['cahill']:.2f}. The hot-pressed samples are close to the
amorphous limit, which agrees independently with the strong extrinsic scattering needed in the L3 fit (L ≈ tens of nm).
These limits are isotropic; Para lying below the Cahill value is common for layered materials.</li>
</ul>"""


def html(figs):
    L = HTML_LABEL
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>MnIn₂Se₄ / Mn₂In₂Se₅ calculation vs experiment</title>
<style>{CSS}</style></head><body><div class="wrap">

<div class="eyebrow">Lattice thermal conductivity · generated by report/make_report_html.py</div>
<h1>{L['MnIn2Se4']} and {L['Mn2In2Se5']}: calculation vs experiment</h1>
<p class="lead">Phonon dispersion, intrinsic κ(T), comparison with the hot-pressed measurements, and the
3-phonon and disorder scattering rates at 300 K. All numbers come from files already in this repository
(see <a href="../PROGRESS_REPORT_20261118.md">PROGRESS_REPORT_20261118.md</a>); only the {L['Mn2In2Se5']}
dispersion is interpolated here from its committed <code>fc2.hdf5</code>.</p>

<h2><span class="n">1</span>Values to quote at 300 K vs experiment</h2>
<p class="lead">W m⁻¹ K⁻¹. Experiment: Perp ≈ in-plane, Para ≈ cross-plane (the platelet c axes tend to align
with the pressing direction). Ratio = calculation / experiment.</p>
<div class="tbl"><table>
<thead><tr><th></th><th class="r">{L['MnIn2Se4']}</th><th class="b">{L['Mn2In2Se5']}</th></tr></thead>
<tbody>
<tr><td>calc. in-plane κ (LBTE + coherence)</td><td><b>4.67</b> <span class="x">cutoff ≈ −10 %</span></td><td><b>3.37</b> <span class="x">upper bound</span></td></tr>
<tr><td>calc. κ_zz (mesh-converged + coherence)</td><td><b>0.71</b> <span class="x">cutoff ≈ −10 %</span></td><td><b>1.09</b> <span class="x">upper bound</span></td></tr>
<tr><td>experiment κ_tot Perp / Para (303 K, hot-pressed)</td><td>0.933 / 0.485</td><td>0.741 / 0.320</td></tr>
<tr><td>texture-independent (2·xx + zz)/3 — experiment</td><td>0.784</td><td>0.601</td></tr>
<tr><td>&nbsp;&nbsp;calc., quoted values</td><td>3.35 <span class="x">4.3×</span></td><td>2.61 <span class="x">4.3×</span></td></tr>
<tr><td>&nbsp;&nbsp;calc., ordered, RTA 20³</td><td>2.67 <span class="x">3.4×</span></td><td>2.33 <span class="x">3.9×</span></td></tr>
<tr><td>&nbsp;&nbsp;calc., VCA antisite disorder, RTA 20³</td><td>1.14 <span class="x">1.5×</span></td><td>1.35 <span class="x">2.3×</span></td></tr>
<tr><td>&nbsp;&nbsp;calc., L3 fit (fitted to these data, §7)</td><td>0.77</td><td>0.62</td></tr>
</tbody></table></div>
<p class="cap">The intrinsic single-crystal κ is 3–4× above the hot-pressed samples. Antisite mass disorder (VCA)
closes part of the gap; the rest needs a nanoscale effective mean free path (L_eff ≈ 10–50 nm).</p>

<h2><span class="n">2</span>Phonon dispersion</h2>
<p class="lead">R-3m primitive cells, path Γ–T–H₂|H₀–L–Γ–S₀|S₂–F–Γ. No imaginary modes in either material; NAC off.</p>
<figure><img alt="Phonon dispersion of both materials" src="{figs['bands']}"></figure>
<p class="cap"><b>{L['MnIn2Se4']}</b>: 4×4×2 FC2, <code>test/MnIn2Se4/fc2_442_phonon/band-dft-cpu.yaml</code>.
<b>{L['Mn2In2Se5']}</b>: 4×4×1 FC2, <code>Mn2In2Se5/fc/fc2.hdf5</code> (phonopy interpolation).
The 4×4×1 FC2 of {L['MnIn2Se4']} gave a too-soft interlayer TA (322.9 m/s) and is not used; the 4×4×2 FC2 gives
1348.3 m/s. The {L['Mn2In2Se5']} 4×4×1 FC2 was checked the same way and is fine (1386.8 m/s).</p>

<h2><span class="n">3</span>Intrinsic κ(T)</h2>
<p class="lead">Ordered crystal, 3-phonon + natural isotopes, RTA. κ_zz is the mesh-converged term (Γ–Z quadrature
error removed). Stars: the 300 K values of §1.</p>
<figure><img alt="Intrinsic kappa versus temperature" src="{figs['intrinsic']}"></figure>
<p class="cap">Meshes 14³ ({L['MnIn2Se4']}) and 16³ ({L['Mn2In2Se5']}). The stars sit above the lines because
in-plane LBTE exceeds RTA (24 % for {L['MnIn2Se4']}) and because they include coherence. Bar: FC3-cutoff
estimate (≈ −10 %). Arrow: {L['Mn2In2Se5']} cutoff not converged, value is an upper bound.</p>

<h2><span class="n">4</span>Calculation vs hot-pressed experiment</h2>
<p class="lead">One panel per material and direction. Dashed: ordered single crystal. Solid: VCA antisite-disordered
single crystal. Dotted: L3 fit (§7). Dash-dot ({L['MnIn2Se4']} only): minimum κ from the measured sound velocities (§6). Filled symbols: 300–623 K (used in the fit). Open symbols: 55–145 K (not fitted).</p>
<figure><img alt="Calculated kappa against the hot-pressed experiment" src="{figs['calc_vs_exp']}"></figure>
<p class="cap">Ordered and VCA: RTA, 20³ (same mesh for each pair), Perp ↔ in-plane, Para ↔ κ_zz, no grain boundary
and no texture, so they are upper references. VCA occupancy: {L['MnIn2Se4']} i = 0.44 (Döll, source still to be
supplied), {L['Mn2In2Se5']} x = 0.13 (taken from the L3 fit, no structural source). VCA lowers κ strongly but stays
well above experiment and cannot produce the low-temperature maximum near 100 K; that needs a nanoscale mean free path.</p>

<h2><span class="n">5</span>Scattering rates at 300 K</h2>
<p class="lead">Largest FC3 cutoff, 20³ mesh, τ⁻¹ = 2Γ. Baseline: real masses + natural isotopes. VCA: averaged
site masses + mass-variance scattering. Top row {L['MnIn2Se4']} (i = 0.44), bottom row {L['Mn2In2Se5']} (x = 0.13).</p>
<figure><img alt="Three-phonon and disorder scattering rates, both materials" src="{figs['scattering']}"></figure>
<p class="cap">Left: the 3-phonon rates of the ordered and VCA crystals agree within a few tens of percent, because VCA
changes the masses, not the anharmonicity. Right: VCA mass-variance scattering is 2–3 orders of magnitude stronger
than natural-isotope scattering, rises as ω² at low frequency and flattens at high frequency. Effect on κ at 300 K:
{L['MnIn2Se4']} in-plane 3.69 → 1.49 (2.48×), κ_zz 0.638 → 0.445 (1.43×); {L['Mn2In2Se5']} in-plane
3.00 → 1.63 (1.84×), κ_zz 0.984 → 0.792 (1.24×). The weaker effect in {L['Mn2In2Se5']} comes from its smaller x.</p>

{checks_html()}

<h2><span class="n">7</span>What the L3 fit is</h2>
<p class="lead">L3 is the extrinsic model of <code>analysis/extrinsic/</code>. <b>Yes, it is fitted to the
experimental values</b>, so its agreement with experiment is by construction and is not a prediction.</p>
<div class="eq">Γ_λ = Γ_3ph,λ(T)  +  Σ g_elem(x)·B_elem,λ  +  |v_λ| / (4π L)
       3-phonon           Mn↔In antisite mass disorder   grain boundary (Casimir)

κ_Para = s·κ_xx + (1−s)·κ_zz          κ_Perp = (1−s/2)·κ_xx + (s/2)·κ_zz</div>
<ul>
<li><b>Inputs that are not fitted</b>: the first-principles 3-phonon rates, phonon velocities and heat capacities
(ordered masses), and the mass-disorder basis rates B, computed once per material on a 24³ mesh.</li>
<li><b>Three fitted parameters</b>: x, the fraction of Mn on In sites; L, an effective boundary mean free path;
s = ⟨sin²θ⟩, the texture (0 = fully textured, 2/3 = random).</li>
<li><b>Data used</b>: only the 300–623 K points of <code>experiment/kappa_tot_digitized.csv</code>, both directions.
The fit minimises the rms of ln(model/experiment) over an (x, L) grid, with s optimised at each point.
The 55–145 K points are not used; they only test the model.</li>
<li><b>Result</b>: {L['MnIn2Se4']} x = 0.45, L = 30.7 nm, s = 0.35, rms 1.2 %; {L['Mn2In2Se5']} x = 0.13,
L = 15.9 nm, s = 0.25, rms 3.2 %.</li>
<li><b>What it does determine</b>: only s (acceptable 0.33–0.38 and 0.15–0.28). x and L both push κ down and are
degenerate: acceptable x spans 0–0.66 and 0.01–0.50. The fitted x is therefore not evidence for any occupancy;
the closeness of 0.45 to Döll's 0.44 is a coincidence of the fit. XRD (texture), SEM (grain size) and a Rietveld
occupancy refinement would break the degeneracy.</li>
<li><b>Not to be combined with VCA</b>: L3 already contains the antisite scattering through x, so the VCA and L3
curves are alternatives, not additive corrections.</li>
<li><b>κ_e</b> is taken as 0. If the electronic part is not negligible, the fitted extrinsic scattering is too strong.</li>
</ul>

<h2><span class="n">8</span>Caveats</h2>
<ul>
<li>The experimental points are digitised from a figure (±0.005 W m⁻¹ K⁻¹ at high T, ±0.01 at low T) and are κ_tot
with κ_e not subtracted. Replace <code>experiment/kappa_tot_digitized.csv</code> once the raw data are available.</li>
<li>{L['Mn2In2Se5']}: the FC3-cutoff steps do not shrink (in-plane −3.9 % → −4.3 %, κ_zz −11.9 % → −9.9 %), so its
values are upper bounds and the κ_zz ordering of the two materials is not robust. The in-plane ordering
({L['MnIn2Se4']} 37 % higher) is.</li>
<li>Between the 145 K and 300 K experimental points the Para values drop much more than Perp ({L['Mn2In2Se5']} −32 %
vs −6 %), possibly two instruments or two samples; this is why the fit uses the high-temperature range only.</li>
</ul>
</div></body></html>
"""


def main():
    figs = {"bands": fig_bands(), "intrinsic": fig_intrinsic(),
            "calc_vs_exp": fig_calc_vs_exp(), "scattering": fig_scattering()}
    OUT.write_text(html(figs), encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
