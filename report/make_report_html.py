#!/usr/bin/env python3
"""Static HTML report: calculation vs experiment for MnIn2Se4 and Mn2In2Se5.

Writes report/calc_vs_experiment_report.html, a single self-contained file
(figures embedded as PNG, no JavaScript, opens offline). Every number comes from
files already in the repository; nothing is recomputed except the Mn2In2Se5
phonon dispersion, which is interpolated from its committed fc2.hdf5.

Sections
  1  300 K values to quote vs experiment (PROGRESS_REPORT_20261118.md §2)
  2  calculation setup and convergence: FC2 / FC3 sources, FC3-cutoff and q-mesh
     tables (10 % criterion), RTA / LBTE / Wigner; data in plot_data/convergence_*.csv
  2  phonon dispersion, both materials
  3  intrinsic kappa(T)            (= paper fig. 1 data)
  4  calculation vs hot-pressed experiment, 2 x 2 (= paper fig. 2 data)
  5  scattering rates at 300 K, both materials (3-phonon and disorder)
  6  independent checks against experiment (lattice, Raman, sound velocity,
     minimum kappa); numbers from report/independent_checks.json
  7  what the L3 fit is
  8  caveats

Pipeline (one source for the HTML and the figure files):
  export_plot_data.py  repository -> report/plot_data/*.csv, *.json
  plot_figures.py      report/plot_data/ -> report/figures/fig*.png (300 dpi) + .pdf
  make_report_html.py  runs both, embeds the figures (150 dpi) in the HTML

Needs: numpy, matplotlib, h5py, pyyaml, phonopy (+ seekpath) for the export;
only numpy + matplotlib with --no-export (reuses report/plot_data/).

  python3 report/make_report_html.py [--no-export]
"""
import base64
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import plot_figures              # noqa: E402  plot_data/ -> figures

# In the repository this file sits in report/ next to plot_data/. In the exported
# bundle it sits in scripts/, with ../data, ../figures and ../report.
BUNDLE = not (HERE / "plot_data").is_dir() and (HERE.parent / "data").is_dir()
DATA = HERE.parent / "data" if BUNDLE else HERE / "plot_data"
FIGDIR = HERE.parent / "figures" if BUNDLE else HERE / "figures"
OUT = (HERE.parent / "report" if BUNDLE else HERE) / "calc_vs_experiment_report.html"
HTML_LABEL = {"MnIn2Se4": "MnIn₂Se₄", "Mn2In2Se5": "Mn₂In₂Se₅"}
CHK = json.loads((DATA / "independent_checks.json").read_text())


CSS = """
:root{--fg:#1b2129;--muted:#5b6573;--rule:#dde2e8;--bg:#f6f7f9;--panel:#fff;--accent:#2b5d8a;--red:#c8372d;--blue:#1f6fb4}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.6 "IBM Plex Sans","Segoe UI",Helvetica,Arial,sans-serif}
.wrap{max-width:1040px;margin:0 auto;padding:32px 16px 64px}
h1{font-size:28px;line-height:1.25;margin:0 0 8px}
h2{font-size:21px;margin:44px 0 6px}
h3{font-size:16px;margin:26px 0 4px}
.lt td,.lt th{text-align:left;white-space:normal;vertical-align:top}
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


def _rows(name):
    import csv as _csv
    with open(DATA / name, encoding="utf-8") as f:
        return list(_csv.DictReader(l for l in f if not l.startswith("#")))


def table1_values():
    """300 K numbers for §1, read from plot_data/: per direction (in-plane, kappa_zz) and the
    texture-independent average (2 in-plane + kappa_zz)/3, for calculation levels and experiment."""
    q = {(r["material"], r["component"]): float(r["kappa_WmK"]) for r in _rows("kappa_quoted_300K.csv")}
    out = {}
    for m in ("MnIn2Se4", "Mn2In2Se5"):
        v = next(r for r in _rows(f"kappa_ordered_vca_{m}.csv") if float(r["T_K"]) == 300)
        e = {r["direction"]: float(r["kappa_tot_WmK"]) for r in _rows("experiment_kappa_tot_digitized.csv")
             if r["material"] == m and r["panel"] == "highT" and float(r["T_K"]) == 303}
        rows = {"quoted": (q[(m, "inplane")], q[(m, "kzz")]),
                "ordered": (float(v["ordered_inplane"]), float(v["ordered_kzz"])),
                "vca": (float(v["vca_inplane"]), float(v["vca_kzz"])),
                "exp": (e["Perp"], e["Para"])}
        out[m] = {k: (a, z, (2 * a + z) / 3) for k, (a, z) in rows.items()}
    return out


def table1_body(labels, notes):
    """Rows of the two §1 tables. labels: row label per level; notes: quoted-value tags per material."""
    V = table1_values()
    ms = ("MnIn2Se4", "Mn2In2Se5")
    r = lambda a, b: f' <span class="x">{a / b:.1f}×</span>'
    t1, t2 = [], []
    for k in ("quoted", "ordered", "vca"):
        c1 = []
        for m in ms:
            a, z, avg = V[m][k]
            ea, ez, eavg = V[m]["exp"]
            b = ("<b>", "</b>") if k == "quoted" else ("", "")
            tag = f' <span class="x">{notes[m]}</span>' if k == "quoted" else ""
            c1 += [f"<td>{b[0]}{a:.2f}{b[1]}{r(a, ea)}{tag}</td>", f"<td>{b[0]}{z:.2f}{b[1]}{r(z, ez)}</td>"]
        t1.append(f"<tr><td>{labels[k]}</td>{''.join(c1)}</tr>")
        t2.append(f"<tr><td>{labels[k]}</td>" + "".join(
            f"<td>{V[m][k][2]:.2f}{r(V[m][k][2], V[m]['exp'][2])}</td>" for m in ms) + "</tr>")
    t1.append(f"<tr><td>{labels['exp']}</td>" + "".join(
        f"<td><b>{V[m]['exp'][0]:.3f}</b></td><td><b>{V[m]['exp'][1]:.3f}</b></td>" for m in ms) + "</tr>")
    t2.append(f"<tr><td>{labels['exp_avg']}</td>" + "".join(f"<td><b>{V[m]['exp'][2]:.3f}</b></td>" for m in ms) + "</tr>")
    return "\n".join(t1), "\n".join(t2)


def ordering_values():
    """Mn2In2Se5 relative to MnIn2Se4 at 300 K, per direction, for experiment and calculation,
    plus the texture mixing with the fitted s, the crossover s and the disorder sweep."""
    V = table1_values()
    A, B = "MnIn2Se4", "Mn2In2Se5"
    s = {m: float(_rows(f"kappa_L3_{m}.csv")[0]["s"]) for m in (A, B)}

    def mix(m, a, z):
        return ((1 - s[m] / 2) * a + s[m] / 2 * z, s[m] * a + (1 - s[m]) * z, (2 * a + z) / 3)

    rel = lambda b, a: [y / x - 1 for x, y in zip(a, b)]
    out = {"s": s,
           "exp": rel(V[B]["exp"], V[A]["exp"]),
           "quoted": rel(V[B]["quoted"], V[A]["quoted"]),
           "ordered": rel(V[B]["ordered"], V[A]["ordered"]),
           "mixed": rel(mix(B, *V[B]["quoted"][:2]), mix(A, *V[A]["quoted"][:2])),
           "mixed_abs": {m: mix(m, *V[m]["quoted"][:2]) for m in (A, B)}}
    (a1, z1), (a2, z2) = V[A]["quoted"][:2], V[B]["quoted"][:2]
    out["s_cross"] = (z2 - z1) / ((a1 - z1) - (a2 - z2))
    sw = {}
    for r in _rows("disorder_sweep_300K.csv"):
        sw.setdefault(r["material"], {})[float(r["x"])] = (float(r["inplane"]), float(r["kzz_total"]))
    out["sweep"] = sw
    return out


def ordering_rows(labels):
    O = ordering_values()
    pc = lambda v: f"{100 * v:+.0f} %"

    def cells(vals, ref=None):
        out = []
        for i, v in enumerate(vals):
            mark = "" if ref is None else (" ✓" if (v < 0) == (ref[i] < 0) else " ✗")
            out.append(f"<td>{pc(v)}{mark}</td>")
        return "".join(out)
    rows = [f"<tr><td>{labels['exp']}</td>{cells(O['exp'])}</tr>"]
    for k in ("quoted", "ordered", "mixed"):
        rows.append(f"<tr><td>{labels[k]}</td>{cells(O[k], O['exp'])}</tr>")
    return "\n".join(rows), O


def _step(new, old):
    return f"{100 * (new / old - 1):+.1f} %"


def _ok(new, old):
    return "✓" if abs(new / old - 1) < 0.10 else "✗"


def methods_html():
    L = HTML_LABEL
    cut = _rows("convergence_fc3_cutoff.csv")
    msh = _rows("convergence_qmesh.csv")
    sv = {(r["material"], r["quantity"]): float(r["value"]) for r in _rows("solver_terms.csv")}

    def cut_rows(m):
        rr = [r for r in cut if r["material"] == m]
        out = []
        for i, r in enumerate(rr):
            a, z = float(r["inplane_population"]), float(r["kzz_rest"])
            if i:
                a0, z0 = float(rr[i - 1]["inplane_population"]), float(rr[i - 1]["kzz_rest"])
                sa, sz = f"{_step(a, a0)} {_ok(a, a0)}", f"{_step(z, z0)} {_ok(z, z0)}"
            else:
                sa = sz = "—"
            b = "<b>" if i == len(rr) - 1 else ""
            e = "</b>" if b else ""
            out.append(f"<tr><td>{L[m]}</td><td>{b}{float(r['cutoff_A']):.2f}{e}</td><td>{a:.3f}</td><td>{sa}</td>"
                       f"<td>{z:.3f}</td><td>{sz}</td></tr>")
        return "\n".join(out)

    def mesh_rows(m):
        rr = [r for r in msh if r["material"] == m]
        out = []
        for i, r in enumerate(rr):
            v = [float(r[k]) for k in ("inplane", "kzz_total", "kzz_rest")]
            if i:
                v0 = [float(rr[i - 1][k]) for k in ("inplane", "kzz_total", "kzz_rest")]
                st = [f"{_step(x, x0)} {_ok(x, x0)}" for x, x0 in zip(v, v0)]
            else:
                st = ["—"] * 3
            n = r["N"]
            out.append(f"<tr><td>{L[m]}</td><td>{n}³</td><td>{v[0]:.4f}</td><td>{st[0]}</td><td>{v[1]:.4f}</td>"
                       f"<td>{st[1]}</td><td>{v[2]:.4f}</td><td>{st[2]}</td></tr>")
        return "\n".join(out)

    m4, m5 = "MnIn2Se4", "Mn2In2Se5"
    r4 = sv[(m4, "inplane_LBTE_14")] / sv[(m4, "inplane_RTA_14")]
    r5 = sv[(m5, "inplane_LBTE_10")] / sv[(m5, "inplane_RTA_10")]
    q4 = sv[(m4, "inplane_LBTE_14")] * sv[(m4, "mesh_correction_inplane")] + sv[(m4, "coherence_inplane")]
    q5 = sv[(m5, "inplane_RTA_32")] * r5 + sv[(m5, "coherence_inplane")]
    z4 = sv[(m4, "population_kzz_rest")] + sv[(m4, "coherence_kzz")]
    z5 = sv[(m5, "population_kzz_rest")] + sv[(m5, "coherence_kzz")]
    return f"""<h2><span class="n">2</span>Calculation setup and convergence</h2>
<p class="lead">How the force constants were obtained, how far κ is converged in the FC3 cutoff and the q-mesh,
and which solver each quoted number comes from. <b>Convergence criterion: a parameter counts as converged when its
last refinement step changes κ by less than 10 %</b> (✓ in the tables), the accuracy target of this project.
Data: <code>plot_data/convergence_fc3_cutoff.csv</code>, <code>convergence_qmesh.csv</code>, <code>solver_terms.csv</code>.</p>

<h3>2.1 Force constants: where FC2 and FC3 come from</h3>
<div class="tbl"><table class="lt">
<thead><tr><th></th><th class="r">{L[m4]}</th><th class="b">{L[m5]}</th></tr></thead>
<tbody>
<tr><td>DFT</td><td colspan="2" style="text-align:center">VASP 6.4.3, PBE + D3 (IVDW = 12), U(Mn 3d) = 4 eV (Dudarev), ferromagnetic, ENCUT 410 eV,
PREC = Accurate, EDIFF = 10⁻⁷, LREAL = .FALSE., KSPACING 0.038 Å⁻¹; NAC off</td></tr>
<tr><td>relaxed structure</td><td>one relaxed 7-atom R-3m cell, used for both FC2 and FC3</td><td>one relaxed 9-atom R-3m cell</td></tr>
<tr><td><b>FC2</b> (harmonic)</td><td><b>4×4×2 supercell</b>, 224 atoms; 7 symmetry-reduced displacements of 0.02 Å (phonopy)</td>
<td rowspan="2"><b>one 4×4×1 supercell for both</b>, 144 atoms; 1 equilibrium + 54 displaced frames, fitted together with pheasy (ALASSO)</td></tr>
<tr><td><b>FC3</b> (anharmonic)</td><td><b>4×4×1 supercell</b>, 112 atoms; 45 random displacements + 1 equilibrium frame, pheasy (RFE-OLS); cutoff 6.45 Å</td></tr>
<tr><td>FC3 cutoff used</td><td><b>6.45 Å</b></td><td><b>6.51 Å</b></td></tr>
<tr><td>check of FC2</td><td>interlayer TA velocity 1348 m/s vs 1338 m/s from the DFT elastic constants (+0.75 %)</td><td>interlayer TA 1387 m/s, C₄₄ 9.5–10.3 GPa</td></tr>
<tr><td>κ code</td><td colspan="2" style="text-align:center">phono3py (3.24.0 and 4.4.0; same population term in both)</td></tr>
</tbody></table></div>
<ul>
<li><b>For {L[m4]}, FC2 and FC3 come from two different supercells of the same crystal.</b> Both use the identical relaxed
cell (same lattice vectors and atomic positions) and the identical DFT input; the INCAR files differ only in the MAGMOM
count, and the k-point density is the same (2×2×1 on 4×4×2 and 2×2×2 on 4×4×1 sample the same k-spacing). Only the size
of the supercell along c differs.</li>
<li><b>Why two supercells.</b> In the 4×4×1 cell (one unit cell, 13.2 Å along c) each atom meets its own periodic
image across the van der Waals gap, which made the interlayer FC2 far too soft: interlayer TA 323 m/s instead of 1338 m/s.
The 4×4×2 FC2 fixes this (1348 m/s). Repeating the anharmonic dataset in 4×4×2 would cost about two days of DFT and does
not raise the usable FC3 cutoff, so FC3 stays on 4×4×1.</li>
<li><b>How they are combined.</b> phono3py takes FC2 and FC3 on different supercells
(<code>phonon_supercell_matrix</code> = 4×4×2 for FC2, <code>supercell_matrix</code> = 4×4×1 for FC3). ShengBTE cannot,
so all {L[m4]} κ values are from phono3py.</li>
<li>{L[m5]} needed no second cell: its 4×4×1 FC2 already gives the correct interlayer TA velocity.</li>
</ul>

<h3>2.2 FC3 cutoff</h3>
<div class="tbl"><table>
<thead><tr><th>material</th><th>cutoff (Å)</th><th>in-plane κ</th><th>step</th><th>κ_zz</th><th>step</th></tr></thead>
<tbody>
{cut_rows(m4)}
{cut_rows(m5)}
</tbody></table></div>
<p class="cap">300 K, RTA population term, natural isotopes, κ_zz = mesh-converged term; meshes 14³ ({L[m4]}) and
16³ ({L[m5]}). <b>Chosen: 6.45 Å and 6.51 Å</b>, the largest cutoffs the 4×4×1 FC3 cells allow under the conservative
supercell rule. By the 10 % criterion all last steps pass. For {L[m4]} the steps also shrink (in-plane −22.9 % → −8.7 %,
κ_zz −50.5 % → −8.2 %), so the remaining error is estimated at about −10 %. For {L[m5]} (in-plane −3.9 % → −4.3 %,
κ_zz −11.9 % → −9.9 %) refits at 6.67–7.43 Å on the same 54 frames were tested by hold-out validation: they are
over-parameterised from 6.90 Å on, so their lower κ is not physics. Fitting noise lowers κ, so the 6.90 Å fit bounds the
next step: <b>no more negative than −2.1 % in-plane and −3.3 % in κ_zz</b> (<code>test/Mn2In2Se5/cutoff_cv/</code>).
Larger cutoffs cannot be tested with the existing frames.</p>

<h3>2.3 q-mesh</h3>
<div class="tbl"><table>
<thead><tr><th>material</th><th>mesh</th><th>in-plane κ</th><th>step</th><th>κ_zz total</th><th>step</th><th>κ_zz converged term</th><th>step</th></tr></thead>
<tbody>
{mesh_rows(m4)}
{mesh_rows(m5)}
</tbody></table></div>
<p class="cap">300 K, RTA, isotropic N³ meshes ({L[m4]} isotopes on; {L[m5]} isotopes off, which changes κ by −1.2 %
in-plane and −0.4 % κ_zz). <b>In-plane κ is converged</b> (steps below 3 % from 14³ on). <b>The total κ_zz keeps
falling</b> by 12–15 % per step: that is a numerical artifact, because the single q-point next to Γ on the Γ–Z line
carries 15–50 % of κ_zz on these meshes and its weight shrinks as N⁻³. Leaving that point out gives the
<b>converged term</b>, which changes by less than 1 % per step from 14³ ({L[m4]}) and 24³ ({L[m5]}).
<b>Chosen</b>: κ_zz = converged term, plateau over 14³–30³ ({L[m4]}, 0.648 ± 0.003) and 24³–32³
({L[m5]}, 0.977 ± 0.008 with isotopes); in-plane {L[m4]} from LBTE at 14³ with the +1.3 % mesh correction measured on
the RTA series, {L[m5]} from RTA at 32³. The κ(T) curves use 14³ and 16³, and the VCA comparison uses 20³ for both
(within 3 % of the converged values).</p>

<h3>2.4 Solvers: RTA, iterative (LBTE) and Wigner</h3>
<ul>
<li><b>RTA</b> (relaxation-time approximation) gives each phonon its own lifetime and treats every scattering event as
resistive. It is cheap and gives the mode-resolved data used for the curves, the mesh and cutoff scans and the
disorder models.</li>
<li><b>LBTE</b> (phono3py's direct solution of the linearised Boltzmann equation, equivalent to the converged iterative
solution of ShengBTE) solves the full collision matrix, so momentum-conserving Normal processes do not destroy the heat
current. In-plane it is higher than RTA by <b>{100 * (r4 - 1):.0f} %</b> for {L[m4]} (14³: {sv[(m4, 'inplane_LBTE_14')]:.4f} vs
{sv[(m4, 'inplane_RTA_14')]:.4f}) and <b>{100 * (r5 - 1):.1f} %</b> for {L[m5]} (10³: {sv[(m5, 'inplane_LBTE_10')]:.4f} vs
{sv[(m5, 'inplane_RTA_10')]:.4f}). For κ_zz, LBTE and RTA agree to 4·10⁻⁵. <b>The quoted in-plane values are therefore
LBTE-based; κ_zz uses RTA.</b></li>
<li><b>Wigner</b> (Simoncelli, Marzari, Mauri, Nat. Phys. 15, 809 (2019)) splits κ into the usual population
(particle-like) term and a <b>coherence</b> (wave-like) term between modes whose frequencies lie closer than their
linewidths. It matters in low-κ, many-branch crystals like these. The coherence term was computed in the RTA (Wigner-RTA)
and is added to the population value. Its size at 300 K: {L[m4]} in-plane {sv[(m4, 'coherence_inplane')]:.3f},
κ_zz {sv[(m4, 'coherence_kzz')]:.3f}; {L[m5]} in-plane {sv[(m5, 'coherence_inplane')]:.3f}, κ_zz {sv[(m5, 'coherence_kzz')]:.3f}
W m⁻¹ K⁻¹, i.e. +2–5 % in-plane and +10–12 % κ_zz.</li>
<li><b>Version rule for the coherence term.</b> phono3py 4.4.0 with <code>SMM19</code> gives a coherence term 12–17× too
large for {L[m4]}, because unregularised soft modes near Γ make it diverge; 3.24.0 and the 4.x <code>IBDB19</code> option
agree with each other. {L[m4]} therefore uses the 3.24.0 value. The {L[m5]} 4.4.0 value passes the same soft-mode test
(−0.3 % in κ_zz when modes below 0.3 THz are removed) and is used.</li>
</ul>
<div class="tbl"><table class="lt">
<thead><tr><th>300 K value</th><th class="r">{L[m4]}</th><th class="b">{L[m5]}</th></tr></thead>
<tbody>
<tr><td>in-plane κ</td><td>LBTE {sv[(m4, 'inplane_LBTE_14')]:.3f} × {sv[(m4, 'mesh_correction_inplane')]:.3f} (mesh) + coherence {sv[(m4, 'coherence_inplane')]:.3f} = <b>{q4:.2f}</b></td>
<td>RTA {sv[(m5, 'inplane_RTA_32')]:.3f} × {r5:.3f} (LBTE/RTA) + coherence {sv[(m5, 'coherence_inplane')]:.3f} = <b>{q5:.2f}</b></td></tr>
<tr><td>κ_zz</td><td>converged term {sv[(m4, 'population_kzz_rest')]:.3f} + coherence {sv[(m4, 'coherence_kzz')]:.3f} = <b>{z4:.2f}</b></td>
<td>converged term {sv[(m5, 'population_kzz_rest')]:.3f} + coherence {sv[(m5, 'coherence_kzz')]:.3f} = <b>{z5:.2f}</b></td></tr>
</tbody></table></div>
<p class="cap">These are the values of §1. The curves in §4–§5 are RTA population only, which is why the §1 values sit
above them.</p>"""


def checks_html():
    L, c = HTML_LABEL, CHK
    lc, le = c["lattice"]["calc"], c["lattice"]["exp"]
    vc = c["velocity_calc"]
    e1, e2 = c["velocity_exp"]
    pct = lambda a, b: f"{100 * (a / b - 1):+.1f} %"
    raman = ", ".join(f"{m['irrep']} {m['cm-1']:.0f}" for m in c["raman_calc"])
    km, kc = c["kmin_exp_velocities"], c["kmin_calc_velocities"]
    return f"""<h2><span class="n">7</span>Independent checks against experiment ({L['MnIn2Se4']})</h2>
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


def scattering_x():
    """Antisite fraction of the VCA scattering data (plot_data/scattering_vca_x.json) and, for
    Mn2In2Se5, whether it is still the x = 0.13 archive instead of the x = 0.19 of the curves."""
    sx = json.loads((DATA / "scattering_vca_x.json").read_text())
    return sx, abs(sx["Mn2In2Se5"] - 0.19) > 1e-9


def html(figs):
    L = HTML_LABEL
    SX, SX_OLD = scattering_x()
    T1, T2 = table1_body({"quoted": "calculation, value to quote", "ordered": "calculation, ordered (RTA 20³)",
                          "vca": "calculation, VCA antisite disorder (RTA 20³)",
                          "exp": "experiment, Perp / Para (303 K)",
                          "exp_avg": "experiment, (2·Perp + Para)/3 (303 K)"},
                         {"MnIn2Se4": "cutoff ≈ −10 %", "Mn2In2Se5": "cutoff ≥ −3 %"})
    ORD, O = ordering_rows({"exp": "experiment, Perp / Para (303 K)",
                            "quoted": "calculation, value to quote, in-plane / κ_zz",
                            "ordered": "calculation, ordered RTA 20³, in-plane / κ_zz",
                            "mixed": "calculation, value to quote mixed with the texture s → Perp / Para"})
    sw, sA, sB = O["sweep"], O["s"]["MnIn2Se4"], O["s"]["Mn2In2Se5"]
    d = lambda m, x0, x1, i: abs(100 * (sw[m][x1][i] / sw[m][x0][i] - 1))
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

<h2><span class="n">1</span>Values at 300 K vs experiment</h2>
<p class="lead">W m⁻¹ K⁻¹, calculation at 300 K, experiment at 303 K (hot-pressed, κ_tot). Each direction is compared
separately; the number after each value is calculation / experiment. Experiment: Perp is measured perpendicular to the
pressing axis and is mostly in-plane, Para is measured along it and is mostly cross-plane (the platelet c axes tend to
align with the pressing direction).</p>
<div class="tbl"><table>
<thead><tr><th></th><th class="r">{L['MnIn2Se4']} in-plane</th><th class="r">{L['MnIn2Se4']} κ_zz</th>
<th class="b">{L['Mn2In2Se5']} in-plane</th><th class="b">{L['Mn2In2Se5']} κ_zz</th></tr></thead>
<tbody>
{T1}
</tbody></table></div>
<p class="cap"><b>Calculation, value to quote</b>: the most accurate single-crystal values (§2.4: in-plane LBTE + coherence,
κ_zz mesh-converged + coherence). {L['MnIn2Se4']}: FC3 cutoff error about −10 %. {L['Mn2In2Se5']}: next cutoff step no more negative
than −2.1 % / −3.3 % (in-plane / κ_zz). <b>Ordered / VCA</b>: RTA at 20³, without / with Mn–In antisite mass disorder (§5); they show how much
disorder alone lowers κ. All calculated values are for a defect-free single crystal without grain boundaries.</p>

<p class="lead"><b>Why a second table.</b> The pellets are only partly textured, so Perp is not pure in-plane and Para is
not pure κ_zz: part of the in-plane conduction leaks into Para (this is why Para / κ_zz ratios above look smaller).
The average (2·Perp + Para)/3 equals (2·κ_in-plane + κ_zz)/3 for <i>any</i> texture, so it compares calculation and
experiment without knowing the texture.</p>
<div class="tbl"><table>
<thead><tr><th>direction average (2·in-plane + κ_zz)/3</th><th class="r">{L['MnIn2Se4']}</th><th class="b">{L['Mn2In2Se5']}</th></tr></thead>
<tbody>
{T2}
</tbody></table></div>
<p class="cap">The single-crystal κ is 3–4× above the hot-pressed samples. Antisite disorder (VCA) closes part of the
gap; the rest needs a nanoscale effective mean free path (L_eff ≈ 10–50 nm, §8); the possible reasons are discussed in §9. The L3 model of §8 is not listed here:
it is fitted to these experimental points, so at 300 K it reproduces them by construction.</p>

<h3>Which material is lower: {L['Mn2In2Se5']} relative to {L['MnIn2Se4']}</h3>
<p class="lead">In the experiment {L['Mn2In2Se5']} is lower in both directions. The calculation reproduces this
in-plane and in the direction average, but not for κ_zz. ✓ / ✗: same / opposite sign as the experiment.</p>
<div class="tbl"><table>
<thead><tr><th>300 K, {L['Mn2In2Se5']} vs {L['MnIn2Se4']}</th><th>in-plane / Perp</th><th>κ_zz / Para</th><th>direction average</th></tr></thead>
<tbody>
{ORD}
</tbody></table></div>
<ul>
<li><b>κ_zz should not be compared with Para directly.</b> Para contains part of the in-plane conduction,
Para = s·κ_in-plane + (1−s)·κ_zz. With the same s for both materials, {L['Mn2In2Se5']} has the lower Para as soon as
s &gt; {O['s_cross']:.2f}. The fitted texture is s = {sA:.2f} ({L['MnIn2Se4']}) and {sB:.2f} ({L['Mn2In2Se5']}); with it
the calculated values give the experimental order in both directions (last row: Para
{O['mixed_abs']['Mn2In2Se5'][1]:.2f} vs {O['mixed_abs']['MnIn2Se4'][1]:.2f}). XRD of both pellet faces would measure s independently.</li>
<li><b>The {L['Mn2In2Se5']} κ_zz is bounded.</b> Its FC3-cutoff steps (−11.9 %, −9.9 %) are followed by a step of at
most −3.3 % (§2.2, hold-out test), so its κ_zz stays about 50 % above {L['MnIn2Se4']}'s within what the data can test.</li>
<li><b>Point defects hardly lower κ_zz.</b> Antisite mass disorder from the ordered to the random limit lowers the
in-plane κ by {d('Mn2In2Se5', 0.0, 0.5, 0):.0f} % but κ_zz by only {d('Mn2In2Se5', 0.0, 0.5, 1):.0f} % in {L['Mn2In2Se5']}
({L['MnIn2Se4']}, x = 0 → 0.44: by {d('MnIn2Se4', 0.0, 0.44, 0):.0f} % and {d('MnIn2Se4', 0.0, 0.44, 1):.0f} %;
<code>plot_data/disorder_sweep_300K.csv</code>). κ_zz is carried by low-frequency interlayer acoustic modes, which mass
disorder barely scatters. A lower κ_zz in the {L['Mn2In2Se5']} pellets therefore points to defects along c, such as
stacking faults, polytype intergrowths or secondary phases from its decomposition between 400 and 800 °C (§9, items 1 and 4).</li>
</ul>

{methods_html()}

<h2><span class="n">3</span>Phonon dispersion</h2>
<p class="lead">R-3m primitive cells, path Γ–T–H₂|H₀–L–Γ–S₀|S₂–F–Γ. No imaginary modes in either material; NAC off.</p>
<figure><img alt="Phonon dispersion of both materials" src="{figs['bands']}"></figure>
<p class="cap"><b>{L['MnIn2Se4']}</b>: 4×4×2 FC2, <code>test/MnIn2Se4/fc2_442_phonon/band-dft-cpu.yaml</code>.
<b>{L['Mn2In2Se5']}</b>: 4×4×1 FC2, <code>Mn2In2Se5/fc/fc2.hdf5</code> (phonopy interpolation).
The 4×4×1 FC2 of {L['MnIn2Se4']} gave a too-soft interlayer TA (322.9 m/s) and is not used; the 4×4×2 FC2 gives
1348.3 m/s. The {L['Mn2In2Se5']} 4×4×1 FC2 was checked the same way and is fine (1386.8 m/s).</p>

<h2><span class="n">4</span>Intrinsic κ(T)</h2>
<p class="lead">Ordered crystal, 3-phonon + natural isotopes, RTA. κ_zz is the mesh-converged term (Γ–Z quadrature
error removed). Stars: the 300 K values of §1.</p>
<figure><img alt="Intrinsic kappa versus temperature" src="{figs['intrinsic']}"></figure>
<p class="cap"><b>Lines</b>: RTA population term, computed at every temperature (meshes 14³ for {L['MnIn2Se4']} and 16³ for
{L['Mn2In2Se5']}). <b>Stars</b>: the 300 K values of §1, which use the more accurate treatment of §2.4 (in-plane LBTE
instead of RTA, plus the Wigner coherence term for both components). That treatment was run at 300 K only, so it
appears as a single point. The stars lie above the lines at 300 K by +28 % (in-plane) and +11 % (κ_zz) for
{L['MnIn2Se4']}, and by +14 % and +8 % for {L['Mn2In2Se5']}; the same correction applies to the curves at other
temperatures in approximate form. Quote numbers from §1, not read off the lines. Bar: FC3-cutoff estimate (≈ −10 %).
{L['Mn2In2Se5']}: bound on the next cutoff step (−2.1 % / −3.3 %) from the hold-out test.</p>

<h2><span class="n">5</span>Calculation vs hot-pressed experiment</h2>
<p class="lead">One panel per material and direction. Dashed: ordered single crystal. Solid: VCA antisite-disordered
single crystal. Dotted: L3 fit (§8). Dash-dot ({L['MnIn2Se4']} only): minimum κ from the measured sound velocities (§7). Filled symbols: 300–623 K (used in the fit). Open symbols: 55–145 K (not fitted).</p>
<figure><img alt="Calculated kappa against the hot-pressed experiment" src="{figs['calc_vs_exp']}"></figure>
<p class="cap">Ordered and VCA: RTA, 20³ (same mesh for each pair), Perp ↔ in-plane, Para ↔ κ_zz, no grain boundary
and no texture, so they are upper references. VCA occupancy: {L['MnIn2Se4']} i = 0.44 (single-crystal refinement, Range et al., Z. Naturforsch.
46b, 1122 (1991)), {L['Mn2In2Se5']} x = 0.19(2) (powder Rietveld of the R-3m phase, Williams et al., J. Mater. Chem. C 12, 1753 (2024), ESI Table S7; see §8). VCA lowers κ strongly but stays
well above experiment and cannot produce the low-temperature maximum near 100 K; that needs a nanoscale mean free path.</p>

<h2><span class="n">6</span>Scattering rates at 300 K</h2>
<p class="lead">Largest FC3 cutoff, 20³ mesh, τ⁻¹ = 2Γ. Baseline: real masses + natural isotopes. VCA: averaged
site masses + mass-variance scattering. Top row {L['MnIn2Se4']} (i = {SX['MnIn2Se4']:g}), bottom row {L['Mn2In2Se5']} (x = {SX['Mn2In2Se5']:g}{
': the per-mode file of the x = 0.19 run is not yet in the repository; the κ values below are at x = 0.19' if SX_OLD else ''}).</p>
<figure><img alt="Three-phonon and disorder scattering rates, both materials" src="{figs['scattering']}"></figure>
<p class="cap">Left: the 3-phonon rates of the ordered and VCA crystals agree within a few tens of percent, because VCA
changes the masses, not the anharmonicity. Right: VCA mass-variance scattering is 2–3 orders of magnitude stronger
than natural-isotope scattering, rises as ω² at low frequency and flattens at high frequency. Effect on κ at 300 K:
{L['MnIn2Se4']} in-plane 3.69 → 1.49 (2.48×), κ_zz 0.638 → 0.445 (1.43×); {L['Mn2In2Se5']} in-plane
3.00 → 1.46 (2.05×), κ_zz 0.984 → 0.753 (1.31×) at x = 0.19 (x = 0.13: 1.63 and 0.792). The weaker effect in {L['Mn2In2Se5']}
comes from its less mixed sites (81 % pure at x = 0.19, against 56/44 on the {L['MnIn2Se4']} octahedral site).</p>

{checks_html()}

<h2><span class="n">8</span>What the L3 fit is</h2>
<p class="lead">L3 is the extrinsic model of <code>analysis/extrinsic/</code>. <b>Yes, it is fitted to the
experimental values</b>, so its agreement with experiment is by construction and is not a prediction.
At 300 K it gives Perp / Para = 0.92 / 0.49 ({L['MnIn2Se4']}, measured 0.93 / 0.49) and 0.78 / 0.34 ({L['Mn2In2Se5']},
measured 0.74 / 0.32). <b>Its purpose is not to supply κ</b> (the measurement already does) but to answer how much extra
scattering, and of what kind, turns the calculated single crystal into the measured pellet: the fitted parameters
(texture s and boundary length L; the antisite fraction x is taken from structure refinements) are the result, and the 55–145 K points, which were not fitted, test it.</p>
<div class="eq">Γ_λ = Γ_3ph,λ(T)  +  Σ g_elem(x)·B_elem,λ  +  |v_λ| / (4π L)
       3-phonon           Mn↔In antisite mass disorder   grain boundary (Casimir)

κ_Para = s·κ_xx + (1−s)·κ_zz          κ_Perp = (1−s/2)·κ_xx + (s/2)·κ_zz</div>
<ul>
<li><b>Inputs that are not fitted</b>: the first-principles 3-phonon rates, phonon velocities and heat capacities
(ordered masses), and the mass-disorder basis rates B, computed once per material on a 24³ mesh.</li>
<li><b>Parameters</b>: x, the fraction of Mn on In sites, <b>fixed</b> at the literature value (below); L, an effective
boundary mean free path, and s = ⟨sin²θ⟩, the texture (0 = fully textured, 2/3 = random), <b>fitted</b>.</li>
<li><b>Data used</b>: only the 300–623 K points of <code>experiment/kappa_tot_digitized.csv</code>, both directions.
The fit minimises the rms of ln(model/experiment) over L, with s optimised at each point.
The 55–145 K points are not used; they only test the model.</li>
<li><b>Result</b> (analysis/disorder/README.md §5.1): {L['MnIn2Se4']} with x = 0.44 (single-crystal XRD, Range et al.,
Z. Naturforsch. 46b, 1122 (1991)): L = 30.7 nm, s = 0.35, rms 1.2 %. {L['Mn2In2Se5']} with x = 0.19(2) (powder Rietveld of the
R-3m phase, Williams et al., J. Mater. Chem. C 12, 1753 (2024), ESI Table S7): L = 22.1 nm, s = 0.23, rms 3.4 %. So L ≈ 20–30 nm
in both materials. (The 0.704/0.296 occupancy of Chernoukhov et al., Molecules 30, 1904 (2025) belongs to a different,
quenched polymorph, P3̄m1, and is not used.)</li>
<li><b>Why x is fixed</b>: with x free as well, the fit determines only s (acceptable 0.33–0.38 and 0.15–0.28). x and L
both push κ down and are degenerate: acceptable x spans 0–0.66 and 0.01–0.50, and the free optimum (x = 0.45, L = 30.7 nm,
rms 1.2 %; x = 0.13, L = 15.9 nm, rms 3.2 %) is not evidence for any occupancy. XRD (texture), SEM (grain size) and a Rietveld
occupancy refinement of our own pellets would break the degeneracy.</li>
<li><b>Not to be combined with VCA</b>: L3 already contains the antisite scattering through x, so the VCA and L3
curves are alternatives, not additive corrections.</li>
<li><b>κ_e</b> is taken as 0. If the electronic part is not negligible, the fitted extrinsic scattering is too strong.</li>
</ul>

<h2><span class="n">9</span>Why the calculation lies above the experiment</h2>
<p class="lead">At 300 K the calculated single-crystal κ is 4.3× above the hot-pressed κ_tot for both materials
(direction average, values to quote, §1). Antisite disorder (VCA, at the literature occupancies) brings this down to 1.5× ({L['MnIn2Se4']}) and
2.0× ({L['Mn2In2Se5']}; RTA, 20³). The candidate reasons for the rest are listed
below, roughly in order of how much each can explain, with the evidence and the test. <b>These are hypotheses.</b>
Only the first two are quantified here.</p>
<div class="tbl"><table class="lt">
<thead><tr><th>#</th><th>mechanism</th><th>possible effect on κ</th><th>evidence</th><th>how to test</th></tr></thead>
<tbody>
<tr><td>1</td><td><b>Planar defects and grain boundaries</b>: a nanoscale effective mean free path</td>
<td>largest; closes the remaining gap with L ≈ 20–30 nm (L3, x fixed from the literature, §8)</td>
<td>below 150 K the measured κ is flat, with no 1/T rise, i.e. a temperature-independent scattering dominates;
the required L is far below typical hot-pressed grain sizes (0.5–5 μm), which points to stacking faults, twins or
polytype intergrowths, common in these layered "124"/"225" compounds; κ_tot lies between the Cahill and diffuson limits (§7)</td>
<td>SEM grain size; TEM/HAADF for stacking faults</td></tr>
<tr><td>2</td><td><b>Mn/In antisite disorder</b></td>
<td>VCA at the literature x: in-plane 2.5× ({L['MnIn2Se4']}) / 2.1× ({L['Mn2In2Se5']}), κ_zz 1.4× / 1.3× at 300 K</td>
<td>single-crystal refinement i = 0.44 (Range 1991); Rietveld x = 0.19(2) (Williams 2024); broad Raman modes in
{L['Mn2In2Se5']} attributed to site disorder (Williams 2024)</td>
<td>Rietveld occupancy of our own pellets. VCA counts only the mass difference, so it is a lower bound: bond-strength and
local-strain disorder are missing (a disordered supercell would include them). LBTE + coherence for VCA is planned
(<code>test/figures/kappa_vs_T/jobs_lbte_wigner/</code>)</td></tr>
<tr><td>3</td><td><b>Paramagnetic spin disorder</b> (Mn²⁺, S = 5/2)</td>
<td>not known; could be tens of per cent</td>
<td>the calculation assumes a static ferromagnetic order; the real crystals are paramagnetic at 300 K and order
(antiferromagnetic / spin glass) only below ≈ 6–7 K, with strong antiferromagnetic exchange (Weiss θ ≈ −198 K for
{L['Mn2In2Se5']}, Williams 2024). Fluctuating spins make the exchange-dependent force constants random and add
spin–phonon scattering</td>
<td>magnetic susceptibility of the pellets; force constants for antiferromagnetic or random spin configurations</td></tr>
<tr><td>4</td><td><b>Off-stoichiometry, vacancies, secondary phases</b></td>
<td>not known; cation vacancies scatter far more strongly than a mass substitution</td>
<td>{L['MnIn2Se4']} has a 5–6 mol % homogeneity range towards In₂Se₃ (Mammadov 2025), i.e. cation vacancies;
{L['Mn2In2Se5']} decomposes into MnIn₂Se₄ + MnSe between about 400 and 800 °C and MnSe admixtures are common
(Chernoukhov 2025), so hot pressing in that range may leave secondary phases and interfaces</td>
<td>XRD of the pressed pellet (phase purity), EDS composition, hot-pressing temperature and time</td></tr>
<tr><td>5</td><td><b>Porosity, pellet softening</b></td>
<td>≈ 10–20 % for 5–10 % porosity</td>
<td>measured sound velocity 20–30 % below the calculation in every direction (§7)</td>
<td>Archimedes density</td></tr>
<tr><td>6</td><td><b>Higher-order anharmonicity, temperature effects</b></td>
<td>typically −10 to −30 % at 300 K in low-κ layered selenides (not computed here)</td>
<td>low-lying optical branches (§3) and strong 3-phonon rates; the calculation uses 0 K harmonic phonons and 3-phonon scattering only</td>
<td>four-phonon scattering, temperature-dependent force constants (expensive)</td></tr>
<tr><td>7</td><td><b>FC3 cutoff</b></td>
<td>{L['MnIn2Se4']} ≈ −10 %; {L['Mn2In2Se5']} next step ≥ −3.3 %</td>
<td>§2.2</td><td>{L['Mn2In2Se5']} done (refit to 7.43 Å + hold-out test); {L['MnIn2Se4']} the same test on its 7.74 Å refits</td></tr>
</tbody></table></div>
<p class="lead"><b>Ruled out as explanations</b></p>
<ul>
<li><b>Texture</b> moves conduction between Perp and Para, but the direction average (2·Perp + Para)/3 does not depend on
it and is still 4.3× too low (§1).</li>
<li><b>Electronic κ_e</b>: the experiment is κ_tot. Subtracting κ_e lowers the experimental lattice κ and widens the gap.</li>
<li><b>Lattice constants</b>: the calculated cell is 0.6–0.7 % larger than measured (§7); a larger cell tends to lower,
not raise, the calculated κ.</li>
<li><b>LBTE and coherence</b> raise the calculated κ, so they widen the gap. <b>Digitising errors</b> (±0.005–0.01 W m⁻¹ K⁻¹) are negligible.</li>
</ul>
<p class="cap">Mechanisms 1 and 2 alone already reproduce the data (that is the L3 fit). Every further mechanism
(3–6) that is real would reduce the boundary scattering needed. The fitted L ≈ 20–30 nm lumps together all scattering
the model leaves out, so the actual spacing of planar defects is at least that large. The cheapest decisive
measurements are the pellet density, SEM/TEM of the microstructure and XRD of the pressed pellet.</p>

<h2><span class="n">10</span>Caveats</h2>
<ul>
<li>The experimental points are digitised from a figure (±0.005 W m⁻¹ K⁻¹ at high T, ±0.01 at low T) and are κ_tot
with κ_e not subtracted. Replace <code>experiment/kappa_tot_digitized.csv</code> once the raw data are available.</li>
<li>{L['Mn2In2Se5']}: the FC3 cutoff is tested only up to 6.90 Å (next step ≥ −2.1 % / −3.3 %); the existing frames
cannot support larger cutoffs. Within that, both orderings hold: in-plane {L['MnIn2Se4']} 37 % higher, κ_zz
{L['Mn2In2Se5']} about 50 % higher.</li>
<li>Between the 145 K and 300 K experimental points the Para values drop much more than Perp ({L['Mn2In2Se5']} −32 %
vs −6 %), possibly two instruments or two samples; this is why the fit uses the high-temperature range only.</li>
</ul>
</div></body></html>
"""


def main():
    if not BUNDLE and "--no-export" not in sys.argv:   # repository: refresh plot_data/ first
        import export_plot_data
        export_plot_data.main()
    png = plot_figures.make_all(DATA, FIGDIR, png_for_html=True)
    figs = {k: "data:image/png;base64," + base64.b64encode(v).decode() for k, v in png.items()}
    OUT.write_text(html(figs), encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e6:.2f} MB) and {FIGDIR}/ (.png 300 dpi + .pdf)")


if __name__ == "__main__":
    main()
