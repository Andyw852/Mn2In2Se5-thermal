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
κ_zz −50.5 % → −8.2 %), so the remaining error is estimated at about −10 %. For {L[m5]} the steps do not shrink (in-plane
−3.9 % → −4.3 %, κ_zz −11.9 % → −9.9 %), so its values are given as <b>upper bounds</b>.</p>

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
<tr><td>&nbsp;&nbsp;calc., L3 fit (fitted to these data, §8)</td><td>0.77</td><td>0.62</td></tr>
</tbody></table></div>
<p class="cap">The intrinsic single-crystal κ is 3–4× above the hot-pressed samples. Antisite mass disorder (VCA)
closes part of the gap; the rest needs a nanoscale effective mean free path (L_eff ≈ 10–50 nm).</p>

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
Arrow: {L['Mn2In2Se5']} cutoff not converged, value is an upper bound.</p>

<h2><span class="n">5</span>Calculation vs hot-pressed experiment</h2>
<p class="lead">One panel per material and direction. Dashed: ordered single crystal. Solid: VCA antisite-disordered
single crystal. Dotted: L3 fit (§8). Dash-dot ({L['MnIn2Se4']} only): minimum κ from the measured sound velocities (§7). Filled symbols: 300–623 K (used in the fit). Open symbols: 55–145 K (not fitted).</p>
<figure><img alt="Calculated kappa against the hot-pressed experiment" src="{figs['calc_vs_exp']}"></figure>
<p class="cap">Ordered and VCA: RTA, 20³ (same mesh for each pair), Perp ↔ in-plane, Para ↔ κ_zz, no grain boundary
and no texture, so they are upper references. VCA occupancy: {L['MnIn2Se4']} i = 0.44 (single-crystal refinement, Range et al., Z. Naturforsch.
46b, 1122 (1991)), {L['Mn2In2Se5']} x = 0.13 (taken from the L3 fit; the literature value for this R-3m phase is 0.19(2), see §8). VCA lowers κ strongly but stays
well above experiment and cannot produce the low-temperature maximum near 100 K; that needs a nanoscale mean free path.</p>

<h2><span class="n">6</span>Scattering rates at 300 K</h2>
<p class="lead">Largest FC3 cutoff, 20³ mesh, τ⁻¹ = 2Γ. Baseline: real masses + natural isotopes. VCA: averaged
site masses + mass-variance scattering. Top row {L['MnIn2Se4']} (i = 0.44), bottom row {L['Mn2In2Se5']} (x = 0.13).</p>
<figure><img alt="Three-phonon and disorder scattering rates, both materials" src="{figs['scattering']}"></figure>
<p class="cap">Left: the 3-phonon rates of the ordered and VCA crystals agree within a few tens of percent, because VCA
changes the masses, not the anharmonicity. Right: VCA mass-variance scattering is 2–3 orders of magnitude stronger
than natural-isotope scattering, rises as ω² at low frequency and flattens at high frequency. Effect on κ at 300 K:
{L['MnIn2Se4']} in-plane 3.69 → 1.49 (2.48×), κ_zz 0.638 → 0.445 (1.43×); {L['Mn2In2Se5']} in-plane
3.00 → 1.63 (1.84×), κ_zz 0.984 → 0.792 (1.24×). The weaker effect in {L['Mn2In2Se5']} comes from its smaller x.</p>

{checks_html()}

<h2><span class="n">8</span>What the L3 fit is</h2>
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
the closeness of 0.45 to the single-crystal 0.44 is a coincidence of the fit. XRD (texture), SEM (grain size) and a Rietveld
occupancy refinement would break the degeneracy.</li>
<li><b>Fixing x from the literature</b> (analysis/disorder/README.md §5.1): {L['MnIn2Se4']} x = 0.44 (single-crystal XRD,
Range et al., Z. Naturforsch. 46b, 1122 (1991)) gives L = 30.7 nm, rms 1.2 %. {L['Mn2In2Se5']} x = 0.19(2) (powder Rietveld of the
R-3m phase, Williams et al., J. Mater. Chem. C 12, 1753 (2024), ESI Table S7) gives L = 22.1 nm, rms 3.4 %. With x fixed, only L and
s are fitted, and L ≈ 20–30 nm in both materials. (The 0.704/0.296 occupancy of Chernoukhov et al., Molecules 30, 1904
(2025) belongs to a different, quenched polymorph, P3̄m1, and is not used.)</li>
<li><b>Not to be combined with VCA</b>: L3 already contains the antisite scattering through x, so the VCA and L3
curves are alternatives, not additive corrections.</li>
<li><b>κ_e</b> is taken as 0. If the electronic part is not negligible, the fitted extrinsic scattering is too strong.</li>
</ul>

<h2><span class="n">9</span>Caveats</h2>
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
    if not BUNDLE and "--no-export" not in sys.argv:   # repository: refresh plot_data/ first
        import export_plot_data
        export_plot_data.main()
    png = plot_figures.make_all(DATA, FIGDIR, png_for_html=True)
    figs = {k: "data:image/png;base64," + base64.b64encode(v).decode() for k, v in png.items()}
    OUT.write_text(html(figs), encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e6:.2f} MB) and {FIGDIR}/ (.png 300 dpi + .pdf)")


if __name__ == "__main__":
    main()
