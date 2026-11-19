#!/usr/bin/env python3
"""Pack the HTML report, its figures, plot data and plotting scripts into one zip.

  python3 report/make_report_html.py          # refresh plot_data/, figures/, HTML first
  python3 report/make_bundle.py [--out DIR]   # writes DIR/<name>.zip (default: report/dist/)

Layout inside the zip:
  README.md
  report/   calc_vs_experiment_report.html  (self-contained, opens offline)
  figures/  fig1..fig4 as .png (300 dpi) and .pdf
  data/     every plotted dataset as CSV / JSON
  scripts/  plot_figures.py (data/ -> figures/, numpy + matplotlib only),
            make_report_html.py (data/ -> figures/ + HTML), and, for reference,
            export_plot_data.py / independent_checks.py (need the full repository)
"""
import argparse
import datetime
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAME = "MnIn2Se4_Mn2In2Se5_calc_vs_experiment_" + datetime.date.today().strftime("%Y%m%d")

FIGURES = {
    "fig1_phonon_dispersion": "Phonon dispersion of both materials (Γ–T–H₂|H₀–L–Γ–S₀|S₂–F–Γ)",
    "fig2_kappa_intrinsic": "Intrinsic κ(T), in-plane and κ_zz; stars = 300 K values to quote",
    "fig3_calc_vs_experiment": "Ordered / VCA / L3 calculation vs the hot-pressed experiment, 2 × 2",
    "fig4_scattering_300K": "3-phonon and disorder scattering rates at 300 K, both materials",
}
DATA = {
    "bands_<material>.csv": "fig1: one row per q-point; segment, path distance, all branch frequencies (THz)",
    "bands_<material>_ticks.csv": "fig1: high-symmetry point positions and labels",
    "kappa_intrinsic_<material>.csv": "fig2: T, κ in-plane, κ_zz (ordered, RTA, isotopes)",
    "kappa_quoted_300K.csv": "fig2: the starred 300 K values",
    "kappa_ordered_vca_<material>.csv": "fig3: ordered and VCA single-crystal κ, 20³ mesh",
    "kappa_L3_<material>.csv": "fig3: L3 extrinsic model, κ_xx, κ_zz and texture-mixed Perp / Para, with s",
    "experiment_kappa_tot_digitized.csv": "fig3: experimental κ_tot (digitised, κ_e not subtracted)",
    "experiment_sound_velocity_MnIn2Se4.csv": "report §6: measured sound velocities",
    "independent_checks.json": "report §6 and fig3: lattice, Raman, velocities, Cahill / diffuson κ",
    "scattering_<material>_300K.csv": "fig4: per-mode frequency and τ⁻¹ (3-phonon, disorder), baseline and VCA",
    "convergence_fc3_cutoff.csv": "report §2.2: κ vs FC3 cutoff, both materials",
    "convergence_qmesh.csv": "report §2.3: κ vs q-mesh (in-plane, κ_zz total, κ_zz converged term)",
    "solver_terms.csv": "report §2.4: RTA / LBTE / coherence terms behind the quoted 300 K values",
}


def readme():
    fig = "\n".join(f"| `{k}.png` / `.pdf` | {v} |" for k, v in FIGURES.items())
    dat = "\n".join(f"| `{k}` | {v} |" for k, v in DATA.items())
    return f"""# MnIn₂Se₄ / Mn₂In₂Se₅ — calculation vs experiment

Exported {datetime.date.today().isoformat()} from the repository `Mn2In2Se5-thermal`
(`report/make_bundle.py`). Every CSV starts with `#` comment lines that state its
source file and conventions.

## Folders

| folder | contents |
|---|---|
| `report/` | `calc_vs_experiment_report_zh.html` (中文版) and `calc_vs_experiment_report.html` (English): the full report, one file each, opens offline in any browser |
| `figures/` | the four report figures, `.png` at 300 dpi and vector `.pdf`; `*_zh` = Chinese labels |
| `data/` | every dataset plotted in the figures, as CSV / JSON |
| `scripts/` | the plotting scripts |

## Figures

| file | content |
|---|---|
{fig}

## Data (`<material>` = `MnIn2Se4` or `Mn2In2Se5`)

| file | used in / content |
|---|---|
{dat}

Units: κ in W m⁻¹ K⁻¹, frequency in THz, τ⁻¹ in ps⁻¹, T in K.

## Redrawing

Needs Python 3 with numpy and matplotlib only.

```bash
cd scripts
python3 plot_figures.py            # data/ -> figures/*.png, *.pdf (English and *_zh Chinese labels)
python3 make_report_html.py        # data/ -> figures/ and report/calc_vs_experiment_report.html
python3 make_report_html_zh.py     # data/ -> figures/*_zh and report/calc_vs_experiment_report_zh.html
```

To change a figure, edit the matching `fig_*` function in `scripts/plot_figures.py`
(Chinese labels: the `ZH` dictionary there; they need a CJK font such as Microsoft YaHei,
SimHei, PingFang SC or Noto Sans CJK SC).
`export_plot_data.py` and `independent_checks.py` are included to show how `data/`
was produced; they need the full repository (force constants, phono3py output)
and phonopy, so they do not run inside this bundle. Links in the HTML that point to
repository files (e.g. the progress report) only work inside the repository.
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "dist"))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    zpath = out / f"{NAME}.zip"
    items = [("README.md", readme())]
    items += [(f"report/{p.name}", p) for p in [HERE / "calc_vs_experiment_report.html",
                                                 HERE / "calc_vs_experiment_report_zh.html"]]
    items += [(f"figures/{p.name}", p) for p in sorted((HERE / "figures").glob("fig*.*"))]
    items += [(f"data/{p.name}", p) for p in sorted((HERE / "plot_data").iterdir())]
    items += [(f"scripts/{n}", HERE / n) for n in
              ("plot_figures.py", "make_report_html.py", "make_report_html_zh.py", "export_plot_data.py",
               "independent_checks.py")]
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for arc, src in items:
            arc = f"{NAME}/{arc}"
            if isinstance(src, str):
                z.writestr(arc, src)
            else:
                z.write(src, arc)
    print(f"wrote {zpath} ({zpath.stat().st_size / 1e6:.2f} MB, {len(items)} files)")


if __name__ == "__main__":
    main()
