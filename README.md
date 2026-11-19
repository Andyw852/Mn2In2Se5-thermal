# MnIn₂Se₄ / Mn₂In₂Se₅ lattice thermal conductivity: first-principles calculation and comparison with experiment

Lattice thermal conductivity of two R-3m layered selenides (van der Waals
interlayer bonding). The calculations give intrinsic single-crystal values;
the experiments are on hot-pressed polycrystals, so the two differ by texture,
antisite disorder, grain boundaries and electronic contribution.
`analysis/extrinsic/` puts those into a single model and fits it to experiment.

## Latest progress report

**[`PROGRESS_REPORT_20261118.md`](PROGRESS_REPORT_20261118.md)** — values to quote, the three paper figures, the uncertainty budget and the open items.

**[`report/calc_vs_experiment_report.html`](report/calc_vs_experiment_report.html)** — one static HTML page (open it in a browser): phonon dispersion, κ(T) vs experiment, and the 300 K scattering rates for both materials, plus what the L3 fit is. Rebuild with `python3 report/make_report_html.py`; the figures (`report/figures/`, PNG + PDF) and their data (`report/plot_data/`, CSV) are written alongside, and `python3 report/make_bundle.py` packs HTML, figures, data and scripts into one zip.

## Status (2026-11-15)

| | MnIn₂Se₄ | Mn₂In₂Se₅ |
|---|---|---|
| Calculation protocol | AutoZT phonon-dft-cpu: PBE-D3 + U(Mn) = 4 eV + FM; FC2 in a **4×4×2** supercell (224 atoms, symmetry-reduced finite displacement, 7 frames); FC3 reused from the earlier 4×4×1 fit (cutoff 6.45 Å; a 4×4×1 cell caps it at 6.507 Å) | AutoZT kl-dft-cpu: PBE-D3 + U(Mn) = 4 eV + FM, 4×4×1 rebased cell (144 atoms), cutoff 6.51 Å |
| κ_L(300 K) xx / zz (W/mK) | **≈ 4.6 / ≈ 0.71** (in-plane: LBTE + coherence; κ_zz: mesh-converged 0.646 + coherence 0.066; as 30³ total + coherence ≈ 0.83). FC3 cutoff ≈ −10 %, estimated. See [report §9.1](REPORT_MnIn2Se4_20261115.md) | ~~2.059 / 0.209~~ (ShengBTE, FC3 export incomplete — invalid) → **≈ 3.37 / ≈ 1.09** (in-plane: RTA 3.05 at 32³, iterative and coherence corrections; κ_zz: mesh-converged 0.977 + coherence 0.112; as 30³ total + coherence ≈ 1.36). **κ_zz FC3 cutoff not converged** (last steps −11.8 %, −9.9 %). See [report §9.1](REPORT_MnIn2Se4_20261115.md) |
| Experimental κ_tot(303 K) Perp / Para | 0.933 / 0.485 | 0.741 / 0.320 |
| Status | FC2 fixed (interlayer TA 323 → 1348 m/s, now within 0.75 % of DFT elastic); FC3 cutoff and κ_zz mesh convergence still open | Force constants reliable; ShengBTE export bug fixed in AutoZT [0127]; S6 to be rerun on a converged κ_zz mesh (see `analysis/shengbte_crosscheck/`) |

## Layout

| Path | Contents |
|---|---|
| [`Mn2In2Se5/`](Mn2In2Se5/) | Production results: κ(T), phonons, cutoff scan, figures, provenance; `fc/` holds force constants |
| [`MnIn2Se4/`](MnIn2Se4/) | Existing data package (39 fit sets, raw κ output, scripts, reports); `redo_441/` holds the structure, k-points and instructions for the recomputation; `archive_2026-10-28_fit-fc-thermal/` is an earlier fit (including the ShengBTE/phono3py discrepancy diagnosis) |
| [`experiment/`](experiment/) | Experimental figures, digitised CSV, and the list of still-missing experimental data |
| [`analysis/extrinsic/`](analysis/extrinsic/) | Model, scripts and fit results for intrinsic κ + antisite mass disorder + grain boundaries + texture |
| [`analysis/shengbte_crosscheck/`](analysis/shengbte_crosscheck/) | ShengBTE vs phono3py cross-check, evidence for the FC3 export bug, and q-mesh convergence |
| [`protocol/`](protocol/) | The AutoZT kl-dft-cpu 4×4×1 rebased-cell protocol |
| [`test/`](test/) | Reproducible inputs/outputs for the κ calculations of both materials |
| [`REPORT_MnIn2Se4_20261115.md`](REPORT_MnIn2Se4_20261115.md) | **Full report of the 4×4×2 FC2 rework and the convergence audit** |

## Connecting calculation and experiment

Hot-pressed polycrystals are not fully textured. Introduce the texture
parameter s = ⟨sin²θ⟩ (θ is the angle between a grain's c axis and the pressing
direction; 0 = fully textured, 2/3 = random):

```
κ_Para = s·κ_xx + (1−s)·κ_zz
κ_Perp = (1−s/2)·κ_xx + (s/2)·κ_zz
(2·κ_Perp + κ_Para)/3 = (2κ_xx + κ_zz)/3      ← texture-independent, directly comparable
```

Mn₂In₂Se₅ at 300 K: the corrected calculation gives (2κ_xx + κ_zz)/3 ≈ 2.3
(16×16×8 RTA) versus 0.60 experimentally — a factor of ≈3.9 (the invalid
ShengBTE value gave 1.44, a factor of 2.4). The gap shrinks with increasing
temperature, the signature of extrinsic (temperature-independent) scattering
superimposed on Umklapp, which texture alone cannot explain. Method and results
are in [`analysis/extrinsic/README.md`](analysis/extrinsic/README.md).

## Repository conventions

- Private repository: `data/provenance` files contain internal cluster paths.
- Raw DFT output (OUTCAR, vasprun.xml, WAVECAR, …) is not committed; force
  constant hdf5 files are (each < 100 MB).
- This repository supersedes `Mn2In2Se5_441-thermal` (merged into
  `Mn2In2Se5/`) and `MnIn2Se4-thermal` (merged into
  `MnIn2Se4/archive_2026-10-28_fit-fc-thermal/`).
- **All documentation in this repository is in English.**
