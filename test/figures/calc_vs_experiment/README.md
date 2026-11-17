# Calculated κ(T) over the experimental κ_tot figure

`calc_vs_experiment.png / .pdf` — the experimental points from
`experiment/kappa_tot_digitized.csv` with our calculations on top.

## Which disorder level is shown (the question that was open)

| element of the figure | level | what it contains |
|---|---|---|
| experiment (symbols) | — | hot-pressed polycrystal, two directions |
| **solid lines** | **L1** | intrinsic single crystal + **natural isotopes only** (no antisite disorder) |
| **dashed lines** | **L3** | L1 + Mn/In antisite mass disorder + grain-boundary mean free path + texture, **fitted to the same experimental CSV** |

So the figure does **not** contain the VCA (L2) numbers; L2 is shown separately in
`../scattering/` and in `../../test/MnIn2Se4/vca_m20/`. L2 and L3 must never be
added together (see `../../analysis/disorder/README.md` §2).

## The L3 parameters now in the figure

Both were re-fitted on 24³ with the corrected force constants (2026-11-17):

| | Mn₂In₂Se₅ | MnIn₂Se₄ |
|---|---:|---:|
| x | 0.13 | **0.45** |
| L | 15.9 nm | 30.7 nm |
| s | 0.250 | 0.35 |
| rms of ln(model/exp) | 3.2 % | 1.2 % |

The MnIn₂Se₄ fit is new: the previous one rested on the superseded 4×4×1 FC2 and
could not be drawn. `kappa_extrinsic.py --dim-fc2/--pmat` made it possible.
Its fitted x (0.45) happens to lie near the i = 0.44 used for the VCA, but this
is **not** independent support: the acceptable region of the same fit spans
x = 0.00–0.66, i.e. the whole range, because x and L are degenerate. The κ data
do not determine x (see `../../analysis/disorder/README.md` §3.3).

κ_zz convention: the L1 lines here use the 14³ rest term (`kzz_rest`). If the
paper quotes 30³ totals (see `../kappa_vs_T/jobs30/`), switch these lines to the
30³ totals as well, so that all figures use one convention.

## Style

Drawn following the project's own plotting programs, versioned in
`analysis/plotting/`: `figsize=(10, 4.2)`, dotted grid, `-o` for in-plane and
`-s` with `mfc="white"` for out-of-plane, `frameon=False` legend, log y.

```bash
cd test/figures/calc_vs_experiment && python3 plot_calc_vs_experiment.py
```

## What the figure says (300 K)

| | calc intrinsic (L1), single-crystal component | calc extrinsic (L3), texture-mixed | experiment |
|---|---:|---:|---:|
| Mn₂In₂Se₅ Perp (L1: in-plane) | 2.99 | 0.77 | 0.74 |
| Mn₂In₂Se₅ Para (L1: cross-plane) | 1.01 | 0.335 | 0.32 |
| MnIn₂Se₄ Perp (L1: in-plane) | 3.65 | 0.91 | 0.93 |
| MnIn₂Se₄ Para (L1: cross-plane) | 0.64 | **0.49** | 0.49 |

The two calculated columns are not the same kind of quantity. The L1 column is
the pure single-crystal component (in-plane, or κ_zz). The L3 column is what the
figure draws: the extrinsic κ_xx and κ_zz combined with the fitted texture s,
Perp = (1 − s/2)κ_xx + (s/2)κ_zz and Para = s·κ_xx + (1 − s)κ_zz, which is what
the measurement corresponds to. An earlier version of this table gave 0.35 for
MnIn₂Se₄ Para, which matches neither the texture-mixed value (0.488) nor the
unmixed extrinsic κ_zz (0.173); the correct, texture-mixed value is 0.488
(`extrinsic_mnin2se4.json`, 300 K, κ_xx = 1.0715, κ_zz = 0.1732, s = 0.35).

How well L3 matches: the fit's rms of ln(model/experiment) over the 16
high-temperature points (300–623 K, both directions) is 1.2 % for MnIn₂Se₄ and
3.2 % for Mn₂In₂Se₅. At 300 K the deviations are +3.8 % / +4.7 % (Mn₂In₂Se₅ Perp /
Para) and −2.0 % / +0.5 % (MnIn₂Se₄ Perp / Para). This is a 3-parameter fit to the
same data (x, L, s; x and L degenerate), so the agreement shows the model *can*
describe the data, not that x or L are determined. The intrinsic curves are
3–4× too high in-plane. At 100 K the gap is larger still (Mn₂In₂Se₅ in-plane:
10.2 calculated vs 0.79 measured), i.e. the 20–150 K range is entirely
extrinsic-limited. The low-temperature points were not used in the fit.
