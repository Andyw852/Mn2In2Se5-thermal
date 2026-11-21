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

**Since 2026-11-21, x is fixed at the literature value** and only L and s are
fitted (`make_extrinsic_json.py`, same 24³ bases, same grid and metric as
`kappa_extrinsic.py fit`). The κ data cannot determine x, so an independent
structure refinement is used for it:

| | Mn₂In₂Se₅ | MnIn₂Se₄ |
|---|---:|---:|
| x (fixed) | **0.19** (Williams *et al.* 2024, ESI Table S7) | **0.44** (Range *et al.* 1991) |
| L | 22.1 nm (acceptable 18.7–22.1) | 30.7 nm |
| s | 0.233 (0.207–0.273) | 0.35 (0.33–0.37) |
| rms of ln(model/exp) | 3.4 % | 1.2 % |

The free three-parameter fit of 2026-11-17 gave, for comparison, Mn₂In₂Se₅
x = 0.13, L = 15.9 nm, s = 0.250, rms 3.2 % and MnIn₂Se₄ x = 0.45, L = 30.7 nm,
s = 0.35, rms 1.2 % (`../../analysis/extrinsic/results/*_m24/fit_*.json`). The
literature x costs only 0.2 % in rms for Mn₂In₂Se₅. Its effect on what the figure
shows is small: s changes from 0.25 to 0.23 and L from 16 to 22 nm.

The MnIn₂Se₄ fit rests on the 4×4×2 FC2 (`kappa_extrinsic.py --dim-fc2/--pmat`).
Its free-fit x (0.45) happens to lie near the i = 0.44 used for the VCA, but this
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
| Mn₂In₂Se₅ Perp (L1: in-plane) | 2.99 | 0.78 | 0.74 |
| Mn₂In₂Se₅ Para (L1: cross-plane) | 1.01 | 0.336 | 0.32 |
| MnIn₂Se₄ Perp (L1: in-plane) | 3.65 | 0.92 | 0.93 |
| MnIn₂Se₄ Para (L1: cross-plane) | 0.64 | **0.49** | 0.49 |

The two calculated columns are not the same kind of quantity. The L1 column is
the pure single-crystal component (in-plane, or κ_zz). The L3 column is what the
figure draws: the extrinsic κ_xx and κ_zz combined with the fitted texture s,
Perp = (1 − s/2)κ_xx + (s/2)κ_zz and Para = s·κ_xx + (1 − s)κ_zz, which is what
the measurement corresponds to. An earlier version of this table gave 0.35 for
MnIn₂Se₄ Para, which matches neither the texture-mixed value (0.488) nor the
unmixed extrinsic κ_zz (0.173); the correct, texture-mixed value is 0.488
(free fit, 300 K, κ_xx = 1.0715, κ_zz = 0.1732, s = 0.35; with x fixed at 0.44 the
file now gives κ_xx = 1.0768, κ_zz = 0.1738 and Para 0.490).

How well L3 matches: the fit's rms of ln(model/experiment) over the 16
high-temperature points (300–623 K, both directions) is 1.2 % for MnIn₂Se₄ and
3.4 % for Mn₂In₂Se₅ (x fixed). At 303 K the deviations are +4.5 % / +4.7 % (Mn₂In₂Se₅
Perp / Para) and −1.9 % / +0.6 % (MnIn₂Se₄ Perp / Para). With x fixed this is a
2-parameter fit (L, s) to the same data, so the agreement shows the model *can*
describe the data with the literature antisite fraction; L is then an effective
mean free path, not a measured grain size. The intrinsic curves are
3–4× too high in-plane. At 100 K the gap is larger still (Mn₂In₂Se₅ in-plane:
10.2 calculated vs 0.79 measured), i.e. the 20–150 K range is entirely
extrinsic-limited. The low-temperature points were not used in the fit.
