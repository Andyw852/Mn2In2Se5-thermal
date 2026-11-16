# Calculated κ(T) over the experimental κ_tot figure

`calc_vs_experiment.png / .pdf` — the experimental points from
`experiment/kappa_tot_digitized.csv` with our calculations drawn on top.

- **Perp = in-plane** (κ_xx), **Para = out-of-plane** (κ_zz), matching the
  labelling in the hot-pressing literature: the c axis of the platelets aligns
  with the pressing direction, so Para ≈ κ_zz.
- **Symbols** = experiment. Colour = material, marker = direction
  (filled circle = Perp, open diamond = Para).
- **Thin solid / dash-dot** = our intrinsic single-crystal κ (in-plane / κ_zz),
  from `../kappa_vs_T/data_*.json`.
- **Thick translucent** = the extrinsic model for Mn₂In₂Se₅
  (x = 0.13 Mn/In antisite disorder, L = 16 nm grain boundary, texture
  s = 0.25) — the model that converts the intrinsic curves into the measured
  polycrystal values. It was fitted to the same experimental CSV in
  `analysis/extrinsic/` (rms of ln(model/exp) = 3.2 %).

## Style

The figure is drawn with the AutoZT kl chain's own settings: the script imports
`skill/kl-dft-cpu/kappa_plots.py` and applies its `RC` block (Arial 9 pt, inward
ticks on all four sides, thin black frame, no grid, 300 dpi) and its direction
palette.

```bash
cd test/figures/calc_vs_experiment && python3 plot_calc_vs_experiment.py
```

## What the figure says

| 300 K | calc intrinsic | calc extrinsic | experiment |
|---|---:|---:|---:|
| Mn₂In₂Se₅ in-plane | 2.99 | 0.77 | 0.74 |
| Mn₂In₂Se₅ cross-plane | 1.01 | 0.34 | 0.32 |
| MnIn₂Se₄ in-plane | 3.65 | — | 0.93 |
| MnIn₂Se₄ cross-plane | 0.64 | — | 0.49 |

The intrinsic single-crystal κ is **3–4× the measured polycrystal** in-plane and
2× cross-plane; once antisite disorder and ~16 nm grain boundaries are included,
Mn₂In₂Se₅ matches to a few percent. The low-temperature panel makes the same
point more strongly: at 100 K the intrinsic in-plane κ of Mn₂In₂Se₅ is 10.2 while
the sample shows 0.79 — the 20–150 K range is entirely extrinsic-limited.

**The MnIn₂Se₄ extrinsic fit is not shown**: the only fit in
`analysis/extrinsic/` for that material rests on the superseded 4×4×1 FC2, whose
interlayer force constants are 4.2× too soft. It needs re-fitting on the 4×4×2
FC2 before it can go on this figure.
