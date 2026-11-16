# Calculated κ(T) over the experimental κ_tot figure

`calc_vs_experiment.png / .pdf` — the experimental points from
`experiment/kappa_tot_digitized.csv` with our calculations on top.

- **Perp = in-plane** (κ_xx), **Para = out-of-plane** (κ_zz) — the labelling used
  for hot-pressed samples, where the platelet c axis aligns with the pressing
  direction.
- **Symbols** = experiment (filled circle Perp, open square Para), colour = material.
- **Solid lines** = our intrinsic single-crystal RTA (in-plane: filled marker,
  κ_zz: open marker), from `../kappa_vs_T/data_*.json`.
- **Dashed** = extrinsic model for Mn₂In₂Se₅ (x = 0.13 Mn/In antisite disorder,
  L = 16 nm grain boundary, texture s = 0.25) — the model that converts the
  intrinsic curves into the measured polycrystal values. Fitted to the same CSV
  in `analysis/extrinsic/` (rms of ln(model/exp) = 3.2 %).

## Style

Drawn following the project's own plotting programs, which are now versioned in
`analysis/plotting/` (recovered from the backup drive):

- `figsize=(10, 4.2)`, two panels, `ax.grid(ls=":", alpha=0.4)`;
- `-o` for in-plane, `-s` with `mfc="white"` for out-of-plane;
- `frameon=False` legend, log y where the range needs it.

```bash
cd test/figures/calc_vs_experiment && python3 plot_calc_vs_experiment.py
```

## What the figure says (300 K)

| | calc intrinsic | calc extrinsic | experiment |
|---|---:|---:|---:|
| Mn₂In₂Se₅ in-plane | 2.99 | 0.77 | 0.74 |
| Mn₂In₂Se₅ cross-plane | 1.01 | 0.34 | 0.32 |
| MnIn₂Se₄ in-plane | 3.65 | — | 0.93 |
| MnIn₂Se₄ cross-plane | 0.64 | — | 0.49 |

The intrinsic single-crystal κ is 3–4× the measured polycrystal in-plane and 2×
cross-plane. At 100 K the gap is larger still (10.2 vs 0.79), i.e. the 20–150 K
range is entirely extrinsic-limited.

**The MnIn₂Se₄ extrinsic curve is deliberately absent**: the only fit in
`analysis/extrinsic/` for that material rests on the superseded 4×4×1 FC2 (whose
interlayer force constants are 4.2× too soft). It must be re-fitted on the 4×4×2
FC2 before it belongs on this figure.
