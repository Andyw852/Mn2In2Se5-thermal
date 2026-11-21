# Mn₂In₂Se₅ — VCA (x = 0.13) versus ordered baseline, mesh 20³

> **Superseded for the VCA curve (2026-11-21).** The project now uses the literature
> x = 0.19 (Williams *et al.* 2024, ESI Table S7): see `../vca_m20_x019/`. This
> directory stays as the record of the x = 0.13 run, and its `kappa_baseline.hdf5`
> (ordered masses) is still the baseline of the pair.

Same route as `../../../MnIn2Se4/vca_m20/`, so the two materials are directly
comparable: phono3py RTA, 20×20×20, 300 K, **FC3 = 6.51 Å** (largest available),
FC2 4×4×1, `primitive_matrix="P"`.

## VCA construction

The antisite model is the one in `analysis/extrinsic/kappa_extrinsic.py`, taken at
**x = 0.13** — the value the κ fit gives for this material:

| site | occupancy | VCA mass (amu) | mass variance g |
|---|---|---:|---:|
| Mn | 0.87 Mn + 0.13 In | **62.73** | **0.1031** |
| In | 0.87 In + 0.13 Mn | **107.03** | **0.0354** |
| Se | natural | 78.96 | 4.63×10⁻⁴ |

62.73 = 0.87×54.938 + 0.13×114.818; 107.03 = 0.87×114.818 + 0.13×54.938.
The primitive cell is the 9-atom cell in the order `Mn Mn In In Se Se Se Se Se`,
so the mass-variance list is
`[0.1031, 0.1031, 0.0354, 0.0354, 4.63e-4 ×5]`.

## Runs

HanHai job **424147**, two RTA runs in one job:

1. **baseline** — real masses, `is_isotope=True`
2. **VCA** — averaged masses, `mass_variances`, `is_isotope=False`

`jobs/run.py` and `jobs/submit.sh` are exactly what ran.

## Results (300 K)

| quantity | baseline | **VCA (x = 0.13)** | suppression |
|---|---:|---:|---:|
| in-plane κ | 2.9984 | **1.6335** | **1.84×** |
| κ_zz total (raw) | 1.3730 | 1.2246 | 1.12× |
| **κ_zz rest of the zone** | **0.9837** | **0.7920** | **1.24×** |

`result.json` holds the raw κ; the rest-of-the-zone terms come from
`kappa_zz_decompose.py` applied to the two hdf5 files.

## Files

| file | contents |
|---|---|
| `result.json` | κ for both runs |
| `kappa_baseline.hdf5`, `kappa_vca.hdf5` | per-mode phono3py output (2.9 MB each) |
| `jobs/` | the calculation as it ran |
| `../../figures/scattering/mn2in2se5/` | the two-panel scattering figure for this pair |
| `../../../analysis/disorder/vca_summary.{csv,json}` | both materials side by side |
