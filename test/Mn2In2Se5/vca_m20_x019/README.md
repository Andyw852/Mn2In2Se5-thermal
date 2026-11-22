# Mn₂In₂Se₅ — VCA at the literature x = 0.19, mesh 20³

Replaces the x = 0.13 VCA curve (`../vca_m20/`, kept as the archive). The antisite
fraction is the powder Rietveld value of the R-3m phase, x = 0.19(2) (Williams *et al.*,
*J. Mater. Chem. C* **12**, 1753 (2024), ESI Table S7; `analysis/disorder/README.md` §5.1).

| site | occupancy | VCA mass (amu) | mass variance g |
|---|---|---:|---:|
| Mn | 0.81 Mn + 0.19 In | 66.315 | 0.12548 |
| In | 0.81 In + 0.19 Mn | 103.441 | 0.05157 |
| Se | natural | 78.96 | 4.63×10⁻⁴ |

## How the numbers were obtained

- **In-plane κ and κ_zz total**: HanHai job 5
  (`../../figures/kappa_vs_T/jobs_lbte_wigner/run_Mn2In2Se5_vca_m20_rta/`,
  `run_lbte_wigner.py`, phono3py 4.4.0, RTA, 20³, FC3 6.51 Å, FC2 4×4×1, 12 T).
  Only its JSON is in the repository.
- **κ_zz rest of the zone** = total − first Γ–Z grid point
  (`kappa_zz_decompose.py` convention). The first point is one irreducible grid point,
  q = (0, 0, 1/20), BZ index 409, weight 2. `first_gz_point.py` recomputes only that
  point with the job-5 driver. It checks that the masses and mass variances equal
  job 5's. In the RTA each grid point is computed independently, so the split is exact.
- **Check of the single-point route.** It was run with the inputs of the full x = 0.13
  run (job 424147, rounded masses 62.73 / 107.03). The stored mode_kappa, gamma and
  frequencies of that point are reproduced to 3×10⁻¹⁵. With the unrounded x = 0.13
  masses the point differs by 0.44 %, because the near-Γ acoustic modes are that
  sensitive to the masses. This is why the x = 0.19 point uses exactly job 5's masses.
- A full 20³ rerun in a 4-core cloud container was stopped after 2 h without output.
  The full run would take several hours there; job 5 took 4.5 min on 64 threads.

## Results

| T (K) | in-plane | κ_zz total | first Γ–Z point | **κ_zz rest** |
|---:|---:|---:|---:|---:|
| 30 | 7.5528 | 14.0632 | 6.6535 | 7.4097 |
| 100 | 3.1203 | 3.5348 | 1.4244 | 2.1104 |
| 200 | 1.9324 | 1.7769 | 0.6828 | 1.0941 |
| **300** | **1.4611** | 1.2048 | 0.4514 | **0.7534** |
| 500 | 1.0157 | 0.7407 | 0.2697 | 0.4710 |
| 700 | 0.7914 | 0.5373 | 0.1924 | 0.3449 |

All 12 temperatures are in `first_gz_point.json` and `decomposition.json`.

At 300 K, against the ordered baseline (2.9984 / 0.9837, job 424496):

| | x = 0.19 (now) | x = 0.13 (archive) |
|---|---:|---:|
| in-plane | 1.4611 (2.05×) | 1.6335 (1.84×) |
| κ_zz rest | 0.7534 (1.31×) | 0.7920 (1.24×) |
| (2κ_xx + κ_zz)/3 against experiment 0.601 | 1.23 (2.0×) | 1.35 (2.3×) |

## Files

| file | contents |
|---|---|
| `first_gz_point.py` | the single-point run and its check; writes `first_gz_point.json` |
| `postprocess.py` | writes `decomposition.json` and `../../figures/kappa_vs_T/data_mn2in2se5_vca.json` (x = 0.19; the x = 0.13 file is kept as `data_mn2in2se5_vca_x013.json`) |

```bash
python first_gz_point.py && python postprocess.py
```

## Scattering-rate figure

The Mn₂In₂Se₅ scattering-rate figure (`../../figures/scattering/mn2in2se5/`) needs the
per-mode hdf5 of job 5. Its `vca_x.json` says which x the figure shows: 0.13 until the
file below is in, 0.19 after.

1. Copy `kappa_vca_rta_m20.hdf5` (about 16 MB) from the job-5 run directory on HanHai
   to this directory.
2. Run `postprocess.py`. With the file present it also:
   - decomposes the full file and checks it against `first_gz_point.json` (2×10⁻⁵);
   - writes the 300 K slice `kappa_vca_300K.hdf5`;
   - replaces `../../figures/scattering/mn2in2se5/kappa_vca.hdf5` and sets `vca_x.json`
     there to 0.19.
3. Rerun `plot_scattering_figure.py` there and `report/make_report_html.py` (`_zh.py`).
   The HTML caption reads x from the data.

The full instructions are in `protocol/local_AI_next_steps_2026-11-22.md`. The route was
tested end to end with the x = 0.13 12-T file (`../../figures/kappa_vs_T/raw/`), which
agreed with its own decomposition to 8×10⁻⁶.
