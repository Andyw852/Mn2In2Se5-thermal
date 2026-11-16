# Raw phono3py output of the kappa_zz scan

One directory per case. Each holds the `kappa-*.hdf5` that phono3py wrote with
`write_kappa=True` (per-mode kappa, gamma, group velocities, heat capacity),
the `result.json` of that run, and for the HanHai cases the slurm log.

| case | mesh | fc3 | machine | who it is used by |
|---|---|---|---|---|
| `m11` | 11³ | 6.45 Å | workstation, phono3py 3.24.0, 155 s | §6.3, §6.5 |
| `m11_b1um` | 11³ + `boundary_mfp=1.0` | 6.45 Å | workstation, 155 s | validation of the a-posteriori boundary reconstruction in §6.4 |
| `m14` | 14³ | 6.45 Å | workstation, 629 s | §6.3, §6.5, reference row of §6.6 |
| `m18` | 18³ | 6.45 Å | workstation, 2674 s | §6.3, §6.5 |
| `m24` | 24³ | 6.45 Å | **HanHai** job 409609, 548 s | §6.3, §6.5 |
| `m30` | 30³ | 6.45 Å | **HanHai** job 409612, 2017 s | §6.3, §6.5 |
| `c560_m14` | 14³ | **5.60 Å** | **HanHai** job 409606, 24 s | §6.6 |
| `c623_m14` | 14³ | **6.23 Å** | **HanHai** job 409608, 24 s | §6.6 |

`kappa_zz_decompose.py` (one level up) reads these files directly:

```bash
python kappa_zz_decompose.py kappa_zz_scan/m*/kappa-m*.hdf5 \
       --json kappa_zz_gamma_line.json --mfp 1 10
```

Two notes on consistency:

- the 14³, 5.60 Å case was run on both machines and gives
  κ_zz = 2.523638910 (phono3py 3.24.0, 12 threads) vs 2.523638946
  (phono3py 4.4.0, 64 threads) — 1.4·10⁻⁸ relative, so mixed provenance is safe;
  the committed file is the HanHai one;
- the 18³ file has κ_zz = 0.9435 and the 24³ one 0.8052, i.e. the totals keep
  falling, which is the Γ–Z artifact documented in `../KAPPA_MESH_SCAN.md`.
  The quantity that converges is the "rest of the zone" column there.
