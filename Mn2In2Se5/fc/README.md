# Mn₂In₂Se₅ force constants (cutoff 6.51 Å)

The production calculation used the **6.51 Å** case from the S5 cutoff scan (`source_fc = cutoff_scan/cut3_6p51` in `kappa_main_T100-800K.json`).

## Files

| file | size | source (cluster) | purpose |
|---|---|---|---|
| `fc2.hdf5` | 227 KB | 3090 `Mn2In2Se5_441/kl-dft-cpu/step5_fc/cutoff_scan/cut3_6p51/` | phono3py recomputation, antisite/grain-boundary fitting |
| `fc3.hdf5` | 10.9 MB | same as above | same as above |
| `phono3py_disp.yaml` | 616 KB | 3090 `step5_fc/` | primitive cell, supercell matrix, displacements |
| `SPOSCAR` | 8.8 KB | 3090 `step4_disp/` | check the atom order (32 Mn + 32 In + 80 Se = 144 atoms) |
| `FORCE_CONSTANTS_2ND` | 4.4 MB | jzzn `step6_kappa/shengbte_runs/cut3_6p51_m16x16x4/` | comparison with phono3py |
| `FORCE_CONSTANTS_3RD_hiphive_incomplete` | 4.6 MB | same as above (original name `FORCE_CONSTANTS_3RD`) | **incomplete, do not use**: the hiphive export lost 3126/9164 triplets crossing the a3′ boundary, and the production κ was computed with it |
| `FORCE_CONSTANTS_3RD_complete` | 6.8 MB | written in full from `fc3.hdf5` under the minimum-image convention on 2026-11-08 (9164 blocks) | goes with `CONTROL` (`../../analysis/shengbte_crosscheck/run_complete_fc3/`); differs from phono3py by < 1% |
| `BTE.KappaTensorVsT_RTA` | 1.1 KB | same as above | ShengBTE RTA κ(T), compared directly with phono3py RTA |
| `BTE.KappaTensorVsT_CONV` | 1.1 KB | same as above | ShengBTE converged solution, for reference |
| `md5.txt` | — | generated locally | md5 of fc2/fc3 and the three FORCE_CONSTANTS files |

## md5

```
b611f13939c205627a1ba7ed45559150  fc2.hdf5
ac7bf6b2477626a0c95ac3c79aa38d74  fc3.hdf5
9e374e9e3e63425e327ef59de2016f2b  FORCE_CONSTANTS_2ND
6477a498d5f2526e443bfc2f00adce89  FORCE_CONSTANTS_3RD_hiphive_incomplete
001e7ceee7c5d79f3f2f07d7f20e3bce  FORCE_CONSTANTS_3RD_complete
```

## ShengBTE RTA results (`BTE.KappaTensorVsT_RTA`, units W m⁻¹ K⁻¹) — **computed with the incomplete 3RD, superseded**

ShengBTE outputs κ_αβ = 1/(N V) Σ f v v τ at each T; the off-diagonal elements are about 1e-10 (numerically zero), and in-plane xx = yy.

| T (K) | κ_xx = κ_yy | κ_zz |
|---|---|---|
| 100 | 5.1506 | 0.6123 |
| 200 | 2.5721 | 0.3117 |
| 300 | 1.7224 | 0.2098 |
| 400 | 1.2958 | 0.1581 |
| 500 | 1.0389 | 0.1268 |
| 600 | 0.8670 | 0.1059 |
| 700 | 0.7440 | 0.0909 |
| 800 | 0.6516 | 0.0796 |

## Cluster locations

- S5 (force constants): 3090 `/data/wangchaoyue852/taskflow/work/Mn2In2Se5_441/kl-dft-cpu/step5_fc/`
- S6 (ShengBTE): jzzn `/public/home/wangchao/AutoZT/work/Mn2In2Se5_441/kl-dft-cpu/step6_kappa/shengbte_runs/cut3_6p51_m16x16x4/`

The largest single file is 10.9 MB, far below the GitHub 100 MB limit.
