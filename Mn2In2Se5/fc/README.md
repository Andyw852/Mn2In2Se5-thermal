# Mn₂In₂Se₅ 力常数（截断 6.51 Å）

生产计算用的是 S5 截断扫描里 **6.51 Å** 那一档（`kappa_main_T100-800K.json` 的 `source_fc = cutoff_scan/cut3_6p51`）。

## 文件

| 文件 | 大小 | 来源（集群） | 用途 |
|---|---|---|---|
| `fc2.hdf5` | 227 KB | 3090 `Mn2In2Se5_441/kl-dft-cpu/step5_fc/cutoff_scan/cut3_6p51/` | phono3py 复算、反位/晶界拟合 |
| `fc3.hdf5` | 10.9 MB | 同上 | 同上 |
| `phono3py_disp.yaml` | 616 KB | 3090 `step5_fc/` | 单胞、超胞矩阵、位移 |
| `SPOSCAR` | 8.8 KB | 3090 `step4_disp/` | 核对原子顺序（32 Mn + 32 In + 80 Se = 144 原子） |
| `FORCE_CONSTANTS_2ND` | 4.4 MB | jzzn `step6_kappa/shengbte_runs/cut3_6p51_m16x16x4/` | 与 phono3py 对照 |
| `FORCE_CONSTANTS_3RD_hiphive_incomplete` | 4.6 MB | 同上（原名 `FORCE_CONSTANTS_3RD`） | **残缺，勿用**：hiphive 导出丢了 3126/9164 个跨 a3′ 边界的三元组，生产 κ 就是用它算的 |
| `FORCE_CONSTANTS_3RD_complete` | 6.8 MB | 2026-11-08 由 `fc3.hdf5` 按最短像全量写出（9164 块） | 与 `CONTROL`（`../../analysis/shengbte_crosscheck/run_complete_fc3/`）配套；与 phono3py 差 < 1% |
| `BTE.KappaTensorVsT_RTA` | 1.1 KB | 同上 | ShengBTE 的 RTA κ(T)，与 phono3py RTA 直接比 |
| `BTE.KappaTensorVsT_CONV` | 1.1 KB | 同上 | ShengBTE 的收敛解，供参考 |
| `md5.txt` | — | 本地生成 | fc2/fc3 与三个 FORCE_CONSTANTS 文件的 md5 |

## md5

```
b611f13939c205627a1ba7ed45559150  fc2.hdf5
ac7bf6b2477626a0c95ac3c79aa38d74  fc3.hdf5
9e374e9e3e63425e327ef59de2016f2b  FORCE_CONSTANTS_2ND
6477a498d5f2526e443bfc2f00adce89  FORCE_CONSTANTS_3RD_hiphive_incomplete
001e7ceee7c5d79f3f2f07d7f20e3bce  FORCE_CONSTANTS_3RD_complete
```

## ShengBTE RTA 结果（`BTE.KappaTensorVsT_RTA`，单位 W m⁻¹ K⁻¹）——**用残缺 3RD 算的，作废**

ShengBTE 按 κ_αβ = 1/(N V) Σ f v v τ 逐 T 输出；非对角元约 1e-10（数值零），面内 xx = yy。

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

## 集群位置

- S5（力常数）：3090 `/data/wangchaoyue852/taskflow/work/Mn2In2Se5_441/kl-dft-cpu/step5_fc/`
- S6（ShengBTE）：jzzn `/public/home/wangchao/AutoZT/work/Mn2In2Se5_441/kl-dft-cpu/step6_kappa/shengbte_runs/cut3_6p51_m16x16x4/`

单文件最大 10.9 MB，远小于 GitHub 100 MB 上限。
