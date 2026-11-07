# Mn₂In₂Se₅ 力常数（待上传）

生产计算用的是 S5 截断扫描里 **6.51 Å** 那一档（`kappa_main_T100-800K.json` 的 `source_fc = cutoff_scan/cut3_6p51`）。
把下面这些文件从集群拷到本目录：

| 文件 | 集群上的位置（材料 `Mn2In2Se5_441`，技能 `kl-dft-cpu`） | 用途 |
|---|---|---|
| `fc2.hdf5` | `step5_fc/cutoff_scan/cut3_6p51/fc2.hdf5` | phono3py 复算、反位/晶界拟合 |
| `fc3.hdf5` | `step5_fc/cutoff_scan/cut3_6p51/fc3.hdf5` | 同上 |
| `FORCE_CONSTANTS_2ND`、`FORCE_CONSTANTS_3RD` | S6 跑 ShengBTE 的目录 `shengbte_runs/cut3_6p51_m16x16x4/` | 与 phono3py 对照（确认两边吃的是同一套力常数） |
| `phono3py_disp.yaml` | `step5_fc/`（没有就用 `step4_disp/`） | 单胞、超胞矩阵、位移 |
| `SPOSCAR` | `step4_disp/SPOSCAR` | 核对原子顺序 |
| `BTE.KappaTensorVsT_RTA` | 同一个 ShengBTE 目录 | ShengBTE 的 RTA，和 phono3py RTA 直接比 |

单个文件超过 95 MB 先停下说明大小（GitHub 单文件上限 100 MB），不要硬推。
仓库根目录的 `.gitignore` 已对本目录放行 `*.hdf5`。
