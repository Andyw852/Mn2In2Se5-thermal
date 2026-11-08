# Mn₂In₂Se₅ S6 重跑结果（0127 修复后，ShengBTE）

> 本目录是 S6 晶格热导率在 [0127] ShengBTE 三阶力常数导出修复后的最终结果，
> 取代旧 kappa（旧结果因 hiphive 导出 FORCE_CONSTANTS_3RD 丢了 34% 跨层三阶项而作废）。
> 对应作业 hanhai25 359384（2026-11-08 12:23→16:25，退出码 0）。

## 计算配置

| 项 | 值 |
|---|---|
| solver | ShengBTE 2.1.0（hanhai25，OpenMPI 5.0.9-ucx，256 核） |
| bte_method | conv（迭代解，主结果）；同时输出 RTA |
| 网格策略 | MESH=auto + MESH_CONV_MAX_LENGTH=50 |
| 最终网格 | 16×16×4（mesh_converged=False，stop_reason=MESH_CONV_MAX_LENGTH=50 reached） |
| 截断 | 6.51 Å（16 shells，pick=largest；κ 在数据可确定范围内未到平台，取最大截断） |
| 温度 | 50–800 K，步长 50（16 点） |
| 同位素 / NAC | ISOTOPE=true / NAC=off |
| 力常数 | source_fc = cutoff_scan/cut3_6p51（见 ../fc/） |

## 关键结果（300 K）

| 方法 | κ_xx=κ_yy | κ_zz | 单位 |
|---|---|---|---|
| CONV（迭代解） | 3.1412 | 0.84092 | W/m/K |
| RTA | 2.9535 | 0.83455 | W/m/K |

## 网格收敛（conv，300 K，选中截断 6.51 Å）

| mesh | κ_xx=κ_yy | κ_zz | 相对变化(xx/zz) |
|---|---|---|---|
| 12 12 3 (L=40 Å) | 2.8984 | 0.81543 | — |
| 16 16 4 (L=50 Å) | 3.1412 | 0.84092 | +8.38% / +3.13% |

> converged=False：xx/yy 差 8.38% 未达 5% 阈，已到 MESH_CONV_MAX_LENGTH=50 上限；
> 下一档 20×20×5 会 OOM（750 GB 内存，实测 rc=137）。
> κ_zz 的网格加密不在此目录，由 phono3py 补（本征 κ_zz 随网格持续增长未收敛，见 analysis/shengbte_crosscheck/）。

## 截断扫描（RTA，300 K）

| 截断 (Å) | 壳层 | in-plane κ | κ_zz |
|---|---|---|---|
| 5.68 | 12 | 3.0229 | 0.9684 |
| 6.12 | 14 | 2.9422 | 0.8801 |
| 6.51 | 16 | 2.8984 | 0.8154 |

## 0127 重导验证

- cut3_6p51：旧导出缺 3126/9164 个三元组 → 重导后 9164 块 / 9164 三元组，回折自检通过。
- FORCE_CONSTANTS_3RD md5 = af86b14cb8013e6ac49426b5a36aefbc（7887357 字节）。
- 旧残缺版在集群上已改名 stale_0127（未传）。

## 文件清单

| 文件 | 说明 |
|---|---|
| kappa_summary.json | 最终 κ 汇总（16×16×4，50–800 K 全表，conv 为主） |
| mesh_convergence.json | 网格收敛记录（12 12 3 → 16 16 4，未收敛） |
| mesh_convergence_reference.json | 网格收敛参考 |
| cutoff_selection.json | 截断选择（chosen_cut=6.51 Å） |
| BTE.KappaTensorVsT_CONV | 迭代解 κ(T)，16×16×4 |
| BTE.KappaTensorVsT_RTA | RTA κ(T)，16×16×4 |
| FORCE_CONSTANTS_3RD | 6.51 Å 三阶力常数（0127 重导后，9164 块） |
| FORCE_CONSTANTS_2ND | 6.51 Å 二阶力常数 |
| shengbte.log | ShengBTE 运行日志（含 0127 重导记录） |

## 说明

- RTA vs CONV：kappa_summary.json 的 source=BTE.KappaTensorVsT_CONV，主结果是 conv；BTE.KappaTensorVsT_RTA 是同一网格的弛豫时间近似。
- 截断未收敛：κ 随截断单调缓降，6.12–6.51 平台跨度 < 0.5 Å 阈，按 pick=largest 取 6.51 Å。
- 上传整理：2026-11-08。
