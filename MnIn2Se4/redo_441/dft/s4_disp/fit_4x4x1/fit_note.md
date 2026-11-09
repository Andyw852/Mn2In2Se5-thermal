# 4x4x1 数据集的双方法力常数拟合（2026-11-10，结论：虚频，弃用）

## 配置

- 技能：fit-fc-thermal（AutoZT），hpc = 3090，单作业内串行两种方法
- `FIT_METHODS = pheasy:RFE-OLS,ALASSO`（RFE-OLS 为主方法、ALASSO 为方法对照）
- `PHEASY_BIN = pheasy-gpu`，`PHEASY_C3_CUTOFF = 6`（gen 逐档扫描 4.30 / 4.71 / 5.00 / 5.60 / 6.23 / 6.45 A）
- `SUPERCELL = 4 4 1`（112 原子），`ENABLE_FC = 3`，`CUT3_SCAN = auto`，`CUT3_BOOTSTRAP = 0`
- 数据集：`dft/s4_disp`（46 帧 x 112 原子），平衡帧残力 0.0005 eV/A 已扣除
- 超胞保对称检查：4 4 1 保有 R-3m 全部 12 个点群操作（通过）
- 完整配置见 `fit_config.json`，方法清单见 `methods.json`

## 结果（两种方法都很好，但都虚频）

| 方法 | 力 RMSE | 相对误差 | 虚频闸 |
|---|---|---|---|
| RFE-OLS（主） | 0.88 meV/A | 0.58 % | **FAIL** min = -0.1388 THz（no-NAC -0.151 THz） |
| ALASSO | 0.84 meV/A | 0.55 % | **FAIL** min = -0.1388 THz |

- 判据：imag_policy `near_gamma_large`（abs(q) < 0.05x1.2 的近 Gamma 声学支、|nu| > 0.10 THz）
- 极小点：branch 0（最低 TA），q = (-0.0357, 0, 0.0357)，|q| ~ 0.0505，145 个负频 q 点
- fc2 对称性 spread = 5.3e-9 THz / 12 个操作（对称性本身没问题）
- 两方法谐波部分一致（HARM 不惩罚、拟合同一批力）→ 虚频与拟合方法无关

## 色散诊断（从 fc2.hdf5 直接算）

- Gamma 点三支声学 ~ 0（ASR 满足，不是平移不变性破缺）
- 最低 TA 支（软方向 (-1,0,1)）：t=0.01 → -0.024，0.05 → -0.109，0.10 → -0.149，0.15 → +0.139 THz
- 面内方向同样：0.01 → -0.012，0.10 → -0.096，0.20 → +0.094 THz
→ 是 |q| ~ 0.03-0.12 一整段的软 TA 支（等效弹性常数偏负），不是单点数值噪声

## 结论与处置（2026-11-10 更正）

**不是超胞大小的问题**：4x4x1 三个方向的垂直宽度都 >= 13.2 A（安全截断 6.50 A），比旧 3x3x3
菱方超胞的约 10.5 A（安全截断 5.17 A）更宽；Mn2In2Se5 用同一套 4x4x1 是稳定的。真正的两个原因：

1. **没加转动不变性 / Huang 平衡条件**（`PHEASY_RASR=auto` 在 bulk 解析成 none）。
   Gamma 点三支声学为 0 只说明平移求和规则满足；|q| ~ 0.03-0.12 的 TA 支为负正是这类约束
   缺失的典型表现。
2. **平衡帧（disp-00000）约 -1.23 kBar 的各向同性残余应力**（对角 -1.246 / -1.246 / -1.207 kBar，
   偏量 < 0.03 kBar）：S1 用 ENCUT=540 弛豫、S4 用 410 取力的 Pulay 应力差，会让长波声学支
   变软（单靠它解释不了 -0.14 THz，但会推一把）。

**处置**：这批 4x4x1 拟合结果作废，但数据保留（本材料唯一与 Mn2In2Se5 同 DFT 设置的力数据集）。
用同一套数据重做拟合，只改 `PHEASY_RASR = H`（Huang 约束），方法和其余参数不变。
若 RASR=H 后仍有 < -0.05 THz 的近 Gamma 负频，下一步先用 ENCUT=540 重算平衡帧（单帧）
再决定是否改 S4 设置。

本目录只保留 fit_config.json / methods.json：集群上的 fc2.hdf5 / fc3.hdf5 / cutoff_scan/（25 MB）/
phonon_summary.json 已被随后同一目录的重拟覆盖，未归档。
（注：`datasets/MnIn2Se4_raw` 那套 3x3x3 / 189 原子数据是 pheasy 自带测试数据集、DFT 设置不明，
不能用于本材料结果，只作参考。）
