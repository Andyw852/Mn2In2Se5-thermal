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

## 结论与处置

4x4x1 沿 c 只有 1 个重复（16.22 A），fc2 力常数截断到 8.1 A，长程相互作用被切掉，
导致近 Gamma TA 支假软。对照同材料 3x3x3 / 189 原子数据集（`MnIn2Se4/02_dataset/`）
RFE-OLS 拟合 min = -0.011 THz（stable），确认是超胞问题而非材料失稳。

**处置**：弃用这批 4x4x1 拟合结果，κ 改用 3x3x3 / 189 原子数据集重拟（见仓库 `MnIn2Se4/` 既有数据 +
后续 redo_441 的 S1_fit 3x3x3 归档）。本目录只保留 fit_config.json / methods.json 作记录：
集群上的 fc2.hdf5 / fc3.hdf5 / cutoff_scan/（25 MB）/ phonon_summary.json 已被随后同一目录的
3x3x3 重拟覆盖，未归档。
