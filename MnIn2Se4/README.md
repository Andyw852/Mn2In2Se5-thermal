# MnIn2Se4 晶格热导率计算数据（力常数拟合 + phono3py RTA）

体系：**MnIn₂Se₄**（菱方 R-3m #166，原胞 7 原子 = Mn₁In₂Se₄；3×3×3 超胞 189 原子，45 个随机位移构型 + 1 个平衡参考帧）。

流程：确定壳层与截断 → pheasy 力常数拟合（LASSO / ALASSO + HARM_DENSE，OLS 对照，symfc 对照）→ phono3py RTA 计算 κ。

> 本目录与 `Mn2In2Se5` 部分并列，同属 MnIn₂Se 体系的晶格热导率工作。

---

## 目录结构

| 目录 | 内容 |
|---|---|
| `00_report/` | **主报告** `MnIn2Se4_报告_20260928.md`、RMSE 与 κ 并表、壳层与截断表、动态稳定性检查、mesh 收敛 |
| `01_structures/` | `POSCAR_primitive`（7 原子原胞）、`SPOSCAR_3x3x3`（189 原子超胞） |
| `02_dataset/` | `dataset_disps.npy` / `dataset_forces.npy`（46 帧）、`disp_matrix.pkl` / `force_matrix.pkl`（45 帧）、`dataset_alignment.json`、`ndata_total.txt` |
| `03_fc_fits/` | **39 个拟合目录**：OLS / LASSO / ALASSO × 12 个截断（2.85–7.00 Å）+ symfc × 3。各含 `POSCAR`、`fc2.hdf5`、`fc3.hdf5`、`run.log` |
| `04_kappa_summary/` | 各截断 κ 汇总、NU 重生成日志、mesh 收敛 |
| `05_kappa_hdf5/` | **39 组 κ 原始输出**：每目录一个 `kappa-m111111.hdf5`（11³ q-mesh、300 K、RTA、含同位素）+ `POSCAR` + `kappa_run.log` |
| `06_scripts/` | 拟合驱动、κ 求解器、symfc 脚本、汇总脚本、官方 `pheasy_fit.sh`(v6c) 与 `prepare_dataset.py` |
| `07_holdout_rmse/` | 整构型 5 折 held-out 力预测 RMSE 日志（`holdout_c3_*.log`） |

---

## 关键参数

- **拟合**：`pheasy-gpu-gh` HEAD `b91da0f`（3090，conda env `wc`，RTX 3090）
  `FIT_ORDER=3  HARM_DENSE=true  C2_CUTOFF=None  C3_CUTOFF=<截断>`
  （HARM_DENSE = FC2 全稠密且不加 L1 罚，只对 FC3 稀疏）
- **截断**：2.85 / 3.70 / 3.95 / 4.30 / 4.70 / 5.00 / 5.20 / 5.50 / 5.80 / 6.10 / 6.50 / 7.00 Å
  （全部落在超胞原子间距壳层的间隙中，见 `00_report/壳层与截断_20260928.md`）
- **独立对照**：`symfc` 1.7.3（2.85 / 4.30 / 7.00 Å）
- **κ**：`phono3py` 4.4.0（hanhai venv mace_cpu），11³ q-mesh，300 K，RTA，
  `is_isotope=True, is_N_U=True, write_kappa=True`，supercell = 3x3x3，primitive = POSCAR

---

## 结果速览

### κ 随三阶截断变化（300 K，11×11×11，RTA + 同位素，单位 W/m/K）

| 截断 (Å) | 2.85 | 3.70 | 3.95 | 4.30 | 4.70 | 5.00 | 5.20 | 5.50 | 5.80 | 6.10 | 6.50 | 7.00 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| κ_xx | 5.2219 | 5.1257 | 3.8368 | 2.7065 | 3.8401 | 3.9725 | 3.9126 | 3.8370 | 3.7306 | 3.3048 | 3.3703 | 3.1752 |
| κ_zz | 4.0167 | 3.5122 | 1.0338 | 0.5038 | 1.4173 | 1.9378 | 2.1302 | 1.8752 | 1.2729 | 1.1477 | 0.8922 | 0.6293 |

（ALASSO + HARM_DENSE，报告 §4.1；κ 对截断敏感，尤其 κ_zz）

### 网格收敛（ALASSO，300 K）

| q-mesh | 11×11×11 | 14×14×14 | 18×18×18 |
|---|---|---|---|
| 长度 (Å) | 40.0 | 50.0 | 62.5 |
| κ_xx = κ_yy | 3.8828 | 3.9221 | 3.8912 |
| κ_zz | 1.9471 | 1.7513 | **1.5324** |
| 相对变化 | — | 10.1% | 12.5% |

**结论**：面内 κ_xx/κ_yy 已收敛（0.8–1.0%）；**κ_zz 尚未收敛**（每级降 10–12.5%，扫描被 `MESH_CONV_MAX_LENGTH=62.5 Å` 上限截断）。层状化合物的层间热输运需要更大 q-mesh。

---

## 数据溯源（重要）

用户提供的四个输入文件与 `pheasy-gpu-gh/test/MnIn2Se4/` 官方测试材料一致：

| 文件 | md5 | 说明 |
|---|---|---|
| `POSCAR` | `cef6af221344890c582da4cfb638f438` | 与官方逐字节一致 |
| `dataset_disps.npy` | `f784694ae79ce1b170cb35987ce76c95` | 与官方逐字节一致（存的是**分数坐标**，不是位移） |
| `dataset_forces.npy` | `22906b11c42f2f32c0964b6b267c7424` | 与官方逐字节一致 |
| `SPOSCAR` | `c702c54dabe657730899e838656be9ba` | 与官方同构，但原子顺序不同 |

用官方 `tools/prepare_dataset.py SPOSCAR dataset_disps.npy dataset_forces.npy --frac --ref -1` 生成的
`disp_matrix.pkl` 与官方同名文件数值最大差 **4e-17**，`force_matrix.pkl` md5 **完全相同** → 预处理与官方一致。

结构：菱方 |a|=|b|=|c|=13.4029 Å、α=17.463°；等价六方描述 **a = 4.069 Å，c = 39.586 Å**（层状，层间范德华）。

---

## 复现

```bash
# 1) 数据集（需要 ase/scipy）
python3 06_scripts/prepare_dataset.py SPOSCAR dataset_disps.npy dataset_forces.npy --frac --align-reference

# 2) 拟合（每个截断一个目录）
bash 06_scripts/pheasy_fit.sh FIT_METHOD=LASSO  HARM_DENSE=true C2_CUTOFF=None C3_CUTOFF=5.20 FIT_ORDER=3 NCPU=4
bash 06_scripts/pheasy_fit.sh FIT_METHOD=ALASSO HARM_DENSE=true C2_CUTOFF=None C3_CUTOFF=5.20 FIT_ORDER=3 NCPU=4

# 3) kappa
python3 06_scripts/kappa_full_solver.py 11
```

---

生成日期：2026-09-28（归档）/ 2026-11-07（上传整理）。
