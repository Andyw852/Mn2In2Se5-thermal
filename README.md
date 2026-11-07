# MnIn2Se 体系晶格热导率数据（MnIn2Se4 / Mn2In2Se5）

本仓库集中存放 **MnIn₂Se₄** 与 **Mn₂In₂Se₅** 两个三方层状材料的力常数拟合与晶格热导率（κ）计算数据，
按「报告 / 结构 / 数据集 / 力常数拟合 / κ 汇总 / κ 原始输出 / 脚本 / 留出集检验」分门别类。

---

## 目录

| 目录 | 体系 | 内容 |
|---|---|---|
| `MnIn2Se4/` | **MnIn₂Se₄**（R-3m #166，7 原子原胞，3×3×3 超胞 189 原子） | 完整数据包，见其 `README.md` |
| `primitive/` | Mn₂In₂Se₅ | 原胞 POSCAR（9 原子） |
| `fcfit_3x3x1/` | Mn₂In₂Se₅ | 3×3×1 超胞（81 原子）拟合：pheasy ALASSO |
| `fcfit_3x3x3/` | Mn₂In₂Se₅ | 3×3×3 超胞（243 原子）拟合：symfc → pheasy RFE-OLS |
| `kappa_3x3x1/` | Mn₂In₂Se₅ | 3×3×1 的 κ 结果 |

---

## MnIn₂Se₄ 摘要

菱方 R-3m #166，原胞 7 原子（Mn₁In₂Se₄）；等价六方 a = 4.069 Å，**c = 39.586 Å**（层状，层间范德华）。
45 个随机位移构型 + 1 个平衡参考帧。

- **拟合**：pheasy LASSO / ALASSO + HARM_DENSE（FC2 全稠密无罚，只对 FC3 稀疏），OLS 与 symfc 对照；
  三阶截断扫描 2.85 / 3.70 / 3.95 / 4.30 / 4.70 / 5.00 / 5.20 / 5.50 / 5.80 / 6.10 / 6.50 / 7.00 Å。
- **κ**：phono3py 4.4.0，11³ q-mesh，300 K，RTA + 同位素。

κ 随截断强烈变化，例如 ALASSO：κ_xx 从 5.22（2.85 Å）降到 3.18（7.00 Å）；
κ_zz 从 4.02 降到 0.63。

网格收敛（300 K）：κ_xx=κ_yy 在 11³→18³ 只变 0.8–1.0%（3.883 / 3.922 / 3.891）；
但 **κ_zz 未收敛**（1.947 → 1.751 → 1.532，每级降 10–12.5%）。

---

## Mn₂In₂Se₅ 两个版本的两个问题（关键）

### ① 3×3×1 超胞丢对称 → κ 乱跳

`fcfit_3x3x1` 的拟合 stable（无虚频，min_freq ≈ 0），但超胞 3×3×1 丢失原胞点群操作，
力常数破坏晶体对称，phono3py 按全对称约化 q 点 → κ 随 q 网格来回跳、不可信：
`kappa_3x3x1/mesh_convergence.json` 里 κ 0.93→0.78→0.96→0.74 乱跳 20–32%。

### ② 3×3×3 超胞保对称，但力常数有虚频

`fcfit_3x3x3` 超胞保对称（supercell_symmetry=strict），但 symfc 与 RFE-OLS 两种方法
拟合结果都虚频 −0.435 / −0.436 THz（min_freq，off_gamma，约 1683 个负频 q 点）。
根因是 **reference frame 未弛豫**：`fc_dataset.json` 里 `equilibrium_max_force = 0.1031 eV/Å`
（平衡位置残余力偏大）。换拟合方法（symfc→RFE-OLS）虚频几乎不变，证明是数据集问题、非方法问题。

**结论**：要得到可信 κ，需先重新弛豫 Mn₂In₂Se₅ 结构（压低平衡位置残余力）→ 重新生成位移-力数据集 → 再拟合 + κ。

---

## 数据溯源（MnIn₂Se₄）

用户提供的四个输入文件与 `pheasy-gpu-gh/test/MnIn2Se4/` 官方测试材料一致（md5 逐字节比对）；
用官方 `prepare_dataset.py` 生成的 `disp_matrix.pkl` / `force_matrix.pkl` 与官方一致
（数值最大差 4e-17 / md5 完全相同）。详见 `MnIn2Se4/README.md`。

---

整理上传：2026-11-07。
