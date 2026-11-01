# Mn2In2Se5 晶格热导率计算（fit-fc-thermal）数据

三方层状材料 Mn₂In₂Se₅（原胞 9 原子，层间距 16.23 Å）。本仓库存放结构与力常数拟合/热导率的关键数据，供后续确定处理方案。

## 目录

- `primitive/` — 原胞 POSCAR（9 原子）
- `fcfit_3x3x1/` — 3×3×1 超胞（81 原子）拟合：pheasy ALASSO
- `fcfit_3x3x3/` — 3×3×3 超胞（243 原子）拟合：phono3py symfc → pheasy RFE-OLS
- `kappa_3x3x1/` — 3×3×1 的 κ 结果

## 两个版本的两个问题（关键）

### ① 3×3×1 超胞丢对称 → κ 乱跳

`fcfit_3x3x1` 的拟合 stable（无虚频，min_freq≈0），但超胞 3×3×1 丢失原胞点群操作，
力常数破坏晶体对称，phono3py 按全对称约化 q 点 → κ 随 q 网格来回跳、不可信：
`kappa_3x3x1/mesh_convergence.json` 里 κ 0.93→0.78→0.96→0.74 乱跳 20–32%。

### ② 3×3×3 超胞保对称，但力常数有虚频

`fcfit_3x3x3` 超胞保对称（supercell_symmetry=strict），但 symfc 与 RFE-OLS 两种方法
拟合结果都虚频 −0.435 / −0.436 THz（min_freq，off_gamma，~1683 个负频 q 点）。
根因是 **reference frame 未弛豫**：`fc_dataset.json` 里 `equilibrium_max_force = 0.1031 eV/Å`
（平衡位置残余力偏大）。换拟合方法（symfc→RFE-OLS）虚频几乎不变，证明是数据集问题、非方法问题。

## 结论

要得到可信 κ，需先重新弛豫 Mn2In2Se5 结构（压低平衡位置残余力）→ 重新生成位移-力数据集 → 再拟合 + κ。
后续方案待定。
