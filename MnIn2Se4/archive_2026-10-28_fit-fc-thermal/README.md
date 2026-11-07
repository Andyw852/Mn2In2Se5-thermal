# MnIn2Se4 晶格热导率计算（fit-fc-thermal）

三方层状范德华材料 MnIn₂Se₄（空间群 R-3m，166；原胞 7 原子，3×3×3 超胞 189 原子，层间距 13.2 Å）的晶格热导率 κ 计算。力常数由 pheasy（ALASSO）从位移-力数据集拟合，phono3py（RTA）算 κ，symfc 做二阶对照。

## 关键结果

| 量 | 值 |
|---|---|
| 面内 κ(300 K) | **4.03 W/m·K**（xx=yy，11→14 网格变化 0.60%，收敛）|
| 面外 κ(300 K) | **1.68 W/m·K**（zz，11→14 网格变化 8.37%，**未收敛**）|
| 拟合相对误差 | 0.646%（ALASSO，OLS 极限 α→0，1303 free IFC）|
| 声子稳定性 | stable，无显著虚频（min_freq −0.019 THz，噪声级）|
| 三阶截断 | 5.0 Å（largest，安全截断 5.17 Å 下可确定的最大）|
| 网格 | 14×14×14（mesh_converged=false，cap 触发）|

κ(T) 全温度（100–800 K）见 `s2_kappa/kappa_summary.json` 与 `kappa_vs_mesh.png` / `kappa_vs_cutoff.png`。

## 主要诊断结论

1. **zz（面外）方向不收敛**：κ_zz 随网格 1.84 → 1.69（8.4%）仍在下行，是层状材料面外慢收敛的固有表现。
2. **根因是截断被安全上限卡死**：`supercell_safe_cutoff = 0.5 × 最短垂直高度 − 0.1 ≈ 5.17 Å`，而层间距 13.2 Å，故**跨层三阶力常数被整段切掉**；换 CUT3_PICK=largest（候选 4.3/4.7/5.0，thirdorder −n 第 2/3/4 近邻）也只是更干净，加不进跨层项。截断扫描到 5.0 Å 仍未到 κ 平台（`cutoff_selection.status=not_converged`）。
3. **ShengBTE 与 phono3py 对本品系差 24×**（ShengBTE 0.159 vs phono3py 3.85，同一套 FORCE_CONSTANTS）：已排除力常数文件格式（FC2 频率一致、FC3 roundtrip 一致）、R 向量缩放、scalebroad、体积、比热、群速度；判定为 ShengBTE 求解器对 7 原子/21 支大体系的三声子散射率虚高（其 RTA 甚至低于自身 small-grain 极限，物理反常）。shengbte-gpu（单卡/多卡）与官方 CPU 版结果一致，排除 GPU fork 的 bug。故 κ 采用 phono3py。

## 目录

- `s1_fit/` — S1 力常数拟合结果（phonon_summary、cutoff_scan、fc_fit_summary、fit_config、fit_manifest、fit_metrics、POSCAR、SPOSCAR）
- `s2_kappa/` — S2 晶格热导率（kappa_summary、mesh_convergence、kappa_vs_mesh、kappa_vs_cutoff）

## 复现要点

- fit-fc-thermal 技能，pheasy ALASSO + phono3py-symfc 对照。
- S2 参数：isotope=false，mesh_conv_max_length=50（cap 14×14×14），cut_scans=[4p30,4p70,5p00]，nominal_cut3_A=5.0。
- 数据：displacement-force 数据集（未含在此仓库，需另行提供）。
