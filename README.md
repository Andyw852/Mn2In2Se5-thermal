# MnIn₂Se₄ / Mn₂In₂Se₅ 晶格热导率：第一性原理计算与实验对比

两个 R-3m 层状硒化物（层间范德华）的晶格热导率。计算给的是单晶本征值，实验是热压多晶，
两者之间差织构、反位无序、晶界和电子热导，`analysis/extrinsic/` 把这几项放进同一个模型里对实验做拟合。

## 现状（2026-11-08）

| | MnIn₂Se₄ | Mn₂In₂Se₅ |
|---|---|---|
| 计算协议 | pheasy 官方测试数据集（DFT 设置不明），3×3×3 菱方 189 原子，截断 ≤ 5.2 Å | AutoZT kl-dft-cpu：PBE-D3 + U(Mn)=4 eV + FM，4×4×1 换基矢 144 原子，截断 6.51 Å |
| κ_L(300 K) xx / zz（W/mK） | 3.9 / 2.1（phono3py RTA，11³，ALASSO 5.2 Å） | ~~2.059 / 0.209~~（ShengBTE 的 3RD 残缺，作废）→ **2.96 / 0.84–1.11**（phono3py RTA，16×16×4 / 16×16×8，κ_zz 未收敛） |
| 实验 κ_tot(303 K) Perp / Para | 0.933 / 0.485 | 0.741 / 0.320 |
| 状态 | **待按 Mn₂In₂Se₅ 同一协议重算** → `MnIn2Se4/redo_441/` | 力常数可靠；ShengBTE 导出 bug 已在 AutoZT [0127] 修复，S6 待按收敛的 κ_zz 网格重跑（见 `analysis/shengbte_crosscheck/`） |

## 目录

| 路径 | 内容 |
|---|---|
| [`Mn2In2Se5/`](Mn2In2Se5/) | 生产计算结果：κ(T)、声子、截断扫描、图、溯源；`fc/` 放力常数 |
| [`MnIn2Se4/`](MnIn2Se4/) | 现有数据包（39 组拟合、κ 原始输出、脚本、报告）；`redo_441/` 是重算用的结构、k 点和指令；`archive_2026-10-28_fit-fc-thermal/` 是更早一次拟合（含 ShengBTE/phono3py 差异的诊断） |
| [`experiment/`](experiment/) | 实验图、读图得到的 CSV、还缺哪些实验数据 |
| [`analysis/extrinsic/`](analysis/extrinsic/) | 本征 κ + 反位质量无序 + 晶界 + 织构的模型、脚本和拟合结果 |
| [`analysis/shengbte_crosscheck/`](analysis/shengbte_crosscheck/) | ShengBTE 与 phono3py 交叉核对、FC3 导出 bug 的证据与 q 网格收敛 |
| [`protocol/`](protocol/) | AutoZT kl-dft-cpu 4×4×1 换基矢协议（Mn₂In₂Se₅ 用过的执行指令，MnIn₂Se₄ 重算照抄） |

## 计算与实验怎么对上

热压多晶里晶粒取向不完全，用织构参数 s = ⟨sin²θ⟩（θ 为晶粒 c 轴与热压方向夹角，0 = 完全织构，2/3 = 随机）：

```
κ_Para = s·κ_xx + (1−s)·κ_zz
κ_Perp = (1−s/2)·κ_xx + (s/2)·κ_zz
(2·κ_Perp + κ_Para)/3 = (2κ_xx + κ_zz)/3      ← 与织构无关，可以直接和计算比
```

Mn₂In₂Se₅ 在 300 K：修正后的计算 (2κ_xx+κ_zz)/3 ≈ 2.3（16×16×8 RTA），实验 0.60，差约 3.9 倍
（作废的 ShengBTE 值是 1.44，差 2.4 倍）。差距随温度升高而缩小，是外在散射（与温度无关的散射率）
叠加在 Umklapp 上的典型表现，单靠织构解释不了。拟合方法和结果见 [`analysis/extrinsic/README.md`](analysis/extrinsic/README.md)。

## 仓库约定

- 私有仓库：`data/provenance` 等文件里有集群内部路径。
- 原始 DFT 输出（OUTCAR、vasprun.xml、WAVECAR…）不入库；力常数 hdf5 入库（单文件 < 100 MB）。
- 本仓库取代了 `Mn2In2Se5_441-thermal`（内容已并入 `Mn2In2Se5/`）和 `MnIn2Se4-thermal`（并入 `MnIn2Se4/archive_2026-10-28_fit-fc-thermal/`）。
