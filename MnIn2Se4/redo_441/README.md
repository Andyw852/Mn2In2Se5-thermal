# MnIn₂Se₄ 按 Mn₂In₂Se₅ 同一协议重算（AutoZT kl-dft-cpu，4×4×1 换基矢超胞）

## 为什么要重算

`MnIn2Se4/` 现有的力常数来自 pheasy 官方测试数据集：DFT 设置不明（是否加 U、是否自旋极化、是否 D3 都不知道），
3×3×3 菱方超胞的安全截断只有约 5.2 Å（跨层三阶力常数被切掉，κ_zz 随截断 0.5–4.0 乱变）。
它和 Mn₂In₂Se₅（PBE-D3 + U(Mn)=4 + FM，144 原子，截断 6.51 Å）不是同一套计算，两个材料没法直接比，
也没法把差别归到材料本身。

## 给执行 AI 的指令

**整体照抄 [`../../protocol/kl-dft-cpu_441_Mn2In2Se5_指令.md`](../../protocol/kl-dft-cpu_441_Mn2In2Se5_指令.md)**，
只改下表这些地方。AutoZT 用当前 master（已含 [0105]–[0107]：自旋/U 继承、保原胞对称化、力解析修复、OOM 自动加内存）。
铁律照旧：`rerun`、`clean`、`-f`、`-y` 逐条请示；不自己写监控循环，用 `autozt ... monitor -d`。

| 项 | Mn₂In₂Se₅（原指令） | **MnIn₂Se₄（本次）** |
|---|---|---|
| 材料名 | `Mn2In2Se5_441` | **`MnIn2Se4_441`** |
| 材料根 POSCAR | 9 原子 | 本目录 `POSCAR`：7 原子，顺序 Mn In In Se×4，已换基矢 |
| 自检输出 | `R-3m (#166)，9 原子，12 个点群操作` | `R-3m (#166)，7 原子，12 个点群操作`，且 `4 4 1 已保全部点群操作，不用换基矢` |
| S1 后的晶格 | \|a\|=\|b\|≈4.03，第三基矢≈18.7 Å | \|a\|=\|b\|≈4.07，第三基矢≈16.2 Å（PBE-D3+U 弛豫后会变一点） |
| 磁序 / U | FM，Mn U=4 | **相同**：`MAGMOM_OVERRIDE="Mn:5 In:0 Se:0"`，`U_OVERRIDE="Mn:4 In:0 Se:0"`（原胞只有 1 个 Mn，FM 是唯一的共线态） |
| S2 总磁矩期望 | ≈ 9.98 μB/胞 | **≈ 5.0 μB/胞**（记为 M₀；明显偏离 5 先停下汇报） |
| KPOINTS | 原胞 8×8×2，超胞 2×2×2 | **相同**（本目录两个 KPOINTS 文件；KSPACING=0.038 下 \|b\|/2π = 0.3012/0.0758 → 8 8 2） |
| S4 超胞 | 4 4 1 = 144 原子 | 4 4 1 = **112 原子**，晶格对照本目录 `SPOSCAR_441_reference` |
| S4 每帧磁矩 | 16 × M₀ ≈ 159.75 | **16 × M₀ ≈ 80** |
| 安全截断 | 6.61 Å | **6.50 Å**（`CUT3_MAX` 保持 7.0，超出的档位 S4 自动剔除） |
| NAC | S3 失败 → S6 `NAC=off` | S3 照跑；**S6 一律 `NAC=off`**，和 Mn₂In₂Se₅ 保持一致，便于两材料直接比 |
| S1 力闸 | `CELL_FORCE_TOL_EV_A=0.002`、`CELL_STAGE_EDIFFG=1E-6` | **相同** |
| S5 | `EQ_FORCE_MAX=0.01`、`FC_IMPORT_DIR=off` | **相同** |

**S6 网格（2026-11-08 定，两个材料相同）**：`MESH = auto` + `MESH_CONV_MAX_LENGTH = 50`
（网格自动加密，封顶在 16×16×4）。**不要设显式 `MESH`**：ShengBTE 迭代解在 20×20×5 上已用掉 750 GB 内存，
16×16×16 只会更大。κ_xx 在 16×16×4 已收敛（20×20×5 只差 0.9%）；κ_zz 的层间网格加密由仓库这边用 phono3py 做，
结果见 `analysis/shengbte_crosscheck/README.md`。

## 汇报与上传

完成后按 `Mn2In2Se5/` 的结构上传到本仓库 `MnIn2Se4/redo_441/`：`data/kappa`、`data/phonon`、`figures`、
`fc/`（选定截断那一档的 fc2.hdf5、fc3.hdf5、phono3py_disp.yaml）。
另外写一个 `dft_summary.txt`：S1 末态最大力、S2 带隙与总磁矩、S4 平衡帧残余力、帧数、每帧磁矩分布。

---

## DFT 阶段归档（2026-11-09，AutoZT kl-dft-cpu）

S1–S3 已完成、S4 输入已生成（**计算未提交**），归档在 [`dft/`](dft/)：

- `dft_summary.txt` —— S1 末态最大力 / S2 带隙与总磁矩 / S3 介电张量与 Born 电荷 / S4 输入核对表
- `s1_std_opt/`、`s2_static/`、`s3_nac/` —— 输入（INCAR/POSCAR/KPOINTS/submit.sh）、结构（CONTCAR）、能量（OSZICAR）
- `s4_disp_inputs/` —— SPOSCAR（112 原子 4×4×1）、位移帧 INCAR/KPOINTS/submit.sh 样例

按仓库 `.gitignore` 约定，原始 `OUTCAR` / `EIGENVAL` 不入库（可从集群重取）。

**关键结果**：S1 max|F| = 0.00064 eV/Å；S2 磁矩 **4.9844 μB/胞**、带隙 **0.4055 eV**（间接）；
S3 介电张量退化为单位阵（LEPSILON 响应为 0，**不影响 NAC=off 的 S6**），Born 有效电荷正常。
