# Mn₂In₂Se₅ 晶格热导率 κ_L(T)（第一性原理，4×4×1 换基矢超胞）

R-3m 层状 Mn₂In₂Se₅（9 原子原胞）的本征晶格热导率：DFT 有限位移 → pheasy 拟合二、三阶力常数 → ShengBTE 迭代解声子玻尔兹曼方程。
由 AutoZT `kl-dft-cpu` 技能链完成（S1 优化 → S2 静态 → S3 NAC → S4 位移 → S5 力常数 → S6 κ），执行协议见
[`../protocol/kl-dft-cpu_441_Mn2In2Se5_指令.md`](../protocol/kl-dft-cpu_441_Mn2In2Se5_指令.md)。

> 这是**单晶、无缺陷、只有三声子 + 天然同位素散射**的本征值。和热压多晶实验对比要另加织构、反位无序、晶界，
> 见 [`../analysis/extrinsic/`](../analysis/extrinsic/)。

> ## ⚠️ 下表作废（2026-11-08）
>
> 集群 S6 喂给 ShengBTE 的 `FORCE_CONSTANTS_3RD` 是 hiphive 导出的，丢了 34% 的三阶项（跨 a₃′ 边界的层间三元组），
> κ_xx 低约 1.7 倍、κ_zz 低约 4–5 倍。同一套 fc2/fc3 用 phono3py、或用补全的 3RD 跑 ShengBTE，结果一致：
> **300 K RTA κ_xx ≈ 2.96、κ_zz ≈ 0.84（16×16×4）/ 1.11（16×16×8，尚未收敛）**。
> 细节见 [`../analysis/shengbte_crosscheck/`](../analysis/shengbte_crosscheck/)；AutoZT [0127] 已修，S6 需重跑。
> 下表只留作记录。

## 原生产结果（作废：ShengBTE 迭代解，q 网格 16×16×4，三阶截断 6.51 Å）

| T (K) | κ_xx = κ_yy | κ_zz |
|---:|---:|---:|
| 30 | 26.862 | 2.1153 |
| 60 | 11.610 | 1.0243 |
| 90 | 7.2043 | 0.6770 |
| 100 | 6.4058 | 0.6099 |
| 120 | 5.2563 | 0.5101 |
| 150 | 4.1569 | 0.4107 |
| 200 | 3.0958 | 0.3106 |
| 300 | **2.0590** | **0.2090** |
| 400 | 1.5452 | 0.1575 |
| 500 | 1.2373 | 0.1264 |
| 600 | 1.0319 | 0.1055 |
| 700 | 0.8852 | 0.0906 |
| 800 | 0.7750 | 0.0793 |

单位 W m⁻¹ K⁻¹。300 K 面内/面外约 9.9 倍。30–150 K 与 100–800 K 两段用同一套力常数与 q 网格。

**核对结果（2026-11-08）**：

1. **求解器交叉核对** → 发现上面那个导出问题，生产 κ 作废（见文首）。
2. **κ_zz 的 q 网格**：自动收敛只把第三个数从 3 加到 4。phono3py RTA 上 16×16×4 → 16×16×8，κ_zz 涨 33%
   （κ_xx 变 0.6%）；12、16 继续涨到 1.304、1.541，本征 κ_zz 不收敛（长自由程层间声子），
   1 μm 晶粒下收敛到 ≈ 0.97。S6 用 `MESH=auto` + `MESH_CONV_MAX_LENGTH=50`（封顶 16×16×4）重跑；
   报告写法见交叉核对页末尾。与实验对比的外在散射模型对网格不敏感。

## 方法

**DFT（VASP 6.4.3）**

| 项 | 设置 |
|---|---|
| 泛函 | PBE + D3 色散（AutoZT 对层状多组分自动选 pbe-d3） |
| 磁性 | **假定铁磁（FM）**，Mn²⁺ 高自旋起点 MAGMOM = 5 μB/Mn；未比较 AFM |
| DFT+U | Mn 3d U_eff = 4 eV（Dudarev），In、Se 不加 |
| 电子态 | S2：带隙 0.254 eV，总磁矩 9.9844 μB/胞（SIGMA 由 0.05 改 0.01 磁矩不变 → 有少量真实部分占据） |
| k 点 | KSPACING = 0.038（2π/Å）：原胞 Γ 8×8×2，超胞 Γ 2×2×2，两者是同一套物理 k 点 |
| 结构优化 | 三段式（ISIF=2 → ISIF=3 ×2），末态最大力闸 0.002 eV/Å，变胞段 EDIFFG = 1E-6 eV；**实际末态最大力 0.000484 eV/Å**；ENCUT 540（1.5×ENMAX，压 Pulay 应力） |
| 取力 | ENCUT 410、PREC=Accurate、EDIFF=1E-7、LREAL=.FALSE.、ISYM=0；**平衡帧残余力 0.000933 eV/Å** （旧 3×3×3 是 0.103） |
| 位移超胞 | 换基矢 (a, b, a₃′) 下 4×4×1 = **144 原子**，保住 R-3m 全部 12 个点群操作；1 平衡帧 + 54 位移帧 |
| NAC | **关**：S3 的 DFPT 失败（ε∞ 得到单位矩阵、Mn Born 电荷 −0.90），结果无效；`data/phonon/BORN_invalid_not_used` 只作记录 |

**力常数（pheasy，GPU 3090）**：ALASSO，二阶全稠密不加罚（HARM_DENSE），三阶截断扫描后取 6.51 Å
（超胞安全截断 6.61 Å）；拟合相对误差 0.415%，2440 个独立参数；声子谱无虚频，fc2 对 12 个操作的对称性破缺 1.5×10⁻¹⁴ THz。

**κ（ShengBTE）**：迭代解（同时输出 RTA），同位素散射开（天然丰度），NAC 关；
截断按 1-SE、平台、数据量三判据选 6.51 Å（`data/kappa/cutoff_scan_S6.json`、`figures/kappa_vs_cutoff.png`）。

## 文件

```
POSCAR                          S1 的输入结构（换基矢后的 9 原子原胞，弛豫前）
POSCAR_relaxed                  S1 弛豫后的原胞（取自 fc/phono3py_disp.yaml，4×4×1 超胞由它生成）
data/dft_inputs/                S1/S2/S4 的 INCAR、KPOINTS 与 dft_summary.txt
config/step6_kappa.step.conf    S6 的参数（求解器、网格、截断判据）
data/kappa/                     κ(T)：CSV、汇总 JSON、ShengBTE 原始 BTE.KappaTensorVsT_CONV、截断扫描
data/phonon/                    phonon_summary.json（稳定性、拟合误差）、fit_manifest.json、BORN_invalid_not_used
data/provenance/                S1–S6 输入文件指纹（含集群内部路径，仓库保持私有）
figures/                        κ(T)、各向异性、截断收敛、声子谱（band_nac 用的是无效 BORN，忽略）
fc/                             力常数（fc2/fc3.hdf5、ShengBTE 文本；3RD 有残缺版与补全版，见 fc/README.md）
```

## 旧版本

旧 3×3×1（81 原子，丢对称，κ 随网格乱跳 0.74–0.96）和 3×3×3（平衡帧残余力 0.103 eV/Å，虚频 −0.44 THz）
已从工作树删除，留在 git 历史里（提交 `5b5c05b`）。
