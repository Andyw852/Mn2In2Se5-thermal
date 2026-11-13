# MnIn2Se4 fc2 (4×4×2 对称有限位移) — 2026-11-13

把二阶力常数从 4×4×1 换到 **4×4×2** 的超胞重算，目的是修 C44 层间剪切模量的问题
（4×4×1 + BHH 下沿 z 的 TA 声速按旧口径只有 317 m/s，DFT 弹性常数给 1338 m/s）。

## 计算设置

| 项 | 值 |
|---|---|
| 方法 | phonopy 对称有限位移（**不是**随机位移 / alm） |
| 超胞 | **4×4×2**（224 原子：32 Mn + 64 In + 128 Se） |
| 对称性 | `strict`（384 个对称操作，保留完整 R-3m 点群） |
| 位移幅度 | **FD_DISTANCE = 0.02 Å** |
| 帧数 | **7**（对称不等价最小完备集） |
| 泛函 | PBE-D3，U(Mn)=4 eV，FM（MAGMOM = 32\*5 192\*0） |
| EDIFF | **1E-7** |
| KPOINTS | 2×2×1（KSPACING 0.038 规则） |
| NAC | 关闭（半金属，DFPT 介电会发散；与 kl 链 NAC=off 口径一致） |
| 集群 | jzzn（cpu192, premium, 64 核/帧） |

## 7 帧的位移（phonopy 约化后）

位移方向全部为 `(0, 0.011600, 0.016292)`（模长 = 0.02 Å），只有正负与位点不同：

| 帧 | 原子 | 元素 | 符号 |
|---|---|---|---|
| disp-00001 | 1 | Mn | + |
| disp-00002 | 33 | In | + |
| disp-00003 | 33 | In | − |
| disp-00004 | 97 | Se | + |
| disp-00005 | 97 | Se | − |
| disp-00006 | 161 | Se | + |
| disp-00007 | 161 | Se | − |

## 结果：声子谱稳定

```
PHONON_DONE = true
stable      = true
verdict     = pass
min_freq    = +0.000719 THz   (无虚频，n_neg_qpoints = 0)
n_qpoint    = 499, n_branch = 21
max_freq    = 7.784 THz
```

对比 kl 链 **4×4×1** 的 min_freq = −0.006 THz（有软模/微虚频），4×4×2 的谱更干净。

每帧墙钟约 4 小时（jzzn 64 核）；hanhai25 128 核对照帧实测 7h46m。

## ✅ 层法向声速与 DFT 弹性常数吻合

**运行 `sound_velocity.py`（本目录）复现下表**：

| 物理量 | 4×4×2 fc2 | DFT 弹性张量 | 偏差 |
|---|---:|---:|---:|
| 横波 TA₁ | **1348.3 m/s** | 1338.3 m/s | **+0.75%** |
| 横波 TA₂ | **1348.3 m/s** | 1338.3 m/s | **+0.75%** |
| 纵波 LA | **3118.6 m/s** | 3253.4 m/s | −4.14% |
| C₄₄ = ρv_T² | **9.543 GPa** | 9.402 GPa | +1.50% |

拟合窗口 t ≤ 0.005 / 0.010 / 0.020 / 0.040 结果稳定（TA 变化 < 0.3 m/s）。
密度统一用声子结构的 5249.3 kg/m³。

### ⚠️ 取声速时的两个陷阱（本脚本已处理）

**陷阱 1：phonopy 的 q 用「所选定原胞」的倒空间坐标，不是 POSCAR 原胞。**

本次 `primitive_matrix="auto"` 把原胞换成了**菱方原胞**，新原胞的
b₁/b₂/b₃ 与层法向夹角**都是 84.91°**（几乎沿层面）。

POSCAR 的 a₁、a₂ 在层面内 → 层法向 = `(0,0,1)`_POSCAR。
同一个层间方向在两套坐标里：

```
(0,0,t)_POSCAR  ≡  (t,t,t)_phonopy_primitive
```

**错误做法**：用 `(0,0,t)` 算频率、再用 POSCAR 的 b₃ 长度换算 |q| →
得到 **4179 / 6256 m/s** 的假值。两套 |b₃| 之比 **3.760** 正好是误差倍数。

**陷阱 2：Γ 点有数值残差（f₀ = 0.000719 THz）。**

直接算 ω/q 时，低 q 会被这个残差污染，q→0 反而**发散**（t=0.002 给 1429，
t=0.010 给 1351）。正确做法是拟合

```
f² = f₀² + (v·q / 2π / 100)²
```

取斜率开方。本脚本即按此实现。

> 上述两点由用户 2026-11-13 独立核查发现（`_fc2_442_direction_audit_20261113/`），
> 本脚本已独立复现其全部数值。

## 遗留问题

1. **纵波 LA 差 4.14%**：弹性计算与声子计算的结构略有不同（层法向晶格分量
   分别为 **13.20318 Å** 与 **13.21409 Å**），应先统一结构与弹性计算口径，
   再判断是超胞收敛性问题还是结构性差异。
2. **4×4×1 的 317/317 m/s 需按同一坐标约定复核**：旧数值很可能也混用了两套基底，
   在复核前**不能**认定它"偏软"、也不能据此断言 4×4×1→4×4×2 修复了 C44
   （虽然 4×4×2 本身与 DFT 吻合到 0.75%）。

## 文件

| 文件 | 说明 |
|---|---|
| `fc2.hdf5` | 二阶力常数（phonopy 格式） |
| `FORCE_CONSTANTS` | 同上，phonopy 文本格式（供 phono3py / ShengBTE 用） |
| `band-dft-cpu.yaml` | 声子谱（499 q 点 × 21 支），**绘图数据源** |
| `phonon_summary.json` | 稳定性判据结果（虚频门禁 / ZA 二次性） |
| `phonon_fit.log` | 拟合日志 |
| `phonopy_disp.yaml` | 位移定义（含 3×3 超胞矩阵与 primitive_matrix） |
| `disp_plan.json` | 位移计划（method / fd_distance / n_disp） |
| `POSCAR` | 原胞结构 |
| `kl_params.txt` | 计算参数快照 |
| **`sound_velocity.py`** | **层法向声速提取（含两个陷阱的处理 + 示意图注释）** |
| `phonon_band_plot/phonon_band_full.png` | 全场声子谱 |
| `phonon_band_plot/phonon_band_lowfreq.png` | 低频放大（看声学支） |
| `phonon_band_plot/phonon_band_summary.json` | 绘图侧汇总 |
| `frames/disp-*/ ` | 7 帧输入（POSCAR / INCAR / KPOINTS / CONTCAR） |

> OUTCAR / vasprun.xml / queue.out 按仓库 .gitignore 未入库（原始 DFT 输出太大，可从集群重取）。
> 需要力数据时重跑 S2_disp 或从 jzzn 取 `step2_disp/disp-*/OUTCAR`。

## 复现

```bash
python sound_velocity.py .
```
