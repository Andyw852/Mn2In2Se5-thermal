# MnIn2Se4 fc2 (4×4×2 对称有限位移) — 2026-11-13

把二阶力常数从 4×4×1 换到 **4×4×2** 的超胞重算，目的是修 C44 层间剪切模量
过软的问题（4×4×1 + BHH 下沿 z 的 TA 声速只有 317 m/s，DFT 弹性常数给 1338 m/s）。

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

## 结果

```
PHONON_DONE = true
stable      = true
verdict     = pass
min_freq    = +0.000719 THz   (无虚频，n_neg_qpoints = 0)
n_qpoint    = 499, n_branch = 21
max_freq    = 7.784 THz
```

对比 kl 链 **4×4×1** 的 min_freq = −0.006 THz（有软模/微虚频），4×4×2 的谱干净。

每帧墙钟约 4 小时（jzzn 64 核）；hanhai25 128 核对照帧实测 7h46m。

## 文件

| 文件 | 说明 |
|---|---|
| `fc2.hdf5` | 二阶力常数（phonopy 格式） |
| `FORCE_CONSTANTS` | 同上，phonopy 文本格式（供 phono3py / ShengBTE 用） |
| `band-dft-cpu.yaml` | 声子谱（499 q 点 × 21 支），**绘图数据源** |
| `phonon_summary.json` | 稳定性判据结果（虚频门禁 / ZA 二次性） |
| `phonon_fit.log` | 拟合日志 |
| `phonopy_disp.yaml` | 位移定义（含 3×3 超胞矩阵） |
| `disp_plan.json` | 位移计划（method / fd_distance / n_disp） |
| `POSCAR` | 原胞结构 |
| `kl_params.txt` | 计算参数快照 |
| `phonon_band_plot/phonon_band_full.png` | 全场声子谱 |
| `phonon_band_plot/phonon_band_lowfreq.png` | 低频放大（看声学支） |
| `phonon_band_plot/phonon_band_summary.json` | 绘图侧汇总 |
| `frames/disp-*/ ` | 7 帧输入（POSCAR / INCAR / KPOINTS / CONTCAR） |

> OUTCAR / vasprun.xml / queue.out 按仓库 .gitignore 未入库（原始 DFT 输出太大，可从集群重取）。
> 需要力数据时重跑 S2_disp 或从 jzzn 取 `step2_disp/disp-*/OUTCAR`。

## ⚠️ 已知问题：层间声速仍不符

从这份 fc2 提出的 Γ 附近声速（q 收敛已检验）：

| 方向 | 本 fc2 (m/s) | DFT 弹性常数 (m/s) | 偏差 |
|---|---|---|---|
| 面内 a1\*/a2\* | 1051 / 1574 / 3355 | 1295 / 1849 / 3540 | 低 19% / 15% / 5% ✅ |
| 层间 b3 | **4179 / 6256 / 13339** | **1338 / 1338 / 3252** | **高 3.1× / 4.7× / 4.1×** ❌ |

**4×4×1 偏软（317）→ 4×4×2 偏硬（4179）**，都没落在 1338 附近。
待查：该 "rebased" 原胞的 a3 = (0, 9.408, 13.214) 相对 Cartesian z 倾斜约 35°，
层法向是否等于 b3 方向尚未确认；也可能是 4×4×2 沿 c 只有 2 个单胞仍不足以刻划
层间长程力常数。
