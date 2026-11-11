# test/ —— κ 交叉核对数据包（MnIn₂Se₄ vs Mn₂In₂Se₅，phono3py vs ShengBTE）

目的：排查 **同一套 fc2/fc3 下，phono3py 与 ShengBTE 给出的 κ 相差 6–10 倍** 的根因。

## 结论性问题（待排查）

MnIn₂Se₄，同一 6.45 Å 截断、同一套力常数、**同一 11×11×11 = 1331 q 点**：

| 300 K (W m⁻¹ K⁻¹) | ShengBTE | phono3py | 比值 |
|---|---:|---:|---:|
| RTA | **0.35647** | **2.24981** | 6.3× |
| 完整解（CONV / LBTE） | **0.40633** | **3.9054** | 9.6× |

已排除：q 网格（同网格仍差 6.3×）、力常数来源（同为 `cutoff_scan/cut3_6p45`）、
min_freq（两者 −3.1488e−05 THz）、NAC（都关）、isotope（都开）。
试过 phono3py 的 σ 展宽 0.05–2.0，结果完全不变（2.24981），说明该路径下 σ 未生效，暂无法排除展宽。

## 目录结构

```
MnIn2Se4/
  fit/                            S1_fit 拟合输入与输出（fit-fc-thermal）
    fit_config.json               pheasy 拟合配置（RFE-OLS + ALASSO, RASR=BHH）
    fit_manifest.json             拟合清单（含 code.optimizer.md5、geometry.ws_extended_pairs）
    phonon_summary.json           拟合后声子稳定性判定（stable, min_freq）
    POSCAR / SPOSCAR              原胞 / 4×4×1 超胞（112 原子）
    dataset_disps.npy             训练集位移（45 帧）
    dataset_forces.npy            训练集力
    fc2.hdf5 / fc3.hdf5           phono3py 格式力常数（下游 κ 的输入）
  phono3py_disp.yaml              phono3py 输入（原胞 + 位移信息）
  shengbte_export/                ShengBTE 格式导出（由上面 fc2/fc3 转出）
    FORCE_CONSTANTS_2ND/3RD
    POSCAR
  kappa_shengbte_20x20x5/         正式口径：ShengBTE 迭代解 @ 20×20×5（2000 q 点）
  kappa_shengbte_11x11x11/        同网格对照：ShengBTE 迭代解 @ 11×11×11（1331 q 点）
  kappa_phono3py_lbte_11x11x11/   phono3py 直接解 @ 11×11×11（1331 q 点）

Mn2In2Se5/
  fit/                            S5_fc 拟合输入与输出（kl-dft-cpu）
    fit_config.json               pheasy 配置（ALASSO, RASR=BHH, 4×4×1, 144 原子）
    fit_manifest.json
    phonon_summary.json           stable, min_freq −2.17e−07 THz
    POSCAR / SPOSCAR
    BORN                          介电/Born 有效电荷（该步 NAC=off，仅供核对）
    dataset_disps.npy / dataset_forces.npy   54 帧
    fc2.hdf5 / fc3.hdf5
  phono3py_disp.yaml
  shengbte_export/                ShengBTE 格式（step6 用）
  kappa_shengbte_16x16x4/         正式结果：ShengBTE 迭代解 @ 16×16×4
```

每个 kappa_* 目录内含：
- `kappa_summary.json`（脚本生成，含 solver / bte_method / mesh / 300 K κ / 稳定性判定）
- `BTE.KappaTensorVsT_RTA`、`BTE.KappaTensorVsT_CONV`（ShengBTE 原始输出，16 个温度）
- `BTE.omega`（声子频率）、`BTE.v`（群速度）、`BTE.qpoints`（q 点列表）
- `CONTROL`（ShengBTE 输入）、`FORCE_CONSTANTS_2ND/3RD`（该次运行实际用的力常数）
- `shengbte.log`
- phono3py 目录另有 `kappa-m111111.hdf5`（phono3py κ 原始输出）、`coleigs-m111111.hdf5`

## 已确认结果

**MnIn₂Se₄ · ShengBTE conv @ 20×20×5（正式口径）**
- 300 K 面内 0.3880、κzz 0.0616；FC3 截断 6.45 Å（not_converged）
- q 网格收敛未验证（MESH_CONV=off，用户指定固定网格）

**MnIn₂Se₄ · ShengBTE conv @ 11×11×11（同网格对照）**
- 300 K 面内 0.40633、κzz 0.058541

**MnIn₂Se₄ · phono3py lbte @ 11×11×11**
- 300 K 面内 3.9054、κzz 0.3629

**Mn₂In₂Se₅ · ShengBTE conv @ 16×16×4**
- 300 K 面内 3.0637、κzz 0.8289；FC3 截断 6.51 Å（not_converged）
- 面内 12×12×3 → 16×16×4 变化 +7.65%（>5% 判据），**未收敛**

## 环境

- 力常数拟合：pheasy-gpu 0.0.2（main = b288268，含 RASR-RANK / WS-IMAGES 两个修复）
- phono3py 4.4.0（hanhai25, conda `atomate2_p_a`）
- ShengBTE 2.1.0（`/home/zyzhuqgroup/wangchao217/software/AutoZT/ShengBTE-2.1.0/ShengBTE`）
- 数据生成：2026-11-11
