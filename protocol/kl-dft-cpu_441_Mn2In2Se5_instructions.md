# Mn2In2Se5 晶格热导率重算指令（AutoZT kl-dft-cpu，4×4×1 换基矢超胞）

> 给执行 AI：先按 AutoZT `AGENTS.md` 铁律 0 读 `AGENTS.md` 与 `skill/kl-dft-cpu/README.md`（含 [0105] 新增两节），
> 再按本文件执行。本文件里列出的 `conf --set`、`register`，用户已通过本文件授权；`rerun`、`clean`、`-f`、`-y`
> 仍须逐条请示。每个检查点不通过就停下来汇报，不要自行绕过。

## 0. 背景（为什么这样算）

- 旧 3×3×1：超胞只保 12 个点群操作中的 4 个，力常数破对称，κ 随 q 网格来回跳（0.93→0.78→0.96→0.74）。
- 旧 3×3×3：平衡帧残余力 0.103 eV/Å，结构不在平衡位置，虚频 −0.44 THz。
- 新方案：用同一个已弛豫结构换一组基矢，(a, b, a3′)。在这组基矢下 **4×4×1 = 144 原子**，保住全部 12 个操作，
  而且仍是对角超胞，pheasy 和 ShengBTE 都能直接用；安全截断 6.6 Å（3×3×3 菱方只有 5.1 Å）。
  已用 phono3py 和 ShengBTE 2.1.0 做过对照：同一套力常数下，4×4×1 与 3×3×3 的 κ 在输出精度内完全相同。
- AutoZT 需要 [0105] 补丁（分支 `claude/jolly-curie-jsvosv`），它修了这条链上会把这种原胞弄坏、或让结果不对称的地方：
  1. S2/S3/S4 以前不带自旋：S1 用 ISPIN=2 + MAGMOM 弛豫，S4 却是非磁取力。这是残余力大、出虚频的一个直接原因。
  2. 对称化以前会把原胞换成 27 原子的惯用胞（就地写回）。
  3. S1 的 KSPACING 以前不能改。
  4. ShengBTE 自动网格以前不保对称。
  5. pheasy 遇到 9 个数的超胞矩阵会静默只取对角元。

### 0.1 磁序与 U（用户已定）

- 磁序固定为 **FM**（`MAGMOM = 2*5 7*0`），不做 AFM 比较，结果里要注明"假定铁磁"。
- Mn 3d 加 **U_eff = 4 eV**（Dudarev，LDAUTYPE=2），In、Se 不加。
- 判据（2026-11-02 修订）：要求的是**绝缘态，即整数占据**，0.3 eV 只是经验门槛。
  U=4 的 S2 结果：带隙 0.254 eV（约 5 倍 SIGMA），磁矩 9.9844。
  磁矩不是整数，用 SIGMA = 0.05 eV 的高斯展宽就能解释：Mn 3d 少数自旋导带很平，展宽尾巴会把少量电子挪过去。
  核对方法：拿 S2 的输入，只把 `SIGMA` 改成 0.01，在 AutoZT 目录之外单独跑一个单点（独立诊断，几分钟）。
  - 磁矩 = 10.00 ± 0.005：判定为绝缘，**U=4 不变**，继续 S3/S4。
    AutoZT 各步保持 SIGMA = 0.05 不改，因为 S1–S4 必须一致。
  - 仍约 9.98：确实有部分占据，停下汇报。
- 不为了凑带隙门槛去调大 U。U 值只按文献或物理理由定，改 U 会让结构和所有声子一起变。
- S4 每帧磁矩的期望值改为 16 × 9.9844 ≈ 159.75（SIGMA = 0.05 下），各帧之间相差应 < 0.2。

## 1. 前置检查（不通过就停）

1. AutoZT 已切到含 [0105] 和 [0106] 的版本（分支 `claude/jolly-curie-jsvosv` 的 `3d13297` 或更新）：
   - `grep -c "total drift" relax_common.py` 的结果 ≥ 1（[0106] 的力解析修复）；
   - `git -C ~/software/AutoZT log --oneline -3` 里有 `[0105]`；
   - `grep -c inherit_spin_u skill/kl-dft-cpu/gen_step4_disp.py` 的结果 ≥ 1；
   - `ls scripts/rebase_for_diag_supercell.py` 存在。
2. 本地自检结构（只读）：`python3 scripts/rebase_for_diag_supercell.py <附件POSCAR> --dim 4 4 1 -o /tmp/x`。
   应输出 `R-3m (#166)，9 原子，12 个点群操作`，以及 `[OK] 输入基矢上的 4 4 1 已保全部点群操作，不用换基矢`。

## 2. 建材料（用新名字，不复用旧目录）

- 材料名：**`Mn2In2Se5_441`**（新名字。旧材料目录里有 3×3×1 的 `step4_disp/disp-*`，还可能有 fit-fc-thermal 的
  `kl_bundle/`；S5 见到 `kl_bundle` 会自动进入"导入旧力常数"模式，那就成了原胞和超胞对不上。）
- 材料根 POSCAR = 附件 `POSCAR`（9 原子，已弛豫结构换基矢；原子顺序 Mn Mn In In Se×5）。
  ```
  autozt -tt kl-dft-cpu -p Mn2In2Se5_441 register --poscar <附件POSCAR路径>
  ```
- **禁止**运行 `reuse_structure.py`：它会用别的技能的 CONTCAR 覆盖材料根 POSCAR，基矢就换回去了。
- `nac` 保持默认开启（极性半导体需要 LO-TO 劈裂）。

## 3. 各步 step.conf（全部用 `conf --set` 写，写完逐步 `conf` 查最终值）

| 步骤 | 设置 | 原因 |
|---|---|---|
| step1_std_opt | `params.KSPACING=0.038` | 原胞 Γ **8×8×2**，保对称（0.03 会给 10×10×2/10×10×3，不保对称） |
| step1_std_opt | `params.CELL_POLICY=none` | 明确不动输入胞（默认 primitive 对已是原胞的输入也不改，写 none 更保险） |
| step1_std_opt | `params.AUTO_MAG=true`、`params.MAGMOM_OVERRIDE="Mn:5 In:0 Se:0"` | Mn²⁺ 高自旋 FM 起点（S2–S4 会自动继承） |
| step1_std_opt | `params.U_OVERRIDE="Mn:4 In:0 Se:0"` | **用户决定：FM + U**。不加 U 时 PBE-D3 得到的是等效金属态（磁矩 9.31、"带隙" 54 meV）。给 Mn 显式写非零 U 会绕过硒化物的阴离子门控；S2–S4 由 [0105] 自动继承 LDAU* |
| step1_std_opt | `params.CELL_FORCE_TOL_EV_A=0.002` | 末态最大力收紧到 2 meV/Å（默认 0.01，不满足 S4 的要求）。**需要 [0106]**：之前的版本读力恒报 0.0000，这道闸从没生效 |
| step1_std_opt | `params.CELL_STAGE_EDIFFG=1E-6` | 变胞段能量判据从 1E-4 收紧到 1E-6 eV。不收紧的话 VASP 在力约 0.009 eV/Å 就停，三遍都过不了 0.002 的力闸 |
| step1_std_opt | `FUNC` 保持 `auto` | 层状多组分会选 pbe-d3；S1 生成后核对 `workflow_method.txt` 的 FUNC；S2–S4 自动继承，不要中途改 |
| step2_static | `params.KSPACING=0.038` | 8×8×2 |
| step3_nac | `params.KSPACING=0.038` | 8×8×2（想更准可用 12×12×3，见第 5 节） |
| step4_disp | `params.SUPERCELL="4 4 1"` | 对角三个数，**不要写 9 个数** |
| step4_disp | `params.KSPACING=0.038` | 超胞 Γ **2×2×2**，与 S1 的 8×8×2 是同一套物理 k 点，参考帧残余力才能真正接近零 |
| step4_disp | 其余保持默认 | `SUPERCELL_SYMMETRY=strict`、`METHOD=alm`、`FORCE_PREC=Accurate`、`FORCE_EDIFF=1E-7`、`FORCE_LREAL=.FALSE.`、`CUT3_MAX=7.0`、`MAX_DISP=500` |
| step5_fc | `params.EQ_FORCE_MAX=0.01` | 平衡帧残余力 > 0.01 eV/Å 直接停（默认 0.2 太松，旧 3×3×3 的 0.103 就是被它放过的） |
| step5_fc | `params.FC_IMPORT_DIR=off` | 强制自己拟合，绝不导入旧力常数 |
| step6_kappa | `params.NAC=off`（2026-11-02 决定） | S3 的 DFPT 失败：ε∞ = 单位矩阵、Mn 的 Born 电荷 −0.9，结果无效。而且 U=4 下磁矩 9.9844 不随 SIGMA 变化，说明有少量真实的部分占据，在这个计算里不是干净的绝缘体，NAC 本身也不适用。只对本材料关闭，不要改项目级的 `nac: false`（那会影响同项目所有 kl 材料） |
| step6_kappa | 其余保持默认 | `SOLVER=shengbte`、`BTE_METHOD=conv`、`MESH=auto`（[0105] 起网格自动上调到保对称）、`ISOTOPE=true` |

命令格式示例：
```
autozt -tt kl-dft-cpu -p Mn2In2Se5_441 -j step1_std_opt conf --set params.KSPACING=0.038
autozt -tt kl-dft-cpu -p Mn2In2Se5_441 -j step4_disp conf --set params.SUPERCELL="4 4 1"
autozt -tt kl-dft-cpu -p Mn2In2Se5_441 -j step4_disp conf
```

**不要**在 S4 的 `[incar]` 里手写 ISPIN/MAGMOM：[0105] 会从 S1 继承，并按超胞展开成 `32*5 112*0`。手写反而可能和 S1 不一致。

## 4. 参考 KPOINTS（AutoZT 用 VASPKIT 按 KSPACING 生成；生成后逐个核对，必须一字不差）

S1 / S2 / S3（9 原子原胞），文件 `KPOINTS_primitive_S1-S3_8x8x2`：
```
Gamma 8x8x2
0
Gamma
   8   8   2
0.0  0.0  0.0
```
S4 每个 `disp-*`（144 原子超胞），文件 `KPOINTS_supercell_S4_2x2x2`：
```
Gamma 2x2x2
0
Gamma
   2   2   2
0.0  0.0  0.0
```
- 已核对：KSPACING=0.038 时，VASPKIT 不论"四舍五入"还是"向上取整"都给出上面两个网格。
- 这组基矢下，原胞 Γ 网格 n n m 保对称的条件是：**前两个数相等，且 m 是 n/gcd(n,4) 的倍数**。
  例：8 8 2、12 12 3、16 16 4 保对称；10 10 2、10 10 3、13 13 3 不保对称。
  超胞（4×4×1 倒格基）n n m 保对称的条件是 m 是 n 的倍数：1 1 1、2 2 2、3 3 3 可以。
- 生成的 KPOINTS 和上面不一样，就停下汇报，不要提交。

## 5. 逐步检查点（每步完成后检查；只读诊断按铁律 4 用 ssh tail/grep）

**S1_opt**
- CONTCAR 仍为 9 原子；|a| = |b| ≈ 4.03 Å，c 方向第三基矢 ≈ 18.7 Å（倾斜）。空间群 R-3m：用 `symmetry_audit.py --mode audit`，1e-5 与 1e-4 都要是 R-3m。
- 收敛判据以 AutoZT 的 `relax_injob` 为准，三段式（relax_common `STAGE_SPEC`）：
  - 段 a（ISIF=2）：力 < 0.05 eV/Å；
  - 段 b/c（ISIF=3）：`EDIFFG = 1E-4` 是**正数，即能量判据**（相邻离子步能量差 < 1e-4 eV），不是 0.1 meV/Å 的力判据；
  - 作业自己另外判：末态最大力 < `CELL_FORCE_TOL_EV_A`（默认 0.01 eV/Å），晶格变化和应力也要稳定，否则再跑一遍。
  - 模板里的 `EDIFFG = -0.001` 会被三段式覆盖，**不生效**。
- S1 结束后看最后一帧的最大力：
  `awk '/TOTAL-FORCE/{f=1;n=0;x=0;next} f&&/^ -+/{if(n>0)f=0;next} f{m=sqrt($4^2+$5^2+$6^2);if(m>x)x=m;n++} END{print x}' OUTCAR`
  - **≤ 0.003 eV/Å**：通过；
  - 0.003–0.01 eV/Å：在 step1_std_opt 设 `CELL_FORCE_TOL_EV_A=0.002`，然后 `retry`（从 CONTCAR 续跑，不是 rerun），再 `start`。
    原因：S4 要求平衡帧 < 0.005 eV/Å，0.01 的余量不够。
  起始残余力 ~1 eV/Å 是正常的：附件结构是在旧的计算设置下弛豫的，这次 S1 用 PBE-D3 + FM，晶格会明显变化。
- INCAR 里有 `ISPIN = 2`、`MAGMOM = 2*5 2*0 5*0`（这是 S1 自己的自动注入，[0105] 管的是 S2–S4）。
- 只看**最终**总磁矩，弛豫中间值不算数。绝缘的 FM Mn²⁺ 应接近整数 10.0 μB/胞：
  - ≈ 10.0：通过；
  - 塌缩到 ≈ 0：停；
  - 明显不是整数（如 9.3）：说明这个泛函下可能是金属或半金属。先报告，等 S2 的带隙结果，再决定是否继续 S3。
- 记下 S1/S2 的最终总磁矩 M₀，后面 S4 的检查要用。

**S2_static**
- gen 日志有 `[OK] S2_static：自旋/DFT+U 继承 S1：ISPIN=2  MAGMOM=2*5 7*0 …`。
- 总磁矩与 S1 相同。若判为金属：先停下汇报，不要直接关 nac（可能是磁态塌缩）。
- gen 日志若出现"检出数值微畸变"：[0105] 起对称化保持原胞。核对 `step2_static/POSCAR` 仍是 9 原子、基矢不变。

**S3_nac**
- 日志有同样的继承行；OUTCAR 含 `MACROSCOPIC STATIC DIELECTRIC TENSOR`，且 ε_xx = ε_yy。
- 可选更密网格（介电更准）：S3 设 `KSPACING=0.024`，生成后看 KPOINTS。若是 `12 12 3` 就保留；若是 `13 13 3`，改成 `0.026` 再生成，应为 `12 12 3`。

**S4_disp（gen 阶段，提交前）**
- 日志有 `[OK] 超胞 4 4 1 保有 R-3m 的全部 12 个点群操作`。若被 `SUPERCELL_SYMMETRY` 拦下，说明原胞基矢被换了，停下汇报。
- `SPOSCAR` 为 144 原子，晶格与附件 `SPOSCAR_441_reference` 相同（逐分量差 < 1e-3 Å；结构经 S1 重新弛豫后会有微小变化）。
- `disp-00000/INCAR` 有 `ISPIN = 2`、`MAGMOM = 32*5 112*0`（日志有 `step4 00000：自旋/DFT+U 继承 S1` 一行）；所有 `disp-*/KPOINTS` 都是 `2 2 2`。
- 帧数 ≤ `MAX_DISP`；超了就把 `CUT3_MAX` 降到 6.0，然后 retry。

**S4_disp（算完）——最关键的一关**
- 平衡帧 `disp-00000` 的最大残余力：
  ```
  awk '/TOTAL-FORCE/{f=1;n=0;next} f&&/^ -+/{if(n>0)f=0;next} f{m=sqrt($4^2+$5^2+$6^2);if(m>x)x=m;n++} END{print x}' disp-00000/OUTCAR
  ```
  **< 0.005 eV/Å 合格**（旧 3×3×1 是 0.004）。0.005–0.01 汇报后继续；> 0.01 停：先查磁矩、KPOINTS、FUNC/IVDW 是否与 S1 一致。
- 每帧总磁矩（取每个 OUTCAR 最后一次）：
  `for f in disp-*/OUTCAR; do grep "number of electron" $f | tail -1 | awk '{printf "%.1f\n",$NF}'; done | sort | uniq -c`，
  应只有一种值，且 ≈ 16 × M₀（M₀ = S2 的最终总磁矩；M₀ = 10 时就是 160）。
  某帧跑到别的磁态，就把该帧列出来汇报。

**S5_fc（3090，pheasy）**
- `fc_dataset.json`：`supercell_matrix = [[4,0,0],[0,4,0],[0,0,1]]`，`natom_super = 144`，`equilibrium_max_force < 0.01`。
- `phonon_summary.json`：`stable = true`、`fc2_symmetry.broken = false`（`n_ops = 12`）、`fit_rmse_relative < 0.02`；
  `min_frequency_THz` 不低于 −0.1。若仍有离开 Γ 点的虚频，停下汇报，附上 q 点和支序号。
- `cutoff_scan.json` 候选截断上限约 6.6 Å（超过安全截断的档位会在 S4 被剔除）。

**S5.1_plot**：只看 `band_nonac*`，0–1 THz 的声学支在 Γ 附近没有负频。
`step3_nac` 目录还在，S5 仍会用那份无效的 vasprun 生成 BORN，所以 `band_nac*` 是错的，忽略。
S6 要确认日志里是 `NAC=False`。

**S6_kappa（ShengBTE conv）**
- 日志中**不能**出现 `rotations are incompatible with the q-point grid`（[0105] 已把自动网格上调到保对称；出现就停下汇报）。
- `kappa_summary.json`：κ_xx = κ_yy，非对角元 ≈ 0；`mesh_converged = true`。若是因为 `MESH_CONV_MAX_LENGTH` 停的，汇报并请示是否调大。
- κ_zz 对层间方向的网格很敏感（对照测试中 16×16×4 与 16×16×16 的 κ_zz 差 60%）。收敛后再做一次确认：
  在收敛网格 n n m 上把第三个数加倍（例如 `MESH="16 16 8"` 加倍到 `"16 16 16"`；n n 2m 也满足保对称条件），
  看 κ_zz 变化是否 < 5%。不满足就继续加密，并汇报。

## 5.5 推进方式（AGENTS.md 铁律 8：不要自己写 bash/ssh 监控循环）

- 巡检只用 `autozt summary --diff`，或 AutoZT 的作用域 monitor：
  `autozt -tt kl-dft-cpu -p Mn2In2Se5_441 monitor -d -i 300`。不要另起 bash 后台循环；已经起了的，按 PID 停掉。
- **S1→S2→S3 和 S4 生成阶段：人工把关，不开自动推进。** 每步 OK 后先核对上面的检查点，再
  `autozt -tt kl-dft-cpu -p Mn2In2Se5_441 start`。
- **S4 先生成、不提交**：`autozt -tt kl-dft-cpu -p Mn2In2Se5_441 -j step4_disp init`，核对 SPOSCAR、
  KPOINTS（2 2 2）、`disp-00000/INCAR`（ISPIN/MAGMOM）、帧数，然后再 `start`。S4 是最贵的一步，生成错了就白算几十帧。
- S4 全部算完、平衡帧残余力和每帧磁矩核对通过之后，可以给这个材料开自动推进：
  `autozt -tt kl-dft-cpu -p Mn2In2Se5_441 auto on`。S5 有 `EQ_FORCE_MAX=0.01` 硬闸兜底。

## 6. "原胞和超胞不匹配"的风险点与对应措施

| 可能出错的地方 | 后果 | 本方案的措施 |
|---|---|---|
| 材料根 POSCAR 被 `reuse_structure.py` 覆盖 | 基矢回到旧的，4 4 1 不再保对称 | 禁止运行；S4 有对称闸兜底 |
| S1 标准化原胞 | 基矢被换 | `CELL_POLICY=none` |
| S2/S3/S4 微畸变对称化 | 以前会换成 27 原子惯用胞 | [0105] 已修（保持原胞）；检查点核对原子数 |
| S4 写 9 个数的 SUPERCELL | 以前 pheasy 会静默只取对角元 | 只写 `4 4 1`；[0105] 起会直接报错 |
| S5 导入旧 kl_bundle | 力常数对应别的胞 | 新材料名 + `FC_IMPORT_DIR=off` |
| S1 磁性、S4 非磁 | 参考结构不在平衡位置，出虚频 | [0105] 自动继承；检查点核对 MAGMOM 和磁矩 |
| k 网格不对称 / S1 与 S4 采样不一致 | 力带对称破缺偏差，残余力变大 | KSPACING=0.038（8×8×2 ↔ 2×2×2，同一套 k 点） |
| ShengBTE 网格不对称 | κ_xx ≠ κ_yy | [0105] 自动上调；检查日志 |

## 7. 汇报内容

每完成一步就按 AGENTS.md 第六节的格式汇报。最终汇报应包含：
- S4 平衡帧残余力、帧数、每帧磁矩；
- S5 的最低频率、拟合误差、选定的三阶截断；
- S6 的 κ(300 K) xx/yy/zz、收敛网格、κ_zz 加密确认的结果；
- 与旧 3×3×1 结果（0.74–0.96 W/mK 来回跳）的对比。
