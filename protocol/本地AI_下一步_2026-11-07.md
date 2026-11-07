# 给本地 AI：2026-11-07 之后的任务

按顺序做。每一步做完汇报结果再进行下一步。

## 1. 合并整理好的仓库

整理在 `Andyw852/Mn2In2Se5-thermal` 的分支 `claude/jolly-curie-jsvosv` 上，先快进合并到 master：

```bash
git clone https://github.com/Andyw852/Mn2In2Se5-thermal.git   # 已有本地副本就 cd 进去
cd Mn2In2Se5-thermal
git fetch origin
git checkout master
git merge --ff-only origin/claude/jolly-curie-jsvosv
git push origin master
```

`--ff-only` 失败说明 master 上有人推过新东西，停下汇报，不要强推。

合并后核对：根目录只有 `Mn2In2Se5/ MnIn2Se4/ experiment/ analysis/ protocol/ README.md .gitignore`；
`Mn2In2Se5/data/kappa/kappa_main_T100-800K.csv` 与原 `Mn2In2Se5_441-thermal` 里的同名文件 md5 相同。

## 2. 上传 Mn₂In₂Se₅ 的力常数和 DFT 输入（只读拷贝，集群上什么都不改）

清单见 `Mn2In2Se5/fc/README.md`。S5 在 3090 上跑、S6 在 hanhai25 上跑，文件可能分在两台机器上，用 `ssh … ls` 找。

```bash
cd Mn2In2Se5-thermal/Mn2In2Se5/fc
# 例：<S5目录> = 3090 上 Mn2In2Se5_441/kl-dft-cpu/step5_fc
scp 3090:<S5目录>/cutoff_scan/cut3_6p51/fc2.hdf5 .
scp 3090:<S5目录>/cutoff_scan/cut3_6p51/fc3.hdf5 .
scp 3090:<S5目录>/phono3py_disp.yaml .          # 没有就从 step4_disp 拿
# ShengBTE 那次运行：<S6目录>/shengbte_runs/cut3_6p51_m16x16x4/
scp hanhai25:<S6目录>/shengbte_runs/cut3_6p51_m16x16x4/FORCE_CONSTANTS_2ND .
scp hanhai25:<S6目录>/shengbte_runs/cut3_6p51_m16x16x4/FORCE_CONSTANTS_3RD .
scp hanhai25:<S6目录>/shengbte_runs/cut3_6p51_m16x16x4/BTE.KappaTensorVsT_RTA .
scp jzzn:<S4目录>/SPOSCAR .
ls -la                                         # 单个文件 > 95 MB 就停下汇报，不要推
md5sum fc2.hdf5 fc3.hdf5 > md5.txt
```

DFT 输入放到 `Mn2In2Se5/data/dft_inputs/`（都是小文本）：
`step1_std_opt/{INCAR,KPOINTS,workflow_method.txt}`、`step2_static/{INCAR,KPOINTS}`、`step4_disp/disp-00000/{INCAR,KPOINTS}`，
再写一个 `dft_summary.txt`：S1 末态最大力、S2 带隙与总磁矩、S4 平衡帧残余力、S4 帧数、每帧磁矩分布（`sort | uniq -c` 的输出）。
POTCAR 不传（有版权）。

```bash
git add Mn2In2Se5/fc Mn2In2Se5/data/dft_inputs
git commit -m "Mn2In2Se5: force constants (cut3 6.51 A) + DFT inputs"
git push origin master
```

## 3. 删除两个多余的仓库（内容已并入本仓库）

| 仓库 | 去向 |
|---|---|
| `Andyw852/Mn2In2Se5_441-thermal` | 全部并入 `Mn2In2Se5/`（`config/` 只保留了 step6 的 step.conf，hpc.yaml 等集群配置没搬） |
| `Andyw852/MnIn2Se4-thermal` | 全部并入 `MnIn2Se4/archive_2026-10-28_fit-fc-thermal/` |

第 1 步核对通过后再删，删除不可恢复：

```bash
gh auth refresh -h github.com -s delete_repo
gh repo delete Andyw852/Mn2In2Se5_441-thermal --yes
gh repo delete Andyw852/MnIn2Se4-thermal --yes
```

别的仓库（taskflow、taskflow-v2.0 等）不在本次范围内，不要动。

## 4. 集群：MnIn₂Se₄ 按同一协议重算

照 `MnIn2Se4/redo_441/README.md` 做（材料名 `MnIn2Se4_441`，与 Mn₂In₂Se₅ 的协议只差表里那几项）。
S1→S2→S3、S4 生成阶段人工把关；S4 提交前核对 SPOSCAR（112 原子）、KPOINTS（2 2 2）、`disp-00000/INCAR` 的 ISPIN/MAGMOM。

## 5. 不用做的

- Mn₂In₂Se₅ 的 S6 不要 rerun：κ_zz 网格加密和 phono3py 交叉核对由仓库这边用上传的力常数做；
  要是交叉核对发现问题，再定集群上怎么补。
