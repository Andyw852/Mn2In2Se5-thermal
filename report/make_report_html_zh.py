#!/usr/bin/env python3
"""Chinese version of the HTML report (calc_vs_experiment_report_zh.html).

Same data, same numbers and same section order as make_report_html.py; only the
text and the figure labels are Chinese (figures fig*_zh.png / .pdf). Run
make_report_html.py first in the repository so that plot_data/ is current.

  python3 report/make_report_html_zh.py
"""
import base64

import make_report_html as en          # paths, CSS, CHK, table helpers
import plot_figures

L = en.HTML_LABEL
M4, M5 = "MnIn2Se4", "Mn2In2Se5"
OUT = en.OUT.with_name("calc_vs_experiment_report_zh.html")
CSS = en.CSS.replace('font:15px/1.6 "IBM Plex Sans","Segoe UI",Helvetica,Arial,sans-serif',
                     'font:15px/1.7 "PingFang SC","Microsoft YaHei","Noto Sans CJK SC","Segoe UI",sans-serif')


def methods_zh():
    cut = en._rows("convergence_fc3_cutoff.csv")
    msh = en._rows("convergence_qmesh.csv")
    sv = {(r["material"], r["quantity"]): float(r["value"]) for r in en._rows("solver_terms.csv")}
    st, ok = en._step, en._ok

    def cut_rows(m):
        rr = [r for r in cut if r["material"] == m]
        out = []
        for i, r in enumerate(rr):
            a, z = float(r["inplane_population"]), float(r["kzz_rest"])
            if i:
                a0, z0 = float(rr[i - 1]["inplane_population"]), float(rr[i - 1]["kzz_rest"])
                sa, sz = f"{st(a, a0)} {ok(a, a0)}", f"{st(z, z0)} {ok(z, z0)}"
            else:
                sa = sz = "—"
            b, e = ("<b>", "</b>") if i == len(rr) - 1 else ("", "")
            out.append(f"<tr><td>{L[m]}</td><td>{b}{float(r['cutoff_A']):.2f}{e}</td><td>{a:.3f}</td><td>{sa}</td>"
                       f"<td>{z:.3f}</td><td>{sz}</td></tr>")
        return "\n".join(out)

    def mesh_rows(m):
        rr = [r for r in msh if r["material"] == m]
        out = []
        for i, r in enumerate(rr):
            v = [float(r[k]) for k in ("inplane", "kzz_total", "kzz_rest")]
            if i:
                v0 = [float(rr[i - 1][k]) for k in ("inplane", "kzz_total", "kzz_rest")]
                s = [f"{st(x, x0)} {ok(x, x0)}" for x, x0 in zip(v, v0)]
            else:
                s = ["—"] * 3
            out.append(f"<tr><td>{L[m]}</td><td>{r['N']}³</td><td>{v[0]:.4f}</td><td>{s[0]}</td><td>{v[1]:.4f}</td>"
                       f"<td>{s[1]}</td><td>{v[2]:.4f}</td><td>{s[2]}</td></tr>")
        return "\n".join(out)

    r4 = sv[(M4, "inplane_LBTE_14")] / sv[(M4, "inplane_RTA_14")]
    r5 = sv[(M5, "inplane_LBTE_10")] / sv[(M5, "inplane_RTA_10")]
    q4 = sv[(M4, "inplane_LBTE_14")] * sv[(M4, "mesh_correction_inplane")] + sv[(M4, "coherence_inplane")]
    q5 = sv[(M5, "inplane_RTA_32")] * r5 + sv[(M5, "coherence_inplane")]
    z4 = sv[(M4, "population_kzz_rest")] + sv[(M4, "coherence_kzz")]
    z5 = sv[(M5, "population_kzz_rest")] + sv[(M5, "coherence_kzz")]
    return f"""<h2><span class="n">2</span>计算设置与收敛性</h2>
<p class="lead">本节说明力常数怎么得到、κ 对三阶力常数截断和 q 网格收敛到什么程度，以及每个推荐值用的是哪种求解方法。
<b>收敛判据：最后一次加密使 κ 变化小于 10 %，即认为收敛</b>（表中 ✓），这也是本项目的精度目标。
数据：<code>plot_data/convergence_fc3_cutoff.csv</code>、<code>convergence_qmesh.csv</code>、<code>solver_terms.csv</code>。</p>

<h3>2.1 力常数：二阶（FC2）和三阶（FC3）的来源</h3>
<div class="tbl"><table class="lt">
<thead><tr><th></th><th class="r">{L[M4]}</th><th class="b">{L[M5]}</th></tr></thead>
<tbody>
<tr><td>DFT</td><td colspan="2">VASP 6.4.3，PBE + D3（IVDW = 12），Mn 3d 加 U = 4 eV（Dudarev），铁磁，ENCUT 410 eV，
PREC = Accurate，EDIFF = 10⁻⁷，LREAL = .FALSE.，KSPACING 0.038 Å⁻¹；不加 NAC</td></tr>
<tr><td>弛豫结构</td><td>同一个弛豫后的 7 原子 R-3m 原胞，FC2 和 FC3 都用它</td><td>同一个弛豫后的 9 原子 R-3m 原胞</td></tr>
<tr><td><b>FC2</b>（简谐）</td><td><b>4×4×2 超胞</b>，224 原子；7 个对称性约化的 0.02 Å 位移（phonopy）</td>
<td rowspan="2"><b>二阶、三阶用同一个 4×4×1 超胞</b>，144 原子；1 个平衡帧 + 54 个位移帧，用 pheasy（ALASSO）一起拟合</td></tr>
<tr><td><b>FC3</b>（非谐）</td><td><b>4×4×1 超胞</b>，112 原子；45 个随机位移帧 + 1 个平衡帧，pheasy（RFE-OLS）拟合</td></tr>
<tr><td>采用的 FC3 截断</td><td><b>6.45 Å</b></td><td><b>6.51 Å</b></td></tr>
<tr><td>FC2 检验</td><td>层间横波声速 1348 m/s，DFT 弹性常数给出 1338 m/s（+0.75 %）</td><td>层间横波 1387 m/s，C₄₄ 9.5–10.3 GPa</td></tr>
<tr><td>热导率程序</td><td colspan="2">phono3py（3.24.0 和 4.4.0，两者粒子项一致）</td></tr>
</tbody></table></div>
<ul>
<li><b>{L[M4]} 的二阶和三阶力常数来自同一晶体的两种超胞，不是两个不同的晶体。</b>两者用的是完全相同的弛豫原胞
（晶格矢量和原子坐标一致）和完全相同的 DFT 输入：INCAR 只有 MAGMOM 的原子数不同，k 点密度也相同
（4×4×2 上 2×2×1，4×4×1 上 2×2×2，对应同一个 k 间距）。唯一的区别是超胞沿 c 方向的层数。</li>
<li><b>为什么要两种超胞。</b>4×4×1 超胞沿 c 只有一个晶胞（13.2 Å），原子会和自己跨过范德华间隙的周期像相互作用，
使层间二阶力常数严重偏软：层间横波声速只有 323 m/s，而弹性常数给出 1338 m/s。改用 4×4×2 后得到 1348 m/s。
三阶数据集如果也在 4×4×2 上重算，要多花约两天 DFT，而且不会提高可用的截断半径，所以 FC3 仍用 4×4×1。</li>
<li><b>如何组合。</b>phono3py 允许 FC2 和 FC3 在不同超胞上（FC2 用 <code>phonon_supercell_matrix</code> = 4×4×2，
FC3 用 <code>supercell_matrix</code> = 4×4×1）。ShengBTE 要求两者在同一超胞上，做不到，所以 {L[M4]} 的热导率全部用 phono3py 计算。</li>
<li>{L[M5]} 不需要第二种超胞：它的 4×4×1 二阶力常数本身就给出正确的层间声速。</li>
</ul>

<h3>2.2 三阶力常数截断</h3>
<div class="tbl"><table>
<thead><tr><th>材料</th><th>截断 (Å)</th><th>面内 κ</th><th>变化</th><th>κ_zz</th><th>变化</th></tr></thead>
<tbody>
{cut_rows(M4)}
{cut_rows(M5)}
</tbody></table></div>
<p class="cap">300 K，RTA 粒子项，含天然同位素，κ_zz 取网格收敛项；网格 14³（{L[M4]}）和 16³（{L[M5]}）。
<b>选取 6.45 Å 和 6.51 Å</b>，即 4×4×1 三阶超胞在保守规则下允许的最大截断。按 10 % 判据，两种材料最后一步都满足。
{L[M4]} 的步长还在缩小（面内 −22.9 % → −8.7 %，κ_zz −50.5 % → −8.2 %），剩余误差估计约 −10 %。
{L[M5]}（面内 −3.9 % → −4.3 %，κ_zz −11.9 % → −9.9 %）用同样 54 帧在 6.67–7.43 Å 重新拟合，并做了留出验证：
从 6.90 Å 起都是过拟合，κ 的下降不是物理效应。拟合噪声平均会把 κ 压低（帧数越少 κ 越低），所以 6.90 Å 的结果给出下一步的界：
<b>面内不低于 −2.1 %，κ_zz 不低于 −3.3 %</b>（<code>test/Mn2In2Se5/cutoff_cv/</code>）。更大的截断用现有数据测不了。</p>

<h3>2.3 q 网格</h3>
<div class="tbl"><table>
<thead><tr><th>材料</th><th>网格</th><th>面内 κ</th><th>变化</th><th>κ_zz 总值</th><th>变化</th><th>κ_zz 收敛项</th><th>变化</th></tr></thead>
<tbody>
{mesh_rows(M4)}
{mesh_rows(M5)}
</tbody></table></div>
<p class="cap">300 K，RTA，各向同性 N³ 网格（{L[M4]} 含同位素；{L[M5]} 不含同位素，同位素只使面内降低 1.2 %、κ_zz 降低 0.4 %）。
<b>面内 κ 已收敛</b>（从 14³ 起每步变化小于 3 %）。<b>κ_zz 总值每步还在下降 12–15 %</b>，这是数值积分误差：
Γ–Z 线上紧挨 Γ 的那一个 q 点在这些网格上占 κ_zz 的 15–50 %，而它的权重按 N⁻³ 减小。去掉这个点后得到<b>收敛项</b>，
从 14³（{L[M4]}）和 24³（{L[M5]}）起每步变化小于 1 %。
<b>选取</b>：κ_zz 取收敛项的平台值，{L[M4]} 14³–30³ 为 0.648 ± 0.003，{L[M5]} 24³–32³ 为 0.977 ± 0.008（含同位素）；
面内 {L[M4]} 用 14³ 的 LBTE 加上由 RTA 序列测得的 +1.3 % 网格修正，{L[M5]} 用 32³ 的 RTA。
κ(T) 曲线用 14³ 和 16³，VCA 对比两种材料都用 20³（与收敛值相差在 3 % 以内）。</p>

<h3>2.4 求解方法：RTA、迭代解（LBTE）和 Wigner</h3>
<ul>
<li><b>RTA</b>（弛豫时间近似）给每个声子一个寿命，把所有散射都当作阻碍热流的过程。计算便宜，并给出逐模式数据，
κ(T) 曲线、网格和截断扫描以及无序模型都用它。</li>
<li><b>LBTE</b>（phono3py 对线性化玻尔兹曼方程的直接解，与 ShengBTE 迭代到收敛的解等价）求解完整的碰撞矩阵，
动量守恒的 Normal 过程不会被当作热阻。面内比 RTA 高：{L[M4]} 高 <b>{100 * (r4 - 1):.0f} %</b>
（14³：{sv[(M4, 'inplane_LBTE_14')]:.4f} 对 {sv[(M4, 'inplane_RTA_14')]:.4f}），{L[M5]} 高 <b>{100 * (r5 - 1):.1f} %</b>
（10³：{sv[(M5, 'inplane_LBTE_10')]:.4f} 对 {sv[(M5, 'inplane_RTA_10')]:.4f}）。κ_zz 两者相差仅 4·10⁻⁵。
<b>因此面内推荐值以 LBTE 为基础，κ_zz 用 RTA。</b></li>
<li><b>Wigner 方法</b>（Simoncelli、Marzari、Mauri，Nat. Phys. 15, 809 (2019)）把 κ 分成通常的粒子项（population）和
<b>相干项</b>（coherence，频率间距小于线宽的模式之间的波动性输运）。对于这类低热导、支数多的晶体，相干项不可忽略。
相干项在 RTA 框架下计算（Wigner-RTA），然后加到粒子项上。300 K 下的大小：{L[M4]} 面内 {sv[(M4, 'coherence_inplane')]:.3f}、
κ_zz {sv[(M4, 'coherence_kzz')]:.3f}；{L[M5]} 面内 {sv[(M5, 'coherence_inplane')]:.3f}、κ_zz {sv[(M5, 'coherence_kzz')]:.3f}
W m⁻¹ K⁻¹，即面内 +2–5 %，κ_zz +10–12 %。</li>
<li><b>相干项的版本规则。</b>phono3py 4.4.0 的 <code>SMM19</code> 算出的 {L[M4]} 相干项偏大 12–17 倍，原因是 Γ 附近的软模
没有正则化导致发散；3.24.0 和 4.x 的 <code>IBDB19</code> 选项结果一致。所以 {L[M4]} 用 3.24.0 的值。
{L[M5]} 的 4.4.0 结果通过了同样的软模检验（去掉 0.3 THz 以下的模式后 κ_zz 只变 −0.3 %），可以使用。</li>
</ul>
<div class="tbl"><table class="lt">
<thead><tr><th>300 K 推荐值</th><th class="r">{L[M4]}</th><th class="b">{L[M5]}</th></tr></thead>
<tbody>
<tr><td>面内 κ</td><td>LBTE {sv[(M4, 'inplane_LBTE_14')]:.3f} × {sv[(M4, 'mesh_correction_inplane')]:.3f}（网格修正）+ 相干项 {sv[(M4, 'coherence_inplane')]:.3f} = <b>{q4:.2f}</b></td>
<td>RTA {sv[(M5, 'inplane_RTA_32')]:.3f} × {r5:.3f}（LBTE/RTA）+ 相干项 {sv[(M5, 'coherence_inplane')]:.3f} = <b>{q5:.2f}</b></td></tr>
<tr><td>κ_zz</td><td>收敛项 {sv[(M4, 'population_kzz_rest')]:.3f} + 相干项 {sv[(M4, 'coherence_kzz')]:.3f} = <b>{z4:.2f}</b></td>
<td>收敛项 {sv[(M5, 'population_kzz_rest')]:.3f} + 相干项 {sv[(M5, 'coherence_kzz')]:.3f} = <b>{z5:.2f}</b></td></tr>
</tbody></table></div>
<p class="cap">这就是第 1 节的数值。第 4–5 节的曲线只含 RTA 粒子项，所以第 1 节的数值落在曲线上方。</p>"""


def checks_zh():
    c = en.CHK
    lc, le = c["lattice"]["calc"], c["lattice"]["exp"]
    vc = c["velocity_calc"]
    e1, e2 = c["velocity_exp"]
    pct = lambda a, b: f"{100 * (a / b - 1):+.1f} %"
    raman = "，".join(f"{m['irrep']} {m['cm-1']:.0f}" for m in c["raman_calc"])
    km, kc = c["kmin_exp_velocities"], c["kmin_calc_velocities"]
    return f"""<h2><span class="n">7</span>与实验的独立对比（{L[M4]}）</h2>
<p class="lead">这些实验数据没有参与任何拟合，用来检验力常数，并给出 κ 的下限参考。由 <code>report/independent_checks.py</code>
基于 4×4×2 二阶力常数计算。</p>
<div class="tbl"><table>
<thead><tr><th>物理量</th><th>实验</th><th>计算</th><th>计算 / 实验</th></tr></thead>
<tbody>
<tr><td>晶格常数 a（Å，六方设定）</td><td>{le['a']:.2f}</td><td>{lc['a']:.3f}</td><td>{pct(lc['a'], le['a'])}</td></tr>
<tr><td>晶格常数 c（Å）</td><td>{le['c']:.2f}</td><td>{lc['c']:.2f}</td><td>{pct(lc['c'], le['c'])}</td></tr>
<tr><td>Raman 峰（cm⁻¹）</td><td>≈ 78、≈ 177（从幻灯片图上估读）</td><td>A1g 75、A1g 177</td><td>在估读误差内</td></tr>
<tr><td>纵波声速 v_L（m/s）</td><td>{e1['vL_ms']:.0f} / {e2['vL_ms']:.0f}</td><td>{vc['avg_vL']}（c 轴 {vc['c_axis'][2]:.0f}，面内 {vc['in_plane'][2]:.0f}）</td><td>{pct(vc['avg_vL'], (e1['vL_ms'] + e2['vL_ms']) / 2)}</td></tr>
<tr><td>横波声速 v_T（m/s）</td><td>{e1['vT_ms']:.0f} / {e2['vT_ms']:.0f}</td><td>{vc['avg_vT']}（范围 {min(vc['in_plane'][0], vc['c_axis'][0]):.0f}–{max(vc['in_plane'][1], vc['c_axis'][1]):.0f}）</td><td>{pct(vc['avg_vT'], (e1['vT_ms'] + e2['vT_ms']) / 2)}</td></tr>
<tr><td>Debye 平均声速 v_m（m/s）</td><td>{e1['v_mean_ms']:.0f} / {e2['v_mean_ms']:.0f}</td><td>{vc['debye_mean']}</td><td>{pct(vc['debye_mean'], (e1['v_mean_ms'] + e2['v_mean_ms']) / 2)}</td></tr>
<tr><td>Cahill 最小热导率 κ_min（W m⁻¹ K⁻¹）</td><td>{km['cahill']:.2f}（用实测声速）</td><td>{kc['cahill']:.2f}（用计算声速）</td><td></td></tr>
<tr><td>diffuson 极限（W m⁻¹ K⁻¹）</td><td>{km['diffuson']:.2f}（用实测声速）</td><td>{kc['diffuson']:.2f}（用计算声速）</td><td></td></tr>
</tbody></table></div>
<ul>
<li><b>晶格和 Raman</b>：弛豫结构与单晶 XRD 相差不到 0.7 %，两个强 Raman 峰对应计算的 A1g 模
（全部 Raman 活性模：{raman} cm⁻¹）。两者都支持力常数是可靠的。</li>
<li><b>声速</b>：超声法，热压块体，两个样品（相差 2 % 以内）。实测 v_L 低于<i>所有</i>方向的计算 v_L
（最低的 c 轴方向也有 {vc['c_axis'][2]:.0f} m/s），所以织构解释不了；更可能是块体的孔隙和缺陷软化。
块体密度还不知道（理论密度 {c['theoretical_density_g_cm3']:.2f} g/cm³），两者暂时无法区分。声波传播方向相对压力方向也还不清楚。</li>
<li><b>最小热导率</b>（Cahill、Watson、Pohl，PRB 46, 6131 (1992)；diffuson 模型，Agne 等，Energy Environ. Sci. 11, 609 (2018)），
原子数密度 n = {c['number_density_A-3']:.4f} Å⁻³：实测 κ_tot（303 K 下 Perp 0.93、Para 0.49）落在 diffuson 值 {km['diffuson']:.2f}
和 Cahill 值 {km['cahill']:.2f} 之间。热压样品已接近非晶极限，这独立地印证了 L3 拟合所需的强外禀散射（L 约几十 nm）。
这两个极限是按各向同性推导的；层状材料的 Para 低于 Cahill 值很常见。</li>
</ul>"""


def html_zh(figs):
    T1, T2 = en.table1_body({"quoted": "计算，推荐值", "ordered": "计算，有序（RTA 20³）",
                             "vca": "计算，VCA 反位无序（RTA 20³）", "exp": "实验，Perp / Para（303 K）",
                             "exp_avg": "实验，(2·Perp + Para)/3（303 K）"},
                            {M4: "截断约 −10 %", M5: "截断 ≥ −3 %"})
    ORD, O = en.ordering_rows({"exp": "实验，Perp / Para（303 K）",
                               "quoted": "计算，推荐值，面内 / κ_zz",
                               "ordered": "计算，有序 RTA 20³，面内 / κ_zz",
                               "mixed": "计算，推荐值按织构 s 混合 → Perp / Para"})
    sw, sA, sB = O["sweep"], O["s"][M4], O["s"][M5]
    d = lambda m, x0, x1, i: abs(100 * (sw[m][x1][i] / sw[m][x0][i] - 1))
    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>MnIn₂Se₄ / Mn₂In₂Se₅ 计算与实验对比</title>
<style>{CSS}</style></head><body><div class="wrap">

<div class="eyebrow">晶格热导率 · 由 report/make_report_html_zh.py 生成</div>
<h1>{L[M4]} 与 {L[M5]}：计算与实验对比</h1>
<p class="lead">内容包括声子谱、本征 κ(T)、与热压样品实验的对比，以及 300 K 下的三声子和无序散射率。所有数值都来自仓库中已有的文件
（见 <a href="../PROGRESS_REPORT_20261118.md">PROGRESS_REPORT_20261118.md</a>）；只有 {L[M5]} 的声子谱是在这里由仓库中的
<code>fc2.hdf5</code> 插值得到的。</p>

<h2><span class="n">1</span>300 K 计算值与实验</h2>
<p class="lead">单位 W m⁻¹ K⁻¹，计算为 300 K，实验为 303 K（热压样品，κ_tot）。两个方向分开对比，每个数后面是 计算/实验 的倍数。
实验中 Perp 垂直于热压方向测量，主要反映面内；Para 沿热压方向测量，主要反映层间（热压时片状晶粒的 c 轴倾向于沿压力方向）。</p>
<div class="tbl"><table>
<thead><tr><th></th><th class="r">{L[M4]} 面内</th><th class="r">{L[M4]} 层间 κ_zz</th>
<th class="b">{L[M5]} 面内</th><th class="b">{L[M5]} 层间 κ_zz</th></tr></thead>
<tbody>
{T1}
</tbody></table></div>
<p class="cap"><b>计算推荐值</b>：最精确的单晶值（第 2.4 节：面内 LBTE + 相干项，κ_zz 网格收敛项 + 相干项）。{L[M4]} 的三阶截断误差约 −10 %；
{L[M5]} 下一步截断不低于 −2.1 % / −3.3 %（面内 / κ_zz）。<b>有序 / VCA</b>：20³ 网格的 RTA，分别不含 / 含 Mn–In 反位质量无序（第 5 节），用来看无序本身能使 κ 降低多少。
所有计算值都对应无缺陷、无晶界的单晶。</p>

<p class="lead"><b>为什么还有第二张表。</b>热压样品只是部分织构，Perp 不完全是面内，Para 也不完全是 κ_zz：有一部分面内导热混进了 Para
（所以上表中 Para 对应的倍数偏小）。而平均值 (2·Perp + Para)/3 对<i>任意</i>织构都等于 (2·κ面内 + κ_zz)/3，
因此不需要知道织构就能直接比较计算和实验。</p>
<div class="tbl"><table>
<thead><tr><th>方向平均 (2·面内 + κ_zz)/3</th><th class="r">{L[M4]}</th><th class="b">{L[M5]}</th></tr></thead>
<tbody>
{T2}
</tbody></table></div>
<p class="cap">单晶 κ 比热压样品高 3–4 倍。反位无序（VCA）缩小了一部分差距，其余需要纳米尺度的有效平均自由程（L_eff ≈ 10–50 nm，第 8 节）；可能的原因见第 9 节。
第 8 节的 L3 模型没有列在这里：它是用这些实验点拟合出来的，300 K 下与实验吻合是拟合的结果，不是独立的计算。</p>

<h3>哪种材料更低：{L[M5]} 相对 {L[M4]}</h3>
<p class="lead">实验中 {L[M5]} 在两个方向上都更低。计算在面内和方向平均上重现了这一点，但 κ_zz 没有。✓ / ✗：与实验同号 / 反号。</p>
<div class="tbl"><table>
<thead><tr><th>300 K，{L[M5]} 相对 {L[M4]}</th><th>面内 / Perp</th><th>κ_zz / Para</th><th>方向平均</th></tr></thead>
<tbody>
{ORD}
</tbody></table></div>
<ul>
<li><b>κ_zz 不应直接与 Para 比较。</b>Para 包含了一部分面内导热：Para = s·κ面内 + (1−s)·κ_zz。两种材料取相同的 s 时，只要
s &gt; {O['s_cross']:.2f}，{L[M5]} 的 Para 就更低。拟合得到的织构为 s = {sA:.2f}（{L[M4]}）和 {sB:.2f}（{L[M5]}）；用它们混合后，
计算值在两个方向上都给出与实验相同的顺序（最后一行：Para {O['mixed_abs'][M5][1]:.2f} 对 {O['mixed_abs'][M4][1]:.2f}）。
测块体两个面的 XRD 可以独立确定 s。</li>
<li><b>{L[M5]} 的 κ_zz 有界。</b>三阶截断步长 −11.9 %、−9.9 % 之后，下一步不超过 −3.3 %（第 2.2 节，留出验证），
所以在数据能检验的范围内，它的 κ_zz 比 {L[M4]} 高约 50 %。</li>
<li><b>点缺陷几乎降不了 κ_zz。</b>反位质量无序从有序到完全无序，{L[M5]} 面内降低 {d(M5, 0.0, 0.5, 0):.0f} %，而 κ_zz 只降低
{d(M5, 0.0, 0.5, 1):.0f} %（{L[M4]}，x = 0 → 0.44：分别降低 {d(M4, 0.0, 0.44, 0):.0f} % 和 {d(M4, 0.0, 0.44, 1):.0f} %；
<code>plot_data/disorder_sweep_300K.csv</code>）。κ_zz 由低频层间声学模承担，质量无序对它们散射很弱。所以 {L[M5]} 块体的 κ_zz 更低，
说明存在沿 c 方向的缺陷，例如层错、多型体共生，或者它在 400–800 °C 分解产生的第二相（第 9 节第 1、4 条）。</li>
</ul>

{methods_zh()}

<h2><span class="n">3</span>声子谱</h2>
<p class="lead">R-3m 原胞，路径 Γ–T–H₂|H₀–L–Γ–S₀|S₂–F–Γ。两种材料都没有虚频；不加 NAC。</p>
<figure><img alt="两种材料的声子谱" src="{figs['bands']}"></figure>
<p class="cap"><b>{L[M4]}</b>：4×4×2 二阶力常数，<code>test/MnIn2Se4/fc2_442_phonon/band-dft-cpu.yaml</code>。
<b>{L[M5]}</b>：4×4×1 二阶力常数，<code>Mn2In2Se5/fc/fc2.hdf5</code>（phonopy 插值）。
{L[M4]} 的 4×4×1 二阶力常数给出的层间横波过软（322.9 m/s），已弃用；4×4×2 给出 1348.3 m/s。
{L[M5]} 的 4×4×1 二阶力常数用同样方法检验过，没有问题（1386.8 m/s）。</p>

<h2><span class="n">4</span>本征 κ(T)</h2>
<p class="lead">有序晶体，三声子 + 天然同位素散射，RTA。κ_zz 取网格收敛项（去掉 Γ–Z 积分误差）。星号为第 1 节的 300 K 推荐值。</p>
<figure><img alt="本征热导率随温度变化" src="{figs['intrinsic']}"></figure>
<p class="cap"><b>曲线</b>：RTA 粒子项，每个温度都计算了（网格 {L[M4]} 14³，{L[M5]} 16³）。
<b>星号</b>：第 1 节的 300 K 推荐值，用的是第 2.4 节中更精确的处理（面内用 LBTE 代替 RTA，两个方向都加上 Wigner 相干项）。
这种处理只在 300 K 计算了，所以只有一个点。300 K 处星号比曲线高：{L[M4]} 面内 +28 %、κ_zz +11 %；{L[M5]} 面内 +14 %、κ_zz +8 %；
其他温度下曲线也大致需要同样的修正。引用数值请用第 1 节的表，不要从曲线上读。
误差棒：{L[M4]} 的三阶截断估计（约 −10 %）。{L[M5]}：留出验证给出的下一步截断的界（−2.1 % / −3.3 %）。</p>

<h2><span class="n">5</span>计算与热压样品实验对比</h2>
<p class="lead">每个面板对应一种材料的一个方向。虚线：有序单晶。实线：VCA 反位无序单晶。点线：L3 拟合（第 8 节）。
点划线（仅 {L[M4]}）：由实测声速得到的最小热导率（第 7 节）。实心点：300–623 K（参与拟合）。空心点：55–145 K（未参与拟合）。</p>
<figure><img alt="计算热导率与热压样品实验对比" src="{figs['calc_vs_exp']}"></figure>
<p class="cap">有序和 VCA：RTA，20³（每对用同一网格），Perp 对应面内、Para 对应 κ_zz，不含晶界和织构，因此是上限参考。
VCA 占位：{L[M4]} i = 0.44（单晶精修，Range 等，Z. Naturforsch. 46b, 1122 (1991)），{L[M5]} x = 0.19(2)（R-3m 相的粉末 Rietveld，
Williams 等，J. Mater. Chem. C 12, 1753 (2024)，ESI 表 S7；见第 8 节）。VCA 显著降低了 κ，但仍远高于实验，而且给不出 100 K 附近的低温极大值，这需要纳米尺度的平均自由程。</p>

<h2><span class="n">6</span>300 K 散射率</h2>
<p class="lead">最大三阶截断，20³ 网格，τ⁻¹ = 2Γ。基线：真实质量 + 天然同位素。VCA：位点平均质量 + 质量方差散射。
上排 {L[M4]}（i = 0.44），下排 {L[M5]}（x = 0.13：x = 0.19 那次计算的逐模式文件还不在仓库里；下面的 κ 值是 x = 0.19 的）。</p>
<figure><img alt="两种材料的三声子与无序散射率" src="{figs['scattering']}"></figure>
<p class="cap">左：有序和 VCA 晶体的三声子散射率相差不超过几十个百分点，因为 VCA 改变的是质量，不是非谐性。
右：VCA 质量方差散射比天然同位素散射强 2–3 个数量级，低频按 ω² 增长，高频趋平。对 300 K 热导率的影响：
{L[M4]} 面内 3.69 → 1.49（2.48 倍），κ_zz 0.638 → 0.445（1.43 倍）；{L[M5]}（x = 0.19）面内 3.00 → 1.46（2.05 倍），
κ_zz 0.984 → 0.753（1.31 倍）（x = 0.13 时为 1.63 和 0.792）。{L[M5]} 的影响较弱，是因为它的位点混合较少（x = 0.19 时 81 % 是纯的，
而 {L[M4]} 八面体位是 56/44）。</p>

{checks_zh()}

<h2><span class="n">8</span>L3 拟合是什么</h2>
<p class="lead">L3 是 <code>analysis/extrinsic/</code> 中的外禀散射模型。<b>它是用实验值拟合的</b>，因此它和实验吻合是拟合出来的结果，不是预测。
300 K 下它给出 Perp / Para = 0.92 / 0.49（{L[M4]}，实测 0.93 / 0.49）和 0.78 / 0.34（{L[M5]}，实测 0.74 / 0.32）。
<b>它的用途不是给出 κ</b>（实验已经测出来了），而是回答"需要多少、什么样的额外散射，才能把计算的单晶变成实测的热压块体"：
拟合出的参数（织构 s、晶界长度 L；反位比例 x 取自结构精修）才是结果，而没有参与拟合的 55–145 K 数据用来检验它。</p>
<div class="eq">Γ_λ = Γ_3ph,λ(T)  +  Σ g_elem(x)·B_elem,λ  +  |v_λ| / (4π L)
       三声子              Mn↔In 反位质量无序             晶界（Casimir）

κ_Para = s·κ_xx + (1−s)·κ_zz          κ_Perp = (1−s/2)·κ_xx + (s/2)·κ_zz</div>
<ul>
<li><b>不参与拟合的输入</b>：第一性原理的三声子散射率、声子群速度和热容（有序质量），以及质量无序的散射率基函数 B，每种材料在 24³ 网格上算一次。</li>
<li><b>参数</b>：x，Mn 占 In 位的比例，<b>固定</b>为文献值（见下）；L，有效晶界平均自由程，和 s = ⟨sin²θ⟩，织构参数（0 = 完全织构，2/3 = 随机取向），<b>参与拟合</b>。</li>
<li><b>使用的数据</b>：只用 <code>experiment/kappa_tot_digitized.csv</code> 中两个方向 300–623 K 的点。对 L 使
ln(模型/实验) 的均方根最小，每个点上 s 单独优化。55–145 K 的点不参与拟合，只用来检验模型。</li>
<li><b>结果</b>（analysis/disorder/README.md §5.1）：{L[M4]} 取 x = 0.44（单晶 XRD，Range 等，Z. Naturforsch. 46b, 1122 (1991)）：
L = 30.7 nm，s = 0.35，误差 1.2 %。{L[M5]} 取 x = 0.19(2)（R-3m 相的粉末 Rietveld，Williams 等，J. Mater. Chem. C 12, 1753 (2024)，ESI 表 S7）：
L = 22.1 nm，s = 0.23，误差 3.4 %。两种材料的 L 都约为 20–30 nm。
（Chernoukhov 等，Molecules 30, 1904 (2025) 的 0.704/0.296 属于另一种淬火得到的多晶型 P3̄m1，不采用。）</li>
<li><b>为什么固定 x</b>：如果 x 也放开拟合，只有 s 能确定（可接受范围 0.33–0.38 和 0.15–0.28）。x 和 L 都压低 κ，彼此可以替代：
可接受的 x 分别为 0–0.66 和 0.01–0.50，自由拟合的最优点（x = 0.45，L = 30.7 nm，误差 1.2 %；x = 0.13，L = 15.9 nm，误差 3.2 %）
不能作为占位的证据。XRD（织构）、SEM（晶粒尺寸）和我们自己样品的 Rietveld 占位精修可以打破这种简并。</li>
<li><b>不能与 VCA 叠加</b>：L3 已经通过 x 包含了反位散射，所以 VCA 曲线和 L3 曲线是二选一，不是可以相加的修正。</li>
<li><b>κ_e</b> 取为 0。如果电子热导不可忽略，拟合出的外禀散射会偏强。</li>
</ul>

<h2><span class="n">9</span>计算值为什么高于实验值</h2>
<p class="lead">300 K 下，两种材料的单晶计算值都比热压样品的 κ_tot 高 4.3 倍（方向平均，推荐值，第 1 节）。考虑反位无序（VCA，取文献占位）后降到
1.5 倍（{L[M4]}）和 2.0 倍（{L[M5]}；RTA，20³）。剩下差距的可能原因列在下面，大致按能解释的程度排序，
并给出依据和检验方法。<b>这些都是猜想</b>，只有前两项在这里做了定量估计。</p>
<div class="tbl"><table class="lt">
<thead><tr><th>#</th><th>机制</th><th>对 κ 的可能影响</th><th>依据</th><th>如何检验</th></tr></thead>
<tbody>
<tr><td>1</td><td><b>面缺陷和晶界</b>：纳米尺度的有效平均自由程</td>
<td>最大；L ≈ 20–30 nm 就能补上剩余差距（L3，x 取文献值，第 8 节）</td>
<td>150 K 以下实测 κ 基本是平的，没有 1/T 上升，说明主导的是与温度无关的散射；所需的 L 远小于热压样品常见的晶粒尺寸（0.5–5 μm），
更像是层错、孪晶或多型体共生，这在 "124"/"225" 这类层状化合物中很常见；κ_tot 已落在 Cahill 和 diffuson 极限之间（第 7 节）</td>
<td>SEM 晶粒尺寸；TEM/HAADF 观察层错</td></tr>
<tr><td>2</td><td><b>Mn/In 反位无序</b></td>
<td>按文献 x 的 VCA：300 K 面内降低 2.5 倍（{L[M4]}）/ 2.1 倍（{L[M5]}），κ_zz 1.4 倍 / 1.3 倍</td>
<td>单晶精修 i = 0.44（Range 1991）；Rietveld x = 0.19(2)（Williams 2024）；{L[M5]} 的 Raman 峰明显展宽，被归因于位点无序（Williams 2024）</td>
<td>对我们自己的块体做 Rietveld 占位精修。VCA 只算了质量差，是下限：键强和局部应变的无序没有包含（无序超胞计算可以包含）。
VCA 的 LBTE + 相干项计算已在计划中（<code>test/figures/kappa_vs_T/jobs_lbte_wigner/</code>）</td></tr>
<tr><td>3</td><td><b>顺磁自旋无序</b>（Mn²⁺，S = 5/2）</td>
<td>未知；可能达到几十个百分点</td>
<td>计算假设了静态铁磁序；实际晶体在 300 K 是顺磁的，只在约 6–7 K 以下才有序（反铁磁 / 自旋玻璃），而且反铁磁交换很强
（{L[M5]} 的 Weiss 温度约 −198 K，Williams 2024）。涨落的自旋使依赖交换作用的力常数变得随机，并带来自旋–声子散射</td>
<td>测块体的磁化率；计算反铁磁或随机自旋构型下的力常数</td></tr>
<tr><td>4</td><td><b>偏离化学计量、空位、第二相</b></td>
<td>未知；阳离子空位的散射远强于质量替代</td>
<td>{L[M4]} 在富 In₂Se₃ 一侧有 5–6 mol% 的均匀区（Mammadov 2025），即存在阳离子空位；{L[M5]} 在约 400–800 °C 之间会分解为
MnIn₂Se₄ + MnSe，且常含 MnSe 杂相（Chernoukhov 2025），在这个温度范围内热压可能留下第二相和界面</td>
<td>热压后块体的 XRD（物相纯度）、EDS 成分、热压温度和时间</td></tr>
<tr><td>5</td><td><b>孔隙、块体软化</b></td>
<td>孔隙率 5–10 % 时约 10–20 %</td>
<td>实测声速在所有方向上都比计算低 20–30 %（第 7 节）</td>
<td>阿基米德法测密度</td></tr>
<tr><td>6</td><td><b>高阶非谐和温度效应</b></td>
<td>在低热导层状硒化物中 300 K 时通常为 −10 到 −30 %（这里没有计算）</td>
<td>存在低频光学支（第 3 节）且三声子散射率高；计算只用了 0 K 简谐声子和三声子散射</td>
<td>四声子散射、温度相关力常数（计算量大）</td></tr>
<tr><td>7</td><td><b>三阶力常数截断</b></td>
<td>{L[M4]} 约 −10 %；{L[M5]} 下一步不低于 −3.3 %</td>
<td>第 2.2 节</td><td>{L[M5]} 已做（重拟合到 7.43 Å + 留出验证）；{L[M4]} 对 7.74 Å 的重拟合做同样的检验</td></tr>
</tbody></table></div>
<p class="lead"><b>可以排除的原因</b></p>
<ul>
<li><b>织构</b>只在 Perp 和 Para 之间重新分配导热，方向平均 (2·Perp + Para)/3 与织构无关，而它仍然低了 4.3 倍（第 1 节）。</li>
<li><b>电子热导 κ_e</b>：实验值是 κ_tot，扣除 κ_e 后实验的晶格热导更低，差距反而更大。</li>
<li><b>晶格常数</b>：计算的晶胞比实测大 0.6–0.7 %（第 7 节），晶胞偏大一般会使计算的 κ 偏低而不是偏高。</li>
<li><b>LBTE 和相干项</b>会使计算值更高，差距更大。<b>数字化误差</b>（±0.005–0.01 W m⁻¹ K⁻¹）可以忽略。</li>
</ul>
<p class="cap">仅机制 1 和 2 就已能复现实验数据（这就是 L3 拟合）。机制 3–6 中只要有一项真实存在，所需的晶界散射就会减少。
拟合出的 L ≈ 20–30 nm 把模型没有包含的所有散射都归到了一起，所以实际面缺陷的间距至少有这么大。
最省事又最有判别力的测量是：块体密度、SEM/TEM 显微结构、热压后块体的 XRD。</p>

<h2><span class="n">10</span>注意事项</h2>
<ul>
<li>实验点是从图上数字化得到的（高温段 ±0.005 W m⁻¹ K⁻¹，低温段 ±0.01），且为 κ_tot，未扣除 κ_e。拿到原始数据后应替换
<code>experiment/kappa_tot_digitized.csv</code>。</li>
<li>{L[M5]}：三阶截断只检验到 6.90 Å（下一步不低于 −2.1 % / −3.3 %），更大的截断现有数据支撑不了。
在这个范围内两个顺序都成立：面内 {L[M4]} 高 37 %，κ_zz {L[M5]} 高约 50 %。</li>
<li>实验在 145 K 和 300 K 之间，Para 下降得比 Perp 多得多（{L[M5]} −32 % 对 −6 %），可能是两台仪器或两个样品，所以拟合只用高温段。</li>
</ul>
</div></body></html>
"""


def main():
    png = plot_figures.make_all(en.DATA, en.FIGDIR, png_for_html=True, lang="zh")
    figs = {k: "data:image/png;base64," + base64.b64encode(v).decode() for k, v in png.items()}
    OUT.write_text(html_zh(figs), encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e6:.2f} MB) and {en.FIGDIR}/ (*_zh.png 300 dpi + .pdf)")


if __name__ == "__main__":
    main()
