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
{L[M5]} 的步长没有缩小（面内 −3.9 % → −4.3 %，κ_zz −11.9 % → −9.9 %），所以它的数值按<b>上限</b>报告。</p>

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

<h2><span class="n">1</span>300 K 推荐值与实验</h2>
<p class="lead">单位 W m⁻¹ K⁻¹。实验中 Perp（垂直于热压方向）≈ 面内，Para（平行于热压方向）≈ 层间（热压时片状晶粒的 c 轴倾向于沿压力方向）。
倍数 = 计算 / 实验。</p>
<div class="tbl"><table>
<thead><tr><th></th><th class="r">{L[M4]}</th><th class="b">{L[M5]}</th></tr></thead>
<tbody>
<tr><td>计算 面内 κ（LBTE + 相干项）</td><td><b>4.67</b> <span class="x">截断约 −10 %</span></td><td><b>3.37</b> <span class="x">上限</span></td></tr>
<tr><td>计算 κ_zz（网格收敛 + 相干项）</td><td><b>0.71</b> <span class="x">截断约 −10 %</span></td><td><b>1.09</b> <span class="x">上限</span></td></tr>
<tr><td>实验 κ_tot Perp / Para（303 K，热压）</td><td>0.933 / 0.485</td><td>0.741 / 0.320</td></tr>
<tr><td>与织构无关的平均 (2·xx + zz)/3：实验</td><td>0.784</td><td>0.601</td></tr>
<tr><td>&nbsp;&nbsp;计算，推荐值</td><td>3.35 <span class="x">4.3×</span></td><td>2.61 <span class="x">4.3×</span></td></tr>
<tr><td>&nbsp;&nbsp;计算，有序，RTA 20³</td><td>2.67 <span class="x">3.4×</span></td><td>2.33 <span class="x">3.9×</span></td></tr>
<tr><td>&nbsp;&nbsp;计算，VCA 反位无序，RTA 20³</td><td>1.14 <span class="x">1.5×</span></td><td>1.35 <span class="x">2.3×</span></td></tr>
<tr><td>&nbsp;&nbsp;计算，L3 拟合（用这些实验数据拟合，见第 8 节）</td><td>0.77</td><td>0.62</td></tr>
</tbody></table></div>
<p class="cap">本征单晶 κ 比热压样品高 3–4 倍。Mn/In 反位质量无序（VCA）缩小了一部分差距，其余需要纳米尺度的有效平均自由程（L_eff ≈ 10–50 nm）。</p>

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
误差棒：{L[M4]} 的三阶截断估计（约 −10 %）。箭头：{L[M5]} 截断未收敛，数值为上限。</p>

<h2><span class="n">5</span>计算与热压样品实验对比</h2>
<p class="lead">每个面板对应一种材料的一个方向。虚线：有序单晶。实线：VCA 反位无序单晶。点线：L3 拟合（第 8 节）。
点划线（仅 {L[M4]}）：由实测声速得到的最小热导率（第 7 节）。实心点：300–623 K（参与拟合）。空心点：55–145 K（未参与拟合）。</p>
<figure><img alt="计算热导率与热压样品实验对比" src="{figs['calc_vs_exp']}"></figure>
<p class="cap">有序和 VCA：RTA，20³（每对用同一网格），Perp 对应面内、Para 对应 κ_zz，不含晶界和织构，因此是上限参考。
VCA 占位：{L[M4]} i = 0.44（单晶精修，Range 等，Z. Naturforsch. 46b, 1122 (1991)），{L[M5]} x = 0.13（取自 L3 拟合；
该 R-3m 相的文献值为 0.19(2)，见第 8 节）。VCA 显著降低了 κ，但仍远高于实验，而且给不出 100 K 附近的低温极大值，这需要纳米尺度的平均自由程。</p>

<h2><span class="n">6</span>300 K 散射率</h2>
<p class="lead">最大三阶截断，20³ 网格，τ⁻¹ = 2Γ。基线：真实质量 + 天然同位素。VCA：位点平均质量 + 质量方差散射。
上排 {L[M4]}（i = 0.44），下排 {L[M5]}（x = 0.13）。</p>
<figure><img alt="两种材料的三声子与无序散射率" src="{figs['scattering']}"></figure>
<p class="cap">左：有序和 VCA 晶体的三声子散射率相差不超过几十个百分点，因为 VCA 改变的是质量，不是非谐性。
右：VCA 质量方差散射比天然同位素散射强 2–3 个数量级，低频按 ω² 增长，高频趋平。对 300 K 热导率的影响：
{L[M4]} 面内 3.69 → 1.49（2.48 倍），κ_zz 0.638 → 0.445（1.43 倍）；{L[M5]} 面内 3.00 → 1.63（1.84 倍），
κ_zz 0.984 → 0.792（1.24 倍）。{L[M5]} 的影响较弱，是因为它的 x 较小。</p>

{checks_zh()}

<h2><span class="n">8</span>L3 拟合是什么</h2>
<p class="lead">L3 是 <code>analysis/extrinsic/</code> 中的外禀散射模型。<b>它是用实验值拟合的</b>，因此它和实验吻合是拟合出来的结果，不是预测。</p>
<div class="eq">Γ_λ = Γ_3ph,λ(T)  +  Σ g_elem(x)·B_elem,λ  +  |v_λ| / (4π L)
       三声子              Mn↔In 反位质量无序             晶界（Casimir）

κ_Para = s·κ_xx + (1−s)·κ_zz          κ_Perp = (1−s/2)·κ_xx + (s/2)·κ_zz</div>
<ul>
<li><b>不参与拟合的输入</b>：第一性原理的三声子散射率、声子群速度和热容（有序质量），以及质量无序的散射率基函数 B，每种材料在 24³ 网格上算一次。</li>
<li><b>三个拟合参数</b>：x，Mn 占 In 位的比例；L，有效晶界平均自由程；s = ⟨sin²θ⟩，织构参数（0 = 完全织构，2/3 = 随机取向）。</li>
<li><b>使用的数据</b>：只用 <code>experiment/kappa_tot_digitized.csv</code> 中两个方向 300–623 K 的点。在 (x, L) 网格上使
ln(模型/实验) 的均方根最小，每个点上 s 单独优化。55–145 K 的点不参与拟合，只用来检验模型。</li>
<li><b>结果</b>：{L[M4]} x = 0.45，L = 30.7 nm，s = 0.35，误差 1.2 %；{L[M5]} x = 0.13，L = 15.9 nm，s = 0.25，误差 3.2 %。</li>
<li><b>它能确定什么</b>：只有 s（可接受范围 0.33–0.38 和 0.15–0.28）。x 和 L 都压低 κ，彼此可以替代：可接受的 x 分别为 0–0.66 和 0.01–0.50。
所以拟合出的 x 不能作为占位的证据，0.45 接近单晶的 0.44 只是巧合。XRD（织构）、SEM（晶粒尺寸）和 Rietveld 占位精修可以打破这种简并。</li>
<li><b>用文献值固定 x</b>（analysis/disorder/README.md §5.1）：{L[M4]} x = 0.44（单晶 XRD，Range 等，Z. Naturforsch. 46b, 1122 (1991)）
给出 L = 30.7 nm，误差 1.2 %。{L[M5]} x = 0.19(2)（R-3m 相的粉末 Rietveld，Williams 等，J. Mater. Chem. C 12, 1753 (2024)，ESI 表 S7）
给出 L = 22.1 nm，误差 3.4 %。固定 x 后只拟合 L 和 s，两种材料的 L 都约为 20–30 nm。
（Chernoukhov 等，Molecules 30, 1904 (2025) 的 0.704/0.296 属于另一种淬火得到的多晶型 P3̄m1，不采用。）</li>
<li><b>不能与 VCA 叠加</b>：L3 已经通过 x 包含了反位散射，所以 VCA 曲线和 L3 曲线是二选一，不是可以相加的修正。</li>
<li><b>κ_e</b> 取为 0。如果电子热导不可忽略，拟合出的外禀散射会偏强。</li>
</ul>

<h2><span class="n">9</span>注意事项</h2>
<ul>
<li>实验点是从图上数字化得到的（高温段 ±0.005 W m⁻¹ K⁻¹，低温段 ±0.01），且为 κ_tot，未扣除 κ_e。拿到原始数据后应替换
<code>experiment/kappa_tot_digitized.csv</code>。</li>
<li>{L[M5]}：三阶截断的步长没有缩小（面内 −3.9 % → −4.3 %，κ_zz −11.9 % → −9.9 %），所以它的数值是上限，两种材料 κ_zz 的大小顺序不可靠；
面内的顺序（{L[M4]} 高 37 %）是可靠的。</li>
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
