# MnIn₂Se₄ / Mn₂In₂Se₅ lattice thermal conductivity — progress report

**Date:** 2026-11-18  
**Scope:** both materials. This report supersedes the status sections of
earlier summaries. The detailed derivations live in
[`REPORT_MnIn2Se4_20261115.md`](REPORT_MnIn2Se4_20261115.md); section numbers (§) below refer to it.

---

## 1. Summary

- **The calculation is complete within the agreed ≈ 10 % target.** The one
  exception is the Mn₂In₂Se₅ FC3-cutoff convergence. The q-mesh, isotopes,
  LBTE versus RTA, and the coherence term have all been settled.
- **The remaining uncertainty is the FC3 cutoff**, and further calculation has
  been ruled out (decision 2026-11-18):
  - MnIn₂Se₄: the cutoff steps shrink, and the remaining error is ≲ 10 %.
  - Mn₂In₂Se₅: the steps do not shrink, so its values are upper bounds.
- **Intrinsic single-crystal κ is 3–4× above the hot-pressed measurements.**
  Mn/In antisite mass disorder (VCA) closes part of the gap. A nanoscale
  phonon mean free path (L_eff ≈ 10–50 nm) is needed to reach the measured
  values. The low-temperature data show this most clearly.

## 2. Values to quote (300 K, W m⁻¹ K⁻¹)

| | MnIn₂Se₄ | Mn₂In₂Se₅ |
|---|---:|---:|
| **in-plane κ** (LBTE + coherence) | **≈ 4.6** (4.67) | **≈ 3.4** (3.37) |
| **κ_zz** (mesh-converged + coherence) | **≈ 0.71** (0.648 + 0.066) | **≈ 1.09** (0.977 + 0.112) |
| κ_zz at 10 μm / 1 μm grain (no coherence) | 0.638 / 0.58 | 0.944 / 0.798 |
| FC3 cutoff, in-plane | ≈ −10 %, estimated (steps −22.9 % → −8.7 %) | upper bound (steps −3.9 % → −4.3 %) |
| FC3 cutoff, κ_zz | ≈ −10 %, estimated (steps −50.5 % → −8.2 %) | **upper bound** (steps −11.9 % → −9.9 %) |
| experiment κ_tot, Perp / Para (303 K, hot-pressed) | 0.933 / 0.485 | 0.741 / 0.320 |

What each number is made of:

- **MnIn₂Se₄ in-plane**: LBTE 4.5219 at 14³, +1.3 % mesh correction from
  the RTA series, + coherence 0.09 (phono3py 3.24.0).
- **Mn₂In₂Se₅ in-plane**: RTA 3.048 (32³), with the iterative and coherence
  corrections (`analysis/disorder/README.md` §4.4).
- **κ_zz**: mesh-converged rest term, i.e. the Γ–Z quadrature error removed
  (§6.3–6.5).
  - MnIn₂Se₄: 0.648 ± 0.003 (14³–30³).
  - Mn₂In₂Se₅: 0.977 ± 0.008 (plateau over 24³–32³, isotopes on).

Two comparisons to keep in mind:

- **In-plane ordering is robust.** MnIn₂Se₄ > Mn₂In₂Se₅ by 37 %, well above
  the cutoff effects.
- **κ_zz ordering is not robust.** Mn₂In₂Se₅ is an upper bound; a further
  20–30 % drop would bring it to the MnIn₂Se₄ value.

**Texture-independent comparison with experiment**, (2κ_xx + κ_zz)/3 against
(2·Perp + Para)/3 at 300 K:

| | MnIn₂Se₄ | Mn₂In₂Se₅ |
|---|---:|---:|
| experiment | 0.784 | 0.601 |
| calc., quoted values (LBTE + coherence) | 3.35 (4.3×) | 2.61 (4.3×) |
| calc., ordered, RTA 20³ | 2.67 (3.4×) | 2.33 (3.9×) |
| calc., VCA, RTA 20³ | 1.14 (1.5×) | 1.35 (2.3×) |
| calc., L3 extrinsic fit | 0.77 (fitted to these data) | 0.62 (fitted to these data) |

## 3. Figures

All three are produced by
[`test/figures/paper/make_paper_figures.py`](test/figures/paper/make_paper_figures.py)
from data in the repository. The κ_zz convention is the same as in §2.

### Fig. 1 — intrinsic κ(T)

![fig1](test/figures/paper/fig1_kappa_vs_T.png)

- **Lines:** RTA population, ordered crystal, natural isotopes; 14³ for
  MnIn₂Se₄, 16³ for Mn₂In₂Se₅.
- **Stars:** the 300 K values of §2.
  - The stars sit above the RTA lines because LBTE exceeds RTA in-plane (24 %
    for MnIn₂Se₄) and because they include the coherence term.
  - Bar on MnIn₂Se₄: the cutoff estimate. Arrow on Mn₂In₂Se₅: an upper bound.

### Fig. 2 — calculation against the hot-pressed experiment

![fig2](test/figures/paper/fig2_calc_vs_exp.png)

- **Dashed / solid:** ordered / VCA single-crystal curves.
  - RTA, 20³ pairs, Perp ↔ in-plane, Para ↔ κ_zz.
  - No boundary or texture, so these are an upper reference.
- **Dotted:** L3 extrinsic model, i.e. ordered-mass κ + antisite mass
  disorder + boundary mean free path + texture, fitted to 300–623 K.
- **Symbols:** experiment; filled is high T, open is low T. It is κ_tot with
  κ_e not subtracted.
- **Main reading:** VCA lowers κ strongly but stays well above experiment, and
  it cannot produce the low-temperature maximum near 100 K. That maximum needs
  a nanoscale mean free path.

### Fig. 3 — the remaining uncertainty: FC3 cutoff

![fig3](test/figures/paper/fig3_cutoff.png)

300 K, normalised to the largest cutoff used (6.45 / 6.51 Å). The MnIn₂Se₄
steps shrink; the Mn₂In₂Se₅ steps do not.

## 4. Disorder and extrinsic parameters

| | MnIn₂Se₄ | Mn₂In₂Se₅ |
|---|---|---|
| VCA (L2) occupancy | i = 0.44 (single-crystal XRD, K.-J. Range, U. Klement, G. Döll, E. Bucher, J. R. Baumann, *Z. Naturforsch.* **46b**, 1122 (1991); verified) | x = 0.13 (taken from the L3 fit). Literature for the R-3m phase: **0.19(2)** (Williams *et al.*, *J. Mater. Chem. C* **12**, 1753 (2024), ESI Table S7); see `analysis/disorder/README.md` §5.1 |
| VCA suppression at 300 K, in-plane / κ_zz | 2.48× / 1.43× | 1.84× / 1.24× |
| L3 fit: texture s | **0.35** (0.33–0.38) | **0.25** (0.15–0.28) |
| L3 fit: x, L_eff | 0.45, 30.7 nm; x **not determined** (0.00–0.66 acceptable) | 0.13, 15.9 nm; x **not determined** (0.01–0.50) |
| L3 fit: rms of ln(model/exp) | 1.2 % | 3.2 % |

- L2 (VCA) and L3 must not be added together. L3 already contains the
  antisite scattering through x (`analysis/disorder/README.md` §2).
- **x and L are degenerate in L3.** The texture s is the only quantity the κ
  data determine.

## 5. Uncertainty budget

| source | status | effect |
|---|---|---|
| q-mesh, in-plane | converged | ±0.2 % |
| q-mesh, κ_zz | converged on the rest term | ±0.5 % (MnIn₂Se₄), ±0.8 % (Mn₂In₂Se₅). The raw totals are 15–25 % high and are not used |
| natural isotopes | included | −1.2 % in-plane, −0.4 % κ_zz at 300 K; −17 % in-plane at 30 K |
| LBTE vs RTA | in-plane quoted from LBTE | +24 % in-plane (MnIn₂Se₄); κ_zz equal to 4·10⁻⁵ |
| coherence | included (rule below) | +2 % in-plane, +10–12 % κ_zz |
| **FC3 cutoff** | **not converged; not pursued** | MnIn₂Se₄ ≲ −10 %; Mn₂In₂Se₅ upper bound |
| electronic κ_e (experiment side) | not subtracted | unknown; needs the electrical conductivity |

## 6. Findings and corrections this round

1. **FC2 fix (MnIn₂Se₄).** The 4×4×1 FC2 gave an interlayer TA of 322.9 m/s.
   The 4×4×2 FC2 gives 1348.3 m/s (DFT elastic 1338.3; C₄₄ 9.543 vs 9.402 GPa).
   Mn₂In₂Se₅'s 4×4×1 FC2 was checked the same way and is fine (1386.8 m/s,
   C₄₄ 9.5–10.3 GPa).
2. **κ_zz q-mesh drift is a quadrature artifact, not physics.** The single
   Γ–Z point nearest Γ carries 50 % → 15 % of κ_zz from 11³ to 30³. Its
   contribution falls as N^-1.8. The rest of the zone is flat (0.6444 →
   0.6473). This withdraws the earlier claim that "intrinsic κ_zz is not
   well-defined (MFP > 1 μm)".
3. **The Mn₂In₂Se₅ rise under z-only refinement is the same artifact.**
   Isotropic 16³ → 24³ lowers κ_zz by 15.6 %.
4. **phono3py 4.4.0 SMM19/NJC23 coherence is wrong for MnIn₂Se₄.** It is
   12–17× too large and mesh-unstable. The cause is unregularised soft modes
   (C/ω diverges): cutting the modes below 0.3 THz restores agreement with
   3.24.0 and IBDB19.
   - **Rule:** coherence from 3.24.0 SMM19 or from 4.x IBDB19 only.
   - Mn₂In₂Se₅'s 4.4.0 value passes the soft-mode test (−0.3 %).
5. **FC3 cutoff ceiling.** The 6.45 / 6.51 Å limits are AutoZT's conservative
   parallelepiped rule. By the minimum-image rule, the existing 4×4×1 cells
   allow ≈ 8 Å (7.96 / 8.01 Å, checked pair by pair). 4×4×2 gains nothing
   under either rule.
6. **Withdrawn numbers:**
   - Mn₂In₂Se₅ κ_zz 0.951 ± 0.032 (an unjustified extrapolation) → 0.977 ± 0.008;
   - the "κ_zz grain suppression −27 %" (compared contaminated totals) → −10 % at 1 μm;
   - "Mn₂In₂Se₅ in-plane cutoff converged to 2 %" (came from ShengBTE 12×12×3) → −3.9 % / −4.3 % per step;
   - "the lower κ_zz, the larger the coherence fraction" (contradicted by the two materials);
   - the MnIn₂Se₄ L3 Para entry 0.35 in `calc_vs_experiment/README.md` → 0.49 (texture-mixed).
7. **ShengBTE results published before the fixes are invalid.** The FC3
   export dropped 41.4 % (MnIn₂Se₄) / 34 % (Mn₂In₂Se₅) of the triplets, and
   the MnIn₂Se₄ CONTROL file misplaced Mn.

## 7. Open items

| # | item | status |
|---|---|---|
| 1 | FC3 cutoff beyond 6.45 / 6.51 Å | **scheduled 2026-11-19** (revising the 11-18 decision): the κ_zz ordering of the two materials depends on it. Steps 0–1 of [`dft/fc3_cutoff_extend/`](dft/fc3_cutoff_extend/) refit the existing frames, no DFT; new frames only if the data/parameter ratio drops below 3 (`protocol/local_AI_next_steps_2026-11-19.md`, task 2) |
| 2 | literature source for Döll i = 0.44 | **closed 2026-11-19**: K.-J. Range, U. Klement, G. Döll, E. Bucher, J. R. Baumann, *Z. Naturforsch.* **46b**, 1122–1124 (1991), atomic-position table: M_oct = 0.56 Mn + 0.44 In, M_tet = 0.78 In + 0.22 Mn. PDF in `references/`. The 1990 *J. Cryst. Growth* paper is not the source |
| 3 | Mn₂In₂Se₅ antisite occupancy | **closed 2026-11-19 for the literature value**: x = 0.19(2) in the R-3m phase (Williams *et al.* 2024, ESI Table S7, powder Rietveld), inside the L3 acceptable range; L3 with x fixed gives L = 22 nm. The 0.704/0.296 of Chernoukhov *et al.* 2025 belongs to a different polymorph (P3̄m1). Optional: VCA rerun at x = 0.19 (estimated −9 % in-plane vs x = 0.13). A refinement of our own samples would still be better |
| 4 | κ_e of the measured samples | open. It needs the electrical conductivity (Wiedemann–Franz) |
| 5 | Tamura and phono3py references | identified but not checked against the papers (`analysis/disorder/README.md` §5) |

## 8. Where things are

| what | path |
|---|---|
| full MnIn₂Se₄ report (derivations, all tables) | `REPORT_MnIn2Se4_20261115.md` |
| paper figures + script | `test/figures/paper/` |
| figure index (which older figures are superseded) | `test/figures/README.md` |
| κ_zz decomposition script | `test/MnIn2Se4/kappa_zz_decompose.py` |
| raw per-mode phono3py output | `test/MnIn2Se4/kappa_zz_scan/`, `test/Mn2In2Se5/kappa_zz_decomposition/scan/`, `test/figures/kappa_vs_T/raw/` |
| disorder treatment, parameters, references | `analysis/disorder/README.md` |
| extrinsic model and fits | `analysis/extrinsic/`, `test/figures/calc_vs_experiment/` |
| coherence version check | `test/figures/kappa_vs_T/jobs30/wig.py` |
| cutoff-extension plan (not pursued) | `dft/fc3_cutoff_extend/` |
| experiment (digitised) | `experiment/kappa_tot_digitized.csv` |
