# Disorder treatment for MnIn₂Se₄ and Mn₂In₂Se₅

This note records **what disorder is included in which number**, the parameters
used, where they come from, and the literature basis. It is the document to review
before any of these κ values is quoted in a paper.

---

## 1. What kind of disorder these materials have

Both are layered chalcogenides with two cation sites:

| material | cation sites per formula unit | disorder |
|---|---|---|
| **MnIn₂Se₄** | 1 octahedral (Mn) + 2 tetrahedral (In) | Mn/In **antisite exchange** between the two sublattices |
| **Mn₂In₂Se₅** | 2 cation sites (2 Mn + 2 In in the 9-atom cell) | same family of antisite exchange; **no site-occupancy model has been fixed in this project yet** |

The antisite exchange keeps the stoichiometry but scrambles the masses: Mn
(54.94 amu) and In (114.82 amu) differ by a factor 2.1, so even a modest exchange
produces strong mass-disorder scattering.

---

## 2. Three levels of treatment, and what each is used for

| level | physics | how it is computed | where it appears |
|---|---|---|---|
| **L0 — none** | ideal ordered crystal | phono3py RTA, `is_isotope=False` | the original Mn₂In₂Se₅ scans (superseded) |
| **L1 — natural isotopes** | mass disorder from the natural isotopic mixture only | `is_isotope=True` (phono3py) | **all MnIn₂Se₄ numbers**; **all Mn₂In₂Se₅ numbers after the 2026-11-17 re-run** |
| **L2 — VCA antisite disorder** | masses averaged over the site occupancy, disorder returned as a mass variance | `mass_variances=[g_i]` with averaged masses (phono3py) | MnIn₂Se₄ at the largest cutoff (this round) |
| **L3 — extrinsic fit** | L1 + antisite mass disorder + grain boundary + texture, fitted to experiment | `analysis/extrinsic/kappa_extrinsic.py` (x, L, s scan) | `analysis/extrinsic/results/`, and the calc-vs-experiment figure |

L1 is the baseline for all reported κ. L2 and L3 are the two ways of adding the
antisite disorder; they are **not** added on top of each other in any number in
this repository.

**L2 and L3 do not treat the same disorder in the same way. Keep this in mind
before comparing them:**

| | L2 (VCA) | L3 (extrinsic fit) |
|---|---|---|
| site masses in the dynamical matrix | **averaged** over the occupancy (81.29 / 101.65 amu) | **ordered** masses kept (54.94 / 114.82 amu) |
| antisite scattering | mass variance g relative to the averaged mass | mass variance g relative to the site's mean mass, same formula |
| amount of antisite disorder | fixed: i = 0.44 (Range et al. 1991, single crystal) | x is a scan parameter; reference value 0.10 |
| grain boundary / texture | none | L and s fitted |

At the same antisite fraction the two still give different κ, because L2 also
changes the harmonic dispersion through the averaged masses.

**The antisite fraction must be the same in both.** With L3's definition
(x = fraction of all Mn on In sites), the L2 occupancy i = 0.44 *is* x = 0.44.
L3's reference value of x = 0.10 is therefore inconsistent with L2 for MnIn₂Se₄.
**Decision (2026-11-19):** i = 0.44 is adopted. It is the single-crystal
refinement of Range et al. 1991 (§5.1), so L3 for MnIn₂Se₄ uses x = 0.44 when x is
fixed from the literature.

---

## 3. Parameters actually used

### 3.1 Natural isotopes (L1)

phono3py's built-in natural-abundance masses and mass variances. No user input.
Effect at 300 K on Mn₂In₂Se₅: in-plane −1.16 %, κ_zz rest −0.39 % (measured
directly, 16³; see `test/Mn2In2Se5/iso_rerun/`).

### 3.2 VCA antisite disorder (L2), MnIn₂Se₄ only

Source: the user's `run_phono3py_kappa_vca.sh`, which states these occupancies are
**the single-crystal values of Range, Klement, Döll et al. 1991 (§5.1), i = 0.44**:

| site | occupancy | VCA mass (amu) | mass variance g |
|---|---|---:|---:|
| octahedral | 0.56 Mn + 0.44 In | **81.29** | **0.1337** |
| tetrahedral | 0.22 Mn + 0.78 In | **101.65** | **0.0595** |
| Se | natural | 78.96 | 4.63×10⁻⁴ |

The Se value 4.63×10⁻⁴ is the natural-abundance g from phonopy's isotope data,
the same one the L1 baseline uses through `is_isotope=True`. The first VCA run
used 0.0003 for Se. The `result.json` in §4.1 is still from that run; the
expected change is far below 0.1 %, but rerun before quoting.

Derivation (reproduced from the script, and verified here):
81.29 = 0.56×54.938 + 0.44×114.818;
g_oct = 0.56(1−54.938/81.29)² + 0.44(1−114.818/81.29)² = 0.1337;
101.65 = 0.22×54.938 + 0.78×114.818;
g_tet = 0.22(1−54.938/101.65)² + 0.78(1−114.818/101.65)² = 0.0595.

⚠️ **This parameterization is MnIn₂Se₄, not Mn₂In₂Se₅**: the occupancy sums
(0.56 + 2×0.22 = 1 Mn, 0.44 + 2×0.78 = 2 In) match the MnIn₂Se₄ formula, and the
mass-variance list has 7 entries = the 7-atom primitive cell
`Mn In In Se Se Se Se`. **No equivalent occupancy model exists in this repository
for Mn₂In₂Se₅**, so its numbers are L1 only.

### 3.3 Extrinsic fit (L3)

`analysis/extrinsic/`: three parameters

| parameter | meaning | what the κ data determine |
|---|---|---|
| x | fraction of all Mn on In sites (antisite) | **not determined in either material**: Mn₂In₂Se₅ accepts x = 0.01–0.50, MnIn₂Se₄ x = 0–0.66, provided L moves with it |
| L | effective phonon mean free path ("grain boundary") | **only together with x**: Mn₂In₂Se₅ 6–60 nm, MnIn₂Se₄ 6–43 nm |
| s | ⟨sin²θ⟩ texture, 0 = fully textured, 2/3 = random | **determined**: Mn₂In₂Se₅ 0.153–0.283, MnIn₂Se₄ 0.330–0.377 |

(x, L) form a degenerate band, because mass disorder and boundary scattering
lower κ almost interchangeably at high temperature. The grid optimum is one point
on that band and must not be quoted as "the fitted x" or "the grain size". x has
to come from a structural refinement (XRD / Rietveld occupancy). L_eff of
tens of nanometres is far below hot-pressed grain sizes, so it stands for
nanoscale scattering (stacking faults, intergrowth, spin disorder), not for grain
boundaries.

### 3.4 Latest fits (2026-11-17): recomputed on 24³ with the corrected force constants

Both materials were re-fitted from scratch: the earlier MnIn₂Se₄ fit used the
superseded 4×4×1 FC2, and the Mn₂In₂Se₅ fit used a z-only-refined mesh.
`kappa_extrinsic.py` now takes `--dim-fc2` and `--pmat` so that a mixed
FC2/FC3 supercell is possible (this was the missing capability).

| | **MnIn₂Se₄** | **Mn₂In₂Se₅** |
|---|---:|---:|
| force constants | FC2 4×4×2, FC3 6.45 Å | `Mn2In2Se5/fc/` (4×4×1), 6.51 Å |
| mesh | 24³ | 24³ (isotropic) |
| rms of ln(model/exp) | **1.2 %** | **3.2 %** |
| intrinsic (no disorder) rms | 87 % | 114 % |
| best x | **0.45** | 0.13 |
| best L | 30.7 nm | 15.9 nm |
| best s | 0.35 | 0.250 |
| fitted g_Mn / g_In | 0.1324 / 0.0609 | 0.1031 / 0.0354 |

**Two consistency checks fell out of this:**

1. The MnIn₂Se₄ fit lands on **x = 0.45**, i.e. on the Döll **i = 0.44**, and its
   fitted mass variances (g_Mn 0.1324, g_In 0.0609) agree with the VCA values
   (0.1337, 0.0595) to 1 % and 2.4 %. The κ data therefore support the VCA
   occupation model independently — that is the argument for using x = 0.44 for
   MnIn₂Se₄ in L3.
2. Re-fitting Mn₂In₂Se₅ on an isotropic 24³ mesh barely moves it
   (s 0.253 → 0.250, x 0.13, L 15.9 nm unchanged). The z-only mesh did not bias
   the high-temperature fit, because that fit is dominated by the extrinsic,
   mesh-insensitive part of κ.

⚠️ **The Mn₂In₂Se₅ fit has to be redone.** The bases it uses
(`results/Mn2In2Se5_m16164`, `_m16168`) are z-only meshes, whose κ_zz carries
the Γ–Z quadrature error. The fit must be repeated on bases from the isotropic
24³ (or 32³) mesh, where κ_zz is converged. At L ≈ 16 nm the boundary
suppresses the long-lifetime Γ–Z modes strongly, so the shift is expected to be
small, but this has not been checked.

MnIn₂Se₄ has a fit as well, but it rests on the **superseded 4×4×1 FC2** and is
therefore not used in any figure. It has to be redone with the 4×4×2 FC2 and an
x consistent with §2.

---

## 4. Results (300 K)

### 4.1 MnIn₂Se₄ — L1 vs L2, at the largest FC3 cutoff

20×20×20, FC3 = 6.45 Å, FC2 = 4×4×2, **RTA**, HanHai job 419636 (Se g = 0.0003 in
this run; see §3.2). The baseline rest term at 20³ (0.6380) is 1.3 % below the
mesh-converged L1 value of 0.648. That is within the tolerance, but the two
numbers are not interchangeable: the L1/L2 ratio below is a same-mesh
comparison.

| quantity | baseline (L1) | **VCA (L2, i = 0.44)** | suppression |
|---|---:|---:|---:|
| in-plane κ | 3.6884 | **1.4879** | **2.48×** |
| κ_zz total (Γ–Z contaminated) | 0.8780 | 0.6954 | 1.26× |
| **κ_zz rest of the zone** | **0.6380** | **0.4453** | **1.43×** |

**These are not comparable to the hot-pressed measurements directly**
(Perp 0.93, Para 0.485). The measured values mix the in-plane and
out-of-plane components through texture, and include boundary scattering. An
L2 value can only be compared with them after the same texture mixing (and,
if wanted, a boundary term). Putting 1.49 next to 0.93, or 0.4453 next to 0.485,
says nothing about whether i = 0.44 is right.

### 4.2 Mn₂In₂Se₅ — L2 (x = 0.13 on 2026-11-17, literature x = 0.19 since 2026-11-21)

Run the same way as MnIn₂Se₄: 20×20×20, FC3 = 6.51 Å, FC2 4×4×1.

**Current: x = 0.19** (Williams *et al.* 2024, §5.1). HanHai job 5
(`../../test/figures/kappa_vs_T/jobs_lbte_wigner/`, phono3py 4.4.0, 12 T) gives the
in-plane κ and the κ_zz total. Its per-mode hdf5 was not in the repository at first,
so the κ_zz rest term was taken as that total minus the first Γ–Z grid point
(q = (0, 0, 1/20), weight 2), recomputed alone with the same driver and the same
masses. That is exact in the RTA, and the single-point run reproduces the full
x = 0.13 file at that point to 10⁻¹⁵ (`../../test/Mn2In2Se5/vca_m20_x019/`). Once the
per-mode file is there, `postprocess.py` decomposes it in full and must agree to 2×10⁻⁵.

| site | occupancy (x = 0.19) | VCA mass (amu) | mass variance g |
|---|---|---:|---:|
| Mn | 0.81 Mn + 0.19 In | **66.315** | **0.12548** |
| In | 0.81 In + 0.19 Mn | **103.441** | **0.05157** |
| Se | natural | 78.96 | 4.63×10⁻⁴ |

| quantity (300 K) | baseline (L1) | **VCA (L2), x = 0.19** | suppression | VCA at x = 0.13 |
|---|---:|---:|---:|---:|
| in-plane κ | 2.9984 | **1.4611** | **2.05×** | 1.6335 (1.84×) |
| κ_zz total (raw) | 1.3730 | 1.2048 | 1.14× | 1.2246 |
| first Γ–Z point | 0.3894 | 0.4514 | | 0.4326 |
| **κ_zz rest of the zone** | **0.9837** | **0.7534** | **1.31×** | 0.7920 (1.24×) |

The first Γ–Z point (the quadrature artifact) rises with x, while the physical
rest term falls; only the rest term is used.

**Archive: x = 0.13** (HanHai job 424147, `../../test/Mn2In2Se5/vca_m20/`), from the
free L3 fit:

| site | occupancy (x = 0.13) | VCA mass (amu) | mass variance g |
|---|---|---:|---:|
| Mn | 0.87 Mn + 0.13 In | **62.73** | **0.1031** |
| In | 0.87 In + 0.13 Mn | **107.03** | **0.0354** |
| Se | natural | 78.96 | 4.63×10⁻⁴ |

| quantity | baseline (L1) | **VCA (L2)** | suppression |
|---|---:|---:|---:|
| in-plane κ | 2.9984 | **1.6335** | **1.84×** |
| κ_zz total (raw) | 1.3730 | 1.2246 | 1.12× |
| **κ_zz rest of the zone** | **0.9837** | **0.7920** | **1.24×** |

### 4.3 Both materials side by side

Same route for both (phono3py RTA, 20³, 300 K, largest FC3 cutoff, baseline with
natural isotopes). Machine-readable copy: `vca_summary.csv` / `vca_summary.json`.

| material | x | case | in-plane κ | κ_zz rest |
|---|---:|---|---:|---:|
| **MnIn₂Se₄** | 0.44 | no VCA | 3.6884 | 0.6380 |
| | | **VCA** | **1.4879** | **0.4453** |
| | | *suppression* | *2.48×* | *1.43×* |
| **Mn₂In₂Se₅** | 0.19 | no VCA | 2.9984 | 0.9837 |
| | | **VCA** | **1.4611** | **0.7534** |
| | | *suppression* | *2.05×* | *1.31×* |
| | 0.13 (archive) | VCA | 1.6335 | 0.7920 |

Two observations:

- The antisite disorder suppresses the **in-plane** component much more strongly
  than the cross-plane one in both materials (2.5× vs 1.4× for MnIn₂Se₄, 2.1× vs
  1.3× for Mn₂In₂Se₅), because the in-plane heat is carried by higher-frequency
  modes where the mass-variance rate (∝ g·ω²) is large.
- MnIn₂Se₄ has the stronger disorder response, because at i = 0.44 each site is
  close to a 50/50 mixture (g = 0.134 / 0.060), while at x = 0.19 the Mn₂In₂Se₅
  sites are still 81 % pure (g = 0.125 / 0.052).

### 4.4 Earlier L1-only statement (superseded by 4.2)

| quantity | value |
|---|---|
| κ_zz intrinsic (isotopes on) | **0.977 ± 0.008** (plateau 24³–32³ × 0.9961; the earlier 0.951 ± 0.032 was an unjustified extrapolation) |
| in-plane κ (RTA 32³, isotopes on) | 3.048 |
| in-plane κ (L1 + LBTE correction + coherence) | ≈ 3.37 |

Note the different basis: the MnIn₂Se₄ table in §4.1 is RTA only, while the
3.37 here includes the LBTE/iterative correction and the coherence term. Compare
the RTA rows only, or apply the same corrections to both.

---

## 5. Literature basis — status of each reference

| item | used for | status |
|---|---|---|
| phono3py `mass_variances` / `is_isotope` | L1 and L2 implementation | **verified**: the keyword and its behaviour were read in the installed phono3py 4.4.0 source (`api_phono3py.py`, `conductivity/`) |
| Tamura-type mass-variance scattering | the formula behind `mass_variances` | **reference identified, not checked against the paper**: S. Tamura, *Phys. Rev. B* **27**, 858 (1983). Check volume and page before submission |
| phono3py itself | every κ in this repository | **reference identified, not checked against the paper**: A. Togo, L. Chaput, I. Tanaka, *Phys. Rev. B* **91**, 094306 (2015); A. Togo, L. Chaput, T. Tadano, I. Tanaka, *J. Phys.: Condens. Matter* **35**, 353001 (2023). Check before submission |
| Simoncelli–Marzari–Mauri (SMM19), Nat. Phys. 15, 809 (2019) | the Wigner/coherence term (`transport_type="SMM19"`) | **verified**: this citation is written in the docstring of phono3py's `conductivity/smm19/__init__.py` |
| single-crystal occupancies at i = 0.44 | the VCA masses and g values; and, by §2, the x of L3 for MnIn₂Se₄ | **verified against the full text (2026-11-19)**: K.-J. Range, U. Klement, G. Döll, E. Bucher, J. R. Baumann, *Z. Naturforsch.* **46b**, 1122–1124 (1991), `references/` (§5.1). Occupancies are in its atomic-position table |

### 5.1 Occupancy literature (checked against the full texts in `../../references/`, 2026-11-19)

**MnIn₂Se₄ — the source of i = 0.44.** K.-J. Range, U. Klement, G. Döll,
E. Bucher, J. R. Baumann, "The Crystal Structure of MnIn₂Se₄, a Ternary Layered
Semiconductor", *Z. Naturforsch.* **46b**, 1122–1124 (1991).

- Single-crystal X-ray diffraction on CVT crystals (AlCl₃ transport). Occupation
  factors of the cations were refined. R-3m (centrosymmetric; refinement in R3m
  was tried and rejected), a = 4.051(1) Å, c = 39.464(2) Å, Z = 3, R = 0.064 for
  609 reflections.
- Atomic-position table: **M_oct = 0.56 Mn + 0.44 In** (3a, 0 0 0) and
  **M_tet = 0.78 In + 0.22 Mn** (6c, z = 0.23404). Both tetrahedral sites of a
  formula unit are one Wyckoff position, so they have the same occupancy. This is
  exactly the two-site model of §3.2.
- Text: "There is some cation disorder on octahedral (M_oct = 0.56 Mn + 0.44 In)
  and tetrahedral (M_tet = 0.78 In + 0.22 Mn) sites."
- Döll et al., *J. Cryst. Growth* **104**, 593 (1990) is the growth paper of the same
  group and is **not** the source of the occupancies. Do not cite it for i = 0.44.
- The crystals are CVT single crystals, like the ones characterised by the
  experimental side (a = 4.05, c = 39.35 Å). The hot-pressed pellets have a
  different thermal history, so their antisite fraction may differ.

**Mn₂In₂Se₅.** Our calculation is the 3R polytype, R-3m, a = 4.030 Å,
c = 48.72 Å (hexagonal setting of `Mn2In2Se5/fc/phono3py_disp.yaml`).

| source | phase | method | reported mixing | x for our model |
|---|---|---|---|---|
| A. J. Williams *et al.*, *J. Mater. Chem. C* **12**, 1753–1762 (2024), and its ESI | **3R, R-3m** (our phase), CVT crystals; a = 4.0189(2), c = 48.830(4) Å | **powder** Rietveld, Rwp = 0.0798 (the crystals were too pliable for a single-crystal solution) | ESI Table S7: Mn site (z = 0.6998) **0.81(2) Mn + 0.19 In**, In site (z = 0.55435) **0.81(2) In + 0.19 Mn**. The main text rounds this to "~21% site mixing"; their Raman supercell swaps 20 % | **0.19(2)** |
| I. V. Chernoukhov *et al.*, *Molecules* **30**, 1904 (2025) | **P-phase, P3̄m1**: a different polymorph, stabilised by quenching from 900 °C | powder XRD | room temperature: Mn1 = 0.704(2) Mn + 0.296 In, In1 = 0.704(2) In + 0.296 Mn; at 850 °C 0.5 / 0.5 | **not applicable** to the R-phase; shows that quenching from high T increases mixing |

So the literature value for our Mn₂In₂Se₅ is **x = 0.19(2)** (Williams 2024, ESI
Table S7). Because Mn and In are equal in number, the fraction of each site taken
by the other cation equals the fraction of all Mn on In sites, i.e. our x.

With x fixed, the existing 24³ L3 fit (`../extrinsic/results/*_m24/fit_*.json`,
`profile_x`, s re-optimised) still describes the 300–623 K data:

| | x | best L | rms of ln(model/exp) |
|---|---:|---:|---:|
| MnIn₂Se₄ | 0.44 (Range 1991) | 30.7 nm | 1.22 % |
| Mn₂In₂Se₅ | 0.19 (Williams 2024) | 22.1 nm | 3.4 % |
| Mn₂In₂Se₅ | 0.13 (fit optimum, for comparison) | 15.9 nm | 3.2 % |

So x is fixed from the literature and only L and s are fitted. The Mn₂In₂Se₅ VCA
(L2) curve now uses x = 0.19 as well (§4.2): 300 K in-plane 1.6335 → 1.4611
(−10.6 %), κ_zz rest 0.7920 → 0.7534 (−4.9 %). The estimate made before that run
from the 24³ L3 bases (ordered masses, Tamura only) was −9 % in-plane and −2 % in
the κ_zz total; the VCA total falls by 1.6 %.

Other literature, from abstracts only (not in `references/`):
- MnIn₂(1−z)Ga₂zSe₄, *Mater. Res. Bull.* (2001): Mn and In randomly distributed in MnIn₂Se₄.
- Mammadov *et al.*, *Condensed Matter and Interphases* **27**(1), 57–66 (2025),
  doi:10.17308/kcmf.2025.27/12487: MnIn₂Se₄ has a 5–6 mol % homogeneity range
  towards In₂Se₃, so off-stoichiometry and vacancies are possible.
- MnIn₂S₄ (sulfide): inversion depends on synthesis temperature
  (*Z. Anorg. Allg. Chem.* 2020, doi:10.1002/zaac.202000014).

---

## 6. Limitations to keep in mind

1. **VCA ignores force-constant disorder.** It assumes the averaged (ordered) FCs
   are valid and puts all the disorder into the masses. Bond-strength and strain
   disorder are therefore missing, so the VCA κ is **expected** to be an upper bound
   relative to a full disordered-supercell (SQS) treatment. This is not guaranteed:
   g = 0.13 is strong disorder, where the second-order (Born) mass-scattering rate is
   not exact, and VCA cannot produce resonant or local modes.
2. **The VCA occupancies were fixed, not fitted.** i = 0.44 is taken from the
   single-crystal refinement; the κ data were not used to adjust it.
3. **The two materials are not on the same disorder footing**: MnIn₂Se₄ has an L2
   number, Mn₂In₂Se₅ only L1. A like-for-like comparison of the two materials
   therefore relies on L1 for both, or needs an occupancy model for Mn₂In₂Se₅.
4. **L2 and L3 are alternatives, not additive.** L3 fits grain boundaries and
   texture to experiment; adding L2 on top would double-count the antisite
   scattering that L3 already contains through its x parameter.
