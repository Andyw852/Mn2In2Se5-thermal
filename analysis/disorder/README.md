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
| amount of antisite disorder | fixed: i = 0.44 (Döll, single crystal) | x is a scan parameter; reference value 0.10 |
| grain boundary / texture | none | L and s fitted |

At the same antisite fraction the two still give different κ, because L2 also
changes the harmonic dispersion through the averaged masses.

**The antisite fraction must be the same in both.** With L3's definition
(x = fraction of all Mn on In sites), the L2 occupancy i = 0.44 *is* x = 0.44.
L3's reference value of x = 0.10 is therefore inconsistent with L2 for MnIn₂Se₄.
**Open decision:** if Döll's i = 0.44 is adopted, L3 for MnIn₂Se₄ must use
x = 0.44. If it is not, L2 must not use i = 0.44 either.

---

## 3. Parameters actually used

### 3.1 Natural isotopes (L1)

phono3py's built-in natural-abundance masses and mass variances. No user input.
Effect at 300 K on Mn₂In₂Se₅: in-plane −1.16 %, κ_zz rest −0.39 % (measured
directly, 16³; see `test/Mn2In2Se5/iso_rerun/`).

### 3.2 VCA antisite disorder (L2), MnIn₂Se₄ only

Source: the user's `run_phono3py_kappa_vca.sh`, which states these occupancies are
**Döll's single-crystal values at i = 0.44**:

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

### 4.2 Mn₂In₂Se₅ — L2 added (2026-11-17)

Run the same way as MnIn₂Se₄: 20×20×20, 300 K, FC3 = 6.51 Å, FC2 4×4×1, HanHai job
424147. The VCA construction uses this material's own fitted x = 0.13:

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
| **Mn₂In₂Se₅** | 0.13 | no VCA | 2.9984 | 0.9837 |
| | | **VCA** | **1.6335** | **0.7920** |
| | | *suppression* | *1.84×* | *1.24×* |

Two observations:

- The antisite disorder suppresses the **in-plane** component much more strongly
  than the cross-plane one in both materials (2.5× vs 1.4× for MnIn₂Se₄, 1.8× vs
  1.2× for Mn₂In₂Se₅), because the in-plane heat is carried by higher-frequency
  modes where the mass-variance rate (∝ g·ω²) is large.
- MnIn₂Se₄ has the stronger disorder response despite the smaller x, because at
  x = 0.44 each site is close to a 50/50 mixture (maximum mass variance), while at
  x = 0.13 the sites stay nearly pure.

### 4.4 Earlier L1-only statement (superseded by 4.2)

| quantity | value |
|---|---|
| κ_zz intrinsic (isotopes on) | **0.951 ± 0.032** |
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
| Döll single-crystal occupancies at i = 0.44 | the VCA masses and g values; and, by §2, the x of L3 for MnIn₂Se₄ | **Bibliographic reference NOT verified here**: the number comes from the user's `run_phono3py_kappa_vca.sh`, which attributes it to a "Döll single-crystal occupancy" without a paper reference. **The exact citation must be supplied before this can go into a manuscript.** Note however that the independent L3 fit now lands on **x = 0.45**, with mass variances matching the VCA values to 1-2 % (section 3.4) -- the kappa data are consistent with i = 0.44 even though they cannot supply the reference themselves |

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
