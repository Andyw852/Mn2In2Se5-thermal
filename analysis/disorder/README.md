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
| Se | natural | 78.96 | 0.0003 |

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

| parameter | meaning | fitted value |
|---|---|---|
| x | fraction of Mn on In sites (antisite) | Mn₂In₂Se₅ 0.13 (acceptable 0.01–0.5) |
| L | effective grain-boundary mean free path | Mn₂In₂Se₅ 15.9 nm |
| s | ⟨sin²θ⟩ texture, 0 = fully textured, 2/3 = random | Mn₂In₂Se₅ 0.253 |

MnIn₂Se₄ has a fit as well, but it rests on the **superseded 4×4×1 FC2** and is
therefore not used in any figure.

---

## 4. Results (300 K)

### 4.1 MnIn₂Se₄ — L1 vs L2, at the largest FC3 cutoff

20×20×20, FC3 = 6.45 Å, FC2 = 4×4×2, HanHai job 419636:

| quantity | baseline (L1) | **VCA (L2, i = 0.44)** | suppression |
|---|---:|---:|---:|
| in-plane κ | 3.6884 | **1.4891** | **2.48×** |
| κ_zz total (Γ–Z contaminated) | 0.8780 | 0.6954 | 1.26× |
| **κ_zz rest of the zone** | **0.6380** | **0.4455** | **1.43×** |

Measured hot-pressed values for comparison: Perp (in-plane) 0.93, Para
(out-of-plane) 0.485.

### 4.2 Mn₂In₂Se₅ — L1 only

| quantity | value |
|---|---|
| κ_zz intrinsic (isotopes on) | **0.951 ± 0.032** |
| in-plane κ (RTA 32³, isotopes on) | 3.048 |
| in-plane κ (L1 + LBTE correction + coherence) | ≈ 3.37 |

---

## 5. Literature basis — status of each reference

| item | used for | status |
|---|---|---|
| phono3py `mass_variances` / `is_isotope` | L1 and L2 implementation | **verified**: the keyword and its behaviour were read in the installed phono3py 4.4.0 source (`api_phono3py.py`, `conductivity/`) |
| Tamura-type mass-variance scattering | the formula behind `mass_variances` | **standard formula, but the exact bibliographic reference was NOT verified in this session** — phono3py's installed source does not carry the citation text. To be filled in from the phono3py paper/documentation before publication |
| Simoncelli–Marzari–Mauri (SMM19), Nat. Phys. 15, 809 (2019) | the Wigner/coherence term (`transport_type="SMM19"`) | **verified**: this citation is written in the docstring of phono3py's `conductivity/smm19/__init__.py` |
| Döll single-crystal occupancies at i = 0.44 | the VCA masses and g values | **NOT verified here**: the number comes from the user's `run_phono3py_kappa_vca.sh`, which attributes it to a "Döll single-crystal occupancy" without a bibliographic entry. **The exact paper must be supplied before this can go into a manuscript** |

---

## 6. Limitations to keep in mind

1. **VCA ignores force-constant disorder.** It assumes the averaged (ordered) FCs
   are valid and puts all the disorder into the masses. Bond-strength and strain
   disorder are therefore missing, which makes the VCA κ an **upper bound** relative
   to a full disordered-supercell (SQS) treatment.
2. **The VCA occupancies were fixed, not fitted.** i = 0.44 is taken from the
   single-crystal refinement; the κ data were not used to adjust it.
3. **The two materials are not on the same disorder footing**: MnIn₂Se₄ has an L2
   number, Mn₂In₂Se₅ only L1. A like-for-like comparison of the two materials
   therefore relies on L1 for both, or needs an occupancy model for Mn₂In₂Se₅.
4. **L2 and L3 are alternatives, not additive.** L3 fits grain boundaries and
   texture to experiment; adding L2 on top would double-count the antisite
   scattering that L3 already contains through its x parameter.
