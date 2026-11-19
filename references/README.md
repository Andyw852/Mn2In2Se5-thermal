# Reference papers

Publisher PDFs of the papers this project relies on. **Read-only material**, not part
of the calculation. See the licence note at the bottom before making the repository
public.

| file | citation | what it is used for here |
|---|---|---|
| `Range_Doll_1991_MnIn2Se4_crystal_structure_ZNatB46_1122.pdf` | K.-J. Range, U. Klement, **G. Döll**, E. Bucher, J. R. Baumann, "The Crystal Structure of MnIn₂Se₄, a Ternary Layered Semiconductor", *Z. Naturforsch.* **46b**, 1122–1124 (1991) | **Source of the cation occupancies used by the VCA** (see below) |
| `Williams_2024_Mn2In2Se5_Mn2Ga2S5_single_crystals_JMaterChemC12_1753.pdf` | A. J. Williams, A. Reifsnyder, B. Yu, C. E. Moore, M. A. Susner, W. Windl, D. W. McComb, J. E. Goldberger, "Single crystal synthesis and properties of the two-dimensional van der Waals frustrated magnets, Mn₂In₂Se₅ and Mn₂Ga₂S₅", *J. Mater. Chem. C* **12**, 1753–1762 (2024), DOI 10.1039/d3tc03550f | Structure of 3R Mn₂In₂Se₅ (R-3m, the phase we calculate): double-octahedral MnSe₆ layers capped by InSe₄ tetrahedra. Powder Rietveld gives **~21 % Mn/In site mixing** → x ≈ 0.21 for the L3 model (occupancies in their ESI Tables S6–S7, not included here) |
| `Chernoukhov_2025_Mn2In2Se5_high_T_layered_Molecules30_1904.pdf` | I. V. Chernoukhov, A. D. Pyreu, A. N. Azarevich, A. N. Samarin, A. V. Bogach, K. O. Znamenkov, A. V. Shevelkov, V. Yu. Verchenko, "High-Temperature Layered Modification of Mn₂In₂Se₅", *Molecules* **30**, 1904 (2025), DOI 10.3390/molecules30091904 | A **different** polymorph, the P-phase (P3̄m1), stabilised by quenching from 900 °C, with Mn1 = 0.704(2) Mn + 0.296 In at room temperature. It is **not** the phase we model (ours is the R-phase, R-3m, c ≈ 48.7 Å); cited only to show that mixing grows when quenched from high temperature |

## The MnIn₂Se₄ occupancy that the VCA uses

Range *et al.* 1991, atomic-position table and discussion:

- octahedral site `M_oct = 0.56 Mn + 0.44 In`
- tetrahedral site `M_tet = 0.78 In + 0.22 Mn`
- `M_oct–Se = 6 × 2.721(1) Å`, `M_tet–Se = 1 × 2.527(2)` and `3 × 2.593(1) Å`
- space group `R3̄m`, a = 4.051(1) Å, c = 39.464(2) Å, Z = 3

This is the "i = 0.44" used in `test/MnIn2Se4/vca_m20/` and in `analysis/disorder/README.md`
(§3.2): with i = 0.44 the averaged site masses come out as 81.29 amu (octahedral) and
101.65 amu (tetrahedral), and the mass variances as 0.1337 and 0.0595.

The paper also states explicitly that *"there is some cation disorder on octahedral and
tetrahedral sites"*, and that the space-group symmetry forces equal occupancy of the
tetrahedral sites — which is what makes a two-site VCA description legitimate.

## Licence / redistribution

| paper | licence |
|---|---|
| *Z. Naturforsch.* 1991 (De Gruyter) | publisher copyright, **not** open access |
| *J. Mater. Chem. C* 2024 (RSC) | publisher copyright, **not** open access unless a CC-BY option was bought |
| *Molecules* 2025 (MDPI) | **CC BY 4.0**, redistribution allowed with attribution |

Only the MDPI paper may be redistributed freely. If this repository is public, these
PDFs are for private reading only; either keep the repository private or replace the
first two entries with their DOIs.
