# Raw kappa(T) files (per-mode phono3py output behind the figures)

Every temperature point of every curve in `../kappa_vs_T.png` and
`../vca_vs_T.png` traces back to one of these files.

## Ordered crystal with natural isotopes

| file | material | temperatures (K) | job |
|---|---|---|---|
| `MnIn2Se4_T100-700.hdf5` | MnIn₂Se₄ | 100–700 step 100 | 416235 (14³, 6.45 Å, isotopes on) |
| `MnIn2Se4_T1-150.hdf5` | MnIn₂Se₄ | **1**, 30, 60, 90, 120, 150 | 416310 |
| `Mn2In2Se5_T100-700.hdf5` | Mn₂In₂Se₅ | 100–700 step 100 | 416236, re-run with isotopes on |
| `Mn2In2Se5_T30-150.hdf5` | Mn₂In₂Se₅ | 30, 60, 90, 120, 150 | 419202 (isotopes on) |

**1 K is present in the MnIn₂Se₄ file but excluded from the figures**: the Bose
occupation diverges there and κ_zz comes out at 1.9×10⁴ W m⁻¹ K⁻¹, which is
numerical, not physics. 30 K is the lowest usable point.

## VCA pair (with and without antisite disorder), same mesh

| file | material | contents | job |
|---|---|---|---|
| `MnIn2Se4_VCA_baseline_T30-700.hdf5` | MnIn₂Se₄ | ordered masses, natural isotopes | 424495 |
| `MnIn2Se4_VCA_disorder_T30-700.hdf5` | MnIn₂Se₄ | averaged masses + mass variance, x = 0.44 | 424495 |
| `Mn2In2Se5_VCA_baseline_T30-700.hdf5` | Mn₂In₂Se₅ | ordered masses, natural isotopes | 424496 |
| `Mn2In2Se5_VCA_disorder_T30-700.hdf5` | Mn₂In₂Se₅ | averaged masses + mass variance, x = 0.13 (archive; the curves now use x = 0.19, `../../../Mn2In2Se5/vca_m20_x019/`) | 424496 |

Each of these carries all **twelve** temperatures (30, 60, 90, 100, 120, 150, 200,
300, 400, 500, 600, 700 K), so the with/without-disorder pair is complete at every
point. They are on a 20³ mesh (both members of a pair share it, so only the
disorder differs); the plain `kappa_vs_T` curves above are on the converged
14³ / 16³ meshes.

Every temperature is decomposed with
`kappa_zz_decompose.py --temp-index <i>`, so κ_zz is the rest-of-the-zone term at
every T, not the Γ–Z-contaminated total.
