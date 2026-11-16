# Mn₂In₂Se₅ — re-run with isotope scattering ON (the project standard)

## Why this directory exists

`Mn2In2Se5/config/step6_kappa.step.conf` sets **`ISOTOPE = true`**, and every
MnIn₂Se₄ run in this study used `is_isotope=True`. All earlier Mn₂In₂Se₅ runs
were made with `is_isotope=False`, which put the two materials on different
footings. Everything was therefore re-run with isotopes on, same force constants
and `primitive_matrix="P"`; only the isotope flag changed, so old and new pair up
one to one. HanHai jobs 419197–419202, 419208.

## The isotope correction is mesh-independent, so only 16³ and 24³ were re-run

The 28³ and 32³ re-runs were started and then cancelled: the isotope effect is a
nearly uniform multiplicative factor that barely moves with the mesh, so the
ratio measured at 16³/24³ can be applied to the already-converged isotope-off
values instead of spending another 4.5 h of cluster time.

| quantity (rest of the zone) | ratio on/off at 16³ | at 24³ |
|---|---:|---:|
| intrinsic κ_zz | 0.99604 | 0.99621 |
| 10 μm grain | 0.99599 | 0.99604 |
| 1 μm grain | 0.99551 | 0.99549 |

The two meshes agree to **0.02 %**, which is what justifies the transfer.

| quantity (300 K) | isotopes off (converged) | **isotopes on (corrected)** |
|---|---:|---:|
| **κ_zz intrinsic** | 0.9547 ± 0.0323 | **0.951 ± 0.032** |
| κ_zz, 10 μm grain | 0.9281 | **0.924** |
| κ_zz, 1 μm grain | 0.7879 | **0.784** |
| in-plane κ (RTA, 32³) | 3.0825 | **3.048** |

In-plane: the ratio is 0.98841 at 16³ and 0.98899 at 24³ (−1.1 %), also
mesh-independent. With the LBTE correction (×1.059, measured at 10³) and the
Wigner coherence term (0.1466 at 16³, isotopes on), the in-plane value becomes
**3.048 × 1.059 + 0.1466 = 3.37**.

## Isotope effect at 16³, 300 K (directly measured)

| quantity | isotopes off | isotopes **on** | change |
|---|---:|---:|---:|
| in-plane κ | 2.9946 | **2.9599** | −1.16 % |
| κ_zz total | 1.5472 | **1.5432** | −0.26 % |
| κ_zz rest of the zone | 1.0090 | **1.0050** | −0.39 % |

Small at 300 K (Mn is monoisotopic; only In and Se contribute, mass variance
≈2×10⁻³) but it grows at low temperature: the in-plane κ at 30 K drops from
63.58 to 52.88 (**−17 %**).

## Files

`cutoff_scan_iso/` — the three FC3 cutoffs at 16³, decomposed:

| r_c | in-plane | κ_zz total | κ_zz rest |
|---|---:|---:|---:|
| 5.68 Å | 3.1440 | 2.0423 | 1.2408 |
| 6.12 Å | 3.0215 | 1.6687 | 1.0934 |
| 6.51 Å | 2.8925 | 1.5187 | 0.9850 |

`wigner_cutoff_iso/` — Wigner-RTA at the same three cutoffs:

| r_c | in-plane P | in-plane C | in-plane TOT | κ_zz P | κ_zz C | κ_zz TOT |
|---|---:|---:|---:|---:|---:|---:|
| 5.68 Å | 3.1440 | 0.1558 | 3.2998 | 2.0423 | 0.1141 | 2.1563 |
| 6.12 Å | 3.0215 | 0.1492 | 3.1708 | 1.6687 | 0.1149 | 1.7836 |
| 6.51 Å | 2.8925 | 0.1466 | 3.0391 | 1.5187 | 0.1144 | 1.6331 |

The coherence term is again essentially cutoff-independent (in-plane −5.9 % over
the range, κ_zz +0.3 %).

## Superseded

Every isotope-off Mn₂In₂Se₅ number earlier in this repository. They are kept for
the record (`../kappa_zz_decomposition/scan/`, `cutoff_scan/`,
`wigner_cutoff/`, `lowtemp/`, `multitemp/`) but must not be quoted next to the
MnIn₂Se₄ values.
