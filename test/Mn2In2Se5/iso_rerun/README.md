# Mn₂In₂Se₅ — re-run with isotope scattering ON (the project standard)

## Why this directory exists

`Mn2In2Se5/config/step6_kappa.step.conf` sets **`ISOTOPE = true`**, and every
MnIn₂Se₄ run in this study used `is_isotope=True`. All earlier Mn₂In₂Se₅ runs
were made with `is_isotope=False`, which put the two materials on different
footings. Everything was therefore re-run with isotopes on, same force constants
and `primitive_matrix="P"`; only the isotope flag changed, so old and new pair up
one to one.

HanHai jobs 419197–419202 and 419208.

## Isotope effect at 16³, 300 K

| quantity | isotopes off | isotopes **on** | change |
|---|---:|---:|---:|
| in-plane κ | 2.9946 | **2.9599** | −1.16 % |
| κ_zz total | 1.5472 | **1.5432** | −0.26 % |
| κ_zz rest of the zone | 1.0090 | **1.0050** | −0.39 % |

Small at 300 K (Mn is monoisotopic; only In and Se contribute, and their mass
variance is ≈2×10⁻³), but it grows at low temperature: the in-plane κ at 30 K
drops from 63.58 to 52.88 (**−17 %**), because isotope scattering takes over as
the anharmonic rate falls.

## Files

`cutoff_scan_iso/` — the three FC3 cutoffs at 16³, decomposed
(`c568/c612/c651_m161616_kappa.hdf5` + the decompose JSON):

| r_c | in-plane | κ_zz total | κ_zz rest |
|---|---:|---:|---:|
| 5.68 Å | 3.1440 | 2.0423 | 1.2408 |
| 6.12 Å | 3.0215 | 1.6687 | 1.0934 |
| 6.51 Å | 2.8925 | 1.5187 | 0.9850 |

`../wigner_cutoff_iso/` — Wigner-RTA at the same three cutoffs, isotopes on:

| r_c | in-plane P | in-plane C | in-plane TOT | κ_zz P | κ_zz C | κ_zz TOT |
|---|---:|---:|---:|---:|---:|---:|
| 5.68 Å | 3.1440 | 0.1558 | 3.2998 | 2.0423 | 0.1141 | 2.1563 |
| 6.12 Å | 3.0215 | 0.1492 | 3.1708 | 1.6687 | 0.1149 | 1.7836 |
| 6.51 Å | 2.8925 | 0.1466 | 3.0391 | 1.5187 | 0.1144 | 1.6331 |

The coherence term is again essentially cutoff-independent (in-plane −5.9 % over
the range, κ_zz +0.3 %).

## Superseded

Every isotope-off Mn₂In₂Se₅ number earlier in this repository. They are kept for
the record (`scan/`, `cutoff_scan/`, `wigner_cutoff/`, `lowtemp/`,
`multitemp/`) but must not be quoted next to the MnIn₂Se₄ values.
