# MnIn2Se4 kappa_zz convergence on the q-mesh

Force constants: FC3 = 6.45 A (the 4x4x1 ceiling), FC2 = **4x4x2**, T = 300 K,
isotope scattering on, phono3py RTA. This file supersedes the earlier version of
the mesh scan, which reported only the total kappa_zz and therefore interpreted
a q-point quadrature artifact as physics.

## 1. Totals

| mesh | in-plane | kappa_zz | wall time |
|---|---:|---:|---:|
| 11x11x11 | 3.5375 | 1.2495 | 107 s |
| 14x14x14 | 3.6474 | 1.0934 | 566 s |
| 18x18x18 | 3.6945 | 0.9435 | 2674 s |
| 24x24x24 | 3.7020 | 0.8052 | 548 s (HanHai, 64 threads) |
| **30x30x30** | **3.6946** | **0.7645** | 2017 s (HanHai, 64 threads) |

In-plane: +3.11 % then +1.28 % per refinement -> **converged**.
kappa_zz: -12.50 % then -13.70 % -> **not converged, and it falls**.

## 2. Why kappa_zz falls: the first Gamma-Z q-point

The drift is localised. Computing the contribution of every grid point from
`mode_kappa` in the raw phono3py output (`kappa_zz_decompose.py`) shows that the
single q-point closest to Gamma on the Gamma-Z line -- the layer normal,
q = (1/N, 1/N, 1/N) in the primitive reciprocal basis -- carries half of the
total on a coarse mesh:

| mesh | kappa_zz total | first Gamma-Z point (+-q) | rest of the zone | first / total |
|---|---:|---:|---:|---:|
| 11^3 | 1.2495 | 0.6244 | 0.6251 | 50.0 % |
| 14^3 | 1.0934 | 0.4490 | 0.6444 | 41.1 % |
| 18^3 | 0.9435 | 0.2974 | 0.6461 | 31.5 % |
| 24^3 | 0.8052 | 0.1594 | 0.6457 | 19.8 % |
| **30^3** | **0.7645** | **0.1172** | **0.6473** | **15.3 %** |

Mechanism:

- the weight of that point in the Brillouin-zone average is 2/N^3;
- its lifetime grows as Gamma is approached: single-point runs up to 60^3 give
  tau ~ omega^-1.3 to omega^-1.9 for the interlayer acoustic branch (TA:
  650 ps at 11^3 -> 10 800 ps at 60^3);
- so its contribution falls as **N^-1.8** (fit 14^3-30^3; step exponents 1.64,
  2.17, 1.38) -- not 1/N. It is a **quadrature error of the q-point integration
  near Gamma along c***. Any exponent of tau below 3 keeps the 3D integral
  finite, so the intrinsic kappa_zz is **finite**.

The **rest of the zone** column is converged from 14^3 on and is the quantity
to quote.

**Mn2In2Se5 is not ruled out.** There kappa_zz *rises* under refinement along c
only (16x16x4 -> 16x16x16: 0.84 -> 1.11 -> 1.30 -> 1.54). With only N_z refined,
the first Gamma-Z point's weight falls as 1/N_z while its lifetime grows as
N_z^p, so its contribution grows as N_z^(p-1) -- a *rising* total for
p ~ 1.3-1.9. The same quadrature error therefore predicts a falling trend under
N^3 refinement and a rising one under 16x16xN_z refinement. Mn2In2Se5 has to be
decomposed (`kappa_zz_decompose.py --poscar ... --pmat ...`) before any of its
kappa_zz values are quoted.

## 3. Finite-grain values from the same output files

RTA kappa is proportional to 1/gamma mode by mode, so a Casimir boundary is
applied a posteriori (`kappa_zz_decompose.py --mfp`) and needs no extra runs.
Validation: at 11^3 the reconstruction gives kappa_zz = 0.8976 at 1 um and the
explicit `boundary_mfp=1.0` run (`kappa_zz_scan/m11_b1um`) gives 0.8976 too.

A boundary shortens the first Gamma-Z point's lifetime but does not remove the
point from the sum, so **the finite-grain values must use the rest term too**.

| mesh | 10 um, rest | 1 um, rest | 10 um, total | 1 um, total |
|---|---:|---:|---:|---:|
| 14^3 | 0.6385 | 0.5957 | 1.0282 | 0.7758 |
| 18^3 | 0.6392 | 0.5911 | 0.8849 | 0.6877 |
| 24^3 | 0.6376 | 0.5846 | 0.7640 | 0.6287 |
| **30^3** | **0.6381** | **0.5823** | 0.7233 | 0.6070 |
| **mesh-converged** | **0.638 +- 0.002** | **0.58 +- 0.01** | -- | -- |

The totals (right) are contaminated: the 30^3 10 um total is above the intrinsic
kappa_zz, which a boundary cannot do. On the rest term the grain effect is
-1.5 % at 10 um and -10 % at 1 um. The 1 um value is the least certain: its rest
term still drifts down by 0.013 from 14^3 to 30^3, while the total falls towards
it from above.

## 4. Intrinsic kappa_zz: estimate

The "rest of the zone" is converged: **0.6444 -> 0.6461 -> 0.6457 -> 0.6473**
for 14^3 -> 30^3 (full spread 0.0029, non-monotonic), while the total still moves
by -5.1 % at 30^3. Dropping the first point also drops the true contribution
c(N) of its cell. With c ~ N^-1.8, the +0.0029 drift from 14^3 to 30^3 gives
c(30^3) ~ +0.001. Result: **0.648 +- 0.003** (the larger of the bias and the
spread).

Independent check: split the rest term at a fixed radius q_c around Gamma. At
q_c = 0.02 1/A the outside part is 0.4486 / 0.4424 / 0.4425 and the inside part
(without the first point) is 0.1974 / 0.2033 / 0.2048 at 18^3 / 24^3 / 30^3.
Both are mesh-independent, so only the overweighted single point is removed, not
real small-q heat current.

| quantity | value (W/m.K) |
|---|---:|
| 18^3 RTA total (**firm upper bound**) | **0.9435** |
| **intrinsic kappa_zz (mesh-converged)** | **0.648 +- 0.003** |
| 10 um grain (rest term, mesh-converged) | 0.638 +- 0.002 |
| 1 um grain (rest term, mesh-converged) | 0.58 +- 0.01 |
| FC3 cutoff (section 5 / report 6.6) | **not converged**: last step -8 %, remaining change unbounded (not included) |

## 5. The cutoff comparison has to be redone

The committed FC3 cutoff scan (`KAPPA_CUTOFF_SCAN.md`) was run on the
**11x11x11** mesh, where the first Gamma-Z point alone is 50 % of kappa_zz. Its
kappa_zz column ("-55 % then +2 %") is therefore the quadrature artifact and says
nothing about the cutoff. The in-plane column is unaffected and remains valid.
The cutoff comparison for kappa_zz is repeated on the converged "rest" term at
14^3 (`kappa_zz_scan/c560_m14`, `c623_m14`) -- see report section 6.6.

## 6. Reproduction

```bash
bash kappa_zz_scan.sh          # the whole scan, CPU only, resumable, ~4 h
python kappa_zz_decompose.py kappa_zz_scan/*/kappa-m*.hdf5 \
       --json kappa_zz_gamma_line.json --mfp 1 10
```

The raw per-mode output (`kappa_zz_scan/m{11,14,18,24,30}/kappa-m*.hdf5`) is
committed (and `kappa_zz_decompose.py` prints every table above, including the
estimate and the fixed-sphere check), so every number above can be regenerated without repeating a
calculation. `kappa_zz_scan.sh` skips any case whose hdf5 already exists.
