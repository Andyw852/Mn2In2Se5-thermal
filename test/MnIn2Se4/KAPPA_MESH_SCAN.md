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
| **24x24x24** | **3.7020** | **0.8052** | 548 s (HanHai, 64 threads) |

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
| 14^3 | 1.0934 | 0.4490 | 0.6444 | 58.9 % |
| 18^3 | 0.9435 | 0.2974 | 0.6461 | 31.5 % |
| **24^3** | **0.8052** | **0.1594** | **0.6457** | **19.8 %** |

Mechanism:

- the weight of that point in the Brillouin-zone average is 2/N^3;
- its mode kappa grows as ~N^2, because the three-phonon scattering rate of the
  interlayer acoustic branch falls as omega^2 near Gamma (tau ~ 1/omega^2);
- so its contribution scales as **1/N** -- it is a **quadrature error of the
  q-point integration near Gamma along c***, not a divergent integral.

The **rest of the zone** column is already converged at ~0.62-0.65 and is the
quantity to quote. The intrinsic kappa_zz is therefore **finite**.

Note the direction of the trend. Here kappa_zz falls with refinement; in
Mn2In2Se5 the same quantity *rises* (16x16x4 -> 16x16x16: 0.84 -> 1.11 -> 1.30
-> 1.54). A rising trend is the signature of unresolved long-wavelength modes;
a falling one is the signature of an overweighted small-q point. The two cannot
have the same cause, so the earlier statement that this is "the same behaviour as
Mn2In2Se5" was wrong and has been removed.

## 3. Finite-grain values from the same output files

RTA kappa is proportional to 1/gamma mode by mode, so a Casimir boundary is
applied a posteriori (`kappa_zz_decompose.py --mfp`) and needs no extra runs.
Validation: at 11^3 the reconstruction gives kappa_zz = 0.8976 at 1 um and the
explicit `boundary_mfp=1.0` run (`kappa_zz_scan/m11_b1um`) gives 0.8976 too.

| mesh | 10 um | 1 um |
|---|---:|---:|
| 11^3 | 1.1878 | 0.8976 |
| 14^3 | 1.0282 | 0.7758 |
| 18^3 | 0.8849 | 0.6877 |
| **24^3** | **0.7640** | **0.6287** |

At 11^3, between infinity and 1 um, the rest term moves 0.6251 -> 0.5857 (-6 %)
while the first Gamma-Z point falls 0.6244 -> 0.3119 (-50 %). The boundary
therefore removes the artifact, not the physical heat current.

## 4. Intrinsic kappa_zz: estimate

The "rest of the zone" is converged and the last step proves it:
**0.6251 -> 0.6444 -> 0.6461 -> 0.6457** for 11^3 -> 14^3 -> 18^3 -> 24^3, i.e.
+3.1 % over the whole range and **-0.06 %** in the last step, while the total
still moves by -14.7 %. A naive linear fit in 1/N gives A = 0.667 with residuals
up to +-0.007 and is dominated by the coarse 11^3 point; the measured sequence
settles at **0.646**.

| quantity | value (W/m.K) |
|---|---:|
| 18^3 RTA total (**firm upper bound**) | **0.9435** |
| **intrinsic kappa_zz (mesh-converged)** | **0.646** |
| statistical (mesh) spread | +-0.003 |
| FC3-cutoff systematic (section 5 / report 6.6) | about -8 % (not included) |
| 1 um grain, 18^3 RTA | 0.6877 |

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

The raw per-mode output (`kappa_zz_scan/m{11,14,18,24}/kappa-m*.hdf5`) is
committed, so every number above can be regenerated without repeating a
calculation. `kappa_zz_scan.sh` skips any case whose hdf5 already exists.
