# MnIn2Se4 kappa vs q-mesh (11x11x11 -> 14x14x14 -> 18x18x18)

Force constants: FC3 = 6.45 A (the 4x4x1 ceiling), FC2 = **4x4x2**, T = 300 K,
isotope scattering on. Written 2026-11-14, extended 2026-11-16 with the 18x18x18
RTA round.

## Results

| mesh | method | in-plane (W/m.K) | kappa_zz (W/m.K) | kappa_zz/in-plane | wall time |
|---|---|---:|---:|---:|---:|
| 11x11x11 | RTA | 3.5375 | 1.2495 | 0.353 | 107 s |
| 11x11x11 | LBTE | 4.4075 | 1.2519 | 0.284 | 939 s |
| **14x14x14** | RTA | 3.6474 | 1.0934 | 0.300 | 566 s |
| **14x14x14** | **LBTE** | **4.5219** | **1.0934** | **0.242** | 6829 s |
| **18x18x18** | **RTA** | **3.6945** | **0.9435** | **0.255** | - |

The 18x18x18 LBTE refinement did not complete and is not included (see
`../../REPORT_MnIn2Se4_20261115.md` section 10).

## Mesh convergence, step by step

| step | method | in-plane | kappa_zz |
|---|---|---:|---:|
| 11 -> 14 | RTA | +3.11 % | -12.50 % |
| 11 -> 14 | LBTE | +2.60 % | -12.66 % |
| 14 -> 18 | RTA | +1.29 % | -13.70 % |

## Verdict

- **In-plane kappa is converged.** The step size falls from +3.1 % to +1.3 % and
  the remaining error is well below the FC3-cutoff uncertainty (-7.9 % in the
  last available cutoff step, see `KAPPA_CUTOFF_SCAN.md`).
- **The intrinsic kappa_zz does not converge.** It falls by 13-14 % for every
  mesh refinement, with no sign of saturation at 18x18x18. This is not a
  numerical defect: the interlayer acoustic branches have mean free paths above
  one micrometre, so the "intrinsic" value of an infinite crystal is set by the
  sample size. A finite-grain value is the physically meaningful quantity; see
  section 7 of `README.md` and section 6.4 of the report
  (1 um grain: 3.5729 / 0.6877).

## A noteworthy observation

At 14x14x14 the **kappa_zz from RTA and LBTE are identical** (1.09334 vs
1.09344, a difference of 8x10^-5), whereas the in-plane value differs by 24 %
(3.6474 -> 4.5219).

Physical explanation: the interlayer coupling is very weak, so the LBTE
iterative correction to kappa_zz (which relies on off-diagonal phonon
transitions) is almost zero along z, hence kappa_zz(LBTE) ~ kappa_zz(RTA); the
in-plane iterative correction is significant. This is consistent with the
picture in which kappa_zz is dominated by weak interlayer coupling.

## Most reliable kappa (mesh 18, RTA, 300 K)

- in-plane **3.6945** W/m.K
- kappa_zz, 1 um grain **0.6877** W/m.K (intrinsic limit 0.9435, not converged)
- anisotropy kappa_zz/in-plane **0.192** at a 1 um grain

Caveats that must be carried along:
- the **FC3 cutoff is not converged** (6.45 A is already the 4x4x1 ceiling, see
  `KAPPA_CUTOFF_SCAN.md`);
- the **intrinsic kappa_zz does not converge** with mesh (13-14 % per step);
- the **coherence term** is not in these numbers (Wigner-RTA adds +2.5 % / +6.0 %
  at mesh 14, see section 8 of `README.md`).

## Reproduction

```bash
python kappa_driver_442.py --mesh 11 11 11 --rta      # about 2 minutes
python kappa_driver_442.py --mesh 11 11 11 --lbte     # about 16 minutes
python kappa_driver_442.py --mesh 14 14 14 --rta      # about 9.5 minutes
python kappa_driver_442.py --mesh 14 14 14 --lbte     # about 114 minutes
python kappa_driver_442.py --mesh 18 18 18 --rta      # about 35 minutes
```
