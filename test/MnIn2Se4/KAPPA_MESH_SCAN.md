# MnIn2Se4 κ mesh convergence (14×14×14)

**2026-11-14** | fc3 = 6.45 Å (the 4×4×1 limit), fc2 = **4×4×2**, 300 K, isotopes on.
This adds one round of **14×14×14** on top of the existing 11×11×11.

## Results

| mesh | method | in-plane (W/m·K) | κzz (W/m·K) | κzz/in-plane | wall time |
|---|---|---:|---:|---:|---:|
| 11×11×11 | RTA | 3.5375 | 1.2495 | 0.353 | 107 s |
| 11×11×11 | LBTE | 4.4075 | 1.2519 | 0.284 | 939 s |
| **14×14×14** | RTA | 3.6474 | 1.0934 | 0.300 | 566 s |
| **14×14×14** | **LBTE** | **4.5219** | **1.0934** | **0.242** | 6829 s |

## Mesh convergence (11 → 14)

| method | in-plane | κzz |
|---|---:|---:|
| RTA | **+3.11%** | **−12.50%** |
| LBTE | **+2.60%** | **−12.66%** |

**Verdict**:
- **The in-plane value is converged**: the two meshes differ by only 2.6–3.1%, and further refinement is expected to change it by <2%.
- **κzz is still not converged**: it drops by 12.7% from 11→14, and 18×18×18 or denser is needed to confirm it.

## A noteworthy observation

At 14×14×14 the **κzz from RTA and LBTE are exactly the same** (1.09334 vs 1.09344, a difference of 8×10⁻⁵),
whereas the in-plane value differs by 24% (3.6474 → 4.5219).

Physical explanation: the interlayer coupling is very weak, so the LBTE iterative correction to κzz (which relies on off-diagonal phonon transitions)
is almost zero along z, hence κzz_LBTE ≈ κzz_RTA; the in-plane iterative correction is significant.
This is consistent with the physical picture in which κzz is dominated by weak interlayer coupling.

## Currently most reliable κ (mesh 14, LBTE, 300 K)

- in-plane **4.5219** W/m·K
- κzz **1.0934** W/m·K
- anisotropy κzz/in-plane **0.242**

⚠️ It must still be noted that the **fc3 cutoff is not converged** (see `KAPPA_CUTOFF_SCAN.md`; 6.45 Å is already the 4×4×1 limit)
and that **the κzz mesh is not converged** (it still changes by −12.7% from 11→14).

## Reproduction

```bash
python kappa_driver_442.py --mesh 14 14 14 --lbte    # about 114 minutes
python kappa_driver_442.py --mesh 14 14 14 --rta     # about 9.5 minutes
```
