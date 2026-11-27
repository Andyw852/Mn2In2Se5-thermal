# Spin-disordered FC2: κ result and diagnosis of the G2 failure (2026-11-26)

Runs: `analysis/extrinsic/kappa_extrinsic.py bases` with each committed FC2 and the
production fc3 (`Mn2In2Se5/fc/fc3.hdf5`), 16×16×4, then `kappa_summary.py`
(output: `kappa_summary.json`). phono3py 4.4.0, about 3.5 min per set on 4 cores. The
`bases.npz` files are not committed; rerun the command to regenerate them.

## κ at 300 K (intrinsic RTA, W m⁻¹ K⁻¹)

| FC2 | κ_xx | κ_zz |
|---|---:|---:|
| production | 3.009 | 0.840 |
| FM_ctl (this round, method control) | 2.979 (**0.990** × prod.) | 1.024 (**1.218** × prod.) |
| PM | 2.553 | 0.609 |
| AFM | 2.557 | 0.827 |
| jackknife PM, 8 sets | 2.52 – 2.59 | **0.47 – 0.88** |

| ratio | in-plane | κ_zz |
|---|---:|---:|
| R = PM / FM_ctl | **0.857 ± 0.018** | 0.595 ± 0.347 |
| R_AFM = AFM / FM_ctl | 0.858 | 0.807 |

The ratios are the same from 100 to 600 K, to within 0.01.

## Verdict

**In-plane: resolved.**
- The method control reproduces the production in-plane κ to 1 %.
- The jackknife error is ±1.8 %.
- PM lowers κ_xx by **14 %**, beyond the 5 % threshold. By the pre-set rule, the
  in-plane L1 moves to the PM FC2.
- PM and AFM give the same ratio (0.857 vs 0.858), so short-range spin order does
  not matter for the in-plane effect.

**κ_zz: not resolved by this dataset.**
- The method control fails: FM_ctl gives +22 % against production.
- The PM jackknife spreads 0.47–0.88 (±58 %).
- The `kappa_summary.py` decision line ("rebase") applies to κ_xx only. No
  statement about the spin effect on κ_zz can be made from these fits.

## Why κ_zz and G2 fail: the interlayer force constants are below the noise

The two failures have the same cause. **The fits do not resolve the weak
interlayer (van der Waals) force constants.**

- **FM_ctl** (all frames VASP 6.4.3, no spin change) already misses the
  production values for exactly the interlayer-sensitive quantities.
  - Γ optical modes: 0.678 / 0.787 / 1.822 THz, against 0.704 / 0.743 /
    1.737 THz in production (−3.6 %, +5.9 %, +4.9 %).
  - LA velocity along z: 3777 m/s, against 3512 m/s in production (+7.5 %).
  - Since κ_zz ∝ v², this alone accounts for most of the +22 %.
- **G2.** At the stiffer Γ modes the AFM−FM shifts agree with the 9-atom result:
  2.46 THz −8.5 % vs −8.1 %; 1.53 THz −2.2 % vs −2.9 %; 3.24 THz +0.3 % vs +0.2 %.
  The 10 pp failures sit on the three soft modes, where FM_ctl itself is off.
  The two AFM frames run with VASP 6.6.0 are therefore **not** needed to explain
  G2. FM_ctl shows the same deviations without them.
- **Why the design under-delivered.** The displacements are 0.01 Å; production
  used σ = 0.017 Å per component, |u| ≈ 0.03 Å.
  - The DFT force noise, judged by the FM_ctl even/odd ratio of 2.2 % at 0.01 Å,
    is about 10⁻³ eV/Å. The self-test assumed 5 × 10⁻⁵ eV/Å, 20× less, so it
    over-predicted the precision.
  - The in-plane force constants are stiff enough to survive this noise. The
    interlayer ones are not.

## Options for κ_zz

1. **Report it as unresolved.** No DFT. The paper states the in-plane effect
   (−14 %) and that the effect on κ_zz could not be resolved. Mn₂In₂Se₅ κ_zz is
   already an upper bound because of the FC3 cutoff.
2. **Paired-difference rerun.** About 3–4 k core-h.
   - Compute PM configurations on the *same displaced geometries* as the
     production 4×4×1 S4 dataset, with the production INCAR verbatim.
   - Fit ΔFC2 = FC2(PM) − FC2(FM) from the force differences frame by frame,
     then use FC2(PM) = FC2(production) + ΔFC2.
   - Geometry-dependent noise largely cancels in the difference, and the
     production amplitude is 3× larger than the 0.01 Å used here.
   - A 2-frame pilot should first measure the noise floor of ΔF.

## Decision (2026-11-27, user)

**Option 1 — report κ_zz as unresolved.** No paired-difference rerun is done. The
paper states the in-plane result (−14 %, R = 0.857 ± 0.018) and that the effect of
spin disorder on κ_zz could not be resolved; Mn₂In₂Se₅ κ_zz stays an upper bound
anyway because of the FC3 cutoff.
