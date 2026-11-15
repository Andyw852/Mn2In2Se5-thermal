# MnIn2Se4 fc3 cutoff convergence scan (κ)

**2026-11-14** | using the existing 46-frame random-displacement dataset, three fc3 cases (**5.6 / 6.23 / 6.45 Å**) were re-fit with pheasy
and κ was computed for each. fc2 is fixed to the new **4×4×2** (`../fc2_442_phonon/fc2.hdf5`).

## Why stop at 6.45 Å

The thickness of the 4×4×1 supercell along the layer normal is **13.214 Å**, so the safe-cutoff limit is h/2 − 0.1 = **6.507 Å**.
The largest shell that can be included is **6.3557 Å**; the next shell at 6.5311 Å already exceeds the limit.

**So 6.45 Å is already the physical limit of 4×4×1 and cannot be pushed any higher.**
Testing a larger cutoff requires a 4×4×2 fc3 supercell (safe limit 8.048 Å) and hence a new DFT dataset.

## Data source

The three fc3 cases are taken directly from `cutoff_scan/` in the original AutoZT fit directory (pheasy RFE-OLS fit, 2026-11-10):

```
$HOME/taskflow/work/MnIn2Se4_441/fit-fc-thermal/step1_fit/cutoff_scan/cut3_5p60/fc3.hdf5
                                                          .../cut3_6p23/fc3.hdf5
                                                          .../cut3_6p45/fc3.hdf5
```

**Consistency check**: the `cut3_6p45` scan copy and the production copy have different sha256 values (the scan is from 11-10, the production run from 11-11 includes fixes),
but the κ computed with the new 4×4×2 fc2 is **completely identical** (in-plane 3.5375 / κzz 1.2495) → the scan copy is trustworthy.

> ⚠️ Refitting locally with **symfc** instead of pheasy was attempted: the result is unusable
> (at the same 6.45 Å cutoff symfc gives an in-plane value of 0.644 versus pheasy 3.5375, a factor of 5.5;
> the three cases differ by a factor of 7 and are non-monotonic). The AutoZT driver comments also note a symfc error of 1.95% versus pheasy 0.578%.
> **This dataset must use pheasy's RFE-OLS sparse fit.**

## Results (300 K, mesh 11×11×11, isotopes on, fc2 = 4×4×2)

### RTA

| fc3 cutoff (Å) | in-plane (W/m·K) | κzz (W/m·K) | κzz/in-plane |
|---|---:|---:|---:|
| 5.60 | 4.9415 | 2.7672 | 0.560 |
| 6.23 | 3.8922 | 1.2332 | 0.317 |
| **6.45** | **3.5375** | **1.2495** | **0.353** |

### LBTE (production setting)

| fc3 cutoff (Å) | in-plane (W/m·K) | κzz (W/m·K) | κzz/in-plane |
|---|---:|---:|---:|
| 5.60 | 5.6473 | 2.7388 | 0.485 |
| 6.23 | 4.7842 | 1.2241 | 0.256 |
| **6.45** | **4.4075** | **1.2519** | **0.284** |

### Change between adjacent cases (LBTE)

| interval | in-plane change | κzz change |
|---|---:|---:|
| 5.60 → 6.23 | **−15.3%** | **−55.3%** |
| 6.23 → 6.45 | **−7.9%** | **+2.3%** |
| 5.60 → 6.45 | −22.0% | −54.3% |

## Conclusions

1. **κzz is essentially flat from 6.23→6.45 (+2.3%)**, but it drops by 55% from 5.60→6.23 — which means
   that the region around 6.45 Å is the **critical range** for κzz, and it sits right at the 4×4×1 ceiling.
2. **The in-plane value still changes by −7.9% in the last step** (−15.3% in the previous one) — convergence is improving but **has not been reached**.
3. **Verdict: the fc3 cutoff is not converged.** The cutoff can no longer be increased within 4×4×1, so a citable κ requires:
   - a **4×4×2 fc3 dataset** (224 atoms, safe cutoff 8.048 Å), scanning 7.0 / 7.6 / 8.05 Å;
   - a mesh convergence study starting from 11×11×11 as well.

**The most reliable numbers currently available (4×4×2 fc2 + 6.45 Å fc3, LBTE, 300 K)**:
in-plane **4.4075** W/m·K, κzz **1.2519** W/m·K — but it must be noted that neither the fc3 cutoff nor the mesh is converged.

## Reproduction

```bash
# after fetching the three fc3 files from 3090
python kappa_cutscan.py     # three cases with RTA (about 5 minutes)
python kappa_lbte_scan.py   # three cases with LBTE (about 47 minutes)
```

Raw data: `kappa_cutoff_scan_results.json`
