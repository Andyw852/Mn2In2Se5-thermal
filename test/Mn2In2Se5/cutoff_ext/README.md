# FC3 cutoff extension for Mn₂In₂Se₅ — refit on the existing 4×4×1 frames

Five **pheasy-gpu ALASSO** fits on the existing 54-frame dataset (144-atom 4×4×1
supercell), 6.51 Å kept as the anchor, in a **separate directory** so the
production `step5_fc/` config is untouched. Plan and the two GPU traps:
`../../../dft/fc3_cutoff_extend/README.md`.

## Fit quality (job 7701, 3 h 48 min for the whole scan)

| FC3 cutoff | free IFCs | data/parameter ratio | train rel. error | train rel. s.e. |
|---|---:|---:|---:|---:|
| 6.51 (anchor) | 3813 | 6.12 | 0.408 % | 0.014 % |
| 6.67 | 3872 | 6.02 | 0.407 % | 0.014 % |
| 6.90 | 4715 | 4.95 | 0.398 % | 0.015 % |
| 7.11 | 6580 | 3.55 | 0.375 % | 0.022 % |
| **7.43** | 8665 | **2.69** | 0.351 % | 0.026 % |

The ratio crosses below 3 at 7.43 Å; that point is kept by user decision
(2026-11-20, exception granted), so no new DFT frames were added.

## κ (16³, 300 K, isotopes **on**, production FC2 `b611f139…` fixed)

phono3py 4.4.0 RTA, `primitive_matrix="P"`; κ_zz split with
`test/MnIn2Se4/kappa_zz_decompose.py` at every cutoff.

| FC3 cutoff | in-plane κ | κ_zz total | **κ_zz rest of zone** | first Γ–Z share |
|---|---:|---:|---:|---:|
| **6.51 (anchor)** | **2.9348** | 1.5161 | **0.9801** | 35.4 % |
| 6.67 | 2.9347 | 1.5935 | 1.0064 | 36.8 % |
| 6.90 | 2.8731 | 1.4666 | 0.9478 | 35.4 % |
| 7.11 | 2.7767 | 1.4021 | 0.8745 | 37.6 % |
| **7.43** | **2.5548** | 1.1868 | **0.7386** | 37.8 % |

Steps relative to the anchor: in-plane −0.0 %, −2.1 %, −5.4 %, **−13.0 %**;
κ_zz rest +2.7 %, −3.3 %, −10.8 %, **−24.6 %**.

**Anchor check.** In-plane 2.9348 is inside the expected 2.89–3.00. κ_zz rest
0.9801 is 0.5 % **below** the existing 6.51 Å value 0.9850 — outside the narrow
0.985–1.009 band quoted for the refit, but well inside the ≈2 % fit-to-fit
scatter that band was derived from.

**Verdict: the FC3 cutoff is not converged for Mn₂In₂Se₅.** Between 6.9 and
7.43 Å both components fall, κ_zz rest by a quarter. Two readings cannot be
separated with this dataset: real long-range interlayer terms that the 4×4×1
cell was missing, or the 7.11–7.43 Å fits absorbing noise (the data/parameter
ratio is 3.55 and 2.69 there, and the s.e. grows). Either way the available
frames do not support a converged value above 7 Å, and a **larger cell (5×5×1,
9.26 Å) or more frames** would be needed; that is a decision for the user, not
something to start without it.

## Files

| file | contents |
|---|---|
| `kappa_vs_cutoff.json` | κ for each cutoff, phono3py 4.4.0 |
| `kappa_cut<m>_m16.hdf5` | per-mode phono3py output of each case (1.6 MB each) |
| `cutoff_scan.json` | the fit log: free IFCs, ratio, train error and s.e. per cutoff |

The fitted `fc3.hdf5` files (9.7–18.4 MB each, 64 MB total) are **not** in git;
they are on the backup drive at
`E:\backup\joint_research_project\xia\raw_hdf5\cutoff_ext\`.
