# Mn₂In₂Se₅ — LBTE/RTA ratio from a 10³ run

The in-plane κ of the reports uses an **iterative (LBTE) correction** that was
borrowed from ShengBTE (×1.06). The 16³ phono3py LBTE (job 414720) stalls in
phono3py's single-threaded loop, so this 10³ run (job 415925, minutes) measures
the ratio on the same force constants instead.

| | in-plane | κ_zz |
|---|---:|---:|
| LBTE | 2.9879 | 1.8972 |
| RTA (same run) | 2.8220 | 1.9275 |
| **ratio** | **1.0588** | **0.9843** |

The in-plane ratio **1.059 agrees with the borrowed ShengBTE 1.06**, which is
what the reports use: in-plane κ = 3.0825 (RTA, 32³ converged) × 1.059 = 3.265,
+ 0.1458 coherence = **3.41**.

At 10³ the κ_zz totals carry a large Γ–Z quadrature artifact (1.90 vs the
converged 0.96), so the κ_zz ratio is not used.

File: `result_lbte_m10.json`.
