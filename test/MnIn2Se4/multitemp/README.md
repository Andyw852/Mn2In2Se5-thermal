# MnIn₂Se₄ — κ(T), 14³ mesh, FC3 6.45 Å, isotopes on

HanHai job 416235. The 300 K row reproduces the committed values exactly
(3.6474 in-plane / 1.0934 κ_zz), which is the self-check that this scan uses the
same force constants as the rest of the report.

| T (K) | in-plane κ | κ_zz |
|---|---:|---:|
| 100 | 11.567 | 3.362 |
| 200 | 5.511 | 1.644 |
| **300** | **3.647** | **1.093** |
| 400 | 2.732 | 0.820 |
| 500 | 2.185 | 0.656 |
| 600 | 1.821 | 0.547 |
| 700 | 1.562 | 0.469 |

Same caveats as everywhere: the FC3 cutoff (6.45 Å) is not converged, and the
κ_zz column is the raw total, which carries the Γ–Z artifact (the converged
quantity is the "rest of the zone" term, 0.6444 at 300 K).

File: `result_multitemp.json`.
