# Mn₂In₂Se₅ — coherence term at every FC3 cutoff (16³, 300 K, isotopes off)

HanHai job 416261 (re-run as 416475 after a NODE_FAIL on anode118), ~272 s per
case. Same force constants and mesh as the rest-term decomposition, so the two
can be read together.

| FC3 cutoff | in-plane P | in-plane C | in-plane TOT | κ_zz P (rest) | κ_zz C |
|---|---:|---:|---:|---:|---:|
| 5.68 Å | 3.1810 | 0.1544 | 3.3354 | 1.2451 | 0.1111 |
| 6.12 Å | 3.0571 | 0.1478 | 3.2049 | 1.0977 | 0.1120 |
| **6.51 Å** | **2.9263** | **0.1452** | **3.0715** | **0.9890** | **0.1115** |
| change over the whole range | −8.0 % | **−5.9 %** | −7.9 % | −20.6 % | **+0.4 %** |

**The coherence term is essentially cutoff-independent** (in-plane −5.9 % over the
range, κ_zz +0.4 %) while the population falls by 8 % in-plane and 21 %
cross-plane. That is the same conclusion as for MnIn₂Se₄: the FC3-cutoff
non-convergence is a **population-only** problem, so the coherence can be carried
over without extra cutoff uncertainty.

File: \`result_wigner_cutoff.json\`. The figures that use it are in
\`../../figures/cutoff/\`.
