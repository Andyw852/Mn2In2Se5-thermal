# MnIn₂Se₄ — κ vs FC3 cutoff, with the Wigner decomposition

14³ mesh, 300 K, isotopes on, FC2 = 4×4×2, phono3py 4.4.0 (`transport_type="SMM19"`,
Wigner-RTA). HanHai job 416260, ~48 s per cutoff.

| FC3 cutoff | in-plane P | in-plane C | **in-plane TOT** | κ_zz P | κ_zz C | κ_zz TOT |
|---|---:|---:|---:|---:|---:|---:|
| 5.60 Å | 5.1834 | 0.6388 | 5.8222 | 2.5236 | 0.5738 | 3.0974 |
| 6.23 Å | 3.9944 | 0.6392 | 4.6337 | 1.1297 | 0.6010 | 1.7306 |
| **6.45 Å** | **3.6474** | **0.6539** | **4.3013** | **1.0934** | **0.5982** | **1.6916** |
| step 5.60→6.23 | −22.9 % | +0.1 % | −20.4 % | −55.2 % | +4.7 % | −44.1 % |
| step 6.23→6.45 | −8.7 % | +2.3 % | −7.2 % | −3.2 % | −0.5 % | −2.3 % |

## The finding that matters

**The coherence term is essentially independent of the FC3 cutoff**: in-plane
0.639 → 0.639 → 0.654 (+2.4 % over the whole range) while the population falls
by 30 %. So the cutoff non-convergence is a **population-only** problem, and the
coherence can be added to whichever population value is used without extra
cutoff uncertainty.

⚠️ The κ_zz columns are the raw totals and carry the Γ–Z quadrature artifact
(half of the 11³ total comes from the single first Γ–Z point). The coherence
term for κ_zz is also flat (0.574 → 0.601 → 0.598), but a rest-term decomposition
of the Wigner κ_zz is not available — `kappa_zz_decompose.py` only handles the
population per-mode output.

Files: `result_wigner_cutoff.json`.
