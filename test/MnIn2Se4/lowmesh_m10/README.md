# MnIn₂Se₄ — low-mesh (10³) cross-checks for the LBTE and coherence terms

The 14³ Wigner-LBTE job sat in phono3py's **single-threaded Python loop** for over
an hour (1 of 64 cores busy, BLAS threads idle — the module has no MPI support in
4.4.0), so these 10³ runs were added as a cheap way to get the *ratios*.

| quantity | 10³ | 14³ (converged run) |
|---|---:|---:|
| Wigner-RTA, population | 3.5027 | 3.6474 |
| Wigner-RTA, coherence | 1.5230 | 0.6539 |
| Wigner-RTA, TOT | 5.0257 | 4.3013 |
| Wigner-LBTE, TOT | 4.3808 | **4.5219** |

## What the 10³ runs establish

1. **The coherence term is not converged at 10³** (1.52 vs 0.65 at 14³), so the
   10³ TOT cannot be used directly — only as a shape check.
2. **The 4.4.0 SMM19-LBTE does not add the coherence to its TOT.** At 14³ the
   LBTE returns `TOT = population_LBTE = 4.5219` while the coherence attribute is
   0.6539; the naive sum would be 5.176. The 10³ run shows the same pattern.
   So the in-plane value is quoted as **≈4.5 with the coherence giving an upper
   bound of 5.2**, rather than as a single 5.18.
3. The population-only LBTE/RTA ratio (1.24 at 14³) is much larger than the
   Wigner-TOT ratio at 10³ (0.87), so the two families must not be mixed.

Files: `result_wigner_rta_m10.json`, `result_wigner_lbte_m10.json`
(HanHai job 415924, a few minutes).
