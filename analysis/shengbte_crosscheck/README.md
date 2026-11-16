# ShengBTE versus phono3py cross-check (Mn₂In₂Se₅, cutoff 6.51 Å, 2026-11-08)

## Values to quote at ≈10 % accuracy (updated 2026-11-16)

*This section supersedes the mesh and finite-grain statements further down,
which were written before the Γ–Z decomposition of κ_zz
(see `test/Mn2In2Se5/kappa_zz_decomposition/`).*

| quantity | value | composition |
|---|---:|---|
| **in-plane κ** | **≈ 3.4 (3.3–3.4)** | RTA **3.0825** at 32³ (converged, last step +0.3 %) **× 1.0590** — the iterative (LBTE) correction now **measured with this same phono3py stack** (LBTE/RTA = 2.9879/2.8220 on a 10³ mesh, job 415925; the ShengBTE value 1.06 agrees) = 3.265, **+ coherence 0.1458** (Wigner-RTA at 16³, job 414720) = **3.41** |
| in-plane FC3 cutoff | **−4.3 % in the last step (−8.0 % over 5.68 → 6.51 Å)** | redone on an isotropic 16³ mesh; inside the 10 % tolerance but not flat. The earlier "−2.7 % then −1.5 %" came from the 12×12×3 mesh and does not survive |
| **κ_zz** | **0.955 ± 0.032** | rest term converged over four isotropic meshes (16³/24³/28³/32³: 1.0090 / 0.9766 / 0.9841 / 0.9801, last step −0.41 %); in-plane RTA converged at 3.0825 (32³). The z-only totals 0.84 → 1.54 are the Γ–Z quadrature artifact and must not be quoted. **Still cutoff-limited**, not mesh-limited |
| κ_zz FC3 cutoff | **not converged: −9.9 % in the last step, −20.6 % over 5.68 → 6.51 Å** | measured on an isotropic 16³ mesh with the rest-term decomposition (job 412029). **So Mn₂In₂Se₅ κ_zz needs a 4×4×2 FC3 to be citable at 10 %**, as MnIn₂Se₄ does. Fit-to-fit noise is ≈2 %, well below the −10 % step. See `test/Mn2In2Se5/kappa_zz_decomposition/README.md` |

**Quote in-plane κ from the iterative solution, not RTA.** The RTA and LBTE
in-plane values differ by ≈6 % here (and by 24 % for MnIn₂Se₄), so the RTA
number alone is not accurate to 10 %.

## Conclusion

**The ShengBTE κ from S6 on the cluster (`Mn2In2Se5/data/kappa/`) is superseded.** The cause lies neither in ShengBTE nor in the force-constant fit,
but in the format conversion `fc3.hdf5 → FORCE_CONSTANTS_3RD`: AutoZT calls hiphive's `write_to_shengBTE`, which only accepts
atom pairs that are already the shortest image inside the supercell without wrapping; the 4×4×1 supercell has only one layer along a₃′, so every triplet crossing that boundary is dropped.

| | number of triplets (first atom in the central cell) |
|---|---:|
| non-zero triplets in `fc3.hdf5` | 9164 |
| `FORCE_CONSTANTS_3RD` on the cluster (hiphive export) | 6038 (missing 3126, 34%) |
| the missing ones: interatomic distance | 2.8–6.35 Å |
| the missing ones: \|Φ\| median / maximum | 0.035 / 4.4 eV/Å³ (retained: 0.042 / 32.7) |

## Same fc2/fc3, 300 K, RTA, q-mesh 16×16×4 (W m⁻¹ K⁻¹)

| calculation | κ_xx = κ_yy | κ_zz |
|---|---:|---:|
| cluster ShengBTE (hiphive-exported 3RD) | 1.7224 | 0.2098 |
| local ShengBTE + the **same** cluster 3RD (`run_original_fc3/`) | 1.7224 | 0.2098 |
| local ShengBTE + the **completed** 3RD (`run_complete_fc3/`) | **2.9535** | **0.8346** |
| phono3py RTA (`analysis/extrinsic/results/Mn2In2Se5_m16164/`) | **2.9743** | **0.8364** |

Local ShengBTE reproduces the cluster results bit for bit with the original file, and agrees with phono3py once the completed file is used (0.7% / 0.2% difference).
The completed file is written from `fc3.hdf5` by `write_fc3_sheng.py`: atoms 2 and 3 of each triplet are taken as the shortest images relative to atom 1,
and the R vectors are relative to the CONTROL coordinates; after writing, everything is folded back by the supercell and compared element by element with `fc3.hdf5`, giving a difference of 5×10⁻¹⁵.

## q-mesh (phono3py RTA; κ_zz is sensitive to the interlayer-direction mesh)

| q-mesh | κ_xx(300 K) | κ_zz(300 K) |
|---|---:|---:|
| 16×16×4 | 2.974 | 0.836 |
| 16×16×8 | 2.958 | 1.109 |
| 20×20×5 | 3.001 | 0.836 |
| 16×16×12 | 2.963 | 1.304 |
| 16×16×16 | 2.959 | 1.541 |

In-plane mesh: 16×16×4 → 20×20×5 changes κ_xx by only +0.9%, so 16×16×4 is already sufficient for κ_xx (the cluster ShengBTE run at 20×20×5
blew up to 750 GB of memory in the iterative solve, so it does not need to be rerun).

κ_zz after adding grain boundaries (Casimir), versus the third number 4 / 8 / 12 / 16:

| grain L | m = 4 | 8 | 12 | 16 |
|---|---:|---:|---:|---:|
| ∞ (intrinsic) | 0.836 | 1.109 | 1.304 | 1.541 |
| 10 μm | 0.826 | 1.071 | 1.229 | 1.392 |
| 1 μm | 0.754 | 0.877 | 0.937 | 0.973 |
| fitted region (x = 0.13, L = 16 nm) | 0.158 | 0.163 | 0.165 | 0.167 |

> **Superseded (2026-11-16).** The explanation below — that the rise is
> long-wavelength physics and that the intrinsic κ_zz "is set by the sample
> size" — is **withdrawn**. The Γ–Z decomposition shows that the rise under
> z-only refinement is the same q-point quadrature artifact that dominates
> MnIn₂Se₄'s κ_zz, and that isotropic refinement *lowers* κ_zz (16³ → 24³:
> 1.547 → 1.305). See `test/Mn2In2Se5/kappa_zz_decomposition/README.md` and the
> "Values to quote" section at the top. The numbers in this table are correct as
> totals; they are simply not values of κ_zz.

κ_xx is converged (< 1%); the intrinsic κ_zz rises by 33% → 18% from 4 → 8 → 12 and is not yet converged (the long-wavelength, low-frequency acoustic modes along the interlayer direction
have very long mean free paths, and a denser mesh captures more of them).

**This is insensitive when comparing with experiment**: adding the extrinsic scattering of the fitted region (e.g. x = 0.13, L = 16 nm), κ_zz at 300 K is
0.158 / 0.163 / 0.165 on the three meshes, and κ_xx is 0.853 in all cases — the modes being cut off are exactly the long-range ones. The texture/extrinsic-scattering conclusions are therefore unaffected by the mesh;
only the "intrinsic single-crystal κ_zz" number needs a dense mesh.

## The AutoZT fix

**[0127]** on the AutoZT branch `claude/jolly-curie-jsvosv`: all ShengBTE 3RD exports now use `fc_common.write_shengbte_fc3`
(full write-out + fold-back self-check); S5/S6 compare against `fc3.hdf5` before reusing an old file and re-export automatically on a mismatch. After merging, a `retry` of S6 is all that is needed.

The note in the old `MnIn2Se4/archive_2026-10-28_fit-fc-thermal/` that "ShengBTE is 24× lower than phono3py, judged to be a ShengBTE solver problem"
may have the same cause, but the provenance of that 3RD (hiphive or the native pheasy export) was not archived, so a conclusion requires checking it first.

## Suggested wording for the paper (Mn₂In₂Se₅, 300 K)

- κ_xx = κ_yy: **≈ 3.3** (RTA 3.020 at 24³, in-plane mesh converged to ≈1 %, corrected by ×1.06 for the iterative solution and + 0.135 for the coherence term; the converged ShengBTE iterative value is 3.14). Do **not** quote the raw RTA 2.96–3.02 as the final number.
- κ_zz: **0.955 ± 0.032** — the rest term is converged over four isotropic meshes (16³/24³/28³/32³: 1.0090 / 0.9766 / 0.9841 / 0.9801, last step −0.41 %); at a 10 μm grain 0.928 ± 0.039, at 1 μm 0.788 ± 0.044. The z-only scan 0.84 → 1.54 must **not** be quoted: it is the Γ–Z quadrature artifact (isotropic 16³ → 24³ lowers κ_zz by 15.6 %). The earlier "finite-grain ≈0.97–1.0 at 1 μm" recommendation is withdrawn. The value is **cutoff-limited**: the last available FC3 step still moves the rest term by −9.9 % (see `test/Mn2In2Se5/kappa_zz_decomposition/README.md`).
- The FC3 cutoff: κ_zz is **not converged** (−9.9 % in the last available step, −20.6 % over 5.68 → 6.51 Å, on an isotropic 16³ mesh with the rest-term decomposition); the in-plane component moves −4.3 % in the last step (−8.0 % over the range). Both are at the 4×4×1 ceiling (6.61 Å), so a 4×4×2 FC3 is the only way to settle them.
- The comparison with hot-pressed samples relies only on the extrinsic-scattering model, which is insensitive to the mesh (last row of the table above).
