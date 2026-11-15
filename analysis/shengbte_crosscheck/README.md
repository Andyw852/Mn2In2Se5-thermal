# ShengBTE versus phono3py cross-check (Mn₂In₂Se₅, cutoff 6.51 Å, 2026-11-08)

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

The intrinsic κ_zz rises by 0.2–0.27 with every step and shows no convergence: the scattering rate of low-frequency acoustic phonons along the interlayer direction falls very rapidly with frequency and their lifetime is extremely long, so under RTA a denser mesh captures more of them, and the "κ_zz of an infinite perfect single crystal" is set by the sample size rather than being a well-defined number. It does converge at finite size: for a 1 μm grain, 16 → extrapolates to about 1.0 (the last step changes by only 3.8%).
The non-converged part comes from long-wavelength phonons with mean free paths above a micrometre along the interlayer direction; the "intrinsic κ_zz" of an infinite single crystal must be reported as a lower bound with a note,
and it is about 0.95 for a finite grain (1 μm).

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

- κ_xx = κ_yy: **3.14** (ShengBTE iterative solution, 16×16×4, third-order cutoff 6.51 Å); RTA 2.96; the in-plane mesh is converged (20×20×5 differs by 0.9%).
- κ_zz: the intrinsic value grows monotonically with the interlayer mesh and does not converge (RTA already reaches 1.54 at 16×16×16), and it is dominated by interlayer acoustic phonons with long mean free paths;
  report the **finite-grain value**: **≈ 0.97–1.0** at L = 1 μm (RTA, 16×16×16, mesh error < 4%), together with the cumulative κ_zz versus mean free path curve.
- The comparison with hot-pressed samples relies only on the extrinsic-scattering model, which is insensitive to the mesh (last row of the table above).
