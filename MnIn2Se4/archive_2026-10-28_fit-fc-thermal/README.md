# MnIn2Se4 lattice thermal conductivity calculation (fit-fc-thermal)

Lattice thermal conductivity κ calculation for the trigonal layered van der Waals material MnIn₂Se₄ (space group R-3m, 166; 7-atom primitive cell, 3×3×3 supercell with 189 atoms, interlayer spacing 13.2 Å). The force constants are fitted by pheasy (ALASSO) from the displacement-force dataset, κ is computed by phono3py (RTA), and symfc provides the second-order control.

## Key results

> Note: the κ values below are superseded — see `../../REPORT_MnIn2Se4_20261115.md`.

| quantity | value |
|---|---|
| in-plane κ(300 K) | **4.03 W/m·K** (xx=yy, 0.60% change from 11→14 mesh, converged) |
| out-of-plane κ(300 K) | **1.68 W/m·K** (zz, 8.37% change from 11→14 mesh, **not converged**) |
| fit relative error | 0.646% (ALASSO, OLS limit α→0, 1303 free IFC) |
| phonon stability | stable, no significant imaginary frequency (min_freq −0.019 THz, noise level) |
| third-order cutoff | 5.0 Å (largest, the largest that can be determined under the safe cutoff 5.17 Å) |
| mesh | 14×14×14 (mesh_converged=false, cap triggered) |

The full-temperature κ(T) (100–800 K) is given in `s2_kappa/kappa_summary.json` and in `kappa_vs_mesh.png` / `kappa_vs_cutoff.png`.

## Main diagnostic conclusions

1. **The zz (out-of-plane) direction is not converged**: κ_zz keeps decreasing with the mesh (1.84 → 1.69, 8.4%), which is the intrinsic slow convergence of out-of-plane transport in layered materials.
2. **The root cause is that the cutoff is pinned by the safety limit**: `supercell_safe_cutoff = 0.5 × shortest perpendicular height − 0.1 ≈ 5.17 Å`, while the interlayer spacing is 13.2 Å, so **the interlayer third-order force constants are cut off entirely**; switching to CUT3_PICK=largest (candidates 4.3/4.7/5.0, thirdorder −n 2nd/3rd/4th neighbours) only makes it cleaner — it cannot add interlayer terms. The cutoff scan up to 5.0 Å still has not reached a κ plateau (`cutoff_selection.status=not_converged`).
3. **ShengBTE and phono3py differ by 24× for this compound** (ShengBTE 0.159 vs phono3py 3.85, same FORCE_CONSTANTS): the force-constant file format (FC2 frequencies agree, FC3 roundtrip agrees), R-vector scaling, scalebroad, volume, heat capacity and group velocity have all been ruled out; the conclusion is that the ShengBTE solver overestimates the three-phonon scattering rate for this large 7-atom/21-branch system (its RTA is even below its own small-grain limit, which is physically anomalous). shengbte-gpu (single- and multi-GPU) agrees with the official CPU version, ruling out a bug in the GPU fork. κ is therefore taken from phono3py.

## Directories

- `s1_fit/` — S1 force-constant fitting results (phonon_summary, cutoff_scan, fc_fit_summary, fit_config, fit_manifest, fit_metrics, POSCAR, SPOSCAR)
- `s2_kappa/` — S2 lattice thermal conductivity (kappa_summary, mesh_convergence, kappa_vs_mesh, kappa_vs_cutoff)

## Key points for reproduction

- fit-fc-thermal skill, pheasy ALASSO + phono3py-symfc control.
- S2 parameters: isotope=false, mesh_conv_max_length=50 (cap 14×14×14), cut_scans=[4p30,4p70,5p00], nominal_cut3_A=5.0.
- Data: displacement-force dataset (not included in this repository; must be provided separately).
