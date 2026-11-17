# Raw kappa(T) files (per-mode phono3py output behind the figures)

These are the hdf5 files that `data_mnin2se4.json` / `data_mnin2se5.json` were
built from, so every temperature point in the figures can be traced to raw data.

| file | material | temperatures (K) | job |
|---|---|---|---|
| `MnIn2Se4_T100-700.hdf5` | MnIn₂Se₄ | 100, 200, 300, 400, 500, 600, 700 | 416235 (14³, 6.45 Å, isotopes on) |
| `MnIn2Se4_T1-150.hdf5` | MnIn₂Se₄ | **1**, 30, 60, 90, 120, 150 | 416310 |
| `Mn2In2Se5_T100-700.hdf5` | Mn₂In₂Se₅ | 100, 200, 300, 400, 500, 600, 700 | 416236 + isotope re-run |
| `Mn2In2Se5_T30-150.hdf5` | Mn₂In₂Se₅ | 30, 60, 90, 120, 150 | 419202 (isotopes on) |

**1 K is present in the MnIn₂Se₄ file but is excluded from the figures**: the Bose
occupation diverges there and κ_zz comes out at 1.9×10⁴ W m⁻¹ K⁻¹, which is
numerical, not physics. 30 K is the lowest usable point.

Every temperature is decomposed with
`kappa_zz_decompose.py --temp-index <i>`, so κ_zz in the figures is the
rest-of-the-zone term at every T, not the Γ–Z-contaminated total.
