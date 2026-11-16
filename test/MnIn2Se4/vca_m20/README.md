# MnIn₂Se₄ — VCA (i = 0.44) versus ordered baseline, mesh 20³

VCA = virtual crystal: the alloy is replaced by a crystal with **averaged masses**
and the (averaged) force constants, and the disorder is put back as **mass-variance
scattering**.

## Inputs

| item | value |
|---|---|
| structure | `Test/MnIn2Se4/fit/POSCAR` equivalent, 7 atoms, order `Mn In In Se Se Se Se` |
| FC2 | 4×4×2, `../fc2_442_phonon/fc2.hdf5` |
| FC3 | 4×4×1 at **6.45 Å** (the largest cutoff this cell allows), `../fit_cut3_6p45/fc3.hdf5` |
| primitive matrix | `[[1,1,2],[-1,-2,-1],[1,1,1]]` |
| mesh | 20×20×20 |
| T | 300 K |
| averaged masses | Mn (octahedral) **81.29** amu, In (tetrahedral) **101.65** amu |
| mass variances | `[0.1337, 0.0595, 0.0595, 0.0003, 0.0003, 0.0003, 0.0003]` |

The occupancies behind these numbers (octahedral 0.56 Mn + 0.44 In, tetrahedral
0.22 Mn + 0.78 In, i = 0.44) come from the user's
`run_phono3py_kappa_vca.sh`; **the bibliographic source for them is not yet
verified** — see `analysis/disorder/README.md` §5.

## Runs

HanHai job **419636**, one job, two phono3py RTA runs:

1. **baseline** — real masses, `is_isotope=True` (natural isotopes)
2. **VCA** — averaged masses, `mass_variances` as above, `is_isotope=False`

`jobs/run.py` and `jobs/submit.sh` are exactly what ran.

## Results

| quantity | baseline | **VCA** | suppression |
|---|---:|---:|---:|
| in-plane κ | 3.6884 | **1.4891** | **2.48×** |
| κ_zz total (raw) | 0.8780 | 0.6954 | 1.26× |
| **κ_zz rest of the zone** | **0.6380** | **0.4455** | **1.43×** |

`result.json` holds the raw κ values; the rest-of-the-zone terms come from
`kappa_zz_decompose.py` applied to the two hdf5 files.

## Files

| file | contents |
|---|---|
| `result.json` | κ for both runs |
| `jobs/run.py`, `jobs/submit.sh` | the calculation as it ran |
| `../../test/figures/scattering/` | the two-panel scattering figure, and the two `kappa_*.hdf5` per-mode outputs |
| `../../analysis/disorder/README.md` | how this fits together with the isotope (L1) and extrinsic-fit (L3) treatments |
