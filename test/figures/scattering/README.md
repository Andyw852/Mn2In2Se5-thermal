# VCA (virtual crystal) versus ordered baseline — MnIn₂Se₄ at the largest FC3 cutoff

**Recipe** (the user's `run_phono3py_kappa_vca.sh`, reproduced here):

> VCA treats the alloy as a "virtual crystal" with averaged masses and averaged
> force constants, and puts the disorder back in as **mass-variance scattering**.

| site | occupancy (i = 0.44) | mass (amu) | mass variance g |
|---|---|---:|---:|
| octahedral | 0.56 Mn + 0.44 In | **81.29** | **0.1337** |
| tetrahedral | 0.22 Mn + 0.78 In | **101.65** | **0.0595** |
| Se | natural | 78.96 | 0.0003 |

The primitive cell is 7 atoms in the order `Mn In In Se Se Se Se` (verified against
`test/MnIn2Se4/phono3py_disp.yaml`), so the mass-variance list is
`[0.1337, 0.0595, 0.0595, 0.0003, 0.0003, 0.0003, 0.0003]`.

## What was run

HanHai job **419636**, mesh **20×20×20**, T = 300 K, **FC3 = 6.45 Å** (the largest
cutoff this cell allows), FC2 = 4×4×2, primitive matrix as in the production runs.
Two phono3py RTA runs in one job:

| run | masses | scattering |
|---|---|---|
| **baseline** | real | natural isotopes (`is_isotope=True`) |
| **VCA** | averaged (81.29 / 101.65) | mass variances (`mass_variances`) |

## κ at 300 K

| quantity | baseline | **VCA (i = 0.44)** | suppression |
|---|---:|---:|---:|
| in-plane κ | 3.6884 | **1.4879** | **2.48×** |
| κ_zz total (raw, Γ–Z contaminated) | 0.8780 | 0.6954 | 1.26× |
| **κ_zz rest of the zone** | **0.6380** | **0.4453** | **1.43×** |

The in-plane component is suppressed much more strongly than the cross-plane one,
because it is carried by higher-frequency modes where the mass-variance rate
(∝ g·ω²) is large.

**Against the hot-pressed experiment** (MnIn₂Se₄: Perp 0.93, Para 0.485 at 300 K):
the VCA cross-plane value 0.446 is already close to the measured 0.485, while the
in-plane value 1.489 is still 1.6× above the measured 0.93 — i.e. mass disorder
alone explains most of the cross-plane suppression but additional scattering
(grain boundaries, strain) is needed in-plane. That is consistent with the fit in
`analysis/extrinsic/`.

## Files

| file | contents |
|---|---|
| `scattering_compare.png / .pdf` | the two-panel figure (see below) |
| `kappa_baseline.hdf5`, `kappa_vca.hdf5` | per-mode phono3py output for both runs |
| `vca_result.json` | the κ values |
| `plot_scattering_figure.py` | `analysis/plotting/plot_scattering_compare.py` with y-axis limits added (the raw data contain near-Γ modes at 10⁻³⁸) |

## Reading the figure

- **Left**: 3-phonon rates of the ordered and the VCA crystal. They agree to within
  a few tens of percent — VCA changes the masses, not the anharmonic strength, so
  the three-phonon channel is essentially preserved.
- **Right**: the disorder channel. VCA mass-variance scattering is **2–3 orders of
  magnitude stronger** than natural isotope scattering, follows ω² at low
  frequency and flattens at high frequency where the mode density saturates.
  This is why adding natural isotopes alone (the step done earlier) barely moves κ,
  while the antisite disorder does.
