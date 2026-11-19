# Extrinsic scattering + texture: converting single-crystal intrinsic κ into hot-pressed polycrystalline experimental values

## Model (settled)

Scattering rates add up (RTA):

```
Γ_λ = Γ_ph-ph,λ(T)  +  Σ_elem g_elem(x) · B_elem,λ  +  |v_λ| / (4π L)
       three-phonon (intrinsic)    mass disorder of Mn↔In antisites (Tamura)   grain boundary (Casimir)
```

The texture parameter s is then used to combine κ_xx and κ_zz into Para / Perp (see the repository-root README). There are three parameters:

| parameter | meaning | range (scan) | note |
|---|---|---|---|
| **x** | fraction of all Mn that sits on In sites (Mn↔In pairwise exchange, stoichiometry unchanged) | 0 → complete-disorder limit, step 0.01 | the Mn site holds (1−x)Mn + x In; the Mn fraction on In sites = x·N_Mn/N_In. Upper bound Mn₂In₂Se₅ = 0.5, MnIn₂Se₄ = 2/3 |
| **L** | grain-boundary mean free path (isotropic) | 5 nm – 100 μm, 61 log points, plus ∞ | hot-pressed grains are typically 0.5–5 μm; stacking faults and twins in layered materials make the effective L far smaller than the grain |
| **s** | ⟨sin²θ⟩, where θ is the angle between the grain c axis and the pressing direction | 0 – 2/3 | 0 = fully textured, 2/3 = random orientation |

**Why only Mn↔In antisites**: the two cations differ in mass by more than a factor of two (54.9 vs 114.8), swapping them costs little within a given structure, and this is the most common
disorder in layered Mn–In–Se compounds; for the Se sublattice only the natural isotopes are counted. Counting the mass difference alone is a **lower bound** — differences in bond strength and local strain make the real scattering stronger,
so the fitted x is an "effective x" that also absorbs the point-defect-like part contributed by vacancies and stacking faults.

**How the fit is done**: only the high-temperature range is used (300–623 K, both directions), minimizing the rms of ln(model/experiment);
a grid scan over (x, L) is performed, with s optimized at every (x, L). The optimum and the "acceptable region" (rms ≤ max(1.25×optimum, optimum+0.02)) are reported;
x and L usually form a degenerate band (both push κ down), in which case the acceptable region governs and XRD / SEM / occupancy refinement are used to break the degeneracy.
The low-temperature range (55–145 K) is not used in the fit; it only checks whether the model matches there.

**Electronic thermal conductivity**: currently treated as 0. If κ_e is not negligible, the fitted extrinsic scattering will be too strong.

## Scripts

`kappa_extrinsic.py`:

1. `bases`: run phono3py once for the intrinsic RTA (three-phonon only), then set g = 1 on each element site to compute only the mass-disorder scattering (read_gamma, without recomputing the three-phonon terms),
   and store `bases.npz`. Afterwards any (x, L) is just a sum over modes — a few milliseconds per point.
2. `eval --x --L`: prints κ_xx and κ_zz(T).
3. `fit --exp ../../experiment/kappa_tot_digitized.csv --material ...`: performs the fit and writes `fit_<material>.json`.

Validation: with the same force constants, the κ computed directly by phono3py with `mass_variances` + `boundary_mfp` is compared against the linear-combination result of this script;
see the "Validation" section below.

## Mn₂In₂Se₅: to be run after the force constants are uploaded (done on the repository side)

```bash
cd analysis/extrinsic
F=../../Mn2In2Se5/fc
# 1) intrinsic + scattering bases, same mesh as the ShengBTE production calculation
python3 kappa_extrinsic.py bases --disp-yaml $F/phono3py_disp.yaml --fc2 $F/fc2.hdf5 --fc3 $F/fc3.hdf5 \
        --mesh 16 16 4 --out work_Mn2In2Se5_m16164
# 2) phono3py RTA vs ShengBTE RTA (same force constants; κ(300 K) should agree within a few percent)
python3 kappa_extrinsic.py eval --bases work_Mn2In2Se5_m16164/bases.npz
# 3) κ_zz mesh: double the third number, then double again (16 16 8, 16 16 16, both symmetry-preserving)
python3 kappa_extrinsic.py bases ... --mesh 16 16 8 --out work_Mn2In2Se5_m16168
# 4) fit on the converged mesh
python3 kappa_extrinsic.py fit --bases work_Mn2In2Se5_m16168/bases.npz \
        --exp ../../experiment/kappa_tot_digitized.csv --material Mn2In2Se5
```

## MnIn₂Se₄ with the 4×4×2 FC2 and the redo for Mn₂In₂Se₅ (2026-11, pending)

`bases` now takes `--dim-fc2` (FC2 on a larger supercell than FC3) and `--pmat`
(primitive matrix; MnIn₂Se₄ needs it so that an N³ mesh is symmetric, as in all
its κ runs). Without these options the behaviour is unchanged.

Check (11³, 300 K, MnIn₂Se₄ 4×4×2 FC2 + 6.45 Å FC3): `eval --x 0` gives
3.5375 / 1.2495, identical to the L1 RTA of report §6.3. At x = 0.44 (ordered
masses, the L3 treatment), in-plane is 1.48. That is close to the VCA (L2)
in-plane of 1.49 (20³), so the mass averaging of L2 matters little in-plane.

```bash
cd analysis/extrinsic
M=../../test/MnIn2Se4
python3 kappa_extrinsic.py bases --poscar $M/fit/POSCAR --dim 4 4 1 --dim-fc2 4 4 2 \
        --pmat 1 1 2 -1 -2 -1 1 1 1 \
        --fc2 $M/fc2_442_phonon/fc2.hdf5 --fc3 $M/fit_cut3_6p45/fc3.hdf5 \
        --mesh 24 24 24 --out work_MnIn2Se4_442_m24
# Mn2In2Se5: same FC set as the kappa_zz decomposition, isotropic mesh
python3 kappa_extrinsic.py bases --disp-yaml ../../Mn2In2Se5/fc/phono3py_disp.yaml \
        --fc2 ../../Mn2In2Se5/fc/fc2.hdf5 --fc3 ../../Mn2In2Se5/fc/fc3.hdf5 \
        --mesh 24 24 24 --out work_Mn2In2Se5_m24
```

See `../disorder/README.md` §2 for which x to use for MnIn₂Se₄ (i = 0.44 from
Range et al. 1991, adopted; Mn₂In₂Se₅ x ≈ 0.21 from Williams et al. 2024, §5.1 there). The fit itself leaves x undetermined.

## Results

### 2026-11-17 — refit on 24³ with the corrected force constants

Both materials were refitted from scratch. MnIn₂Se₄ had to wait for
`--dim-fc2/--pmat` (its FC2 lives on 4×4×2 while FC3 is 4×4×1).

Bases (24³, RTA, three-phonon only):

| | command | intrinsic κ(300 K) xx / zz |
|---|---|---:|
| MnIn₂Se₄ | `bases --dim 4 4 1 --dim-fc2 4 4 2 --pmat 1 1 2 -1 -2 -1 1 1 1` | 3.7576 / 0.8107 |
| Mn₂In₂Se₅ | `bases --disp-yaml Mn2In2Se5/fc/phono3py_disp.yaml` | 3.0204 / 1.3054 |

(The Mn₂In₂Se₅ value is exactly the committed 24³ isotope-off RTA, which confirms
the pipeline; the MnIn₂Se₄ value differs from the committed 24³ RTA by +1.5 %
because the bases are three-phonon-only while the committed number has isotopes.)

Fits:

| | **MnIn₂Se₄** | **Mn₂In₂Se₅** |
|---|---:|---:|
| rms of ln(model/exp) | **1.2 %** | **3.2 %** |
| intrinsic (no disorder) rms | 87 % | 114 % |
| best x | **0.45** | 0.13 |
| best L | 30.7 nm | 15.9 nm |
| best s | 0.35 | 0.250 |
| acceptable x | 0–0.66 | 0.01–0.50 |
| acceptable L | 5.9–42.7 nm | 5.9–59.5 nm |
| **acceptable s** | **0.330–0.377** | **0.153–0.283** |
| fitted g_Mn / g_In | 0.1324 / 0.0609 | 0.1031 / 0.0354 |

Two things to keep in mind when reading this:

1. **x and L remain degenerate** in both materials; only s is constrained. Report
   the bands, not the single optimum.
2. In **MnIn₂Se₄** the fitted x (0.45) coincides with the single-crystal
   i = 0.44 of Range et al. 1991, and the fitted mass variances (0.1324 / 0.0609) coincide with the VCA
   values (0.1337 / 0.0595) to 1 % / 2.4 %. That is independent support for the
   VCA occupation model — see `../disorder/README.md` §3.4.
3. In **Mn₂In₂Se₅** the isotropic 24³ refit barely moves the earlier z-only result
   (s 0.253 → 0.250, x and L unchanged), so the high-temperature fit was not biased
   by the mesh.

Tracked results: `results/MnIn2Se4_442_m24/` and `results/Mn2In2Se5_m24/` (bases.npz + bases_meta.json + fit JSON). Generated in the scratch directories
(`bases.npz` + `bases_meta.json` + `fit_*.json`). They can be regenerated with the
`bases` commands above (12–15 min on 64 cores each).

## Results (earlier, superseded)

The experimental values used in the fit are digitized from the figure (`experiment/kappa_tot_digitized.csv`); once the raw data are available, simply rerun `fit`.

### Validation: linear combination = direct calculation

With the same set of MnIn₂Se₄ force constants (11³ mesh) and 5 combinations of (x, L) (x = 0.05/0.10/0.20, x=0.10 + L=1 μm, L=0.1 μm),
the κ computed directly by phono3py with `mass_variances` + `boundary_mfp` and the linear-combination result of this script differ by at most **0.01%**
(`results/MnIn2Se4_oldfc_ALASSO5.20_m11/validation_vs_direct_phono3py.txt`).

### MnIn₂Se₄ (old force constants; method demonstration only)

The force constants are the ALASSO 5.20 Å case from the pheasy test dataset; the DFT settings are unknown, so the numbers cannot be taken as conclusions — only the trends matter.
Intrinsic κ(300 K) = 3.97 / 2.12 (xx / zz), (2xx+zz)/3 = 3.35, experiment 0.78, a factor of 4.3. Fit to the 16 high-temperature points:

| scenario | best rms(ln) | meaning |
|---|---|---|
| intrinsic (x=0, no grain boundary, s optimal) | 1.25 | 3.5× too high |
| antisites only, x up to complete disorder 2/3 | 0.76 | still about 2.1× too high |
| any antisites + grain boundary L ≥ 0.5 μm | ≥ 0.34 | still ≥ 40% too high |
| **any antisites + L ≈ 6–36 nm** | **0.01–0.03** | matches in both directions over the whole temperature range |

- **Texture s = 0.33–0.37**, stable across the whole acceptable region (fully random is 0.67); it is the best-determined quantity in this dataset.
- **x cannot be pinned down**: any x from 0.01 to 0.66 can be fitted well with a matching L (larger x goes with larger L); at high temperature the effects of mass disorder and grain boundaries are almost interchangeable.
- The low-temperature range (55–140 K, excluded from the fit) serves as a check: along the acceptable region, model/experiment = 0.98–1.39, closest around x ≈ 0.1–0.2 and L ≈ 10–14 nm.

**Interpretation**: matching the experiment requires pinning the phonon mean free path at **about 10–30 nm**, two orders of magnitude smaller than the hot-pressed grains (micrometre scale),
so the "grain boundary" term actually stands for finer structure: stacking faults / stacking domains, intergrowth of the two layered phases, or spin-phonon scattering of paramagnetic Mn²⁺;
or the intrinsic κ may simply be overestimated (old force constants with unknown settings, cutoff ≤ 5.2 Å). A judgement must wait for force constants of MnIn₂Se₄ and Mn₂In₂Se₅ recalculated with the same protocol.

### Parameter choices (settled)

1. **s**: take it from the fit and report the acceptable region. Once XRD is available, cross-check it independently with the Lotgering factor.
2. **x**: the κ data cannot determine it, so it is fixed by structure. Until a Rietveld occupancy refinement is available, use **x = 0.10** uniformly as the reference value and report the sensitivity at x = 0 and 0.20;
   complete disorder (Mn₂In₂Se₅ 0.5, MnIn₂Se₄ 2/3) serves as the upper bound.
3. **L**: fitted at fixed x and called the "effective mean free path L_eff", not the grain size. When L_eff is much smaller than the SEM grain size,
   there are additional nanoscale scattering sources, and only then is it worth doing dedicated calculations of stacking faults / spin disorder.

### Mn₂In₂Se₅ (provisional: the q-mesh for the intrinsic κ_zz is not yet converged)

The force constants are reliable (S4 equilibrium-frame residual force 0.0009 eV/Å); the intrinsic κ comes from phono3py RTA (the ShengBTE production values are superseded, see `../shengbte_crosscheck/`).
At 300 K the intrinsic (2xx+zz)/3 ≈ 2.3 versus experiment 0.60. The 16 high-temperature points:

| scenario | 16×16×4 | 16×16×8 |
|---|---|---|
| intrinsic (s optimal) rms(ln) | 0.95 | 1.06 |
| antisites up to complete disorder x = 0.5, no grain boundary | 0.36 | 0.57 |
| any x + grain boundary L ≈ 0.5–0.6 μm | ≥ 0.27 | ≥ 0.35 |
| **any x + L_eff ≈ 6–60 nm** | **0.03–0.05** | **0.03–0.05** |
| texture s (acceptable region) | 0.16–0.28 | 0.14–0.28 |

The picture is the same as for MnIn₂Se₄: it only matches if the mean free path is pinned at a few tens of nanometres; x cannot be pinned down (L_eff ≈ 14 nm at x = 0.1 and ≈ 50 nm at x = 0.5).
The final numbers will be fixed once the κ_zz mesh is converged and S6 has been rerun with the correct 3RD to produce an iterative solution.
