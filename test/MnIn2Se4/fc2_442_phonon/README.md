# MnIn2Se4 fc2 (4×4×2 symmetry-reduced finite displacements) — 2026-11-13

The second-order force constants were recomputed in a **4×4×2** supercell instead of 4×4×1, in order to fix the C44 interlayer shear modulus problem
(with 4×4×1 + BHH the TA sound velocity along z is only 317 m/s under the old convention, whereas the DFT elastic constants give 1338 m/s).

## Calculation settings

| item | value |
|---|---|
| method | phonopy symmetry-reduced finite displacements (**not** random displacements / alm) |
| supercell | **4×4×2** (224 atoms: 32 Mn + 64 In + 128 Se) |
| symmetry | `strict` (384 symmetry operations, keeping the full R-3m point group) |
| displacement amplitude | **FD_DISTANCE = 0.02 Å** |
| number of frames | **7** (minimal complete set of symmetry-inequivalent displacements) |
| functional | PBE-D3, U(Mn)=4 eV, FM (MAGMOM = 32\*5 192\*0) |
| EDIFF | **1E-7** |
| KPOINTS | 2×2×1 (KSPACING 0.038 rule) |
| NAC | off (semimetal; the DFPT dielectric response diverges; consistent with the NAC=off convention of the kl chain) |
| cluster | jzzn (cpu192, premium, 64 cores/frame) |

## The 7 displacement frames (after phonopy reduction)

All displacements have the direction `(0, 0.011600, 0.016292)` (norm = 0.02 Å); only the sign and the site differ:

| frame | atom | element | sign |
|---|---|---|---|
| disp-00001 | 1 | Mn | + |
| disp-00002 | 33 | In | + |
| disp-00003 | 33 | In | − |
| disp-00004 | 97 | Se | + |
| disp-00005 | 97 | Se | − |
| disp-00006 | 161 | Se | + |
| disp-00007 | 161 | Se | − |

## Results: the phonon spectrum is stable

```
PHONON_DONE = true
stable      = true
verdict     = pass
min_freq    = +0.000719 THz   (no imaginary frequency, n_neg_qpoints = 0)
n_qpoint    = 499, n_branch = 21
max_freq    = 7.784 THz
```

For comparison, the **4×4×1** of the kl chain has min_freq = −0.006 THz (a soft mode / tiny imaginary frequency); the 4×4×2 spectrum is cleaner.

Wall time is about 4 hours per frame (jzzn, 64 cores); a control frame on hanhai25 with 128 cores took 7h46m.

## ✅ The layer-normal sound velocity agrees with the DFT elastic constants

**Run `sound_velocity.py` (in this directory) to reproduce the table below**:

| quantity | 4×4×2 fc2 | DFT elastic tensor | deviation |
|---|---:|---:|---:|
| transverse TA₁ | **1348.3 m/s** | 1338.3 m/s | **+0.75%** |
| transverse TA₂ | **1348.3 m/s** | 1338.3 m/s | **+0.75%** |
| longitudinal LA | **3118.6 m/s** | 3253.4 m/s | −4.14% |
| C₄₄ = ρv_T² | **9.543 GPa** | 9.402 GPa | +1.50% |

The result is stable for fit windows t ≤ 0.005 / 0.010 / 0.020 / 0.040 (TA varies by < 0.3 m/s).
The density is uniformly taken as 5249.3 kg/m³ from the phonon structure.

### ⚠️ Two pitfalls when extracting the sound velocity (both handled in this script)

**Pitfall 1: phonopy's q uses the reciprocal-space coordinates of the "chosen primitive cell", not those of the POSCAR primitive cell.**

Here `primitive_matrix="auto"` replaced the primitive cell with the **rhombohedral primitive cell**; in the new cell
b₁/b₂/b₃ all make an angle of **84.91°** with the layer normal (i.e. they lie almost in the layer plane).

a₁ and a₂ of the POSCAR lie in the layer plane → layer normal = `(0,0,1)`_POSCAR.
The same interlayer direction in the two coordinate systems:

```
(0,0,t)_POSCAR  ≡  (t,t,t)_phonopy_primitive
```

**Wrong approach**: compute frequencies with `(0,0,t)` and then convert |q| using the length of b₃ from the POSCAR →
this gives spurious values of **4179 / 6256 m/s**. The ratio of the two |b₃| values, **3.760**, is exactly the error factor.

**Pitfall 2: there is a numerical residue at the Γ point (f₀ = 0.000719 THz).**

When ω/q is computed directly, low q is contaminated by this residue and actually **diverges** as q→0 (t=0.002 gives 1429,
t=0.010 gives 1351). The correct approach is to fit

```
f² = f₀² + (v·q / 2π / 100)²
```

taking the square root of the slope. This script does exactly that.

> Both points above were found by the user in an independent audit on 2026-11-13 (`_fc2_442_direction_audit_20261113/`),
> and this script independently reproduces all of their numbers.

## ✅ 4×4×1 vs 4×4×2: the overly soft C44 problem is fixed

Both were re-checked with the **same** phonopy primitive-cell reciprocal-space convention (script `../compare_441_vs_442.py`,
results `../compare_441_vs_442_results.json`):

| supercell | TA₁ | TA₂ | LA | C₄₄ (GPa) | TA deviation | LA deviation |
|---|---:|---:|---:|---:|---:|---:|
| **4×4×1** | **322.9** | 322.9 | 2491.5 | 0.547 | **−75.87%** | −23.42% |
| **4×4×2** | **1348.3** | 1348.3 | 3118.5 | 9.543 | **+0.75%** | −4.15% |
| DFT elastic | 1338.3 | 1338.3 | 3253.4 | 9.402 | — | — |

**Conclusions**:
1. **The interlayer TA of 4×4×1 really is only ~323 m/s** — the previously reported 317 is essentially correct and is **not** the spurious value
   caused by mixing basis sets (although the 4179 from the 4×4×2 run was).
2. **4×4×2 rises to 1348 m/s, only 0.75% away from the DFT value of 1338** → the interlayer long-range force constants along c had indeed been
   truncated by a 1-unit-cell thickness; enlarging to 2 unit cells **fixes the overly soft C44**.

## Open issues

**The longitudinal LA differs by about 4.1%** (3118.6 vs 3253.4). The structures used in the two calculations differ slightly:
the layer-normal lattice component is **13.21409 Å** (this phonon structure) and **13.20318 Å** (elastic calculation).
The structure and the elastic-calculation setup should first be unified (same CONTCAR, same ENCUT/functional/k-mesh),
and only then can one decide whether the 4.1% comes from the structural difference, supercell convergence, or the elastic-constant fit itself.

> Note: the LA deviation of 4×4×1 (−23.4%) is much larger than that of 4×4×2 (−4.1%), which shows that the convergence direction along c is correct,
> but whether 2 unit cells are already sufficient can only be settled by comparison with DFT elastic constants computed under a unified setup.

## Files
## Files

| file | description |
|---|---|
| `fc2.hdf5` | second-order force constants (phonopy format) |
| `FORCE_CONSTANTS` | as above, phonopy text format (for phono3py / ShengBTE) |
| `band-dft-cpu.yaml` | phonon spectrum (499 q-points × 21 branches), **the plotting data source** |
| `phonon_summary.json` | stability-criterion results (imaginary-frequency gate / ZA quadratic behaviour) |
| `phonon_fit.log` | fit log |
| `phonopy_disp.yaml` | displacement definition (including the 3×3 supercell matrix and primitive_matrix) |
| `disp_plan.json` | displacement plan (method / fd_distance / n_disp) |
| `POSCAR` | primitive-cell structure |
| `kl_params.txt` | snapshot of the calculation parameters |
| **`sound_velocity.py`** | **layer-normal sound-velocity extraction (handles both pitfalls, with comments)** |
| `phonon_band_plot/phonon_band_full.png` | full phonon spectrum |
| `phonon_band_plot/phonon_band_lowfreq.png` | low-frequency zoom (to see the acoustic branches) |
| `phonon_band_plot/phonon_band_summary.json` | plotting-side summary |
| `frames/disp-*/ ` | inputs for the 7 frames (POSCAR / INCAR / KPOINTS / CONTCAR) |

> OUTCAR / vasprun.xml / queue.out are not committed, per the repository .gitignore (the raw DFT output is too large; it can be re-fetched from the cluster).
> If force data are needed, rerun S2_disp or fetch `step2_disp/disp-*/OUTCAR` from jzzn.

## Reproduction

```bash
python sound_velocity.py .                    # 4×4×2 layer-normal sound velocity vs DFT
cd .. && python compare_441_vs_442.py         # 4×4×1 vs 4×4×2 comparison
```
