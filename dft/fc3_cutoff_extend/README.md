# Pushing the FC3 cutoff past 6.45 / 6.51 Å — first without new DFT

**Status 2026-11-19: scheduled** (`protocol/local_AI_next_steps_2026-11-19.md`, task 2), steps 0 and 1 only. The
2026-11-18 decision not to pursue it is revised because the κ_zz ordering of the two materials depends on it
(report §1, "Which material is lower").

## The limit of the existing 4×4×1 cells is ≈ 8 Å, not 6.5 Å

The 6.51 / 6.45 Å ceiling quoted in the report is AutoZT's `cut3_safe_cutoff`:
half the shortest **interplanar distance of the supercell parallelepiped**,
minus 0.1 Å. For a non-orthogonal cell this depends on the basis chosen. It is
only needed by codes that find periodic images by rounding fractional
coordinates in that basis.

The basis-independent condition is that no atom pair has two periodic images
inside the cutoff, i.e. the cutoff is below half the **shortest supercell
lattice vector**. Checked pair by pair on the Mn₂In₂Se₅ 4×4×1 supercell:

| cutoff | pairs | pairs with two images inside the cutoff (real aliasing) | pairs where fractional rounding in the skewed basis picks the wrong image |
|---|---:|---:|---:|
| 6.61 Å | 3248 | 0 | 0 |
| 7.00 Å | 4320 | 0 | 144 |
| 7.96 Å | 5552 | **0** | 208 |
| 8.50 Å | 7048 | 504 | 468 |

| supercell | parallelepiped rule | minimum-image rule |
|---|---:|---:|
| Mn₂In₂Se₅ 4×4×1 (existing) | 6.61 Å | **7.96 Å** |
| Mn₂In₂Se₅ 4×4×2 | 6.61 Å | 7.96 Å |
| Mn₂In₂Se₅ 5×5×1 (rebased) | 8.02 Å | 9.26 Å |
| MnIn₂Se₄ 4×4×1 (existing) | 6.51 Å | **8.01 Å** |
| MnIn₂Se₄ 4×4×2 | 6.55 Å | 8.05 Å |

So **4×4×2 gains nothing for either material under either rule**. The
existing 4×4×1 frames can support a fit up to ≈ 8 Å, provided the fitting code
uses true minimum-image vectors. The report's "4×4×2 safe limit 8.05 Å" used
the minimum-image rule, while the 4×4×1 limit beside it used the parallelepiped
rule, so the two were not comparable.

## Step 0 — does the fitting code use true minimum-image vectors? (minutes)

This must be settled before trusting any fit above 6.6 Å. In pheasy (the
fit used for both materials: RFE-OLS for MnIn₂Se₄, ALASSO for Mn₂In₂Se₅ — see each `test/*/fit/fit_config.json`), find how pair/cluster vectors are formed
from the supercell:

- **Shortest image over neighbouring cells**, as in phonopy's
  `get_smallest_vectors` or ASE neighbour lists: fine up to ≈ 7.9 Å.
- **Fractional coordinates rounded in the cell basis**, `d − round(d)`: wrong
  for 144–208 pairs between 7.0 and 7.96 Å in this cell. Stop and ask.

An empirical check that needs no code reading: for the 4×4×1 SPOSCAR, compare
the list of pair distances pheasy reports up to 7.9 Å with phonopy's smallest
vectors. They must agree.

## Step 1 — refit the existing frames at larger cutoffs (no DFT)

Both materials, existing 4×4×1 datasets:

- Mn₂In₂Se₅: 55 frames.
- MnIn₂Se₄: 46 frames, FC3 only. Its FC2 stays the 4×4×2 one.

1. pheasy fit with the same method and settings as before (MnIn₂Se₄ RFE-OLS, Mn₂In₂Se₅ ALASSO), FC3 cutoff at the next shells:
   - Mn₂In₂Se₅: 6.67 / 6.90 / 7.11 / 7.43 Å (the shell list in
     `Mn2In2Se5/data/kappa/cutoff_scan_S6.json` goes to 7.43);
   - MnIn₂Se₄: the corresponding shells up to ≈ 7.9 Å.
2. For each cutoff, record the data/parameter ratio and the train/CV error.
   AutoZT requires a ratio ≥ 3. The Mn₂In₂Se₅ ratio was 12.3 / 9.1 / 6.1 at
   5.68 / 6.12 / 6.51 Å, and a rough extrapolation gives ≈ 4 at 7.0 Å and
   ≈ 2–2.5 at 7.5–7.9 Å.
3. κ: phono3py RTA, isotopes on. Use 16³ for Mn₂In₂Se₅ (`primitive_matrix="P"`)
   and 14³ for MnIn₂Se₄ (4×4×2 FC2). Split κ_zz with
   `test/MnIn2Se4/kappa_zz_decompose.py` and compare the rest term with the
   existing values (Mn₂In₂Se₅ 0.9850, MnIn₂Se₄ 0.6444 at their last cutoffs).
4. Decide:
   - **ratio ≥ 3 up to 7.4–7.9 Å and the κ_zz steps shrink below ≈ 5 %**:
     cutoff converged, done, no DFT.
   - **ratio < 3 above ≈ 7 Å**: go to step 2.
   - **κ_zz keeps falling ≈ 10 % per step up to 7.9 Å**: not converging in
     reach. Only then is a larger cell (5×5×1, 9.26 Å) worth discussing.

## Step 2 — only if needed: more frames in the same 4×4×1 cell

`make_frames.py` writes frames that can be appended to the existing datasets.
The supercell, atom order and positions are identical to the existing SPOSCARs
(checked to 1e-15). The displacement statistics (Gaussian, σ = 0.017 Å) and the
S4 INCAR/KPOINTS are also the same.

```bash
cd dft/fc3_cutoff_extend
python make_frames.py --material Mn2In2Se5 --out new_frames --first 55 --count 30
cp <the 4x4x1 S4 POTCAR> new_frames/POTCAR
bash submit_frames.sh new_frames 55 56      # two frames first, check them
python check_frames.py new_frames
bash submit_frames.sh new_frames 57 84
```

Then copy the new `disp-*` into the existing `step4_disp/` and refit (step 1).
For MnIn₂Se₄ use `--material MnIn2Se4 --first <last existing + 1>`.

Cost: 144-atom (Mn₂In₂Se₅) / 112-atom (MnIn₂Se₄) frames at the S4 settings,
≈ 2–3 h each on jzzn with 48 cores (estimated from the 224-atom 4×4×2 frames at
≈ 4 h). 30 frames ≈ 4–6 k core-hours, ≈ 1 day of wall time with a few in
parallel. To reach a ratio of ≈ 3 at 7.9 Å, Mn₂In₂Se₅ needs roughly 85–90 frames
in total, i.e. ≈ 30–35 more.

## Files

| file | what |
|---|---|
| `make_frames.py` | extra frames for either 4×4×1 dataset; prints both cutoff limits for any `--dim`; never overwrites |
| `submit_frames.sh` | one jzzn job per frame (cpu192, 48 cores, regular QOS, the modules and MPI settings of the earlier S4 runs) |
| `check_frames.py` | completion, SCF convergence, max force, force drift per frame |
