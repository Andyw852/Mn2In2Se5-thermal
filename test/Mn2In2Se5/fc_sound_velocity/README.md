# Mn₂In₂Se₅ — is the 4×4×1 FC2 too soft between the layers?

**Question.** Mn₂In₂Se₅'s FC2 and FC3 both come from a 4×4×1 supercell, and the
cell is only 16.24 Å along the layer normal — one layer thick. That is exactly
the defect that made MnIn₂Se₄'s old 4×4×1 FC2 4.2× too soft between the layers
and suppressed its κ_zz by ≈3.5×. If Mn₂In₂Se₅ has the same problem, every κ_zz
it has (including the 24³ run) is wrong and a denser q-mesh would be wasted.

**Answer: it does not.** The interlayer transverse sound velocity from the 4×4×1
FC2 is **1335–1387 m/s**, i.e. C₄₄ = **9.5–10.3 GPa** — the value a layered
selenide should have, and 4× above MnIn₂Se₄'s broken 4×4×1.

## Method

Phonopy reads the FC2 and the layer-normal sound velocity is fitted from small-q
frequencies as f² = f0² + (v q / 2π / 100)², with the Γ-point frequency offset
f0 removed (the two pitfalls documented in
`../MnIn2Se4/fc2_442_phonon/sound_velocity.py`). The layer normal is derived
from a₁×a₂ of the unit cell and converted into the reduced coordinates of the
*primitive* reciprocal basis, so the script works for any primitive matrix. The
two transverse branches are identified by their degeneracy for q ∥ c.

Pure post-processing — no DFT, no new force constants.

```bash
cd test/Mn2In2Se5/fc_sound_velocity
./run_fc_sound_velocity.sh          # ~4 s per case
```

`fc_sound_velocity.json` holds the numbers; the per-case tables are printed to
stdout.

## Validation ladder

The same script applied to the two MnIn₂Se₄ FC2s reproduces the numbers that are
already established for that material (DFT elastic constants: TA 1338.3 m/s,
C₄₄ 9.402 GPa):

| FC2 | interlayer TA | C₄₄ from q ∥ c | vs DFT C₄₄ | verdict |
|---|---:|---:|---:|---|
| MnIn₂Se₄ 4×4×1 (the broken one) | 323.2 m/s | 0.55 GPa | **−94.2 %** | reproduces the known 4.2× softening |
| MnIn₂Se₄ 4×4×2 (the corrected one) | 1348.1 m/s | 9.54 GPa | **+1.5 %** | reproduces the published result |
| **Mn₂In₂Se₅ 4×4×1 (`Mn2In2Se5/fc/`, the set the κ runs used)** | **1386.8 m/s** | **10.27 GPa** | — | **normal** |
| Mn₂In₂Se₅ 4×4×1 (`test/Mn2In2Se5/fit/`, the other revision) | 1334.5 m/s | 9.51 GPa | — | normal |

So the script lands on both established MnIn₂Se₄ values, and Mn₂In₂Se₅'s
one-layer FC2 comes out where a layered selenide should — not in the collapsed
regime.

For reference, the whole elastic picture from the `fc/` set: C₃₃ = 66.1 GPa,
C₄₄ = 10.3 GPa, in-plane acoustic velocities 1335.5 / 1798.3 / 3821.2 m/s.

## Caveats

1. **There is no DFT elastic-constant calculation for Mn₂In₂Se₅**, so the
   comparison above is indirect: it rests on the MnIn₂Se₄ ladder and on the
   magnitude being right for this chemistry. If a quantitative reference is
   wanted, a C₄₄ strain calculation in the primitive cell is small (a few static
   points at ε_xz = ±0.005, ±0.010, C₄₄ = 2ΔE/(V ε²)).
2. **The fitted FC2 breaks the hexagonal C₄₄ = C₄₄' identity.** With q ∥ c the
   `fc/` set gives C₄₄ = 10.27 GPa, while the c-polarised branch at q in the
   layer plane gives 8.48 GPa (−17 %). The `fit/` revision is consistent
   (9.51 vs 9.52 GPa, +0.1 %). So quote C₄₄ ≈ 9.5–10.3 GPa, not a third digit,
   and prefer `fit/` if symmetry-consistency matters. The conclusion — "not
   collapsed" — is unaffected either way.
3. This checks **FC2 only**. Whether the 4×4×1 *FC3* is adequate is a separate
   question, and its cutoff ceiling is 6.61 Å.

## Consequence

The κ_zz mesh scans for Mn₂In₂Se₅ are meaningful, so isotropic **28³ and 32³**
are worth running to see whether the rest term flattens (it is 0.9766 at 24³ and
still falling by 3.2 % per step). No 4×4×2 FC2 is needed for this material.
