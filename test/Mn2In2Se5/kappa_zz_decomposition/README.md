# Mn₂In₂Se₅ — Γ–Z decomposition of the κ_zz mesh scan

**This is the check the report asks for in §10 item 4.** The κ_zz of Mn₂In₂Se₅
*rises* under refinement along c only (16×16×4 → 16×16×16). The question was
whether that rise is the same Γ–Z quadrature artifact that dominates MnIn₂Se₄
κ_zz, or real physics.

## What was run

phono3py RTA, 300 K, **isotopes off** (the setting of the committed reference
scan), fc2/fc3 = `../../Mn2In2Se5/fc/{fc2,fc3}.hdf5` (md5 `b611f139…` /
`ac7bf6b2…`), 4×4×1 supercell, **`primitive_matrix="P"`** and the cell from
`Mn2In2Se5/fc/phono3py_disp.yaml`.

`primitive_matrix="P"` matters: with the `fit_config.json` matrix
`[[2,1,1],[-1,-1,-2],[1,1,1]]` a 16×16×4 mesh breaks the grid symmetry and
phono3py refuses to run. With `"P"` the primitive cell is the 9-atom cell
(identity matrix), which is what `analysis/extrinsic/kappa_extrinsic.py` used
for the committed reference scan.

This reproduces the reference numbers from `Mn2In2Se5/data/kappa_v2/mesh_convergence_reference.json`
and `analysis/shengbte_crosscheck/README.md`:

| mesh | this run κ_zz | committed reference |
|---|---:|---:|
| 16×16×4 | 0.8409 | 0.84092 |
| 16×16×8 | 1.1157 | 1.109 |
| 16×16×12 | 1.3090 | 1.304 |
| 16×16×16 | 1.5472 | 1.541 |

(0.5 % or better everywhere; the residual is phono3py 4.4.0 vs the original run.)

Ran on HanHai, slurm job **410627**, one node / 64 threads; 11–90 s per mesh.
Job scripts in `hanhai_jobs/`, raw per-mode output in `scan/`.

## Result: the rise is **partly** the artifact, and it is not the whole story

`kappa_zz_decompose.py` split, isotropic 16³ included:

| mesh | κ_zz total | first Γ–Z point | rest of the zone | first / total |
|---|---:|---:|---:|---:|
| 16×16×4 | 0.8409 | 0.1899 | 0.6510 | 22.6 % |
| 16×16×8 | 1.1157 | 0.3312 | 0.7845 | 29.7 % |
| 16×16×12 | 1.3090 | 0.4066 | 0.9024 | 31.1 % |
| 16×16×16 | 1.5472 | 0.5382 | 1.0090 | 34.8 % |

Two things happen at once and z-only refinement cannot separate them:

- **The artifact is present and its share grows.** The first Γ–Z point goes
  0.1899 → 0.5382 (×2.8) and its share of the total goes 22.6 % → 34.8 %. That is
  what the quadrature error predicts for refinement along one axis only: the
  weight falls as 1/N_z while the lifetime grows as N_z^p, so the term grows as
  N_z^(p−1). A rising κ_zz under z-only refinement is therefore **not** by itself
  evidence of real long-wavelength physics.
- **And the rest of the zone also rises**, 0.6510 → 1.0090 (×1.55). So z-only
  refinement alone cannot separate the artifact from a genuine under-resolution
  of the interlayer direction.

## The decisive test: isotropic refinement

An **isotropic 24³** mesh (job 410627, 26.8 min) settles it:

| mesh | κ_zz total | first Γ–Z point | rest of the zone | first / total |
|---|---:|---:|---:|---:|
| 16×16×16 | 1.5472 | 0.5382 | 1.0090 | 34.8 % |
| **24×24×24** | **1.3054** | 0.3288 | 0.9766 | 25.2 % |

**Isotropic refinement lowers κ_zz by 15.6 %** (1.5472 → 1.3054), the opposite
direction from the z-only series. That is the signature of the quadrature
artifact, not of real long-wavelength physics: the z-only scan was walking *up*
the artifact, the isotropic scan walks *down* it. In-plane κ goes 2.9946 →
3.0204, i.e. it is converged to +0.9 %.

Mesh-converged estimate on the rest term (two cubic meshes, cell-bias
correction, first-point exponent 1.22):

| quantity | estimate |
|---|---:|
| **intrinsic κ_zz** | **0.926 ± 0.051** |
| κ_zz at a 10 μm grain | 0.912 ± 0.036 |
| κ_zz at a 1 μm grain | 0.786 ± 0.036 |

The rest term is still drifting (1.0090 → 0.9766, −3.2 % per isotropic step) and
there are only two cubic meshes, so **±0.05 is a real uncertainty, not a
formality**: a 32³ mesh is needed before Mn₂In₂Se₅ κ_zz can be quoted as a
number. What can already be said:

- it is near **0.93**, not 1.54;
- the reference values 0.84 → 1.54 from the z-only scan span a factor of 1.8 for
  the *same* force constants, so quoting 1.541 (16×16×16 RTA) would overstate
  κ_zz by ≈65 %;
- the earlier finite-grain claim of ≈0.97–1.0 at 1 μm
  (`analysis/shengbte_crosscheck/README.md`) is too high — on the rest term the
  1 μm value is 0.786 ± 0.036;
- the +4.5 % / +12.4 % coherence ratios at 16×16×4 are measured against a
  population term that carries the artifact, so they must not be quoted either.

**Consequences**

1. The old statement in `REPORT_MnIn2Se4_20261115.md` §6.3 that a rising trend
   proves the Mn₂In₂Se₅ value is physical ("missing long-wavelength modes,
   intrinsic value not well-defined") is **withdrawn**.
2. The MnIn₂Se₄ result must not be transplanted verbatim either: for MnIn₂Se₄ the
   rest term is flat (0.6461 → 0.6473, 18³ → 30³) and κ_zz is settled on the
   mesh; for Mn₂In₂Se₅ the rest term is still falling at 24³ (−3.2 % per step),
   so one more isotropic mesh is needed.
3. Both materials follow the same rule: **quote the rest of the zone, never the
   total, and never a z-only mesh scan.**

## Reproduce

```bash
cd test/Mn2In2Se5/kappa_zz_decomposition
python ../../MnIn2Se4/kappa_zz_decompose.py "scan/*/kappa-m*.hdf5" \
       --poscar ../fit/POSCAR --pmat 1 0 0 0 1 0 0 0 1 \
       --json kappa_zz_mesh_decomposition.json
```

This reproduces `kappa_zz_mesh_decomposition.json` byte for byte.
The estimate uses the two cubic meshes (16³, 24³); the z-only cases only appear
in the decomposition tables, which is why `--no-estimate` must **not** be used
any more.

`primitive_matrix="P"` in phono3py resolves to the **identity** for this cell,
which is why `--pmat 1 0 0 0 1 0 0 0 1` is the right matrix here. `../fit/POSCAR`
is the same 9-atom cell as `Mn2In2Se5/fc/phono3py_disp.yaml`.
