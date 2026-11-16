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
- **But the rest of the zone also rises**, 0.6510 → 1.0090 (×1.55), and it is not
  flat. So the artifact does not explain the whole increase either: the
  interlayer direction is genuinely under-resolved in this material, and its
  κ_zz does not converge under z-only refinement.

**Consequences**

1. The old statement in `REPORT_MnIn2Se4_20261115.md` §6.3 that a rising trend
   proves the Mn₂In₂Se₅ value is physical ("missing long-wavelength modes,
   intrinsic value not well-defined") is **withdrawn**: the artifact contributes.
2. Equally, the MnIn₂Se₄ conclusion cannot be transplanted as "the same
   artifact" — here the rest term does *not* converge, whereas for MnIn₂Se₄ it
   does (0.6461 → 0.6473 from 18³ to 30³).
3. **Mn₂In₂Se₅ κ_zz is not converged and no value of it should be quoted yet**,
   including the +4.5 % / +12.4 % coherence ratios from the RTA runs at 16×16×4,
   whose population term is itself mesh-dependent. Settling it needs isotropic
   meshes at least to 24³ (one 24³ run is queued) plus a finite-grain analysis on
   the rest term, in the same way as §6.4 of the MnIn₂Se₄ report.

## Reproduce

```bash
cd test/Mn2In2Se5/kappa_zz_decomposition
python ../../MnIn2Se4/kappa_zz_decompose.py "scan/*/kappa-m*.hdf5" \
       --poscar ../fit/POSCAR --pmat 1 0 0 0 1 0 0 0 1 \
       --no-estimate --json kappa_zz_mesh_decomposition.json
```

This reproduces `kappa_zz_mesh_decomposition.json` byte for byte.
`--no-estimate` is needed because only **one** cubic mesh (16³) has been run so
far and the estimator needs at least two; add the queued 24³ and drop the flag.

`primitive_matrix="P"` in phono3py resolves to the **identity** for this cell,
which is why `--pmat 1 0 0 0 1 0 0 0 1` is the right matrix here. `../fit/POSCAR`
is the same 9-atom cell as `Mn2In2Se5/fc/phono3py_disp.yaml`.
