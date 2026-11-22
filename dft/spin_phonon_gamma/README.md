# Spin–lattice coupling screen for Mn₂In₂Se₅: FM vs AFM Γ phonons

Every κ in this repository assumes ferromagnetic (FM) Mn order. The samples are
paramagnetic at 300 K, with strong antiferromagnetic exchange (Weiss θ ≈ −198 K,
Williams 2024). If the force constants depend on the spin arrangement, the
paramagnetic crystal has random force constants, which add scattering, and its
mean phonons differ from the FM ones. This test measures that dependence where it
is cheapest, at Γ of the 9-atom primitive cell. It is item 3 of the "Why the
calculation lies above the experiment" table in the report.

Instructions for running it: `../../protocol/local_AI_spin_phonon_2026-11-22.md`.

## What is computed

| set | cell | spins (Mn1, Mn2 = the two Mn layers of a slab) | runs | gives |
|---|---|---|---:|---|
| `FM/` | 9 atoms | ↑ ↑ (production) | 1 + 54 | Γ phonons, energy, moments |
| `AFM/` | 9 atoms | ↑ ↓ (layers antiparallel, each layer FM) | 1 + 54 | Γ phonons, energy, moments |
| `energy_2x1/` | 18 atoms (2a₁, a₂, a₃) | FM; stripe in each layer, two interlayer variants | 3 | in-plane exchange scale |

- Geometry: the production FM structure (`Mn2In2Se5/POSCAR_relaxed`) for every spin
  configuration. The test therefore isolates how the force constants depend on the
  spins at fixed geometry. The AFM residual force is reported as the magnetostriction
  signal.
- DFT: the production S2 settings (PBE + D3(BJ), U(Mn) = 4 eV, ENCUT 410, Γ 8×8×2), plus
  ISYM = 0 and ADDGRID as in the S4 force runs. EDIFF is 1E-8. `prepare.py` checks the
  shared tags against `Mn2In2Se5/data/dft_inputs/step2_static/INCAR`.
- Displacements: ±0.01 Å along a₁, a₂, a₃ for every atom, with no symmetry reduction.
  This avoids phonopy's magnetic symmetry search, which fails for this cell
  (spglib error, phonopy 4.4.0). It also treats FM and AFM identically, and the
  splitting of the FM E doublets then measures the numerical noise.

## Analysis and decision rule (`analyse.py`, fixed before the runs)

- Γ force constants: forces minus the ref forces, central differences, translational
  invariance imposed.
- Each FM optical mode i gets the overlap-weighted AFM frequency
  ω̃ᵢ = Σⱼ |⟨eᵢ(FM)|eⱼ(AFM)⟩|² ωⱼ, and δᵢ = ω̃ᵢ/ωᵢ − 1.
- Checks: every run finished below NELM; Mn moments |m| > 4 μB with the intended signs;
  FM against the production Γ frequencies from `Mn2In2Se5/fc/fc2.hdf5` (expected
  within ≈ 1–2 %).

| verdict | condition | meaning / next step |
|---|---|---|
| weak | max\|δ\| ≤ max(2 %, 3 × noise) and rms ≤ 1 % | spin-dependent force constants negligible for κ; item 3 can be closed for this channel |
| intermediate | otherwise | Γ phonons of the stripe configuration (18 atoms) before deciding |
| strong | max\|δ\| ≥ 5 % or rms ≥ 3 % | paramagnetic force-constant disorder matters: spin-disordered supercells |

If the in-plane stripe energy is more than 2× the interlayer AFM energy (per Mn),
`analyse.py` warns. The 9-atom AFM then probes the weaker exchange path, and a
"weak" verdict is not conclusive.

Why 2 %: a frequency shift δ corresponds to a force-constant change of about 2δ, and
disorder scattering scales with its square. A 2 % shift gives (0.04)² ≈ 0.002, compared
with the VCA mass-variance g = 0.05–0.13. The mean-field effect on κ (through the
velocities) is of order 2–3δ, i.e. ≲ 5 %.

## Tested here (2026-11-22, no VASP)

- `analyse.py --selftest` builds forces from the production Γ force constants. The FM
  set reproduces the production Γ frequencies to 1×10⁻⁸ THz. An AFM-like perturbation
  (a stiffer Mn2–Se bond) matches direct diagonalisation to 2×10⁻⁷ THz.
- `prepare.py` followed by fabricated OUTCARs was run through `analyse.py` end to end:
  parsing, residual-force subtraction, noise estimate, energies and verdict.

## Limits

- Γ only. Phonons are computed for one AFM pattern (interlayer); the in-plane exchange
  enters only through the energies.
- Static force constants at fixed geometry. Dynamic spin–phonon (magnon) scattering is
  not tested.
- U = 4 eV is fixed. The exchange energies depend on U.

## Files

| file | what |
|---|---|
| `prepare.py` | writes the run directories, INCAR/KPOINTS/POSCAR, `displacements.json`, `manifest.json` |
| `run_config.sh` | Slurm job (jzzn cpu192, settings of `../fc3_cutoff_extend/submit_frames.sh`), stages `ref` / `disp` / `all` |
| `analyse.py` | the analysis and the decision rule; `--selftest` checks it without VASP |
