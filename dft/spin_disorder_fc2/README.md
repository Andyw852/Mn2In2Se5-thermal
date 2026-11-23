# Paramagnetic (spin-disordered) FC2 of Mn₂In₂Se₅

Follow-up to `../spin_phonon_gamma/`. That screen found **strong** spin–lattice
coupling (verdict by the pre-set rule): at Γ, interlayer AFM order shifts the FM phonons
by up to 8.1 %, rms 4.3 % below 3 THz, noise 0.36 %. The energies also show that FM,
which every κ in this repository assumes, is not the ground state. Interlayer AFM lies
10.4 meV/Mn lower. The samples are paramagnetic over the whole measured range. This
directory computes the paramagnetic harmonic force constants and passes them to the
existing κ machinery.

Instructions for running it: `../../protocol/local_AI_spin_disorder_fc2_2026-11-23.md`.

## 1. What is known (literature and our data)

| quantity | value | source |
|---|---|---|
| Curie–Weiss | θ_W = −198 K (fit above 50 K), μ_eff 6.0 μB/Mn (Mn²⁺ d⁵ high spin) | Williams *et al.*, *J. Mater. Chem. C* **12**, 1753 (2024), Fig. 5 (`references/`) |
| low-T state | spin-glass freezing at T_f ≈ 7 K (AC χ, Mydosh parameter 0.0066); frustration index > 20 | same |
| exchange paths | 180° Mn–Se–Mn between the two Mn layers of a slab: strongly AFM; 90° paths: weak, FM or AFM; triangular in-plane lattice → frustration | same, Goodenough–Kanamori discussion, Fig. 6 |
| DFT in that paper | AFM "spins alternating across the layers" (PBE+U 5 eV, D2) | same, Methods |
| our 9/18-atom energies | AFM(interlayer) − FM = −10.44, stripe A −7.99, stripe B −3.86 meV/Mn | `../spin_phonon_gamma/results/` |
| our exchange fit (E = Σ J sᵢsⱼ, S² absorbed, J > 0 AFM) | J⊥(3.89 Å, 3 nbrs) = 1.58; J⊥(5.60 Å, 3 nbrs, the 180° path) = 1.90; J∥(4.03 Å, 6 nbrs) = 0.18 meV | 3 energies, 3 parameters, no redundancy |
| implied θ_CW (classical, S = 5/2) | ≈ −60 K | −(S+1)/S · Σ zJ / 3k_B |

Consistency with the literature: the strongest coupling is the 180° interlayer path, and
the ground state is interlayer AFM. Both agree with Williams 2024. The magnitude does
not agree. The DFT exchange gives |θ| about 3× smaller than the measured 198 K. A Weiss
fit starting at 50 K ≈ |θ|/4 tends to overestimate |θ|, and U = 4 eV is fixed here. So
the experimental temperatures are T/|θ_W| ≈ 1.5 (300 K) to 3 (623 K). That is
paramagnetic, but with short-range AFM correlations that are not negligible near
300 K (see §4).

## 2. How the literature treats a paramagnet with spin–lattice coupling

From memory; check volume and pages before citing.

| approach | idea | examples |
|---|---|---|
| ordered proxy | compute FM or AFM and call it the material | common in thermal-transport papers. Not adequate when T ≫ T_N and the phonons depend on the order, which is our case |
| disordered local moments (DLM) | paramagnet = random, uncorrelated local moments | Gyorffy *et al.*, *J. Phys. F* **15**, 1337 (1985) |
| magnetic SQS / spin-space averaging (**this work**) | static random collinear spins in a supercell whose pair correlations vanish shell by shell. Forces are averaged over spin configurations, giving the effective paramagnetic FC2 | Alling, Marten, Abrikosov, *PRB* **82**, 184430 (2010) (CrN, magnetic SQS); Körmann *et al.*, *PRB* **85**, 125104 (2012) (spin-space averaging, paramagnetic bcc Fe phonons); Ikeda, Seko, Togo, Tanaka, *PRB* **90**, 134106 (2014) (bcc Fe, SQS + phonopy); SQS: Zunger *et al.*, *PRL* **65**, 353 (1990) |
| DLM-MD / TDEP | spins re-randomised during AIMD; effective temperature-dependent FC2 (and FC3) fitted from the trajectory | Steneteg, Alling, Abrikosov, *PRB* **85**, 144404 (2012); Shulumba *et al.*, *PRB* **89**, 174108 (2014) (CrN); TDEP: Hellman, Abrikosov, Simak, *PRB* **84**, 180301 (2011) |
| coupled spin + lattice dynamics | spin dynamics with AIMD; gives the extra phonon damping by spin fluctuations | Stockem *et al.*, *PRL* **121**, 125902 (2018) (paramagnetic CrN: anomalous lifetime shortening) |
| order-dependent phonons (ordered phase) | FM vs AFM force constants differ; phonon splitting below T_N | Fennie & Rabe, *PRL* **96**, 205505 (2006) (ZnCr₂O₄); Massidda *et al.*, *PRL* **82**, 430 (1999) (MnO) |
| experiment / phenomenology | extra paramagnetic spin–phonon scattering term in a Callaway model | Slack, *Phys. Rev.* **122**, 1451 (1961) (MnF₂ vs ZnF₂: paramagnetic MnF₂ has the lower κ); Slack & Newman, *PRL* **1**, 359 (1958) (MnO, NiO) |
| fitting algorithm used here | symmetry-adapted least squares for FC from random displacements | Seko & Togo, *PRB* **110**, 214302 (2024) (symfc) |

The usual hierarchy:
1. Mean effect. Replace the FM FC2 by the spin-averaged paramagnetic FC2 (SQS,
   T → ∞ limit of the correlations). This changes frequencies, velocities and the phase
   space for three-phonon scattering.
2. Fluctuation effect. The spread of the FC between spin configurations acts as
   force-constant disorder. It scatters phonons roughly like point defects and is
   nearly independent of T above |θ|.
3. Correlations. Near |θ|, short-range order is included through correlated spin
   configurations or spin–lattice dynamics.

This directory does step 1 fully and step 2 as a measured number (the excess fit
residual). Step 3 is bracketed by the FM and AFM bounds.

## 3. What is computed

Cell: the production 4×4×1 supercell (144 atoms, 32 Mn, atom order of `Mn2In2Se5/fc/SPOSCAR`)
at the production FM-relaxed geometry. DFT is the production S4 force run (PBE + D3(BJ),
U = 4 eV, ENCUT 410, Γ 2×2×2, ISYM = 0, ADDGRID), with EDIFF 1E-8 and NELM 300.
`prepare.py` checks every shared tag against the production INCAR.

| set | spins | runs | role |
|---|---|---:|---|
| `FM_ctl` | all ↑ (production) | 1 + 8 + 1 | method control against production fc2; frame 9 = −frame 1 checks even-order (anharmonic) contamination |
| `AFM` | Mn1 layer ↑, Mn2 layer ↓ | 1 + 8 | ordered bound (strongest correlations); cross-check with the 9-atom Γ result |
| `PM01`…`PM08` | 8 ↑ / 8 ↓ per layer; ⟨sᵢsⱼ⟩ = 0 **exactly** in all five shells < 7.2 Å (3.89, 5.60, 6.90 Å interlayer; 4.03, 6.98 Å intralayer) | 8 × (1 + 2) | paramagnet (uncorrelated limit) |
| total | | **43** | 144-atom single points |

- The spin configurations (`spin_configs.json`) and displacements (`frames.npz`) are
  committed. They were generated once with fixed seeds (`prepare.py --regenerate`). No
  two PM configurations are related by any of the 192 space-group operations of the
  supercell or by a global flip.
- Displacements: every atom moves by 0.01 Å in a random direction. A random frame
  constrains the FC2 much more than a single-atom displacement does, and every frame
  samples all local spin environments at once. That is why the spin average works with
  few frames.
- Each displaced frame restarts from its set's ref WAVECAR, so it stays in the same
  spin state. Forces are taken relative to the set's own ref. This removes the static
  force of the spin arrangement, up to ~0.05 eV/Å in the 9-atom AFM.

**Fit** (`analyse.py`): phonopy + symfc with the full R-3m symmetry and the acoustic
sum rule.
- `FM_ctl` and `AFM` are fitted on their own frames.
- `PM` is one pooled fit over all 16 PM frames, which is the spin-space average.
  Symmetrising a fit over random spin configurations is the same averaging as in the
  SQS / DLM-MD papers above.
- Jackknife: 8 refits, each leaving out one spin configuration, give the error from
  spin sampling.

AFM keeps the full R-3m symmetry for the force constants. Every operation either
preserves all spins or flips all of them (time reversal), and force constants are even
under time reversal. The 9-atom AFM doublets came out degenerate, as they should.

**Tested here without VASP** (`analyse.py --selftest`; `prepare.py` → fake OUTCARs →
full analysis):
- Synthetic model: production FC2 plus spin-dependent Mn–Mn springs, λ = 0.08, 0.03,
  0.05 eV/Å² on the 3.89, 4.03 and 5.60 Å shells. With these values the synthetic
  AFM−FM Γ shifts (max 11 %, rms below 3 THz 3 %) resemble the measured ones.
- Also included: random static forces of 0.02 eV/Å in every non-FM ref and 5×10⁻⁵ eV/Å
  force noise.

| check | result |
|---|---|
| FM_ctl fit vs exact FC2, rms below 3 THz | 0.08 % |
| AFM fit vs exact FC2, rms below 3 THz | 0.12 % |
| PM pooled fit vs exact spin average, rms below 3 THz | 0.32 % |
| PM − FM signal, rms below 3 THz | 2.0 % |
| jackknife error (median) | 0.25 % |
| even/odd ratio | 0.07 % |
| energies vs Heisenberg | exact |
| full fc2 hdf5 read by phono3py (the reader `kappa_extrinsic.py` uses) | identical frequencies |
| `kappa_extrinsic.py bases` with production FORCE_CONSTANTS_2ND written in this full format + production fc3, 16×16×4 | reproduces `Mn2In2Se5_m16164` κ₀ to 0.18 % (all T, xx and zz; 2 min 45 s here) |

Design study, same synthetic model:
- 4 frames per ordered set gave single-mode errors up to 1.3 %; 8 frames give 0.3 %.
- For PM, 8 configurations × 2 frames gave 0.42 % rms below 3 THz; 8 × 4 gave 0.27 %;
  16 × 2 gave 0.34 %.
- 8 × 2 was chosen as the best accuracy for the cost.

## 4. Gates and decision rules (fixed before the runs)

Gates. All must pass, otherwise stop:

| gate | quantity | threshold |
|---|---|---|
| G1 | FM_ctl vs production fc2, 6×6×2 mesh, modes > 0.5 THz | rms ≤ 1 % |
| G2 | AFM vs FM_ctl Γ shifts vs the 9-atom shifts, optical modes < 3 THz | max deviation ≤ 1.5 pp |
| G3 | PM jackknife, modes < 3 THz | median relative error ≤ 1 % |

κ (`kappa_summary.py`): intrinsic RTA κ from `analysis/extrinsic/kappa_extrinsic.py bases`.
- Each FC2 set (FM_ctl, PM, AFM and the 8 jackknife PM sets) is run with the production
  fc3 on the 16×16×4 mesh, exactly like `analysis/extrinsic/results/Mn2In2Se5_m16164`.
- R = κ(PM)/κ(FM_ctl) uses the same method on both sides, so the method error cancels.

| result at 300 K | meaning |
|---|---|
| \|R − 1\| ≤ 5 % (xx and zz) | the mean spin effect is within the other uncertainties; the FM-based L1 stands, and the paper states the bound |
| otherwise | L1 and the L3 fit are redone on the PM FC2 (the bases are already computed) |

Reported in either case:
- AFM/FM_ctl, as the ordered bound;
- the excess fit residual of PM over FM_ctl, as the size of the FC spread between spin
  configurations (step 2 of §2, the input for a spin-disorder scattering term);
- PM energies against the Heisenberg prediction (−5.74 meV/Mn), an out-of-sample test
  of the exchange fit.

Short-range order: at 300 K, T/|θ_W| ≈ 1.5. The true FC2 lies between PM (uncorrelated)
and AFM (fully correlated interlayer), closer to PM. If R and R_AFM differ substantially,
a correlation-weighted interpolation needs no new DFT: the PM, FM and AFM sets give the
linear dependence on the interlayer correlation.

## 5. Limits

- FC3 stays FM (production). Only the harmonic part is spin-averaged. Spin-averaged
  FC3 would need DLM-MD or TDEP-type data.
- Static, collinear Ising-like spins. Dynamic spin–phonon (magnon) scattering is not
  included. The FC spread is measured but not yet turned into a scattering rate.
- Fixed geometry (production FM structure). The AFM static force (magnetostriction) is
  reported, not relaxed.
- U = 4 eV, as in production. The exchange depends on U.
- The 4×4×1 cell resolves pair correlations up to 7.2 Å.

## Files

| file | what |
|---|---|
| `prepare.py` | run directories, INCAR/KPOINTS/POSCAR, `displacements.json`, `manifest.json`; `--list-configs`, `--regenerate` |
| `spin_configs.json`, `frames.npz` | the committed spin configurations and random displacements |
| `run_set.sh` | Slurm job per set (jzzn cpu192, settings of `../fc3_cutoff_extend/submit_frames.sh`), stages `ref` / `disp` |
| `analyse.py` | `--selftest`; `--refs-only` (step 2); full fit, comparisons, gates; writes `results/` and `jackknife/` |
| `kappa_summary.py` | κ ratios, jackknife error, decision (run on the cloud side from the committed fc2 files; fc3 = `Mn2In2Se5/fc/fc3.hdf5`) |
