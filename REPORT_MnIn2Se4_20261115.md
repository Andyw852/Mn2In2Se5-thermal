# MnIn₂Se₄ thermal conductivity — 4×4×2 FC2 rework and full convergence audit

**Date:** 2026-11-15
**Scope:** MnIn₂Se₄ (R-3m layered selenide). Mn₂In₂Se₅ results are referenced for comparison.

---

## 1. Motivation

The previous MnIn₂Se₄ κ calculation used a **4×4×1** supercell for the second-order
force constants (FC2). A direct check of the interlayer transverse sound velocity
from that FC2 gave

```
TA along the layer normal = 322.9 m/s      (DFT elastic constants: 1338.3 m/s)
```

a factor of **4.1 too soft**. Because C₄₄ governs the shear channel that carries
heat across the van der Waals gap, a 4× too soft TA velocity suppresses κ_zz by
roughly an order of magnitude. This was identified as a supercell-truncation
problem: a 4×4×1 cell contains only **one** unit cell along the layer normal
(13.214 Å), so interlayer force constants beyond 6.51 Å are folded back
(aliased) into the cell.

**Remedy:** recompute FC2 in a **4×4×2** supercell (224 atoms, safe FC3 cutoff
raised to 8.05 Å), then re-evaluate κ.

---

## 2. What does *not* work: the FC2/FC3 supercell mismatch in ShengBTE

Mn₂In₂Se₅ uses one supercell (4×4×1) for both FC2 and FC3. For MnIn₂Se₄ the
correct FC2 must live in 4×4×2 while FC3 remains in 4×4×1. **phono3py supports
this** through `phonon_supercell_matrix`; **ShengBTE does not** — it requires
`FORCE_CONSTANTS_2ND` and `FORCE_CONSTANTS_3RD` on the *same* supercell.

Three work-arounds were considered and rejected:

| Option | Verdict |
|---|---|
| Run ShengBTE with the old 4×4×1 FC2 | Physically wrong (this is the bug we are fixing) |
| Downsample 4×4×2 FC2 → 4×4×1 | Re-introduces the aliasing that caused the error |
| Build a 4×4×2 FC3 dataset | ~2 days of DFT; the only rigorous ShengBTE route |

In addition, writing ShengBTE input for the 4×4×1 FC2 is **fundamentally
ill-defined**: with one unit cell along c, atom pairs have two degenerate
minimum images (atom 0–atom 2: both images at 8.1479 Å, difference 1.8×10⁻¹⁵),
so the real-space R-vectors cannot be assigned uniquely.

> **Consequence:** for MnIn₂Se₄ the ShengBTE cross-check of the *corrected*
> force constants is not achievable. phono3py with 4×4×2 FC2 is the only
> valid path.

---

## 3. New 4×4×2 FC2

| Item | Value |
|---|---|
| Method | phonopy symmetry-reduced finite displacement (**not** random displacement / ALM) |
| Supercell | **4×4×2** (224 atoms: 32 Mn + 64 In + 128 Se) |
| Symmetry | `strict` — 384 operations, full R-3m point group preserved |
| Displacement amplitude | **0.02 Å** |
| Frames | **7** (the minimal symmetry-complete set) |
| Functional | PBE-D3, U(Mn) = 4 eV, FM (MAGMOM = 32\*5 192\*0) |
| EDIFF | 1E-7 |
| KPOINTS | 2×2×1 |
| NAC | off (semimetal; matches the kl-chain convention) |
| Hardware | jzzn, 64 cores/frame, ≈4 h/frame; cross-check frame on hanhai25 with 128 cores took 7 h 46 min |

### The 7 displacements

All displacement vectors are `(0, 0.011600, 0.016292)` Å (magnitude exactly
0.02 Å); only the sign and the atom differ.

| Frame | Atom | Element | Sign |
|---|---|---|---|
| disp-00001 | 1 | Mn | + |
| disp-00002 | 33 | In | + |
| disp-00003 | 33 | In | − |
| disp-00004 | 97 | Se | + |
| disp-00005 | 97 | Se | − |
| disp-00006 | 161 | Se | + |
| disp-00007 | 161 | Se | − |

### Result

```
PHONON_DONE = true
stable      = true
verdict     = pass
min_freq    = +0.000719 THz   (0 imaginary q-points)
n_qpoint    = 499, n_branch = 21
max_freq    = 7.784 THz
```

The 4×4×1 FC2 gave min_freq = −0.006 THz (soft modes), so the 4×4×2 spectrum is
cleaner as well.

---

## 4. ⚠️ Two pitfalls when extracting sound velocities from FC2

Both were found during this work; the first one initially produced a wrong
conclusion and is documented here to prevent repetition.

### Pitfall 1 — phono3py uses the *selected primitive* reciprocal basis

`primitive_matrix="auto"` replaced the 7-atom cell by a **rhombohedral**
primitive whose b₁/b₂/b₃ all make **84.91°** with the layer normal (i.e. they lie
almost in the layer plane). The layer normal, expressed in *that* basis, is

```
(0,0,t)_POSCAR  ≡  (t,t,t)_phonopy_primitive
```

Computing frequencies along `(0,0,t)` and then converting with the **POSCAR** b₃
length gives **4179 / 6256 m/s** — a spurious result. The ratio of the two |b₃|
is exactly **3.760**, which is the error factor. (This mirrors a previously
reported Mn₂In₂Se₅ issue of identical origin.)

### Pitfall 2 — the Γ-point residual

The FC2 has a small ASR violation, f₀ = 0.000719 THz. Computing ω/q directly is
then contaminated at small q and *diverges* as q → 0 (t = 0.002 → 1429 m/s,
t = 0.010 → 1351 m/s). The correct procedure is to fit

```
f² = f₀² + (v·q / 2π / 100)²
```

and take the square root of the slope. Results are then stable for fit windows
t ≤ 0.005 … 0.040.

**Corrected values (4×4×2 FC2, layer normal, ρ = 5249.3 kg/m³):**

| Quantity | 4×4×2 FC2 | DFT elastic | Deviation |
|---|---:|---:|---:|
| TA₁ | **1348.3 m/s** | 1338.3 m/s | **+0.75 %** |
| TA₂ | **1348.3 m/s** | 1338.3 m/s | **+0.75 %** |
| LA | **3118.6 m/s** | 3253.4 m/s | −4.14 % |
| C₄₄ = ρv_T² | **9.543 GPa** | 9.402 GPa | +1.50 % |

---

## 5. 4×4×1 vs 4×4×2 under one consistent convention

The 4×4×1 pheasy fit uses the *same* `primitive_matrix`, so the two can be
compared directly:

| Supercell | TA₁ | TA₂ | LA | C₄₄ (GPa) | TA dev. | LA dev. |
|---|---:|---:|---:|---:|---:|---:|
| **4×4×1** | **322.9** | 322.9 | 2491.5 | 0.547 | **−75.87 %** | −23.42 % |
| **4×4×2** | **1348.3** | 1348.3 | 3118.5 | 9.543 | **+0.75 %** | −4.15 % |
| DFT elastic | 1338.3 | 1338.3 | 3253.4 | 9.402 | — | — |

**Conclusions**

1. The old 317–323 m/s value was **real**, not a basis error (unlike the 4179
   value, which was).
2. 4×4×2 raises the interlayer TA to 1348 m/s, within **0.75 %** of DFT — **the
   C₄₄ problem is solved**.

---

## 6. κ with the new FC2 (phono3py, 300 K, isotope scattering on)

### 6.1 Effect of the FC2 change alone (mesh 11×11×11, FC3 = 6.45 Å)

| FC2 | Method | In-plane | κ_zz | κ_zz/in-plane |
|---|---|---:|---:|---:|
| 4×4×1 | RTA | 2.2498 | 0.3508 | 0.156 |
| 4×4×1 | LBTE | 3.9054\* | 0.3629\* | 0.093 |
| **4×4×2** | RTA | 3.5375 | 1.2495 | 0.353 |
| **4×4×2** | LBTE | **4.4075** | **1.2519** | **0.284** |

\* historical values, same mesh and FC3.

LBTE change (4×4×1 → 4×4×2): in-plane **+12.9 %**, κ_zz **+245 % (3.45×)**.
The κ_zz jump is the direct consequence of the interlayer TA velocity moving
from 323 → 1348 m/s.

**Self-check:** running the same driver with the *old* FC2 reproduces the
historical RTA value exactly (2.2498 vs 2.24981).

> **Scope of the κ_zz column here.** These are mesh-11 totals and still contain
> the Γ–Z quadrature error of §6.3 (the first Γ–Z q-point alone is 0.6244 of
> 1.2495). They quantify the *effect of the FC2 change*, which is what this
> table is for, but they are not values of κ_zz. The intrinsic κ_zz is 0.648
> ± 0.003 (q-mesh only, §6.5), and LBTE and RTA agree on it to 4·10⁻⁵.

### 6.2 FC3 cutoff scan (mesh 11×11×11)

A 4×4×1 cell caps the FC3 cutoff at **6.507 Å** (half the layer-normal repeat
minus 0.1 Å); the largest contained shell is 6.3557 Å and the next shell
(6.5311 Å) already exceeds the limit. **6.45 Å is therefore the physical
ceiling.** Three fits were taken from the existing pheasy `cutoff_scan`
(5.60 / 6.23 / 6.45 Å; the 6.45 Å scan file reproduces the production κ
exactly, confirming its validity). The scan is run on the **11×11×11** mesh, so
its 6.45 Å row is identical to the mesh-11 row of §6.1 and §6.3.

| Cutoff (Å) | RTA in-plane | RTA κ_zz | LBTE in-plane | LBTE κ_zz |
|---|---:|---:|---:|---:|
| 5.60 | 4.9415 | 2.7672 | 5.6473 | 2.7388 |
| 6.23 | 3.8922 | 1.2332 | 4.7842 | 1.2241 |
| **6.45** | **3.5375** | **1.2495** | **4.4075** | **1.2519** |

Adjacent-step change (LBTE): 5.60 → 6.23 gives **−15.3 %** in-plane and
**−55.3 %** in κ_zz; 6.23 → 6.45 gives **−7.9 %** and **+2.3 %**.

**Verdict: the in-plane FC3 cutoff is not converged** — it still changes by
−7.9 % in the last step and the 4×4×1 ceiling has been reached. Extending the
cutoff requires a 4×4×2 FC3 dataset (~2 days of DFT).

> **The κ_zz column of this table cannot be read as a cutoff effect.** The scan
> runs on the 11×11×11 mesh, where the single first Γ–Z q-point carries 50 % of
> κ_zz (§6.3); its "−55 % then +2 %" pattern is the mesh quadrature, not the
> cutoff. The cutoff comparison for κ_zz has to be made on the converged
> "rest" term at 14³ or 18³ (`kappa_zz_scan/c560_m14`, `c623_m14`, §6.6).

> A symfc-based refit was attempted and **rejected**: at the same 6.45 Å cutoff
> symfc gives an in-plane κ of 0.644 versus pheasy's 3.5375 (5.5× off), and the
> three cutoffs scatter by a factor of 7 non-monotonically. The 46-frame random
> displacement dataset requires pheasy's RFE-OLS sparse regression.

### 6.3 q-mesh convergence

| Mesh | In-plane | Δ | κ_zz | Δ |
|---|---:|---:|---:|---:|
| 11×11×11 | 3.5375 | — | 1.2495 | — |
| 14×14×14 | 3.6474 | +3.11 % | 1.0934 | −12.50 % |
| 18×18×18 | 3.6945 | +1.28 % | 0.9435 | −13.70 % |
| 24×24×24 | 3.7020 | +0.20 % | 0.8052 | −14.66 % |
| **30×30×30** | **3.6946** | **−0.20 %** | **0.7645** | **−5.05 %** |

- **In-plane is converged**: +3.11 % → +1.28 % → +0.20 % → −0.20 %, i.e. it
  oscillates inside ±0.2 % from 18³ on.
- The **total** κ_zz is not converged and is still falling by 14 % per step — but
  that is entirely the quadrature error of one q-point, as the next table shows.

**The cause is the q-point integration near Γ along c*, not a divergence of the
integral.** The phono3py output resolves this directly: the dominant contribution
to the drift comes from the single q-point closest to Γ on the Γ–Z line — the
layer normal, q = (1/N, 1/N, 1/N) in the primitive reciprocal basis. Its weight
in the Brillouin-zone average is 2/N³, while its mode κ grows with N because the
lifetime of the interlayer acoustic branch grows as Γ is approached (measured
τ ∝ ω^-1.3 to ω^-1.9). Its contribution falls as ≈N^-1.8 (fitted over 14³–30³,
§6.5), and on a coarse mesh it is large enough to carry half of the total:

| Mesh | κ_zz total | first Γ–Z point (±q) | rest of the zone | first / total |
|---|---:|---:|---:|---:|
| 11³ | 1.2495 | 0.6244 | 0.6251 | 50.0 % |
| 14³ | 1.0934 | 0.4490 | 0.6444 | 41.1 % |
| 18³ | 0.9435 | 0.2974 | 0.6461 | 31.5 % |
| 24³ | 0.8052 | 0.1594 | 0.6457 | 19.8 % |
| **30³** | **0.7645** | **0.1172** | **0.6473** | **15.3 %** |

The "rest of the zone" column is the converged part and the "first Γ–Z point"
column is the quadrature error. The numbers make the separation unambiguous:
across 11³ → 30³ the total falls by 39 % while the **rest term moves from 0.6251
to 0.6473 (+3.6 % in total, and +0.25 % in the last step)**. So **the intrinsic
κ_zz is finite** — the three-dimensional integral of mode κ converges — and it
is the rest column that converges to it (§6.5).

Two checks that this is the right reading:

- A dense mesh is not missing long-wavelength modes; a *coarse* mesh
  overweights the few Γ–Z points that stand in for the whole small-q region.
  Under *isotropic* refinement (N³) the first point's weight falls as N^-3,
  faster than its lifetime grows, so the total falls.
- **This does not rule out the same artifact in Mn₂In₂Se₅.** There κ_zz rises
  under refinement along c only (16×16×4 → 16×16×16: 0.84 → 1.11 → 1.30 →
  1.54). With only N_z refined, the first Γ–Z point's weight falls as 1/N_z
  while its lifetime grows as N_z^p, so its contribution grows as N_z^(p−1).
  For p ≈ 1.3–1.9 that is a *rising* total. A rising trend there and a falling
  trend here are therefore what one quadrature error predicts for the two
  refinement schemes. Whether Mn₂In₂Se₅ κ_zz is affected has to be checked with
  `kappa_zz_decompose.py --poscar <its cell> --pmat ...` on its meshes. It is
  open (§10), and the Mn₂In₂Se₅ κ_zz values must not be quoted until it is done.
- The rest term is only mildly sensitive to a Casimir boundary (−10 % at
  1 μm, §6.4), whereas the first Γ–Z point collapses by half at 1 μm and to
  4 % of the 30³ total.

### 6.4 Finite-grain κ (Casimir boundary)

All values are **RTA, 300 K, FC3 = 6.45 Å**. They are rebuilt from the per-mode
output of the §6.3 runs by adding the Casimir rate mode by mode. The rebuild is
validated against an explicit phono3py boundary run (11³, 1 μm: 0.8976 both
ways) and reproduces the direct 18³ boundary runs (3.6810 / 0.8849 and
3.5729 / 0.6877). Every number comes from `kappa_zz_decompose.py` (§11).

**A finite-grain value has to use the same estimator as the intrinsic one.** A
boundary shortens the lifetime of the first Γ–Z point, but it does not remove
that point from the sum. The finite-grain *totals* therefore still carry the
quadrature error of §6.3. At 30³ the 10 μm total (0.7233) is even *above* the
intrinsic κ_zz (0.648), which is impossible because boundary scattering can only
lower κ. The totals are listed for comparison with earlier versions of this
report only.

**κ_zz, rest of the zone (the quantity to quote)**

| Boundary | 14³ | 18³ | 24³ | 30³ | mesh-converged (§6.5) |
|---|---:|---:|---:|---:|---:|
| ∞ (intrinsic) | 0.6444 | 0.6461 | 0.6457 | 0.6473 | **0.648 ± 0.003** |
| 10 μm | 0.6385 | 0.6392 | 0.6376 | 0.6381 | **0.638 ± 0.002** |
| 1 μm | 0.5957 | 0.5911 | 0.5846 | 0.5823 | **0.58 ± 0.01** |

**For comparison only: κ_zz totals (contaminated, not values of κ_zz)**

| Boundary | 14³ | 18³ | 24³ | 30³ |
|---|---:|---:|---:|---:|
| ∞ | 1.0934 | 0.9435 | 0.8052 | 0.7645 |
| 10 μm | 1.0282 | 0.8849 | 0.7640 | 0.7233 |
| 1 μm | 0.7758 | 0.6877 | 0.6287 | 0.6070 |

**In-plane (total; no Γ–Z issue)**: 3.6946 (∞), 3.6790 (10 μm), 3.5618 (1 μm) at
30³, i.e. −0.4 % and −3.6 %.

Reading:

- **The grain effect on κ_zz is small.** It is −1.5 % at 10 μm and −10 % at
  1 μm. The earlier "−5.4 % / −20.6 %" (and "−27 %" before that) compared
  contaminated totals, so most of it was the boundary suppressing the
  quadrature error, not real heat current.
- **The 1 μm value is the least well-determined.** Its rest term still drifts
  down by 0.0134 from 14³ to 30³, while the total falls towards it
  (0.6070 − 0.5823 = 0.025 at 30³, shrinking as N^-2.6). Both bracket the
  converged value from either side, so it is quoted as 0.58 ± 0.01.
- Under a boundary the first-point term decays faster (exponent 2.6 at 1 μm,
  2.0 at 10 μm, against 1.8 without a boundary). That is expected: the
  boundary caps the long lifetimes, so the first point becomes an ordinary
  quadrature sample, whose share scales towards N^-3.
- **Essentially no κ_zz is carried by modes with mean free paths above 10 μm**:
  the 10 μm value is 1.5 % below the intrinsic one at every mesh (0.638 vs
  0.648). The "MFP > 1 μm" fraction of the *totals* (31.5 % at 18³) is not
  long-MFP physics at all — it is the first Γ–Z point (§6.3), i.e. the
  quadrature error.

### 6.5 Intrinsic κ_zz: estimate

The drift in §6.3 is a quadrature error of the q-point sum near Γ, so the
quantity to converge is the **rest of the zone**:

| Mesh | rest of the zone | first Γ–Z point | total |
|---|---:|---:|---:|
| 11³ | 0.6251 | 0.6244 | 1.2495 |
| 14³ | 0.6444 | 0.4490 | 1.0934 |
| 18³ | 0.6461 | 0.2974 | 0.9435 |
| 24³ | 0.6457 | 0.1594 | 0.8052 |
| **30³** | **0.6473** | **0.1172** | **0.7645** |

**Decay of the first-point term.** Fitted over 14³–30³, the first-point term
falls as **N^-1.8**, with step exponents of 1.64, 2.17 and 1.38. It is *not* 1/N,
and "the scattering rate falls as ω²" is not what the data show. A separate
single-point lifetime estimate of the interlayer acoustic branch gives
τ ∝ ω^-p with p ≈ 1.3–1.9, consistent with the measured N^-1.8; that estimate is
not part of this repository, and **the exponent quoted here (N^-1.8) is measured
directly from the committed data**. Any p < 3 leaves the 3D integral finite, so
the intrinsic κ_zz is well defined.

**Bias of the rest term.** Dropping the first point also drops the true
contribution c(N) of the cell around it, so rest(N) = κ_zz − c(N). Taking
c(N) ∝ N^-1.8, the measured drift of the rest term from 14³ to 30³ (+0.0029)
gives c(30³) ≈ +0.001. The uncertainty is the larger of that bias and the full
spread of the rest term over 14³–30³ (0.0029; note the non-monotonic
18³ → 24³ step).

**Independent check (fixed sphere).** Split the rest term at a fixed physical
radius q_c around Γ. Once q_c spans several mesh spacings, both parts must be
mesh-independent. At q_c = 0.02 Å⁻¹ (no 2π) the part outside the sphere gives
0.4486 / 0.4424 / 0.4425 and the part inside it, without the first point,
gives 0.1974 / 0.2033 / 0.2048 at 18³ / 24³ / 30³. At q_c = 0.03 Å⁻¹ the two
parts are 0.3356 / 0.3239 / 0.3279 and 0.3105 / 0.3218 / 0.3194. The converged
total therefore includes a real ≈0.2–0.3 from the small-q region. What is
removed is only the overweighted single point, not the small-q physics.

**Result (RTA, 300 K, FC3 = 6.45 Å):**

> **κ_zz(intrinsic) = 0.648 ± 0.003 W m⁻¹ K⁻¹ — q-mesh uncertainty only.**

This is **not** a converged κ_zz. The FC3 cutoff is not converged (§6.6). On this
same estimator the last available cutoff step (6.23 → 6.45 Å) still lowers κ_zz
by 8.2 %, and in-plane κ by 8.7 %. Three cutoffs with a hard ceiling at 6.45 Å
cannot be extrapolated, so neither the size nor the sign of the remaining
change is known. In particular, "≈8 % downward" is the size of the last step,
not a bound. The next shells (6.53, 6.78, 6.89 Å) need the 4×4×2 FC3 dataset.

Quote: **κ_zz = 0.648 (mesh-converged, ±0.003); FC3 cutoff not converged (last
step −8 %).** The totals 0.9435 (18³), 0.8052 (24³) and 0.7645 (30³) are
quadrature-contaminated and are not values of κ_zz. The coherence term of §7
(+0.066 at 14³) is +10 % of this value, not the +6.0 % quoted there against the
contaminated 14³ population term.

### 6.6 FC3 cutoff comparison on the converged term

The committed cutoff scan runs at 11³, where the first Γ–Z point is 50 % of
κ_zz, so its κ_zz column is the artifact (§6.2). The comparison is repeated at
14³ on the *rest* term:

| fc3 cutoff (Å) | in-plane (14³) | κ_zz total (14³) | κ_zz rest of the zone (14³) |
|---|---:|---:|---:|
| 5.60 | 5.1834 | 2.5236 | 1.4175 |
| 6.23 | 3.9944 | 1.1297 | 0.7020 |
| **6.45** | **3.6474** | **1.0934** | **0.6444** |

| interval | in-plane | κ_zz rest | κ_zz rest, as measured on the contaminated 11³ total |
|---|---:|---:|---:|
| 5.60 → 6.23 | −22.9 % | **−50.5 %** | −55.3 % |
| 6.23 → 6.45 | −8.7 % | **−8.2 %** | +2.3 % |

**Verdict: on the converged term the κ_zz cutoff behaves exactly like the
in-plane one — it drops by half between 5.60 and 6.23 Å and then by a further
8.2 %, which is the same −8.7 % the in-plane component still shows.** So the
earlier "+2.3 % in the last step" was the Γ–Z artifact cancelling part of a real
−8 % cutoff effect, and the conclusion is not "κ_zz is flat at 6.45 Å" but
**"the FC3 cutoff is not converged for either component; the last available
step still moves both by ≈8 %, and the remaining change is unbounded, and the 4×4×1 cell cannot go further"**. That is what makes the
4×4×2 FC3 dataset the first priority of §10: it is needed for κ_zz by the same
8 % as for the in-plane, and, being one layer thick, the 4×4×1 cell is also the
suspect for the interlayer terms that set the κ_zz lifetimes.

---

## 7. Wigner (coherence) contribution

All results above use phono3py's default `conductivity_type=None`, i.e. the
**population** (particle-like) term only. phono3py also implements the Wigner
decomposition, which adds the **coherence** (wave-like, off-diagonal) channel.
That term is not included in any of the earlier numbers.

**Wigner-RTA, MnIn₂Se₄, mesh 14×14×14, FC3 = 6.45 Å, 300 K**

| Component | In-plane | κ_zz |
|---|---:|---:|
| **Total (TOT)** | **3.7492** | **1.1595** |
| Population (P) | 3.6584 | 1.0934 |
| Coherence (C) | 0.0908 | 0.0661 |
| C / P | **+2.5 %** | **+6.0 %** |

The population term reproduces the plain RTA κ_zz exactly (1.0934), confirming
internal consistency.

**For comparison, Mn₂In₂Se₅ (16×16×4, FC3 = 6.51 Å, 300 K):**

| Component | In-plane | κ_zz |
|---|---:|---:|
| Total (TOT) | 3.1102 | 1.0122 |
| Population | 2.9752 | 0.9007 |
| Coherence | 0.1350 | 0.1115 |
| C / P | +4.5 % | **+12.4 %** |

Both materials follow the same rule: **the lower κ_zz, the larger the coherence
fraction.** For Mn₂In₂Se₅ the calculation requires `primitive_matrix="P"`;
with `"auto"` the 16×16×4 mesh breaks grid symmetry (phono3py silently falls
back to a 5×5×5 grid, giving a spurious κ_zz of 2.90 instead of 0.90).

---

## 8. FC3 → ShengBTE export audit (both materials)

Writing `fc3.hdf5` to ShengBTE's `FORCE_CONSTANTS_3RD` via hiphive's
`write_to_shengBTE` silently drops every triplet that is not already the
shortest image inside the supercell. For a 4×4×1 cell this removes all
triplets crossing the layer-normal boundary.

| Material | Non-zero triplets in fc3.hdf5 | hiphive export | Complete writer |
|---|---:|---:|---:|
| Mn₂In₂Se₅ | 9164 | 6038 (**−34 %**) | — |
| **MnIn₂Se₄** | **5586** | **3276 (−41.4 %)** | **5586** |

### MnIn₂Se₄ CONTROL file bug

The 20×20×5 production CONTROL contains

```
positions(:,1)=0.0000000000 0.0000000000 0.0000000000     ← wrong
```

whereas the correct value is `(1, 0, 1)` — Mn sits at fractional
(0.99999999996, 0, 0.99999999999), which was mis-snapped to 0. This is the
AutoZT bug fixed in commits `12115b1` (rebase, snap values near 0/1) and
`a226758` (FORCE_CONSTANTS_2ND ordering). **All ShengBTE-derived κ for
MnIn₂Se₄ published before this fix are invalid.**

The 20×20×5 triplet count is complete (5586), so the failure there is the atom
placement, not the triplet filter — but either defect invalidates the result.

---

## 9. Summary of κ values

All at 300 K, isotope scattering on.

| Configuration | Method | In-plane | κ_zz | Status |
|---|---|---:|---:|---|
| FC2 4×4×1, mesh 11, 6.45 Å | RTA | 2.2498 | 0.3508 | superseded (FC2 too soft) |
| FC2 4×4×1, mesh 11, 6.45 Å | LBTE | 3.9054 | 0.3629 | superseded |
| **FC2 4×4×2, mesh 11, 6.45 Å** | LBTE | 4.4075 | 1.2519 | cutoff not converged |
| **FC2 4×4×2, mesh 14, 6.45 Å** | LBTE | 4.5219 | 1.0934 | cutoff not converged |
| **FC2 4×4×2, mesh 18, 6.45 Å** | RTA | 3.6945 | 0.9435 | **upper bound on the intrinsic κ_zz** |
| **FC2 4×4×2, mesh 24, 6.45 Å** | RTA | 3.7020 | 0.8052 | in-plane converged; κ_zz total still contaminated (§6.5) |
| FC2 4×4×2, mesh 30, 6.45 Å | RTA | 3.6946 | 0.7645 | same |
| FC2 4×4×2, mesh 14, 6.45 Å | Wigner-RTA | 3.7492 | 1.1595 | coherence included |
| FC2 4×4×2, 14³–30³, Γ–Z quadrature removed, 10 μm | RTA | 3.6790 | **0.638 ± 0.002** | finite grain (§6.4) |
| FC2 4×4×2, 14³–30³, Γ–Z quadrature removed, 1 μm | RTA | 3.5618 | **0.58 ± 0.01** | finite grain (§6.4) |
| **FC2 4×4×2, 14³–30³, Γ–Z quadrature removed** | RTA | 3.6946 | **0.648 ± 0.003** | **intrinsic κ_zz, q-mesh converged; FC3 cutoff not converged** (§6.5) |

**Reading this table.** The in-plane component is converged (+0.20 % from 18³ to
24³). The κ_zz column needs care: the totals at 11³–24³ — **including the LBTE
ones** — are contaminated by a q-point quadrature error near Γ along c* (§6.3),
which is why they fall by 14 % per refinement even though the rest of the zone
is converged. Removing that contribution gives the intrinsic value
**0.648 ± 0.003 W m⁻¹ K⁻¹** (§6.5). The remaining uncertainty is the FC3 cutoff
(§6.6: the last step still moves κ_zz by −8 %, with no bound on what follows),
not the mesh, and it applies to RTA and LBTE alike (for
κ_zz the two agree to 4·10⁻⁵).

---

## 10. Open items

Ordered by how much each item can still move κ_zz.

| # | Item | Cost | Why it matters |
|---|---|---|---|
| **1** | **FC3 in a 4×4×2 supercell (cutoff > 6.45 Å)** | ~2 days DFT | **Now the largest uncertainty in κ_zz, not only an in-plane cutoff issue.** The 4×4×1 FC3 cell is one layer thick (13.214 Å along c) — the same geometric defect that destroyed the 4×4×1 FC2 (§4). The interlayer acoustic modes that carry κ_zz take their lifetimes from exactly the interlayer FC3 terms that this cell cannot represent. |
| 2 | 24³ and 30³ q-mesh RTA | 9 min + 34 min on HanHai | **done** (§6.3, §6.5): the rest term is converged, κ_zz = 0.648 ± 0.003 on the q-mesh (the FC3 cutoff, item 1, is not settled) |
| 3 | FC3 cutoff comparison for κ_zz on the converged "rest" term (14³) | ~1 h | **done** — §6.6: on the converged term κ_zz falls −50.5 % then −8.2 %, i.e. like the in-plane −22.9 % / −8.7 %; the old "+2.3 %" was the artifact |
| 4 | **Γ–Z decomposition of the Mn₂In₂Se₅ κ_zz meshes** | minutes (post-processing) | **done** — `test/Mn2In2Se5/kappa_zz_decomposition/`. The artifact is present and its share grows (first point 0.190 → 0.538, share 22.6 % → 34.8 % for 16×16×4 → 16³), **but the rest of the zone also grows 0.651 → 1.009**, unlike here where it settles at 0.648. So the Mn₂In₂Se₅ rise is *not* the artifact alone, and its κ_zz remains unconverged and must not be quoted (the runs reproduce the committed reference scan to 0.5 %). An isotropic 24³ run is queued. |
| 5 | Unify the structure for the LA comparison | small | phonon cell has 13.21409 Å vs elastic 13.20318 Å along c |

**Dropped: 18³ LBTE and the Wigner-LBTE refinement.** At both 11³ and 14³ the
LBTE and RTA values of κ_zz agree to 4·10⁻⁵ (§6.3), so inverting the full
collision matrix on a denser mesh buys nothing for κ_zz; and the coherence term
is already quantified at RTA level (§7). Neither is worth the machine time.

**Note on the two runs that were lost.** Both are single long phono3py calls that
keep the whole collision matrix in memory and cannot be checkpointed. The
workstation was rebooted repeatedly during them, so each attempt died before
completion. They are superseded by items 2–4 above and no conclusion depends on
them.

---

## 11. Data and reproducibility — repository file map

Every calculation in this report has its inputs and outputs committed to the
repository. The full index is **[test/MnIn2Se4/README.md](test/MnIn2Se4/README.md)**;
the table below is the direct lookup.

All paths are relative to the repository root.

### 11.1 Calculation → inputs / outputs

| § | Calculation | Inputs (committed) | Outputs (committed) | Driver |
|---|---|---|---|---|
| 3 | **FC2 in a 4×4×2 supercell** (7 frames, 0.02 Å) | `test/MnIn2Se4/fc2_442_phonon/frames/disp-0000{1..7}/` (INCAR, KPOINTS, POSCAR, CONTCAR), `fc2_442_phonon/disp_plan.json`, `fc2_442_phonon/POSCAR` | `fc2_442_phonon/fc2.hdf5`, `fc2_442_phonon/FORCE_CONSTANTS`, `fc2_442_phonon/band-dft-cpu.yaml`, `fc2_442_phonon/phonon_summary.json`, `fc2_442_phonon/phonon_fit.log`, `fc2_442_phonon/phonon_band_plot/*.png` | phonopy finite displacement (AutoZT phonon-dft-cpu S1→S3) |
| 4 | **Interlayer sound velocity / C₄₄** | `test/MnIn2Se4/fc2_442_phonon/fc2.hdf5`, `fc2_442_phonon/POSCAR` | printed table (values in §4 of this report) | `test/MnIn2Se4/fc2_442_phonon/sound_velocity.py` |
| 5 | **4×4×1 vs 4×4×2** | `test/MnIn2Se4/fit/POSCAR`, `test/MnIn2Se4/fit_cut3_6p45/fc2.hdf5`, `test/MnIn2Se4/fc2_442_phonon/fc2.hdf5` | `test/MnIn2Se4/compare_441_vs_442_results.json` | `test/MnIn2Se4/compare_441_vs_442.py` |
| 6.1 | **κ vs FC2 supercell** | the two FC2 files above + `test/MnIn2Se4/fit_cut3_6p45/fc3.hdf5` | `test/MnIn2Se4/kappa_results_442.json`, `test/MnIn2Se4/kappa_results_mn441.json` | `test/MnIn2Se4/kappa_driver_442.py` |
| 6.2 | **FC3 cutoff scan** (5.60 / 6.23 / 6.45 Å) | `test/MnIn2Se4/fc3_cutoff_variants/fc3_cut5p60.hdf5`, `fc3_cut6p23.hdf5`, `fc3_cut6p45.hdf5` + `fc2_442_phonon/fc2.hdf5` | `test/MnIn2Se4/kappa_cutoff_scan_results.json`, `kappa_cutoff_convergence.png`, `KAPPA_CUTOFF_SCAN.md`; κ_zz cutoff comparison on the converged term: `kappa_zz_scan/c560_m14/`, `c623_m14/` | `test/MnIn2Se4/kappa_cutscan.py` (RTA), `test/MnIn2Se4/kappa_lbte_scan.py` (LBTE), `kappa_zz_scan.sh` |
| 6.3 | **κ vs q-mesh** (11³ / 14³ / 18³ / 24³ / 30³) and the Γ–Z decomposition | `fc2_442_phonon/fc2.hdf5`, `fit_cut3_6p45/fc3.hdf5` | `kappa_zz_scan/m{11,14,18,24,30}/kappa-m*.hdf5` (raw phono3py output, per-mode κ), `kappa_zz_gamma_line.json`, `kappa_mesh_scan_results.json`, `kappa_mesh_convergence.png`, `KAPPA_MESH_SCAN.md` | `test/MnIn2Se4/kappa_zz_scan.sh` (driver of the whole scan), `kappa_driver_442.py --mesh N N N --rta --write-kappa`, `kappa_zz_decompose.py` (Γ–Z split + finite-grain rebuild) |
| 6.4 | **Finite-grain κ** (Casimir 10 μm / 1 μm), on the rest term | same force constants; the values are rebuilt from the per-mode data of the 6.3 runs, and validated against `kappa_zz_scan/m11_b1um/` (a real `boundary_mfp=1.0` run) | `test/MnIn2Se4/kappa_zz_gamma_line.json` (`finite_grain` per case, `estimate` mesh-converged); `kappa_finite_grain.json` is the older 18³ **totals** and is superseded | `test/MnIn2Se4/kappa_zz_decompose.py --mfp 1 10`; `kappa_driver_442.py --boundary-mfp 1.0` for the validation run |
| 7 | **Wigner coherence** | same force constants | `test/MnIn2Se4/kappa_wigner_mn441.json` | `test/MnIn2Se4/kappa_driver_442.py` (`conductivity_type="wigner"`), `test/MnIn2Se4/kappa_run_14_wigner_lbte.py` |
| 8 | **FC3 → ShengBTE export audit** | `test/MnIn2Se4/fit_cut3_6p45/fc3.hdf5`, `test/MnIn2Se4/phono3py_disp.yaml` | `test/MnIn2Se4/shengbte_corrected/CONTROL_fixed`, `shengbte_corrected/FORCE_CONSTANTS_3RD_complete` (5586 triplets), `shengbte_corrected/triplet_audit.json` | `analysis/shengbte_crosscheck/write_fc3_sheng.py`, `analysis/shengbte_crosscheck/write_fc2_sheng.py` |

### 11.2 Legacy (superseded) data

| Path | Contents | Status |
|---|---|---|
| `test/MnIn2Se4/fit/` | pheasy FC2+FC3 fit in the 4×4×1 cell, dataset displacements/forces, fit config and manifest | input provenance |
| `test/MnIn2Se4/fit_cut3_6p45/` | 4×4×1 FC2 + FC3 at cutoff 6.45 Å | still used as the FC3 of record |
| `test/MnIn2Se4/kappa_shengbte_11x11x11/`, `kappa_shengbte_20x20x5/` | ShengBTE runs | **κ invalid** — see §8 |
| `test/MnIn2Se4/kappa_phono3py_lbte_11x11x11/` | phono3py LBTE with the old 4×4×1 FC2 | superseded by §6 |
| `test/MnIn2Se4/kappa_finite_grain.json` | κ_zz **totals** at 18³ with 10 μm / 1 μm / ∞ boundary | **superseded** — the totals contain the Γ–Z quadrature error; §6.4 now quotes the rest term, stored in `kappa_zz_gamma_line.json` |
| `test/MnIn2Se4/shengbte_export/` | the defective `FORCE_CONSTANTS_3RD` / `CONTROL` pair | kept as evidence for §8 |
| `MnIn2Se4/` (repository root) | the original 3×3×3 test dataset (κ ≈ 3.9 / 2.1) | **superseded** by this 4×4×2 work |

### 11.3 Aggregate numbers

- `test/MnIn2Se4/kappa_results_mn441.json` — one JSON with every κ value from
  this report (FC2 comparison, cutoff scan, mesh scan, finite grain, Wigner,
  FC3 triplet audit, and the explicit conclusions block).
- `test/MnIn2Se4/kappa_results_442.json` — the raw output of `kappa_driver_442.py`.
- `test/MnIn2Se4/kappa_zz_gamma_line.json` — the Γ–Z decomposition of every
  mesh of §6.3: total κ_zz, the contribution of the first Γ–Z q-point, the
  remaining "rest" term (the mesh-converged quantity), and the 1 μm / 10 μm
  finite-grain reconstruction from the same file.
- `test/MnIn2Se4/kappa_zz_scan/m{11,14,18,24,30}/kappa-m*.hdf5` — the raw phono3py
  output of the scan, including per-mode κ, so `kappa_zz_decompose.py` can be
  re-run without repeating any calculation.
- `test/MnIn2Se4/k18_lbte.json` / `k14_wigner_lbte.json` — outputs of the two
  LBTE refinement drivers; they are **not** part of any conclusion (see §10).

### 11.4 Execution hosts

Two machines were used, and every number in this report is tagged by which one
produced it:

| Machine | Software | Which § |
|---|---|---|
| workstation (12 cores) | phonopy 2.47.1, phono3py 3.24.0, Python 3.11 | §3–§8: FC2 fit, sound velocity, the 11³/14³/18³ meshes, the cutoff scan at 11³, Wigner, the 1 μm validation run, and the 14³ cutoff cases c560/c623 |
| **HanHai (USTC SCC, `hanhai25`)** | phonopy 4.4.0, phono3py 4.4.0, Python 3.11.16, numpy 2.4.6 | the **24³ and 30³** meshes of §6.3 (1 node, 64 OMP threads, partition `CPU-256C256GB`, qos `qos_cpu-256c256gb`) |

**Cross-check of the two versions and of the two thread counts.** The 14³,
fc3 = 5.60 Å case was run on both machines:

| | in-plane | κ_zz | wall time |
|---|---:|---:|---:|
| workstation, phono3py 3.24.0, 12 threads | 5.18341238 | 2.523638910 | 601 s |
| HanHai, phono3py 4.4.0, 64 threads | 5.18341203 | 2.523638946 | 24 s |

Relative differences 7·10⁻⁸ and 1.4·10⁻⁸ — the two versions are interchangeable
here, and 64 threads give a 25× speed-up. The HanHai job record (scripts,
inputs, md5 of the uploaded force constants, job ids) is kept in §11.6.

### 11.5 Environment

phonopy 2.47.1, phono3py 3.24.0, symfc 1.6.0, Python 3.11.
`sound_velocity.py` and the `kappa_*.py` drivers need only numpy + phonopy +
phono3py; no proprietary code is required to reproduce the numbers in this
report.

### 11.6 HanHai job record

The 24³ / 30³ meshes and the two 14³ cutoff cases ran on HanHai under
`$HOME/work/mn441_kappa_zz` (`wangchao217@hanhai25`): slurm jobs **409606** =
c560_m14, **409608** = c623_m14, **409609** = m24, **409612** = m30. Layout:

```
~/work/mn441_kappa_zz/
├── in/POSCAR  in/fc2_442.hdf5  in/fc3_cut5p60.hdf5  in/fc3_cut6p23.hdf5  in/fc3_cut6p45.hdf5
│     md5 9a2b3d55… 5f37ee8a… 1def9327… d7365ad9…   (identical to the committed files)
├── run_kzz.py          # the phono3py driver, same physics as kappa_driver_442.py
├── submit_case.sh      # sbatch submit_case.sh <case> <mesh> <fc3-file>
├── log/slurm-kzz-<jobid>.{out,err}
└── out/<case>/         # cwd of the run: kappa-m*.hdf5 + result.json
```

The per-case `result.json`, the `kappa-*.hdf5` files and the slurm logs are
committed under `test/MnIn2Se4/kappa_zz_scan/<case>/`, so the repository holds
one raw-output directory per mesh and per cutoff.

### 11.7 Reproducing a single number

```bash
cd <repo-root>/test/MnIn2Se4
python compare_441_vs_442.py            # §5  TA 322.9 vs 1348.3 m/s
python fc2_442_phonon/sound_velocity.py .   # §4  TA/LA/C44 of the 4x4x2 FC2
python kappa_cutscan.py                 # §6.2  RTA cutoff scan (~5 min)
python kappa_lbte_scan.py               # §6.2  LBTE cutoff scan (~47 min)
python kappa_driver_442.py              # §6.1  RTA + LBTE vs FC2 supercell
bash   kappa_zz_scan.sh                 # §6.3  whole mesh + cutoff scan (CPU only, resumable)
python kappa_zz_decompose.py kappa_zz_scan/*/kappa-m*.hdf5 \
       --json kappa_zz_gamma_line.json --mfp 1 10     # §6.3/§6.4 table
```

The two LBTE refinement drivers (`kappa_run_18_lbte.py`,
`kappa_run_14_wigner_lbte.py`) are kept for completeness but are not part of any
conclusion, for the reasons in §10.

> Raw DFT output (OUTCAR, vasprun.xml, WAVECAR) is not committed — see the
> repository README. The per-frame INCAR/KPOINTS/POSCAR/CONTCAR **are**
> committed, so the FC2 fit can be reproduced with one phonopy call once the
> forces are regenerated.
