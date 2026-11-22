# Empirical κ models for Mn₂In₂Se₅, from the calculated phonons

Question (2026-11-22): do empirical formulas based on the phonon spectrum agree with
the hot-pressed experiment? `empirical_mn2in2se5.py` evaluates them with inputs from
the production force constants only. Output: `empirical_mn2in2se5.json`.

## Inputs (all calculated)

| quantity | value |
|---|---|
| acoustic velocities, averaged over 3000 directions | v_L 3536, v_T 1694, Debye mean 1830 m/s |
| θ_D, θ_a = θ_D / n^(1/3) (n = 9 atoms) | 186 K, 89 K |
| mode Grüneisen parameter (from FC3, volumetric, C_v-weighted at 300 K) | all modes 0.82; acoustic mean 0.08 (positive and negative branches cancel), acoustic rms 0.53 |
| number density (experimental lattice) | 0.0395 Å⁻³ |

## Results (W m⁻¹ K⁻¹)

| model | 300 K | 623 K |
|---|---:|---:|
| Slack, γ = 0.82 (all-mode mean) | 6.0 | 2.9 |
| Slack, γ = 0.53 (acoustic rms) | 12.0 | 5.8 |
| first principles, ordered, direction average (LBTE + coherence / RTA 20³) | 2.61 / 2.33 | |
| first principles, VCA x = 0.19, RTA 20³ | 1.23 | |
| **experiment, (2·Perp + Para)/3** | **0.60** | **0.48** |
| Cahill κ_min (full integral) | 0.43 | 0.44 |
| diffuson | 0.28 | 0.28 |

## Reading

- **Slack** estimates the same intrinsic, Umklapp-limited quantity as the BTE
  calculation, only more crudely. With the γ of these phonons it lies above the BTE
  value, and 10–20× above experiment. It reaches the experiment only with γ ≈ 2.4,
  three times the calculated value, and then γ is just a fitted number. It also falls as
  1/T, whereas the experiment falls as about T^−0.3. A match at one temperature would
  therefore be a coincidence.
- **Cahill and diffuson** are lower limits. They are almost temperature-independent
  above θ_D, like the measurement. The measured direction average is only 1.4× the
  Cahill value (Para 0.32 is between the two limits; Perp 0.74 is 1.7× Cahill). So the
  pellets conduct heat close to the amorphous limit. This points to scattering far
  stronger than three-phonon processes, in line with L3 (L ≈ 20 nm). MnIn₂Se₄ from the
  calculated velocities: Cahill 0.40, measured average 0.78, i.e. 2.0×
  (`report/independent_checks.json`).
- Fixed in the report on the same day: the Cahill value in
  `report/independent_checks.py` lacked the factor 1/2 of the high-temperature limit and
  was twice too high (0.62 → 0.31 with the measured MnIn₂Se₄ velocities).
