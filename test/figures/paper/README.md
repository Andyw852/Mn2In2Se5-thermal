# Paper figures

```bash
cd test/figures/paper && python make_paper_figures.py     # seconds, no recomputation
```

**One convention for all three figures**, the same as REPORT §9.1:

- κ_zz is the **mesh-converged term**: the Γ–Z quadrature error is removed
  (report §6.3–6.5). The 30³ totals are not plotted.
- The curves are RTA population, with natural isotopes.
- The stars in fig. 1 are the values quoted in the text at 300 K.

| figure | content | data |
|---|---|---|
| `fig1_kappa_vs_T` | intrinsic κ(T) of the ordered crystal, both materials, in-plane and κ_zz. Stars: 300 K quoted values (in-plane LBTE + coherence 4.67 / 3.37; κ_zz + coherence 0.71 / 1.09). Bars: FC3-cutoff uncertainty, MnIn₂Se₄ ≈ −10 % (estimate), Mn₂In₂Se₅ −2.1 % / −3.3 % (bound on the next step from the hold-out test, `../../Mn2In2Se5/cutoff_cv/`; it used to be drawn as an open upper bound) | `../kappa_vs_T/data_mnin2se4.json` (14³), `data_mnin2se5.json` (16³) |
| `fig2_calc_vs_exp` | 2 × 2 panels, material × direction. Hot-pressed experiment; ordered and VCA single-crystal curves (Perp ↔ in-plane, Para ↔ κ_zz); L3 extrinsic fit, texture-mixed | `../kappa_vs_T/data_*_vca.json` (20³ pairs), `../calc_vs_experiment/extrinsic_*.json`, `experiment/kappa_tot_digitized.csv` |
| `fig3_cutoff` | κ vs FC3 cutoff at 300 K, normalised to the largest cutoff, with the step sizes. MnIn₂Se₄ steps shrink. For Mn₂In₂Se₅ an open triangle marks the 6.90 Å refit on the same frames, which is over-parameterised, so it is a lower bound on the next step (≥ −2.1 % / −3.3 %) | `../cutoff/data_*.json` (population and κ_zz rest term; the coherence columns there are 4.4.0 SMM19 values and are **not** used) |

Points to state in the captions:

- **Fig. 1.** The lines (RTA) and the stars (LBTE + coherence) differ by
  construction. For MnIn₂Se₄ in-plane, LBTE is 24 % above RTA.
- **Fig. 2.**
  - Ordered and VCA are single-crystal components without boundary scattering
    or texture. They are an upper reference, not a fit.
  - The antisite occupancies (VCA and L3) are the literature values:
    MnIn₂Se₄ i = 0.44 (single-crystal refinement of Range et al.,
    Z. Naturforsch. 46b, 1122 (1991)) and Mn₂In₂Se₅ x = 0.19(2) (powder Rietveld
    of the R-3m phase, Williams et al., J. Mater. Chem. C 12, 1753 (2024), ESI
    Table S7; `analysis/disorder/README.md` §5.1). Until 2026-11-21 the
    Mn₂In₂Se₅ curves used x = 0.13 from the free L3 fit.
  - L3 has x fixed at that value and fits L and s to the same 300–623 K points
    (Mn₂In₂Se₅ L = 22 nm, s = 0.23; MnIn₂Se₄ L = 31 nm, s = 0.35).
  - The experiment is κ_tot, with κ_e not subtracted.
- **Fig. 3.** Mn₂In₂Se₅ κ_zz steps are −11.9 % / −9.9 % in these data (the
  report quotes −11.8 % from the 1.2451 run; same conclusion). The 6.90 Å
  point comes from a different refit run (`cutoff_ext`) and is drawn as a ratio
  to that run's own 6.51 Å value.
