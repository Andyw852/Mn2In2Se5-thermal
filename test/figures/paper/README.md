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
| `fig1_kappa_vs_T` | intrinsic κ(T) of the ordered crystal, both materials, in-plane and κ_zz. Stars: 300 K quoted values (in-plane LBTE + coherence 4.67 / 3.37; κ_zz + coherence 0.71 / 1.09). The bar on MnIn₂Se₄ is the ≈ −10 % cutoff estimate; the arrow on Mn₂In₂Se₅ marks an upper bound | `../kappa_vs_T/data_mnin2se4.json` (14³), `data_mnin2se5.json` (16³) |
| `fig2_calc_vs_exp` | 2 × 2 panels, material × direction. Hot-pressed experiment; ordered and VCA single-crystal curves (Perp ↔ in-plane, Para ↔ κ_zz); L3 extrinsic fit, texture-mixed | `../kappa_vs_T/data_*_vca.json` (20³ pairs), `../calc_vs_experiment/extrinsic_*.json`, `experiment/kappa_tot_digitized.csv` |
| `fig3_cutoff` | κ vs FC3 cutoff at 300 K, normalised to the largest cutoff, with the step sizes. It shows the one remaining uncertainty: MnIn₂Se₄ steps shrink, Mn₂In₂Se₅ steps do not | `../cutoff/data_*.json` (population and κ_zz rest term; the coherence columns there are 4.4.0 SMM19 values and are **not** used) |

Points to state in the captions:

- **Fig. 1.** The lines (RTA) and the stars (LBTE + coherence) differ by
  construction. For MnIn₂Se₄ in-plane, LBTE is 24 % above RTA.
- **Fig. 2.**
  - Ordered and VCA are single-crystal components without boundary scattering
    or texture. They are an upper reference, not a fit.
  - The VCA occupancies are MnIn₂Se₄ i = 0.44 (single-crystal refinement of
    Range et al., Z. Naturforsch. 46b, 1122 (1991)) and Mn₂In₂Se₅ x = 0.13
    (from the degenerate L3 fit; the literature value for the R-3m phase is
    0.19(2), Williams et al. 2024, ESI Table S7, see `analysis/disorder/README.md` §5.1).
  - L3 is a 3-parameter fit to the same 300–623 K points. Only the texture s
    is determined; x and L are degenerate.
  - The experiment is κ_tot, with κ_e not subtracted.
- **Fig. 3.** Mn₂In₂Se₅ κ_zz steps are −11.9 % / −9.9 % in these data (the
  report quotes −11.8 % from the 1.2451 run; same conclusion).
