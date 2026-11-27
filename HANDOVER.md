# Handover and decision log

One entry per decision that changes a quoted number or a published conclusion, newest
last. Each entry: what was decided, by whom, and what it changes.

## Decisions

### 2026-11-27 — In2Mn2Se5 in-plane moves to the paramagnetic FC2 (spin disorder)

**Decision (user, 2026-11-27).** The in-plane lattice thermal conductivity of
In2Mn2Se5 is quoted from the spin-disordered (paramagnetic) FC2, not from the
ferromagnetic one. kappa_zz and In2MnSe4 are unchanged.

- Ratio: `R_xx(T) = kappa_xx(PM) / kappa_xx(FM_ctl)`, 16x16x4 mesh, production fc3,
  `R_xx = 0.822 (30 K), 0.837 (60 K), 0.848 (100 K), 0.853 (150 K), 0.855 (200 K),
  0.857 (300 K), 0.858 (>= 400 K)`, jackknife standard error `+/- 0.023` at 30 K and
  `+/- 0.018` from 300 K up. Stored in
  `dft/spin_disorder_fc2/results/R_xx_vs_T.json`, derived with
  `dft/spin_disorder_fc2/build_R_xx_vs_T.py` from
  `analysis/extrinsic/kappa_extrinsic.py bases` runs on the committed FC2 files.
- Effect: spin disorder lowers the in-plane single-crystal kappa by 14 %, above the
  5 % rule fixed before the run, so the in-plane L1 uses the PM FC2.
- 300 K quoted value: `3.37 -> 2.89` (3.37 x 0.857). The RTA ratio is applied to the
  LBTE + coherence value; that the same ratio holds for that composition is assumed,
  not demonstrated.
- Why kappa_zz is not rescaled: the same dataset cannot resolve the spin effect there
  (PM/FM_ctl jackknife 0.47-0.88, and the FM control is itself +22 % off production
  because the 0.01 A displacements do not resolve the interlayer force constants).
  kappa_zz therefore keeps the FM value and stays an upper bound from the FC3 cutoff.
  The paired-difference rerun that could resolve it was declined (option 1 of
  `dft/spin_disorder_fc2/results/KAPPA_VERDICT.md`, ~3-4 k core-h).
- Why In2MnSe4 is not rescaled: its spin-disorder effect was not computed.
- L3 is not refitted: switching only its in-plane basis would mix PM and FM inputs
  (a mixed basis would need a script change). The note added to L3 states that a PM
  in-plane basis would lower the intrinsic in-plane kappa by 14 % and mainly change
  the fitted L_eff.
- Figures: fig1 and fig2 scale the In2Mn2Se5 in-plane curves by R_xx(T) (fig2 also for
  the VCA curve, assuming spin and antisite effects multiply); kappa_zz is not scaled.
- Files touched: `PROGRESS_REPORT_20261118.md`, `REPORT_MnIn2Se4_20261115.md` (§9.1),
  `README.md`, `test/figures/paper/README.md`, `analysis/disorder/README.md`,
  `test/figures/paper/make_paper_figures.py` and the regenerated fig1/fig2.
