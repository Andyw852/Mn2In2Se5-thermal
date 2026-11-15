# Experimental thermal conductivity (hot-pressed polycrystal)

`kappa_tot_figure.png`: experimental total thermal conductivity κ_tot, two samples (labelled In₂Mn₂Se₅ and In₂MnSe₄ in the figure), two directions:

- **Para**: parallel to the pressing direction. During hot pressing the c axis of the platelet grains tends to align with the pressure direction, so Para ≈ out-of-plane (dominated by κ_zz);
- **Perp**: perpendicular to the pressing direction, ≈ in-plane (dominated by κ_xx).

Two panels: left 300–623 K (filled symbols), right 2–145 K (open symbols).

`kappa_tot_digitized.csv`: **numbers read off the figure**, not raw data. For the high-temperature panel the centroid of each symbol colour was taken, with an uncertainty of about ±0.005 W/mK;
in the low-temperature panel the symbols overlap, so only individual symbols that could be separated were taken (55–140 K), with an uncertainty of about ±0.01. Replace this file once the raw data are available, keeping the column names unchanged:

```
material,direction,T_K,kappa_tot_WmK,kappa_e_WmK,panel
```

`kappa_e_WmK` is the electronic thermal conductivity and is currently left empty (treated as 0). Once the electrical conductivity σ(T) (preferably also the Seebeck coefficient S(T)) is available, fill it in via
κ_e = L·σ·T with L = 1.5 + exp(−|S|/116) ×10⁻⁸ W Ω K⁻² (S in μV/K); the fit then subtracts it automatically.

## Two things noticed while digitizing the figure

1. The mismatch between the low-temperature panel at 145 K and the high-temperature panel at 300 K differs greatly between the two directions: for Mn₂In₂Se₅, Perp goes from 0.79 to 0.74 (−6%),
   while Para goes from about 0.47 to 0.32 (−32%). Low temperature is usually PPMS steady-state and high temperature laser flash, so this may be two instruments or two samples.
   We therefore fit only the high-temperature range first and compare the low-temperature range separately, rather than forcing the two ranges into a single fit.
2. Both samples show a plateau in Perp at 100–145 K and lack the κ ∝ 1/T rise of an intrinsic crystal (the calculated in-plane κ of Mn₂In₂Se₅
   is 6.4 W/mK at 100 K, while the experimental Perp is only 0.79). This shows that the low-temperature range is completely dominated by extrinsic scattering (grain boundaries, point defects, possibly also spin-phonon).

## Experimental data still needed

| data | purpose |
|---|---|
| raw data for both panels (CSV/Excel) | replace the digitized values |
| electrical conductivity σ(T), Seebeck S(T) (both directions) | subtract κ_e |
| XRD: spectra of both faces of the hot-pressed sample (for the Lotgering factor or the (00l)/(hk0) intensity ratio) | independent determination of the texture s, compared with the fitted value |
| SEM: fracture-surface grain size (platelet diameter and thickness) | constrain the grain-boundary scattering length L |
| density (relative to theoretical) | porosity correction |
| if available: Rietveld-refined Mn/In occupancy (antisite fraction) | gives the disorder parameter x directly |
| if available: magnetic susceptibility / magnetic transition temperature | determine whether there is spin-phonon scattering in the low-temperature range |
