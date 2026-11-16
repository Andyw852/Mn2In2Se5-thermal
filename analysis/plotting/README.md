# Plotting programs (recovered from the working directories)

These are the plotting scripts used earlier in the project. They were living only
on the backup drive (`E:\backup\joint_research_project\xia\Mn2In2Se5\fix\`) and
are versioned here now. `plot_phonon_band.py` is left out on purpose: its comments are in Chinese and this repository is English-only, so it remains on the backup drive under `Mn2In2Se5/fix/`.

| script | what it plots |
|---|---|
| `plot_kappa_compare.py` | κ_L(T) from two phono3py hdf5 files: in-plane vs out-of-plane (left panel), isotropic average (right) |
| `plot_mfp_compare.py` | mean free path vs frequency, and cumulative κ vs MFP |
| `plot_scattering_compare.py` | 3-phonon and isotope linewidths vs frequency |
| `plot_phonon_band.py` | phonon band structure |

They read phono3py `kappa-*.hdf5` directly and take the file names as arguments,
so they work on any run under `test/`.

## Conventions used by these scripts

- `figsize=(10, 4.2)`, two panels;
- `ax.grid(True, ls=":", alpha=0.4)`;
- `-o` for in-plane and `-s` with `mfc="white"` for out-of-plane;
- `C_BASE = "#1f77b4"`, `C_VCA = "#d62728"`;
- `ax.legend(fontsize=8, frameon=False)`;
- log y-scale when the range spans more than a decade.

The experiment-versus-calculation figure in `test/figures/calc_vs_experiment/`
follows the same conventions.
