# κ(T) at 30³, ordered + VCA, κ_zz as the plain total

Convention for the figures from here on: κ_zz is the **phono3py total on a 30³
mesh**, without the Γ–Z decomposition. At 300 K the 30³ total of MnIn₂Se₄ is
0.7645, about 15–18 % above the mesh-extrapolated 0.646, and the 24 → 30 step
is still −5 %. State it in the text as "30³, may still be high by ~15 %".
In-plane κ is converged and is not affected.

| file | what |
|---|---|
| `run_kappa_vs_T.py` | one material, `ordered` and/or `vca`, 12 temperatures (30–700 K), writes per-mode hdf5 + `result_<mat>_m30.json` |
| `submit.sh` | HanHai job, one per material (`MAT=MnIn2Se4` / `MAT=Mn2In2Se5`) |

```bash
sbatch --export=ALL,MAT=MnIn2Se4  submit.sh
sbatch --export=ALL,MAT=Mn2In2Se5 submit.sh
# afterwards, from test/figures/kappa_vs_T/:
python plot_theory_vs_exp.py            # default --source m30-total
python plot_theory_vs_exp.py --source m20-rest   # the earlier 20^3 rest-term version
```

Expected cost on HanHai (64 threads), from the measured 30³ single-T run
(2017 s for MnIn₂Se₄): about 1–2 h per case for MnIn₂Se₄ and 2–4 h for Mn₂In₂Se₅
(9 atoms); each job runs two cases.

## Results (HanHai jobs 432549 / 432550, 64 threads)

MnIn₂Se₄ (job 432549): ordered 743 s, VCA 741 s. The 300 K ordered point
reproduces the committed 30³ RTA **exactly** — in-plane 3.6946, κ_zz total
0.7645 — which validates the driver against the existing data.

| MnIn₂Se₄, 300 K | ordered 30³ | VCA 30³ | suppression |
|---|---:|---:|---:|
| in-plane | 3.6946 | **1.5340** | 2.41× |
| κ_zz total | 0.7645 | **0.5904** | 1.30× |

Mn₂In₂Se₅ (job 432550): ordered 2326 s, VCA comparable.

The 30³ totals are the "plain total" convention; the rest-of-the-zone values are
lower by ≈0.12 at 300 K (0.6444 at 14³, 0.6473 at 30³) and are what the earlier
20³ VCA tables used.

## Files

`run_<MAT>/result_<MAT>_m30.json` is tracked. The per-mode hdf5 of each case is
41 MB, so they are kept on the backup drive
(`E:\backup\joint_research_project\xia\raw_hdf5\m30\`) rather than in git;
`run_kappa_vs_T.py` regenerates them, and the plotting only needs the JSON.

VCA parameters are those of the existing 20³ runs (MnIn₂Se₄ i = 0.44, Mn₂In₂Se₅
x = 0.13). Their sources are still open: see `analysis/disorder/README.md` §5.
