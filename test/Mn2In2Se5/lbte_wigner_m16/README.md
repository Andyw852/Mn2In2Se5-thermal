# Mn₂In₂Se₅ — LBTE and Wigner-LBTE at the largest FC3 cutoff (6.51 Å)

The Mn₂In₂Se₅ in-plane κ quoted in the reports is **RTA × 1.06 (iterative
correction borrowed from ShengBTE) + 0.135 coherence (from a Wigner-RTA run
whose data was never committed)**. This directory replaces both borrowed numbers
with phono3py runs on the same force constants as the mesh scans.

| item | setting |
|---|---|
| FC2 / FC3 | `Mn2In2Se5/fc/` (4×4×1), FC3 cutoff **6.51 Å** — the ceiling of this cell |
| q-mesh | 16×16×16 |
| temperature | 300 K, isotopes **off** (same as the κ_zz mesh scans) |
| solutions | (1) `is_LBTE=True`, type None; (2) `is_LBTE=True`, `conductivity_type="wigner"` |

HanHai job **412628** (1 node / 64 threads, `--exclude=anode132` after the same
bad-node failure described in `../../MnIn2Se4/lbte_wigner_m14/README.md`).

| file | contents |
|---|---|
| `result.json` | LBTE κ and the Wigner-LBTE TOT / population / coherence split |
| `kappa-m161616.hdf5` | per-mode output of the LBTE run |
| `jobs/run.py`, `jobs/submit.sh` | exactly what ran |

This also closes the reproducibility gap in report §7, where the +4.5 % / +12.4 %
coherence ratios for Mn₂In₂Se₅ were quoted without a committed data file.
