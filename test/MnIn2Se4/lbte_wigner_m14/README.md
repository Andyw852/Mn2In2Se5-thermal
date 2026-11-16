# MnIn₂Se₄ — LBTE and Wigner-LBTE at the largest FC3 cutoff (6.45 Å)

These two runs close the last gap in the quotable in-plane κ. The report quotes
in-plane κ as **LBTE + a coherence term**, but until now the coherence term was
measured only at RTA level and only at 14³, and the LBTE value came from a
different (older) run.

| item | setting |
|---|---|
| FC2 | 4×4×2 (`fc2_442_phonon/fc2.hdf5`) |
| FC3 | 4×4×1, **6.45 Å** (`fit_cut3_6p45/fc3.hdf5`) — the ceiling of this cell |
| q-mesh | 14×14×14 |
| temperature | 300 K, isotope scattering **on** |
| solutions | (1) `is_LBTE=True`, type None; (2) `is_LBTE=True`, `conductivity_type="wigner"` |

HanHai job **412627** (1 node / 64 threads). The script was first submitted
without `--exclude`, failed instantly on a bad node (`anode132`, exit 0:53), and
ran on resubmit with `--exclude=anode132` — noted here because the same failure
mode hit an earlier job.

| file | contents |
|---|---|
| `result.json` | LBTE κ and the Wigner-LBTE TOT / population / coherence split |
| `kappa-m141414.hdf5` | per-mode output of the LBTE run |
| `jobs/run.py`, `jobs/submit.sh` | exactly what ran |

Purpose: replace the mixed-provenance numbers in report §9.1
(4.5219 LBTE at 14³ + 0.091 coherence from Wigner-RTA) with a single LBTE run
that yields both components at the same cutoff and mesh.
