# Coherence-term checks (report §7)

Scripts and raw output of the three jobs that back the §7 statements. Nothing here
is a production run; `run_cutoff.py` in particular removes modes from the
population term as well, so its P is a diagnostic, not a result.

| job | script | output | what it shows |
|---|---|---|---|
| 437142 | `submit_wig.sh` -> `wig.py` | `logs/wig-437142.out` | MnIn₂Se₄ 11³ on phono3py **4.4.0**: SMM19 and NJC23 give C = 1.2421/1.1654, IBDB19 gives C = 0.1141/0.0670. Reproduces the local run exactly, so the blow-up is 4.4.0 itself, not the `phonors` pin. |
| 437030 | `run.py` | `result_mesh.json`, `logs/slurm-coh-chk-437030.out` | Mn₂In₂Se₅ 6.51 Å, 11³/12³/14³, SMM19 vs IBDB19: C in-plane 0.1515/0.1574/0.1492, κ_zz 0.1079/0.1114/0.1102 (SMM19); IBDB19 agrees to 1.3–1.7 % in-plane and 3.6–3.8 % in κ_zz. |
| 440011 | `run_cutoff.py` | `result_cutoff.json`, `logs/cut-440011.out` | Mn₂In₂Se₅ 16³, `cutoff_frequency=0.3`: C changes by −7 % in-plane and −0.3 % in κ_zz, so the +12.4 % coherence fraction is insensitive to the soft-mode cutoff. |

## Files

| file | contents |
|---|---|
| `wig.py` | the MnIn₂Se₄ version-comparison script (`<mesh> [SMM19|NJC23|IBDB19] [cutoff_THz]`), same file that sits one directory up |
| `run.py` | Mn₂In₂Se₅ mesh/variant scan |
| `run_cutoff.py` | Mn₂In₂Se₅ cutoff diagnostic |
| `submit_wig.sh`, `submit_cutoff.sh` | the two submit scripts that ran |
| `result_mesh.json`, `result_cutoff.json` | raw κ_xx / κ_zz for population, coherence and total in every case |
| `logs/` | Slurm stdout of the three jobs |

Reference values used in §7 are the 300 K entries of these JSONs.
