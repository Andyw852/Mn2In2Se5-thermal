# HanHai job record — 24³ / 30³ meshes and the 14³ cutoff cases

The 24×24×24 and 30×30×30 q-meshes of §6.3 and the two 14³ cutoff cases of §6.6
were run on HanHai (USTC SCC, `hanhai25`, user `wangchao217`), not on a
workstation. This directory is the complete, verified record of those runs.

| file | purpose |
|---|---|
| `run_kzz.py` | phono3py driver: reads `in/POSCAR`, `in/fc2_442.hdf5`, `in/<--fc3>`; RTA, isotopes on, 300 K, `write_kappa=True` |
| `submit_case.sh` | `sbatch submit_case.sh <case> <mesh-N> <fc3-file>` — 1 node, 64 OMP threads, partition `CPU-256C256GB`, qos `qos_cpu-256c256gb` |
| `logs/slurm-kzz-<jobid>.out` | stdout of each job (the `.err` files were empty) |
| `MANIFEST.sha256` | sha256 of **every** file that ever existed in the remote run directory, so the local copy can be audited |

Remote layout (`$HOME/work/mn441_kappa_zz/`):

```
in/POSCAR  in/fc2_442.hdf5  in/fc3_cut5p60.hdf5  in/fc3_cut6p23.hdf5  in/fc3_cut6p45.hdf5
run_kzz.py  submit_case.sh  kappa_zz_decompose.py
log/slurm-kzz-<jobid>.{out,err}
out/<case>/     kappa-m*.hdf5 + result.json
```

Slurm jobs: **409606** = c560_m14 (24 s), **409608** = c623_m14 (24 s),
**409609** = m24 (548 s), **409612** = m30 (2017 s).

## Everything was pulled back and verified

All raw output is committed next door under `../kappa_zz_scan/<case>/`. The
transfer was checked with `MANIFEST.sha256`: 18 of 18 files match, i.e. the
committed `kappa-*.hdf5`, `result.json`, both scripts and all five input files
are byte-identical to what ran on HanHai. To re-verify:

```bash
cd test/MnIn2Se4
sha256sum -c hanhai_jobs/MANIFEST.sha256      # paths are relative to the remote tree
```

(That `-c` only works after recreating the remote layout, which is why the
mapping below is spelled out.)

| `MANIFEST.sha256` path | committed path |
|---|---|
| `./run_kzz.py` | `hanhai_jobs/run_kzz.py` |
| `./submit_case.sh` | `hanhai_jobs/submit_case.sh` |
| `./in/POSCAR` | `fit/POSCAR` |
| `./in/fc2_442.hdf5` | `fc2_442_phonon/fc2.hdf5` |
| `./in/fc3_cut5p60.hdf5` | `fc3_cutoff_variants/fc3_cut5p60.hdf5` |
| `./in/fc3_cut6p23.hdf5` | `fc3_cutoff_variants/fc3_cut6p23.hdf5` |
| `./in/fc3_cut6p45.hdf5` | `fit_cut3_6p45/fc3.hdf5` |
| `./out/c560_m14/` | `kappa_zz_scan/c560_m14/` |
| `./out/c623_m14/` | `kappa_zz_scan/c623_m14/` |
| `./out/m24/` | `kappa_zz_scan/m24/` |
| `./out/m30/` | `kappa_zz_scan/m30/` |
| `./out/local_m11,m14,m18/` | `kappa_zz_scan/m11,m14,m18/` (uploaded to HanHai for the joint decomposition) |
| `./kappa_zz_decompose.py` | intermediate version; the authoritative one is `../kappa_zz_decompose.py` |

Environment: phonopy 4.4.0, phono3py 4.4.0, Python 3.11.16, numpy 2.4.6,
`$HOME/software/pybin/python3`. Cross-check against the workstation's
phono3py 3.24.0 on the same 14³ / 5.60 Å case: κ_zz 2.523638946 vs 2.523638910
(1.4·10⁻⁸ relative), so the two versions are interchangeable here.
