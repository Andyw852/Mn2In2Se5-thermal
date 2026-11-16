# HanHai job record (24³ / 30³ meshes and the 14³ cutoff cases)

The 24×24×24 and 30×30×30 q-meshes of §6.3 and the two 14³ cutoff cases of §6.6
were run on HanHai (USTC SCC, `hanhai25`, user `wangchao217`), not on a
workstation. These two files are the exact scripts that ran.

| file | purpose |
|---|---|
| `run_kzz.py` | phono3py driver: reads `in/POSCAR`, `in/fc2_442.hdf5`, `in/<--fc3>`; RTA, isotopes on, 300 K, `write_kappa=True` |
| `submit_case.sh` | `sbatch submit_case.sh <case> <mesh-N> <fc3-file>` — 1 node, 64 OMP threads, partition `CPU-256C256GB`, qos `qos_cpu-256c256gb` |

Remote layout (`$HOME/work/mn441_kappa_zz/`):

```
in/POSCAR  in/fc2_442.hdf5  in/fc3_cut5p60.hdf5  in/fc3_cut6p23.hdf5  in/fc3_cut6p45.hdf5
run_kzz.py  submit_case.sh
log/slurm-kzz-<jobid>.{out,err}
out/<case>/     kappa-m*.hdf5 + result.json
```

The uploaded force constants are md5-identical to the committed ones:

| file | md5 |
|---|---|
| `fc2_442.hdf5` | `9a2b3d55deb306b1313900204cdb613f` |
| `fc3_cut5p60.hdf5` | `5f37ee8a3d95bc57fd7ae1d7156d8344` |
| `fc3_cut6p23.hdf5` | `1def93276574178d13ef0653794ccd3a` |
| `fc3_cut6p45.hdf5` | `d7365ad9b9c99aa537f86adf5264ad83` |

Slurm jobs: **409606** = c560_m14, **409608** = c623_m14, **409609** = m24,
**409612** = m30.

Environment: phonopy 4.4.0, phono3py 4.4.0, Python 3.11.16, numpy 2.4.6,
`$HOME/software/pybin/python3`. Cross-check against the workstation's
phono3py 3.24.0 on the same 14³ / 5.60 Å case: κ_zz 2.523638946 vs 2.523638910
(1.4·10⁻⁸ relative), so the two versions are interchangeable here.

Results are committed under `../kappa_zz_scan/<case>/` (raw `kappa-*.hdf5`,
`result.json`, and the slurm log).
