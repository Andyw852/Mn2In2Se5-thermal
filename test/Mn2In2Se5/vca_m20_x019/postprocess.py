#!/usr/bin/env python3
"""Mn2In2Se5 VCA curve at the literature x = 0.19 (20^3, RTA, 12 T) -> the kappa(T) data file.

Route used (2026-11-21): first_gz_point.py computes the first Gamma-Z grid point at x = 0.19
with the driver of HanHai job 5, and rest(T) = kzz_total(job 5, T) - first_point(T). Its
output first_gz_point.json is the input here. This script

  1. writes decomposition.json (the per-T split, same keys as kappa_zz_decompose.py);
  2. writes test/figures/kappa_vs_T/data_mn2in2se5_vca.json at x = 0.19 (baseline_* taken
     unchanged from the x = 0.13 file of job 424496: same mesh, same temperatures, ordered masses);
  3. only if the per-mode hdf5 of job 5 has been copied here as kappa_vca_rta_m20.hdf5
     (HanHai: run_Mn2In2Se5_vca_m20_rta/, written with WRITE_KAPPA=1): decomposes it at all 12 T
     with test/MnIn2Se4/kappa_zz_decompose.py, checks it against step 1, and writes the 300 K slice
     kappa_vca_300K.hdf5 and test/figures/scattering/mn2in2se5/kappa_vca.hdf5 for the
     scattering-rate figure. Without that file the scattering figure stays at x = 0.13.

  python first_gz_point.py && python postprocess.py
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SPLIT = HERE / "first_gz_point.json"
H5 = HERE / "kappa_vca_rta_m20.hdf5"
JOB5 = REPO / "test/figures/kappa_vs_T/jobs_lbte_wigner/run_Mn2In2Se5_vca_m20_rta/result_Mn2In2Se5_vca_m20.json"
OLD = REPO / "test/figures/kappa_vs_T/data_mn2in2se5_vca_x013.json"
NEW = REPO / "test/figures/kappa_vs_T/data_mn2in2se5_vca.json"


def full_hdf5_check(dec, temps):
    """Decompose the full per-mode file of job 5 and compare with the single-point route."""
    import h5py
    sys.path.insert(0, str(REPO / "test" / "MnIn2Se4"))
    import kappa_zz_decompose as kd

    rec = kd.reciprocal_rows(str(REPO / "test/Mn2In2Se5/fit/POSCAR"), [1, 0, 0, 0, 1, 0, 0, 0, 1])
    with h5py.File(H5) as f:
        assert [float(t) for t in f["temperature"][()]] == temps, "hdf5 temperatures differ"
    full = [kd.decompose(str(H5), [], rec, [], temp_index=i)["intrinsic"] for i in range(len(temps))]
    dev = max(abs(a[k] - b[k]) for a, b in zip(full, dec)
              for k in ("kappa_zz_first_gamma_z_point", "kappa_zz_rest", "inplane_total"))
    print(f"full hdf5 vs single-point route: max abs deviation {dev:.1e} W/mK (JSON totals have 5 decimals)")
    assert dev < 2e-5, "full decomposition disagrees with first_gz_point.json"

    i300 = temps.index(300.0)
    out = HERE / "kappa_vca_300K.hdf5"
    with h5py.File(H5) as f, h5py.File(out, "w") as g:
        nT = len(temps)
        for k in f:
            d = f[k][()]
            if isinstance(d, np.ndarray) and d.ndim >= 1 and d.shape[0] == nT and k in (
                    "gamma", "heat_capacity", "kappa", "mode_kappa", "temperature"):
                d = d[i300:i300 + 1]
            g.create_dataset(k, data=d)
    shutil.copy(out, REPO / "test/figures/scattering/mn2in2se5/kappa_vca.hdf5")
    print(f"wrote {out.name} (300 K slice) and copied it to test/figures/scattering/mn2in2se5/kappa_vca.hdf5")


def main():
    split, job5 = json.loads(SPLIT.read_text()), json.loads(JOB5.read_text())
    temps = [float(t) for t in job5["T_K"]]
    dec = split["rows"]
    assert [r["T_K"] for r in dec] == temps
    (HERE / "decomposition.json").write_text(json.dumps(
        {"source": "first_gz_point.json (job 5 totals minus the first Gamma-Z point computed here)",
         "x": 0.19, "mesh": split["mesh"], "rows": dec}, indent=1))

    old = json.loads(OLD.read_text())
    assert [float(t) for t in old["T_K"]] == temps, "baseline temperatures differ"
    v = job5["vca"]["site"]
    new = {
        "material": "Mn2In2Se5", "mesh": "20x20x20", "fc2": "4x4x1", "fc3_cutoff_A": 6.51,
        "T_K": old["T_K"],
        "note": ("VCA comparison on a single mesh (20^3). baseline_* = ordered masses + natural isotopes "
                 "(HanHai job 424496, unchanged); vca_* = averaged site masses + mass-variance scattering at "
                 "the literature x = 0.19 (Williams et al. 2024, ESI Table S7), phono3py 4.4.0: in-plane and "
                 "kappa_zz total from HanHai job 5, kappa_zz = rest of the zone = total minus the first "
                 "Gamma-Z grid point, computed with the same driver (test/Mn2In2Se5/vca_m20_x019/). "
                 "The x = 0.13 version is kept as data_mn2in2se5_vca_x013.json."),
        "vca": {"x": 0.19, "x_source": job5["vca"]["x_source"],
                "mass": {"Mn": v["Mn"]["mass"], "In": v["In"]["mass"], "Se": 78.96},
                "g": {"Mn": v["Mn"]["g"], "In": v["In"]["g"], "Se": 0.000463}},
        "baseline_inplane": old["baseline_inplane"],
        "baseline_kzz_rest": old["baseline_kzz_rest"],
        "vca_inplane": [round(r["inplane_total"], 4) for r in dec],
        "vca_kzz_rest": [round(r["kappa_zz_rest"], 4) for r in dec],
    }
    NEW.write_text(json.dumps(new, indent=1))
    print(f"wrote {NEW.relative_to(REPO)}")
    i = temps.index(300.0)
    bi, bz = old["baseline_inplane"][i], old["baseline_kzz_rest"][i]
    print(f"300 K: in-plane {new['vca_inplane'][i]} (baseline {bi}, {bi / new['vca_inplane'][i]:.2f}x), "
          f"kzz rest {new['vca_kzz_rest'][i]} (baseline {bz}, {bz / new['vca_kzz_rest'][i]:.2f}x)")

    if H5.exists():
        full_hdf5_check(dec, temps)
    else:
        print(f"{H5.name} not here: scattering figure left at x = 0.13")


if __name__ == "__main__":
    main()
