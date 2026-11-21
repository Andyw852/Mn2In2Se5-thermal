#!/usr/bin/env python3
"""Post-process the Mn2In2Se5 VCA run at the literature x = 0.19 (20^3, RTA, 12 T).

Input  (this directory):
  result_Mn2In2Se5_vca_m20.json   written by run_lbte_wigner.py
  kappa_vca_rta_m20.hdf5          per-mode output of the same run (12 temperatures)

Does three things:
  1. checks the run against HanHai job 5 (same driver, same settings), whose JSON is
     test/figures/kappa_vs_T/jobs_lbte_wigner/run_Mn2In2Se5_vca_m20_rta/;
  2. splits kappa_zz at every temperature with test/MnIn2Se4/kappa_zz_decompose.py
     (Gamma-Z first point removed; Mn2In2Se5 cell: --poscar test/Mn2In2Se5/fit/POSCAR,
     identity primitive matrix) and writes decomposition.json;
  3. writes test/figures/kappa_vs_T/data_mn2in2se5_vca.json at x = 0.19 (baseline_* taken
     unchanged from the x = 0.13 file of job 424496, same mesh, same temperatures) and a
     300 K slice kappa_vca_300K.hdf5 for the scattering-rate figure.

  python postprocess.py
"""
import json
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "test" / "MnIn2Se4"))
import kappa_zz_decompose as kd  # noqa: E402

RUN = HERE / "result_Mn2In2Se5_vca_m20.json"
H5 = HERE / "kappa_vca_rta_m20.hdf5"
JOB5 = REPO / "test/figures/kappa_vs_T/jobs_lbte_wigner/run_Mn2In2Se5_vca_m20_rta/result_Mn2In2Se5_vca_m20.json"
OLD = REPO / "test/figures/kappa_vs_T/data_mn2in2se5_vca_x013.json"
NEW = REPO / "test/figures/kappa_vs_T/data_mn2in2se5_vca.json"


def main():
    run, job5 = json.loads(RUN.read_text()), json.loads(JOB5.read_text())
    assert run["T_K"] == job5["T_K"], "temperature grids differ"
    a, b = run["rta"]["population"], job5["rta"]["population"]
    dev = max(abs(x / y - 1) for c in ("inplane", "kzz_total") for x, y in zip(a[c], b[c]))
    print(f"check against HanHai job 5: max relative deviation over 12 T, in-plane and kzz total = {dev:.2e}")

    rec = kd.reciprocal_rows(str(REPO / "test/Mn2In2Se5/fit/POSCAR"), [1, 0, 0, 0, 1, 0, 0, 0, 1])
    with h5py.File(H5) as f:
        temps = [float(t) for t in f["temperature"][()]]
    rows = [kd.decompose(str(H5), [], rec, [], temp_index=i) for i in range(len(temps))]
    dec = [{"T_K": r["T_K"], **r["intrinsic"]} for r in rows]
    (HERE / "decomposition.json").write_text(json.dumps(
        {"file": H5.name, "x": 0.19, "mesh": rows[0]["mesh"], "rows": dec}, indent=1))
    for r in dec:
        print(f"  {r['T_K']:6.0f} K  in-plane {r['inplane_total']:8.4f}  kzz total {r['kappa_zz_total']:8.4f}"
              f"  first G-Z {r['kappa_zz_first_gamma_z_point']:7.4f}  rest {r['kappa_zz_rest']:8.4f}")

    old = json.loads(OLD.read_text())
    assert [float(t) for t in old["T_K"]] == temps, "baseline temperatures differ"
    v = run["vca"]["site"]
    new = {
        "material": "Mn2In2Se5", "mesh": "20x20x20", "fc2": "4x4x1", "fc3_cutoff_A": 6.51,
        "T_K": old["T_K"],
        "note": ("VCA comparison on a single mesh (20^3). baseline_* = ordered masses + natural isotopes "
                 "(HanHai job 424496, unchanged); vca_* = averaged site masses + mass-variance scattering at "
                 "the literature x = 0.19 (Williams et al. 2024, ESI Table S7), phono3py 4.4.0, same driver and "
                 "settings as HanHai job 5, rerun 2026-11-21 for the per-mode output "
                 "(test/Mn2In2Se5/vca_m20_x019/). kappa_zz = rest of the zone at every T. The x = 0.13 version "
                 "is kept as data_mn2in2se5_vca_x013.json."),
        "vca": {"x": 0.19, "x_source": run["vca"]["x_source"],
                "mass": {"Mn": v["Mn"]["mass"], "In": v["In"]["mass"], "Se": 78.96},
                "g": {"Mn": v["Mn"]["g"], "In": v["In"]["g"], "Se": 0.000463}},
        "baseline_inplane": old["baseline_inplane"],
        "baseline_kzz_rest": old["baseline_kzz_rest"],
        "vca_inplane": [round(r["inplane_total"], 4) for r in dec],
        "vca_kzz_rest": [round(r["kappa_zz_rest"], 4) for r in dec],
    }
    NEW.write_text(json.dumps(new, indent=1))
    print(f"wrote {NEW.relative_to(REPO)}")

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


if __name__ == "__main__":
    main()
