#!/usr/bin/env python3
"""kappa_zz rest term of the 20^3 Mn2In2Se5 VCA curve at x = 0.19, without the full 20^3 run.

kappa_zz(total) = kappa_zz(first Gamma-Z point) + kappa_zz(rest)
(test/MnIn2Se4/kappa_zz_decompose.py). HanHai job 5 already has the 20^3 totals at
x = 0.19 for all 12 temperatures (run_lbte_wigner.py, phono3py 4.4.0), but its per-mode
hdf5 is not in the repository, and the full 20^3 RTA run takes several hours on 4 cores.
The first Gamma-Z point is one irreducible grid point, q = (0, 0, 1/20), BZ grid index 409,
weight 2. This script computes only that grid point with the same driver (same FC, masses,
mass variances, mesh), so

    rest(T) = kzz_total(job 5, T) - first_point(T)

with no approximation: in the RTA every grid point is computed independently and its
mode_kappa does not depend on which other points are run.

Validation: the single-point run with the inputs of the full 20^3 x = 0.13 run
(HanHai job 424147, test/Mn2In2Se5/vca_m20/jobs/run.py: rounded masses 62.73 / 107.03,
g 0.1031 / 0.0354) must reproduce the stored mode_kappa, gamma and frequencies of that grid
point (test/figures/scattering/mn2in2se5/kappa_vca.hdf5, first point 0.432621). It does, to
1e-15. With the unrounded x = 0.13 masses of site_mix() the first point comes out 0.44 %
higher: the near-Gamma acoustic modes are that sensitive to the masses, so the x = 0.19 point
must be (and is) computed with exactly the masses of job 5, which the script checks.

  python first_gz_point.py            # ~ 10 s on 4 cores; writes first_gz_point.json
"""
import importlib.util
import json
import time
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DRIVER = REPO / "test/figures/kappa_vs_T/jobs_lbte_wigner/run_lbte_wigner.py"
JOB5 = DRIVER.parent / "run_Mn2In2Se5_vca_m20_rta/result_Mn2In2Se5_vca_m20.json"
REF013 = REPO / "test/figures/scattering/mn2in2se5/kappa_vca.hdf5"
MESH, GP, N_BZ = 20, 409, 20 ** 3

spec = importlib.util.spec_from_file_location("run_lbte_wigner", DRIVER)
drv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(drv)


def first_point(x, temps, mix=None):
    """mode_kappa_zz and gamma of grid point GP, VCA at antisite fraction x
    (or with explicit site masses / mass variances mix = {element: (mass, g)})."""
    site_mix = drv.site_mix
    if mix is not None:
        drv.site_mix = lambda material, x_: mix
    try:
        ph3, gl, masses = drv.build(drv.CFG["Mn2In2Se5"], "vca", "Mn2In2Se5", x, MESH)
    finally:
        drv.site_mix = site_mix
    addr = np.asarray(ph3.grid.addresses[GP], float) / MESH
    assert np.allclose(addr, [0, 0, 0.05]), f"grid point {GP} is {addr}, not (0, 0, 1/20)"
    t0 = time.time()
    ph3.run_thermal_conductivity(temperatures=temps, mass_variances=gl, grid_points=[GP], log_level=0)
    tc = ph3.thermal_conductivity
    mk = np.asarray(tc.mode_kappa)
    g = np.asarray(tc.gamma)
    if mk.ndim == 5:                       # (sigma, T, gp, band, 6)
        mk, g = mk[0], g[0]
    assert mk.shape[1] == 1 and int(tc.grid_weights[0]) == 2
    kzz_band = mk[:, 0, :, 2] / N_BZ       # (T, band), same normalisation as kappa_zz_decompose.py
    return {"x": x, "masses": masses, "mass_variances": [round(v, 6) for v in gl],
            "seconds": round(time.time() - t0, 1),
            "kzz_band": kzz_band, "gamma": g[:, 0, :], "frequency": np.asarray(tc.frequencies)[0]}


def main():
    # 1. validation against the full x = 0.13 run, with that run's (rounded) inputs
    with h5py.File(REF013) as f:
        i = int(np.where(f["grid_point"][()] == GP)[0][0])
        ref_kzz = f["mode_kappa"][0, i, :, 2] / N_BZ
        ref_g = f["gamma"][0, i]
        ref_f = f["frequency"][i]
        assert float(f["temperature"][0]) == 300.0
    v = first_point(0.13, [300.0], mix={"Mn": (62.73, 0.1031), "In": (107.03, 0.0354)})
    rel = abs(v["kzz_band"][0].sum() / ref_kzz.sum() - 1)
    m = ref_g > 1e-6
    dg = float(np.max(np.abs(v["gamma"][0][m] / ref_g[m] - 1)))
    df = float(np.max(np.abs(v["frequency"] - ref_f)))
    print(f"x = 0.13 (job 424147 inputs), 300 K, gp {GP}: first point {v['kzz_band'][0].sum():.6f} "
          f"(full run {ref_kzz.sum():.6f}, rel. dev. {rel:.1e}); max rel. d gamma {dg:.1e}; "
          f"max d freq {df:.1e} THz  [{v['seconds']} s]")
    assert rel < 1e-10 and dg < 1e-10 and df < 1e-10, "single-point run does not reproduce the full run"

    # 2. x = 0.19, the 12 temperatures of job 5
    job5 = json.loads(JOB5.read_text())
    temps = job5["T_K"]
    r = first_point(0.19, temps)
    assert r["masses"] == job5["masses"] and r["mass_variances"] == job5["mass_variances"], \
        "masses / mass variances differ from job 5"
    tot = job5["rta"]["population"]["kzz_total"]
    first = r["kzz_band"].sum(axis=1)
    rows = [{"T_K": T, "kappa_zz_total": t, "kappa_zz_first_gamma_z_point": round(float(f1), 6),
             "kappa_zz_rest": round(t - float(f1), 5), "inplane_total": ip}
            for T, t, f1, ip in zip(temps, tot, first, job5["rta"]["population"]["inplane"])]
    out = {"what": "kappa_zz split at x = 0.19: totals from HanHai job 5, first Gamma-Z point "
                   f"(BZ grid point {GP}, q = (0, 0, 1/20), weight 2) computed here",
           "x": 0.19, "mesh": [MESH] * 3, "phono3py": job5["phono3py"],
           "validation_x013_300K": {"first_point": round(float(v["kzz_band"][0].sum()), 6),
                                    "full_run": round(float(ref_kzz.sum()), 6),
                                    "rel_dev": float(f"{rel:.2e}"), "max_rel_dgamma": float(f"{dg:.2e}")},
           "seconds": r["seconds"], "rows": rows}
    (HERE / "first_gz_point.json").write_text(json.dumps(out, indent=1))
    for row in rows:
        print(f"  {row['T_K']:6.0f} K  in-plane {row['inplane_total']:8.4f}  kzz total {row['kappa_zz_total']:8.4f}"
              f"  first G-Z {row['kappa_zz_first_gamma_z_point']:7.4f}  rest {row['kappa_zz_rest']:8.4f}")
    print("wrote first_gz_point.json")


if __name__ == "__main__":
    main()
