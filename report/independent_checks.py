#!/usr/bin/env python3
"""Checks of the MnIn2Se4 calculation against experimental data that were NOT used in any fit.

  1  lattice constants (hexagonal setting) of the relaxed structure vs single-crystal XRD
  2  Gamma-point Raman-active modes of the 4x4x2 FC2 vs the measured Raman spectrum
  3  acoustic velocities of the 4x4x2 FC2 (along c, in-plane, averaged over 3000 directions)
     vs the ultrasonic velocities of the hot-pressed samples
  4  Cahill and diffuson minimum kappa (high-temperature limits) from the measured velocities

Writes report/independent_checks.json, read by make_report_html.py.
Needs numpy + phonopy.   python3 report/independent_checks.py
"""
import csv
import json
from pathlib import Path

import numpy as np
import phonopy

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "test" / "MnIn2Se4" / "fc2_442_phonon"
KB = 1.380649e-23
THZ_TO_CM = 33.35641

# single-crystal XRD (experimental slide, 2026-11-19); Raman peaks read off the slide image, +-5 cm-1
EXP_LATTICE = {"a": 4.05, "c": 39.35}
EXP_RAMAN_APPROX = [78, 177]
FORMULA_MASS = 54.938045 + 2 * 114.818 + 4 * 78.96     # amu, 600.41


def main():
    ph = phonopy.load(str(D / "phonopy_disp.yaml"), force_constants_filename=str(D / "fc2.hdf5"),
                      produce_fc=False, is_nac=False, log_level=0)
    lat = ph.primitive.cell                     # rhombohedral primitive, rows = vectors
    # hexagonal setting: a = |a1 - a2|, c = |a1 + a2 + a3|
    a_hex = float(np.linalg.norm(lat[0] - lat[1]))
    c_hex = float(np.linalg.norm(lat.sum(axis=0)))

    # Gamma irreps -> Raman-active (g) modes
    ph.set_irreps([0, 0, 0], degeneracy_tolerance=1e-3)
    ir = ph.irreps
    modes = [{"irrep": lab, "THz": round(float(ir._freqs[d[0]]), 3),
              "cm-1": round(float(ir._freqs[d[0]]) * THZ_TO_CM, 1)}
             for d, lab in zip(ir._degenerate_sets, ir._ir_labels) if lab]
    raman = [m for m in modes if m["irrep"].endswith("g")]

    # acoustic velocities: slope of the three lowest branches at |q| = 0.004 1/A (no 2 pi)
    dq = 0.004

    def vel(dirs):
        dirs = np.asarray(dirs, float)
        dirs = dirs / np.linalg.norm(dirs, axis=1)[:, None]
        ph.run_qpoints((dirs * dq) @ lat.T)     # q_frac = lattice . q_cart
        f = np.sort(ph.get_qpoints_dict()["frequencies"][:, :3], axis=1)
        return f * 1e12 / (dq * 1e10)            # m/s

    v_c = vel([[0, 0, 1]])[0]
    v_ab = vel([[1, 0, 0]])[0]
    rng = np.random.default_rng(0)
    v = vel(rng.normal(size=(3000, 3)))
    vL, vT = float(v[:, 2].mean()), float(v[:, :2].mean())
    v_debye = float((np.mean((1 / v ** 3).sum(axis=1)) / 3) ** (-1 / 3))

    # experiment
    with open(ROOT / "experiment" / "sound_velocity_MnIn2Se4.csv", encoding="utf-8") as f:
        rows = list(csv.DictReader(l for l in f if not l.startswith("#")))
    exp = [{k: float(r[k]) for k in ("vL_ms", "vT_ms", "v_mean_ms")} | {"sample": r["sample"]} for r in rows]

    # number density and theoretical density from the experimental lattice (3 f.u. per hexagonal cell)
    V = np.sqrt(3) / 2 * EXP_LATTICE["a"] ** 2 * EXP_LATTICE["c"] * 1e-30
    n = 21 / V
    rho = 3 * FORMULA_MASS * 1.66053907e-27 / V

    def kmin(vl, vt):
        # high-temperature limits. Cahill, Watson & Pohl (1992): (pi/6)^(1/3) kB n^(2/3) sum_i v_i (T/Th_i)^2
        # int_0^(Th_i/T) x^3 e^x/(e^x-1)^2 dx; the integral tends to (Th_i/T)^2 / 2, hence the factor 1/2
        # (= 1.21/3, the form quoted by Agne et al. 2018). Until 2026-11-22 the 1/2 was missing (kmin 2x too high).
        cahill = 0.5 * (np.pi / 6) ** (1 / 3) * KB * n ** (2 / 3) * (vl + 2 * vt)
        diffuson = 0.76 * KB * n ** (2 / 3) * (vl + 2 * vt) / 3
        return round(float(cahill), 3), round(float(diffuson), 3)

    vl_e = float(np.mean([e["vL_ms"] for e in exp]))
    vt_e = float(np.mean([e["vT_ms"] for e in exp]))
    out = {
        "lattice": {"calc": {"a": round(a_hex, 3), "c": round(c_hex, 2)}, "exp": EXP_LATTICE},
        "raman_calc": raman, "raman_exp_approx_cm-1": EXP_RAMAN_APPROX, "gamma_modes": modes,
        "velocity_calc": {"c_axis": np.round(v_c).tolist(), "in_plane": np.round(v_ab).tolist(),
                          "avg_vL": round(vL), "avg_vT": round(vT), "debye_mean": round(v_debye)},
        "velocity_exp": exp,
        "number_density_A-3": round(n * 1e-30, 5), "theoretical_density_g_cm3": round(rho / 1000, 3),
        "kmin_exp_velocities": dict(zip(("cahill", "diffuson"), kmin(vl_e, vt_e))),
        "kmin_calc_velocities": dict(zip(("cahill", "diffuson"), kmin(vL, vT))),
    }
    (ROOT / "report" / "independent_checks.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
