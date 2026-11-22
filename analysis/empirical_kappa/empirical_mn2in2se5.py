#!/usr/bin/env python3
"""Empirical kappa models for Mn2In2Se5 from the calculated phonons, against the hot-pressed experiment.

Two kinds of model, both fed only with quantities from the production force constants
(Mn2In2Se5/fc/fc2.hdf5, fc3.hdf5) and the experimental lattice (Williams et al. 2024):

  intrinsic estimate  Slack (Morelli & Slack 2006 form):
                        kappa = A(g) M theta_a^3 delta n^(1/3) / (g^2 T),
                        A(g) = 2.43e-6 / (1 - 0.514/g + 0.228/g^2)   [W/mK, amu, K, A]
                      theta_a = theta_D / n^(1/3), theta_D from the Debye mean of the acoustic velocities
                      (3000 directions), g = mode Grueneisen parameter from FC3 (phono3py, volumetric
                      = trace/3), heat-capacity weighted at 300 K on a 10^3 mesh.
  lower limits        Cahill, Watson & Pohl 1992 (full integral at T, and its high-T limit
                      0.5 (pi/6)^(1/3) kB n^(2/3) (vL + 2 vT)); diffuson, Agne et al. 2018,
                      0.76 kB n^(2/3) (vL + 2 vT) / 3.

Writes empirical_mn2in2se5.json.   python empirical_mn2in2se5.py   (needs phono3py)
"""
import json
from pathlib import Path

import numpy as np
import phono3py
import phonopy
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5
from phono3py.phonon3.gruneisen import Gruneisen

KB, HBAR = 1.380649e-23, 1.054571817e-34
ROOT = Path(__file__).resolve().parents[2]
FC = ROOT / "Mn2In2Se5" / "fc"
A_EXP, C_EXP = 4.0189e-10, 48.830e-10          # Williams et al. 2024, 3 f.u. per hexagonal cell
EXP_300 = {"Perp": 0.741, "Para": 0.320}        # 303 K, experiment/kappa_tot_digitized.csv
EXP_623 = {"Perp": 0.590, "Para": 0.264}


def main():
    p3 = phono3py.load(str(FC / "phono3py_disp.yaml"), produce_fc=False, log_level=0)
    fc2 = read_fc2_from_hdf5(str(FC / "fc2.hdf5"))
    fc3 = read_fc3_from_hdf5(str(FC / "fc3.hdf5"))

    # acoustic velocities (q in the basis of the phonopy primitive cell used here)
    ph = phonopy.Phonopy(p3.unitcell, supercell_matrix=p3.supercell_matrix, primitive_matrix="P")
    ph.force_constants = fc2
    lat, dq = np.array(ph.primitive.cell), 0.004

    def vel(dirs):
        dirs = np.asarray(dirs, float)
        dirs = dirs / np.linalg.norm(dirs, axis=1)[:, None]
        ph.run_qpoints((dirs * dq) @ lat.T)
        return np.sort(ph.qpoints.frequencies[:, :3], axis=1) * 1e12 / (dq * 1e10)

    v = vel(np.random.default_rng(0).normal(size=(3000, 3)))
    vL, vT = float(v[:, 2].mean()), float(v[:, :2].mean())
    vD = float((np.mean((1 / v ** 3).sum(axis=1)) / 3) ** (-1 / 3))

    vhex = np.sqrt(3) / 2 * A_EXP ** 2 * C_EXP
    n = 27 / vhex
    nat = len(p3.primitive.masses)
    mbar = float(np.mean(p3.primitive.masses))
    delta = (vhex / 27) ** (1 / 3)

    def cahill(T):
        s = 0.0
        for vi in (vL, vT, vT):
            th = vi * HBAR / KB * (6 * np.pi ** 2 * n) ** (1 / 3)
            x = np.linspace(1e-6, th / T, 4000)
            s += vi * (T / th) ** 2 * np.trapezoid(x ** 3 * np.exp(x) / (np.exp(x) - 1) ** 2, x)
        return float((np.pi / 6) ** (1 / 3) * KB * n ** (2 / 3) * s)

    k_cah_hi = float(0.5 * (np.pi / 6) ** (1 / 3) * KB * n ** (2 / 3) * (vL + 2 * vT))
    k_dif = float(0.76 * KB * n ** (2 / 3) * (vL + 2 * vT) / 3)

    # mode Grueneisen parameters
    g = Gruneisen(fc2, fc3, p3.supercell, p3.primitive)
    g.set_sampling_mesh([10, 10, 10], primitive_symmetry=p3.primitive_symmetry)
    g.run()
    gam = np.trace(np.array(g.gruneisen_parameters), axis1=2, axis2=3) / 3
    fr, w = np.array(g.frequencies), np.array(g._weights, float)
    x = np.clip(fr, 1e-6, None) * 1e12 * 2 * np.pi * HBAR / (KB * 300.0)
    cv = np.where(fr > 0.05, x ** 2 * np.exp(x) / (np.exp(x) - 1) ** 2, 0.0) * w[:, None]
    g_ac = float((cv[:, :3] * gam[:, :3]).sum() / cv[:, :3].sum())
    g_ac_rms = float(np.sqrt((cv[:, :3] * gam[:, :3] ** 2).sum() / cv[:, :3].sum()))
    g_all = float((cv * gam).sum() / cv.sum())

    th_d = HBAR * vD / KB * (6 * np.pi ** 2 * n) ** (1 / 3)
    th_a = th_d / nat ** (1 / 3)

    def slack(gm, T):
        a = 2.43e-6 / (1 - 0.514 / gm + 0.228 / gm ** 2)
        return float(a * mbar * th_a ** 3 * delta * 1e10 * nat ** (1 / 3) / (gm ** 2 * T))

    g_needed = float(next(gm for gm in np.arange(0.5, 6, 0.001)
                          if slack(gm, 300) <= (2 * EXP_300["Perp"] + EXP_300["Para"]) / 3))
    exp_avg = {T: (2 * e["Perp"] + e["Para"]) / 3 for T, e in ((300, EXP_300), (623, EXP_623))}
    out = {
        "velocities_m_s": {"vL_avg": round(vL), "vT_avg": round(vT), "v_debye": round(vD)},
        "number_density_A-3": round(n * 1e-30, 5), "mean_mass_amu": round(mbar, 2), "delta_A": round(delta * 1e10, 3),
        "theta_D_K": round(th_d, 1), "theta_a_K": round(th_a, 1),
        "gruneisen_300K": {"acoustic_mean": round(g_ac, 3), "acoustic_rms": round(g_ac_rms, 3),
                           "all_modes_mean": round(g_all, 3)},
        "slack_W_mK": {f"gamma={gm:.2f} ({lab})": {"300K": round(slack(gm, 300), 2), "623K": round(slack(gm, 623), 2)}
                       for gm, lab in ((g_all, "all-mode mean"), (g_ac_rms, "acoustic rms"))},
        "slack_gamma_needed_to_match_experiment_300K": round(g_needed, 2),
        "cahill_W_mK": {"300K": round(cahill(300), 3), "623K": round(cahill(623), 3), "high_T_limit": round(k_cah_hi, 3)},
        "diffuson_W_mK": round(k_dif, 3),
        "experiment_W_mK": {"303K": EXP_300, "623K": EXP_623,
                            "direction_average_303K": round(exp_avg[300], 3), "direction_average_623K": round(exp_avg[623], 3)},
        "first_principles_direction_average_300K": {"ordered_quoted_LBTE_coherence": 2.61, "ordered_RTA_20": 2.33,
                                                    "VCA_x019_RTA_20": 1.23},
    }
    (Path(__file__).with_name("empirical_mn2in2se5.json")).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
