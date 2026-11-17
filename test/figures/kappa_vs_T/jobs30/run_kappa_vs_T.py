#!/usr/bin/env python3
"""kappa(T) at a chosen mesh, ordered and VCA, kappa_zz as the plain total.

One material per run; writes result_<MAT>_m<mesh>.json (schema consumed by
../plot_theory_vs_exp.py) and the per-mode phono3py hdf5 of each case.

  run_kappa_vs_T.py --material MnIn2Se4 --mesh 30 --cases ordered vca

Convention (see README.md here): kappa_zz is the phono3py total on the mesh, i.e.
no Gamma-Z decomposition. At 300 K the 30^3 total of MnIn2Se4 is 0.765, about
15 % above the mesh-extrapolated 0.646.

VCA: site masses averaged over the Mn/In antisite occupancy, disorder returned as
mass-variance (Tamura) scattering. Parameters are those of the 20^3 runs
(MnIn2Se4 i = 0.44, Mn2In2Se5 x = 0.13); their sources are listed in
analysis/disorder/README.md section 5.
"""
import argparse
import json
import os
import time

import numpy as np
import phono3py
from phonopy.interface.vasp import read_vasp
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))   # repository root
TEMPS = [30.0, 60.0, 90.0, 100.0, 120.0, 150.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]
G_SE = 4.63e-4

CFG = {
    "MnIn2Se4": {
        "poscar": "test/MnIn2Se4/fit/POSCAR",
        "fc2": "test/MnIn2Se4/fc2_442_phonon/fc2.hdf5",
        "fc3": "test/MnIn2Se4/fit_cut3_6p45/fc3.hdf5",
        "sc": [[4, 0, 0], [0, 4, 0], [0, 0, 1]],
        "psc": [[4, 0, 0], [0, 4, 0], [0, 0, 2]],
        "pmat": [[1, 1, 2], [-1, -2, -1], [1, 1, 1]],
        "mass": {"Mn": 81.29, "In": 101.65},
        "g": {"Mn": 0.1337, "In": 0.0595, "Se": G_SE},
        "x": 0.44, "cut": 6.45,
    },
    "Mn2In2Se5": {
        "yaml": "Mn2In2Se5/fc/phono3py_disp.yaml",
        "fc2": "Mn2In2Se5/fc/fc2.hdf5",
        "fc3": "Mn2In2Se5/fc/fc3.hdf5",
        "pmat": "P",
        "mass": {"Mn": 62.73, "In": 107.03},
        "g": {"Mn": 0.1031, "In": 0.0354, "Se": G_SE},
        "x": 0.13, "cut": 6.51,
    },
}


def load_cell(cfg):
    """returns unitcell, supercell_matrix, phonon_supercell_matrix (or None)"""
    if "poscar" in cfg:
        return (read_vasp(os.path.join(ROOT, cfg["poscar"])),
                np.array(cfg["sc"], int), np.array(cfg["psc"], int))
    p = phono3py.load(os.path.join(ROOT, cfg["yaml"]), produce_fc=False, log_level=0)
    return p.unitcell, np.array(p.supercell_matrix, int), None


def build(cfg, unitcell, sc, psc, mesh):
    kw = dict(primitive_matrix=cfg["pmat"], log_level=0)
    kw["supercell_matrix"] = sc
    if psc is not None:
        kw["phonon_supercell_matrix"] = psc
    ph3 = phono3py.Phono3py(unitcell, **kw)
    ph3.fc3 = read_fc3_from_hdf5(os.path.join(ROOT, cfg["fc3"]))
    ph3.fc2 = read_fc2_from_hdf5(os.path.join(ROOT, cfg["fc2"]))
    ph3.mesh_numbers = [mesh, mesh, mesh]
    return ph3


def take(k, T):
    k = np.asarray(k).reshape(len(T), 6)
    return [round(float(0.5 * (r[0] + r[1])), 4) for r in k], [round(float(r[2]), 4) for r in k]


ap = argparse.ArgumentParser()
ap.add_argument("--material", required=True, choices=sorted(CFG))
ap.add_argument("--mesh", type=int, default=30)
ap.add_argument("--cases", nargs="+", default=["ordered", "vca"], choices=["ordered", "vca"])
ap.add_argument("--temps", nargs="+", type=float, default=TEMPS)
a = ap.parse_args()
cfg = CFG[a.material]
out = {"material": a.material, "mesh": a.mesh, "fc3_cutoff_A": cfg["cut"],
       "T_K": a.temps, "convention": "kappa_zz = plain mesh total (no Gamma-Z decomposition)",
       "vca": {"x": cfg["x"], "mass": cfg["mass"], "g": cfg["g"]}}

for case in a.cases:
    uc, sc, psc = load_cell(cfg)
    if case == "vca":
        masses = np.array(uc.masses, dtype=float)
        for i, s in enumerate(uc.symbols):
            if s in cfg["mass"]:
                masses[i] = cfg["mass"][s]
        uc.masses = masses
        gl = [cfg["g"][s] for s in uc.symbols]
        print("VCA primitive order", list(uc.symbols), flush=True)
        print("VCA mass variances ", gl, flush=True)
    print("start %s %s %s^3 %d T" % (a.material, case, a.mesh, len(a.temps)), flush=True)
    ph3 = build(cfg, uc, sc, psc, a.mesh)
    ph3.init_phph_interaction()
    t0 = time.time()
    if case == "ordered":
        ph3.run_thermal_conductivity(temperatures=a.temps, is_isotope=True,
                                     write_kappa=True, log_level=0)
    else:
        ph3.run_thermal_conductivity(temperatures=a.temps, is_isotope=False,
                                     mass_variances=gl, write_kappa=True, log_level=0)
    ip, zz = take(ph3.thermal_conductivity.kappa, a.temps)
    out[case] = {"inplane": ip, "kappa_zz_total": zz, "seconds": round(time.time() - t0, 1)}
    src = "kappa-m{0}{0}{0}.hdf5".format(a.mesh)
    if os.path.exists(src):
        os.replace(src, "kappa_%s_m%d.hdf5" % (case, a.mesh))
    # written into the working directory (jobs30/run_<MAT>/), which is where
    # ../plot_theory_vs_exp.py looks for it
    fn = os.path.join(os.getcwd(), "result_%s_m%d.json" % (a.material, a.mesh))
    json.dump(out, open(fn, "w"), indent=1)
    print("%s done %s s -> %s" % (case, out[case]["seconds"], fn), flush=True)
print("DONE", flush=True)
