#!/usr/bin/env python3
"""Mn2In2Se5 RTA kappa with ISOTOPE SCATTERING ON.

The project configuration (Mn2In2Se5/config/step6_kappa.step.conf) sets
ISOTOPE = true, and every MnIn2Se4 run in this study used isotopes on. All the
Mn2In2Se5 runs made so far used is_isotope=False, which puts the two materials on
different footings. This driver redoes them with isotopes on, same force
constants, same primitive_matrix="P".

  run.py --mesh 16 16 16 --fc3 in/fc3.hdf5 --temps 300
  run.py --mesh 16 16 16 --fc3 in/fc3_cut5p68.hdf5 --temps 300

Writes kappa-m<mesh>.hdf5 and result.json into the current directory.
"""
import argparse
import json
import os
import time

import numpy as np
import phono3py
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

BASE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--mesh", nargs=3, type=int, default=[16, 16, 16])
ap.add_argument("--fc3", default="in/fc3.hdf5")
ap.add_argument("--temps", nargs="+", type=float, default=[300.0])
ap.add_argument("--tag", default="")
a = ap.parse_args()

p = phono3py.load(f"{BASE}/in/phono3py_disp.yaml", produce_fc=False, log_level=0)
ph3 = phono3py.Phono3py(p.unitcell, supercell_matrix=np.array(p.supercell_matrix, int),
                        primitive_matrix="P", log_level=0)
ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/{a.fc3}")
ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2.hdf5")
ph3.mesh_numbers = a.mesh
ph3.init_phph_interaction()
print("start mesh=%s fc3=%s isotope=True T=%s" % (a.mesh, a.fc3, a.temps), flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=a.temps, is_isotope=True,
                             write_kappa=True, log_level=0)
k = np.asarray(ph3.thermal_conductivity.kappa).reshape(-1, 6)
out = {"material": "Mn2In2Se5", "mesh": a.mesh, "fc3": a.fc3,
       "isotope": True, "primitive_matrix": "P", "tag": a.tag,
       "T_K": list(a.temps),
       "kappa_inplane": [round(float(0.5 * (r[0] + r[1])), 4) for r in k],
       "kappa_zz": [round(float(r[2]), 4) for r in k],
       "seconds": round(time.time() - t0, 1)}
json.dump(out, open("result.json", "w"), indent=1)
for T, ip, zz in zip(a.temps, out["kappa_inplane"], out["kappa_zz"]):
    print("  T=%.0f in-plane %.4f  kappa_zz %.4f" % (T, ip, zz), flush=True)
print("DONE", flush=True)
