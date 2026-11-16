#!/usr/bin/env python3
"""phono3py RTA kappa for Mn2In2Se5 at a chosen q-mesh AND a chosen FC3 cutoff.

Same physics as ../mn2in2se5_kzz_fcset/run_kzz.py (primitive_matrix="P",
isotopes off, 300 K, write_kappa) but the FC3 file is an argument, so the
cutoff scan of kappa_zz can be redone on an isotropic mesh.

Writes kappa-m<mesh>.hdf5 in the current directory plus result.json.
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
ap.add_argument("--mesh", nargs=3, type=int, required=True)
ap.add_argument("--fc3", required=True, help="file name inside in/")
ap.add_argument("--out", required=True)
a = ap.parse_args()

p = phono3py.load(f"{BASE}/in/phono3py_disp.yaml", produce_fc=False, log_level=0)
ph3 = phono3py.Phono3py(p.unitcell,
                        supercell_matrix=np.array(p.supercell_matrix, int),
                        primitive_matrix="P", log_level=0)
ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/in/{a.fc3}")
ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2.hdf5")
ph3.mesh_numbers = a.mesh
ph3.init_phph_interaction()
print("start mesh=%s fc3=%s" % (a.mesh, a.fc3), flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_isotope=False,
                             write_kappa=True, log_level=0)
dt = time.time() - t0
k = np.asarray(ph3.thermal_conductivity.kappa)[0, 0]
ip = float(0.5 * (k[0] + k[1]))
out = {"material": "Mn2In2Se5", "method": "RTA", "mesh": " ".join(map(str, a.mesh)),
       "fc3": a.fc3, "T_K": [300.0], "isotope": False, "primitive_matrix": "P",
       "inplane": ip, "kappa_xx": float(k[0]), "kappa_zz": float(k[2]),
       "seconds": round(dt, 1), "phono3py": phono3py.__version__}
print("mesh=%s fc3=%s in-plane=%.4f kzz=%.4f (%.0fs)"
      % (a.mesh, a.fc3, ip, k[2], dt), flush=True)
json.dump(out, open(a.out, "w"), indent=1)
print("DONE")
