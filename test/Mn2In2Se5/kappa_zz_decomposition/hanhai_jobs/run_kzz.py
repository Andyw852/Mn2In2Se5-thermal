#!/usr/bin/env python3
"""phono3py RTA kappa for Mn2In2Se5 at a chosen q-mesh (HanHai job entry point).

Reproduces the setup of the committed reference mesh scan
(analysis/extrinsic/kappa_extrinsic.py, cmd_bases):

  cell + supercell matrix : in/phono3py_disp.yaml (4x4x1)
  primitive_matrix        : "P"    <- NOT the fit_config matrix: with
                             [[2,1,1],[-1,-1,-2],[1,1,1]] a 16x16x4 mesh breaks
                             the grid symmetry and phono3py refuses to run
  is_isotope              : False  (the reference numbers 0.836 / 1.109 / 1.304 /
                             1.541 are without isotope scattering)
  write_kappa             : True   (per-mode kappa, for kappa_zz_decompose.py)

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
ap.add_argument("--out", required=True)
ap.add_argument("--disp-yaml", default=f"{BASE}/in/phono3py_disp.yaml")
ap.add_argument("--isotope", action="store_true")
a = ap.parse_args()

p = phono3py.load(os.path.abspath(a.disp_yaml), produce_fc=False, log_level=0)
smat = np.array(p.supercell_matrix, int)
ph3 = phono3py.Phono3py(p.unitcell, supercell_matrix=smat,
                        primitive_matrix="P", log_level=0)
ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/in/fc3.hdf5")
ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2.hdf5")
ph3.mesh_numbers = a.mesh
ph3.init_phph_interaction()
print("start mesh=%s supercell=%d atoms isotope=%s"
      % (a.mesh, len(ph3.supercell), a.isotope), flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(
    temperatures=[300.0], is_isotope=a.isotope,
    write_kappa=True, log_level=0,
)
dt = time.time() - t0
k = np.asarray(ph3.thermal_conductivity.kappa)[0, 0]
ip = float(0.5 * (k[0] + k[1]))
out = {"material": "Mn2In2Se5", "method": "RTA",
       "mesh": " ".join(map(str, a.mesh)), "T_K": [300.0],
       "isotope": bool(a.isotope), "primitive_matrix": "P",
       "inplane": ip, "kappa_xx": float(k[0]), "kappa_zz": float(k[2]),
       "seconds": round(dt, 1), "phono3py": phono3py.__version__}
print("mesh=%s in-plane=%.4f kzz=%.4f (%.0fs)" % (a.mesh, ip, k[2], dt), flush=True)
json.dump(out, open(a.out, "w"), indent=1)
print("DONE")
