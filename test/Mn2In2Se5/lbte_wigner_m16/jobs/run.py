#!/usr/bin/env python3
"""Mn2In2Se5: LBTE and Wigner-LBTE at the largest FC3 cutoff (6.51 A), 16^3.

Same cell / force constants as the committed RTA mesh runs (Mn2In2Se5/fc/,
primitive_matrix="P", isotopes off, T = 300 K), so the iterative solution and
the coherence term end up on the same footing as the rest-term decomposition.

Runs two phono3py solutions in one job:
  1) is_LBTE=True, conductivity_type=None  -> the iterative (LBTE) kappa
  2) is_LBTE=True, conductivity_type="wigner" -> TOT / population / coherence
"""
import json
import os
import time

import numpy as np
import phono3py
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

BASE = os.path.dirname(os.path.abspath(__file__))
MESH = [16, 16, 16]

p = phono3py.load(f"{BASE}/in/phono3py_disp.yaml", produce_fc=False, log_level=0)
ph3 = phono3py.Phono3py(p.unitcell, supercell_matrix=np.array(p.supercell_matrix, int),
                        primitive_matrix="P", log_level=0)
ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/in/fc3.hdf5")
ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2.hdf5")
ph3.mesh_numbers = MESH
ph3.init_phph_interaction()

def snap(k):
    # see the note in ../../MnIn2Se4/lbte_wigner_m14/jobs/run.py: the kappa
    # arrays can be (1, 6) or (1, 1, 6); this crashed after a successful solve.
    a = np.asarray(k).reshape(-1, 6)[0]
    return {"inplane": float(0.5 * (a[0] + a[1])), "kappa_zz": float(a[2])}

out = {"material": "Mn2In2Se5", "mesh": " ".join(map(str, MESH)),
       "fc2": "4x4x1", "fc3_cutoff_A": 6.51, "T_K": 300.0, "isotope": False,
       "primitive_matrix": "P"}

print("start LBTE", flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=True, is_isotope=False,
                             write_kappa=True, log_level=0)
out["LBTE"] = snap(ph3.thermal_conductivity.kappa)
out["LBTE"]["seconds"] = round(time.time() - t0, 1)
print("LBTE", out["LBTE"], flush=True)

print("start Wigner-LBTE", flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=True, is_isotope=False,
                             conductivity_type="wigner", write_kappa=False, log_level=0)
tc = ph3.thermal_conductivity
tot, pop = snap(tc.kappa_TOT_exact), snap(tc.kappa_P_exact)
coh = {"inplane": tot["inplane"] - pop["inplane"],
       "kappa_zz": tot["kappa_zz"] - pop["kappa_zz"]}
out["Wigner_LBTE"] = {"TOT": tot, "population": pop, "coherence": coh,
                      "seconds": round(time.time() - t0, 1)}
print("Wigner-LBTE", json.dumps(out["Wigner_LBTE"]), flush=True)

json.dump(out, open(f"{BASE}/result.json", "w"), indent=1)
print("DONE")
