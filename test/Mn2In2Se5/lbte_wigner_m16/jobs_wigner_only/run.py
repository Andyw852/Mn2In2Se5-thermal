#!/usr/bin/env python3
"""Mn2In2Se5: Wigner-LBTE only (16^3, FC3 6.51 A, isotopes off, 300 K).

Same fixed indexing as ../mn441_wigner_only/run.py. The LBTE part of job 412628
runs first and writes kappa-m161616.hdf5; this script adds the Wigner-LBTE solve
so the coherence term is on the same footing as the committed RTA and LBTE runs.
"""
import json
import os
import time

import numpy as np
import phono3py
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

BASE = os.path.dirname(os.path.abspath(__file__))
MESH = [16, 16, 16]

def snap(k):
    a = np.asarray(k).reshape(-1, 6)[0]
    return {"inplane": float(0.5 * (a[0] + a[1])), "kappa_zz": float(a[2])}

p = phono3py.load(f"{BASE}/in/phono3py_disp.yaml", produce_fc=False, log_level=0)
ph3 = phono3py.Phono3py(p.unitcell, supercell_matrix=np.array(p.supercell_matrix, int),
                        primitive_matrix="P", log_level=0)
ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/in/fc3.hdf5")
ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2.hdf5")
ph3.mesh_numbers = MESH
ph3.init_phph_interaction()

print("start Wigner-LBTE", flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=True, is_isotope=False,
                             conductivity_type="wigner", write_kappa=False, log_level=0)
tc = ph3.thermal_conductivity
tot, pop = snap(tc.kappa_TOT_exact), snap(tc.kappa_P_exact)
coh = {"inplane": tot["inplane"] - pop["inplane"],
       "kappa_zz": tot["kappa_zz"] - pop["kappa_zz"]}
out = {"material": "Mn2In2Se5", "mesh": " ".join(map(str, MESH)), "fc2": "4x4x1",
       "fc3_cutoff_A": 6.51, "T_K": 300.0, "isotope": False, "primitive_matrix": "P",
       "Wigner_LBTE": {"TOT": tot, "population": pop, "coherence": coh,
                       "seconds": round(time.time() - t0, 1)}}
print("Wigner-LBTE", json.dumps(out["Wigner_LBTE"]), flush=True)
json.dump(out, open(f"{BASE}/result.json", "w"), indent=1)
print("DONE")
