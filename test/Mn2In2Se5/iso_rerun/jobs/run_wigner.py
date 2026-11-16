#!/usr/bin/env python3
"""Mn2In2Se5 Wigner-RTA (coherence) at each FC3 cutoff, ISOTOPES ON.

Same as test/Mn2In2Se5/wigner_cutoff/ but with is_isotope=True, to match the
project configuration (ISOTOPE = true) and the MnIn2Se4 runs.
"""
import json
import os
import time

import numpy as np
import phono3py
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

BASE = os.path.dirname(os.path.abspath(__file__))
CASES = [("5p68", "fc3_cut5p68.hdf5"), ("6p12", "fc3_cut6p12.hdf5"),
         ("6p51", "fc3_cut6p51_scan.hdf5")]

def snap(k):
    a = np.asarray(k).reshape(-1, 6)[0]
    return {"inplane": float(0.5 * (a[0] + a[1])), "kappa_zz": float(a[2])}

out = {"material": "Mn2In2Se5", "mesh": "16 16 16", "isotope": True,
       "primitive_matrix": "P", "T_K": 300.0, "cases": {}}
for tag, f in CASES:
    p = phono3py.load(f"{BASE}/in/phono3py_disp.yaml", produce_fc=False, log_level=0)
    ph3 = phono3py.Phono3py(p.unitcell, supercell_matrix=np.array(p.supercell_matrix, int),
                            primitive_matrix="P", log_level=0)
    ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/in/{f}")
    ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2.hdf5")
    ph3.mesh_numbers = [16, 16, 16]
    ph3.init_phph_interaction()
    print("Wigner-RTA", tag, flush=True)
    t0 = time.time()
    ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=False, is_isotope=True,
                                 transport_type="SMM19", write_kappa=False, log_level=0)
    tc = ph3.thermal_conductivity
    out["cases"][tag] = {"TOT": snap(tc.kappa), "population": snap(tc.kappa_intra),
                         "coherence": snap(tc.kappa_inter),
                         "seconds": round(time.time() - t0, 1)}
    print("  ", out["cases"][tag]["TOT"], flush=True)
    json.dump(out, open(f"{BASE}/result_wigner_cutoff_iso.json", "w"), indent=1)
print("DONE", flush=True)
