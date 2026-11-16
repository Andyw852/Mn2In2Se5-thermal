#!/usr/bin/env python3
"""MnIn2Se4: LBTE and Wigner-LBTE at the largest FC3 cutoff (6.45 A), 14^3.

Runs two phono3py solutions in one job:
  1) is_LBTE=True, conductivity_type=None  -> the iterative (LBTE) kappa
  2) is_LBTE=True, conductivity_type="wigner" -> TOT / population / coherence
Both with isotope scattering ON and T = 300 K, matching the committed RTA runs.
FC2 = 4x4x2 phonon supercell, FC3 = 4x4x1 at 6.45 A (the cell's ceiling).
"""
import json
import os
import time

import numpy as np
import phono3py
from phonopy.interface.vasp import read_vasp
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

BASE = os.path.dirname(os.path.abspath(__file__))
MESH = [14, 14, 14]

ph3 = phono3py.Phono3py(
    read_vasp(f"{BASE}/in/POSCAR"),
    supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 1]],
    phonon_supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 2]],
    primitive_matrix=[[1, 1, 2], [-1, -2, -1], [1, 1, 1]],
    log_level=0,
)
ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/in/fc3_645.hdf5")
ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2_442.hdf5")
ph3.mesh_numbers = MESH
ph3.init_phph_interaction()

def snap(k):
    # phono3py returns these kappa arrays as (1, 6) or (1, 1, 6) depending on the
    # solver; flatten everything except the last (Voigt) axis. The first version
    # of this script indexed only one axis and crashed AFTER a successful LBTE
    # solve, before writing the JSON.
    a = np.asarray(k).reshape(-1, 6)[0]
    return {"inplane": float(0.5 * (a[0] + a[1])), "kappa_zz": float(a[2])}

out = {"material": "MnIn2Se4", "mesh": " ".join(map(str, MESH)),
       "fc2": "4x4x2", "fc3_cutoff_A": 6.45, "T_K": 300.0, "isotope": True}

print("start LBTE", flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=True, is_isotope=True,
                             write_kappa=True, log_level=0)
out["LBTE"] = snap(ph3.thermal_conductivity.kappa)
out["LBTE"]["seconds"] = round(time.time() - t0, 1)
print("LBTE", out["LBTE"], flush=True)

print("start Wigner-LBTE", flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=True, is_isotope=True,
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
