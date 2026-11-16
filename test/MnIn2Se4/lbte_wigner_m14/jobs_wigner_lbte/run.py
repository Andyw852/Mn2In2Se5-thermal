#!/usr/bin/env python3
"""MnIn2Se4 Wigner-LBTE at 14^3 / 6.45 A -- guaranteed TOT capture.

Job 414857 ran the solve successfully in 1 h 06 min but then died on
tc.kappa_intra, which does not exist on LBTECalculator. The 4.4.0 names for the
SMM19 LBTE solver are:

    kappa             = TOT (intra + inter)
    kappa_intra_exact = LBTE population part
    kappa_intra_RTA   = RTA population part
    kappa_inter       = coherence part

This script writes the TOT to JSON *immediately* after the solve (before touching
any of the part attributes), then adds the parts inside try/except. So a wrong
attribute name can no longer throw away an hour of compute.
"""
import json
import os
import time
import traceback

import numpy as np
import phono3py
from phonopy.interface.vasp import read_vasp
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

BASE = os.path.dirname(os.path.abspath(__file__))
MESH = [14, 14, 14]

def snap(k):
    a = np.asarray(k).reshape(-1, 6)[0]
    return {"inplane": float(0.5 * (a[0] + a[1])), "kappa_zz": float(a[2])}

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

print("start Wigner-LBTE (SMM19)", flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=True, is_isotope=True,
                             transport_type="SMM19", write_kappa=True, log_level=0)
dt = round(time.time() - t0, 1)
tc = ph3.thermal_conductivity

out = {"material": "MnIn2Se4", "mesh": " ".join(map(str, MESH)), "fc2": "4x4x2",
       "fc3_cutoff_A": 6.45, "T_K": 300.0, "isotope": True,
       "solver": "phono3py 4.4.0, transport_type=SMM19, is_LBTE=True",
       "seconds": dt, "TOT": snap(tc.kappa)}
json.dump(out, open(f"{BASE}/result_tot.json", "w"), indent=1)
print("saved TOT", out["TOT"], flush=True)

try:
    out["population_LBTE"] = snap(tc.kappa_intra_exact)
    out["population_RTA"] = snap(tc.kappa_intra_RTA)
    out["coherence"] = snap(tc.kappa_inter)
except BaseException:
    with open(f"{BASE}/parts_error.txt", "w") as f:
        f.write(traceback.format_exc())
    print("part attributes failed; TOT is saved", flush=True)

json.dump(out, open(f"{BASE}/result.json", "w"), indent=1)
print("DONE", flush=True)
