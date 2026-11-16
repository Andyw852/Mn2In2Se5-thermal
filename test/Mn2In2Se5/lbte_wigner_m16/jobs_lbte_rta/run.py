#!/usr/bin/env python3
"""Mn2In2Se5: coherence term and iterative solution at 16^3, FC3 6.51 A.

Two failures made this a separate, incremental script:
  - the combined job (412628) was SIGKILLed by the kernel at 11 GB RSS because
    the node was out of memory (node-level pressure, not a cgroup limit), and
    no partial output survived;
  - phono3py 4.4.0 has no conductivity_type keyword -- the coherence solver is
    transport_type="SMM19".

So: run the cheap Wigner-RTA solve FIRST (that is the coherence term the reports
need), write its JSON immediately, then the plain LBTE solve (the iterative
correction), writing its JSON immediately too. A kill now loses only the step in
progress.

Order matters: result_wigner_rta.json is written before the expensive step.
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

def dump(name, obj):
    json.dump(obj, open(f"{BASE}/{name}", "w"), indent=1)
    print("wrote", name, flush=True)

p = phono3py.load(f"{BASE}/in/phono3py_disp.yaml", produce_fc=False, log_level=0)
ph3 = phono3py.Phono3py(p.unitcell, supercell_matrix=np.array(p.supercell_matrix, int),
                        primitive_matrix="P", log_level=0)
ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/in/fc3.hdf5")
ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2.hdf5")
ph3.mesh_numbers = MESH
ph3.init_phph_interaction()

meta = {"material": "Mn2In2Se5", "mesh": " ".join(map(str, MESH)), "fc2": "4x4x1",
        "fc3_cutoff_A": 6.51, "T_K": 300.0, "isotope": False, "primitive_matrix": "P"}

# ---- 1) Wigner-RTA: the coherence term --------------------------------------
print("start Wigner-RTA (SMM19)", flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=False, is_isotope=False,
                             transport_type="SMM19", write_kappa=True, log_level=0)
tc = ph3.thermal_conductivity
# phono3py 4.4.0 interband (SMM19/NJC23/IBDB19) solvers expose
#   kappa        = TOT (intra + inter)
#   kappa_intra  = population part
#   kappa_inter  = coherence part
# There is no kappa_TOT_RTA / kappa_C on this object in 4.4.0.
tot = snap(tc.kappa)
pop = snap(tc.kappa_intra)
coh = snap(tc.kappa_inter)
out = dict(meta)
out["Wigner_RTA"] = {"TOT": tot, "population": pop, "coherence": coh,
                     "seconds": round(time.time() - t0, 1)}
dump("result_wigner_rta.json", out)

# ---- 2) plain LBTE: the iterative correction --------------------------------
print("start LBTE", flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=True, is_isotope=False,
                             write_kappa=True, log_level=0)
tc = ph3.thermal_conductivity
rta = snap(getattr(tc, "kappa_RTA")) if hasattr(tc, "kappa_RTA") else None
out2 = dict(meta)
out2["LBTE"] = snap(tc.kappa)
out2["RTA_same_run"] = rta
out2["seconds"] = round(time.time() - t0, 1)
dump("result_lbte.json", out2)
print("DONE")
