#!/usr/bin/env python3
"""Mn2In2Se5: is the 4.4.0 SMM19 coherence term mesh-stable?

The report's +12.4 % coherence ratio for Mn2In2Se5 came from HanHai job 416475,
which used transport_type="SMM19" on phono3py 4.4.0 -- i.e. the same
implementation that is 1.8x mesh-unstable and ~17x too large for MnIn2Se4.
This script reruns the 6.51 A case at several meshes so the stability can be
judged directly, and also with IBDB19 for reference.

  python run.py 11 12 14      # meshes to run, all with SMM19 + IBDB19
"""
import json
import os
import sys
import time

import numpy as np
import phono3py
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

BASE = os.path.dirname(os.path.abspath(__file__))
MESHES = [int(x) for x in sys.argv[1:]] or [11, 12, 14]
print("phono3py", phono3py.__version__, flush=True)


def snap(k):
    a = np.asarray(k).reshape(-1, 6)[0]
    return {"inplane": float(0.5 * (a[0] + a[1])), "kappa_zz": float(a[2])}


out = {"material": "Mn2In2Se5", "fc2": "4x4x1", "fc3": "6.51 A", "T_K": 300.0,
       "isotope": False, "transport_types": ["SMM19", "IBDB19"], "meshes": {}}
for N in MESHES:
    out["meshes"][str(N)] = {}
    for variant in ("SMM19", "IBDB19"):
        p = phono3py.load(f"{BASE}/in/phono3py_disp.yaml", produce_fc=False, log_level=0)
        ph3 = phono3py.Phono3py(p.unitcell, supercell_matrix=np.array(p.supercell_matrix, int),
                                primitive_matrix="P", log_level=0)
        ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/in/fc3.hdf5")
        ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2.hdf5")
        ph3.mesh_numbers = [N, N, N]
        ph3.init_phph_interaction()
        t0 = time.time()
        ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=False, is_isotope=False,
                                     transport_type=variant, write_kappa=False, log_level=0)
        tc = ph3.thermal_conductivity
        out["meshes"][str(N)][variant] = {
            "TOT": snap(tc.kappa), "population": snap(tc.kappa_intra),
            "coherence": snap(tc.kappa_inter), "seconds": round(time.time() - t0, 1)}
        c = out["meshes"][str(N)][variant]
        print("%2d^3 %-7s P %.4f/%.4f  C %.4f/%.4f  %.0f s" % (
            N, variant, c["population"]["inplane"], c["population"]["kappa_zz"],
            c["coherence"]["inplane"], c["coherence"]["kappa_zz"], c["seconds"]), flush=True)
        json.dump(out, open(f"{BASE}/result_mesh.json", "w"), indent=1)
print("DONE", flush=True)
