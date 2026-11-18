#!/usr/bin/env python3
"""Mn2In2Se5: does the 4.4.0 SMM19 coherence term depend on soft modes?

Same test as for MnIn2Se4 in report section 7. cutoff_frequency = 0.3 THz removes
the soft acoustic modes near Gamma, where the SMM19 C/omega term diverges, and
the coherence term is compared with and without the cutoff, for SMM19 and IBDB19.

Diagnostic only: the cutoff also removes those modes from the population term, so
the P printed here is not a result.

  python run_cutoff.py [mesh] [cutoff_THz]
"""
import json
import os
import sys
import time

import numpy as np
import phono3py
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

BASE = os.path.dirname(os.path.abspath(__file__))
MESH = int(sys.argv[1]) if len(sys.argv) > 1 else 16
CUT = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
print("phono3py", phono3py.__version__, "mesh", MESH, "cutoff", CUT, "THz", flush=True)


def snap(k):
    a = np.asarray(k).reshape(-1, 6)[0]
    return {"inplane": float(0.5 * (a[0] + a[1])), "kappa_zz": float(a[2])}


out = {"material": "Mn2In2Se5", "fc2": "4x4x1", "fc3": "6.51 A", "T_K": 300.0,
       "isotope": False, "mesh": MESH, "cutoff_THz": CUT, "cases": {}}
for variant in ("SMM19", "IBDB19", "SMM19_default"):
    cut = 0.0 if variant.endswith("_default") else CUT
    ttype = variant.replace("_default", "")
    p = phono3py.load(f"{BASE}/in/phono3py_disp.yaml", produce_fc=False, log_level=0)
    kwargs = dict(supercell_matrix=np.array(p.supercell_matrix, int),
                  primitive_matrix="P", log_level=0)
    if cut > 0:
        kwargs["cutoff_frequency"] = cut
    ph3 = phono3py.Phono3py(p.unitcell, **kwargs)
    ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/in/fc3.hdf5")
    ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2.hdf5")
    ph3.mesh_numbers = [MESH, MESH, MESH]
    ph3.init_phph_interaction()
    t0 = time.time()
    ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=False, is_isotope=False,
                                 transport_type=ttype, write_kappa=False, log_level=0)
    tc = ph3.thermal_conductivity
    out["cases"][variant] = {
        "cutoff_THz": cut, "TOT": snap(tc.kappa), "population": snap(tc.kappa_intra),
        "coherence": snap(tc.kappa_inter), "seconds": round(time.time() - t0, 1)}
    c = out["cases"][variant]
    print("%-14s cut=%.1f  P %.4f/%.4f  C %.4f/%.4f  TOT %.4f/%.4f  (%.0f s)" % (
        variant, cut, c["population"]["inplane"], c["population"]["kappa_zz"],
        c["coherence"]["inplane"], c["coherence"]["kappa_zz"],
        c["TOT"]["inplane"], c["TOT"]["kappa_zz"], c["seconds"]), flush=True)
    json.dump(out, open(f"{BASE}/result_cutoff.json", "w"), indent=1)
print("DONE", flush=True)
