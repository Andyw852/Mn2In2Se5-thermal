#!/usr/bin/env python3
"""MnIn2Se4 kappa_zz on the q-mesh -- HanHai job entry point.

FC2 : in/fc2_442.hdf5   (4x4x2 phonon supercell; this is the corrected FC2)
FC3 : in/<--fc3>        (4x4x1, 5.60 / 6.23 / 6.45 A variants)
Method: phono3py RTA, isotope scattering on, T = 300 K, write_kappa so the
        per-mode output is kept for kappa_zz_decompose.py.
The driver is deliberately identical in physics to
test/MnIn2Se4/kappa_driver_442.py in the repository; only the paths differ.
"""
import argparse
import json
import os
import time

import numpy as np
import phono3py
from phonopy.interface.vasp import read_vasp
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

BASE = os.path.dirname(os.path.abspath(__file__))

ap = argparse.ArgumentParser()
ap.add_argument("--mesh", nargs=3, type=int, required=True)
ap.add_argument("--fc3", required=True, help="file name inside in/")
ap.add_argument("--out", required=True)
ap.add_argument("--boundary-mfp", type=float, default=None)
a = ap.parse_args()

ph3 = phono3py.Phono3py(
    read_vasp(f"{BASE}/in/POSCAR"),
    supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 1]],
    phonon_supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 2]],
    primitive_matrix=[[1, 1, 2], [-1, -2, -1], [1, 1, 1]],
    log_level=0,
)
ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/in/{a.fc3}")
ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2_442.hdf5")
ph3.mesh_numbers = a.mesh
ph3.init_phph_interaction()
print(f"start mesh={a.mesh} fc3={a.fc3}", flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(
    temperatures=[300.0], is_LBTE=False, is_isotope=True,
    boundary_mfp=a.boundary_mfp, write_kappa=True, log_level=0,
)
dt = time.time() - t0
k = np.asarray(ph3.thermal_conductivity.kappa)[0, 0]
ip = float(0.5 * (k[0] + k[1]))
out = {"method": "RTA", "mesh": " ".join(map(str, a.mesh)),
       "fc3": a.fc3, "boundary_mfp_um": a.boundary_mfp, "T_K": [300.0],
       "inplane": ip, "kappa_xx": float(k[0]), "kappa_zz": float(k[2]),
       "seconds": round(dt, 1),
       "phono3py": phono3py.__version__}
print(f"mesh={a.mesh} fc3={a.fc3} in-plane={ip:.4f} kzz={k[2]:.4f} ({dt:.0f}s)", flush=True)
json.dump(out, open(a.out, "w"), indent=1)
print("DONE")
