#!/usr/bin/env python3
"""MnIn2Se4: Wigner-LBTE thermal conductivity on a 14x14x14 q-mesh.

Force constants
  FC2 : fc2_442_phonon/fc2.hdf5            (4x4x2 phonon supercell)
  FC3 : fit_cut3_6p45/fc3.hdf5             (4x4x1, 6.45 A cutoff)
Method : phono3py conductivity_type="wigner" with is_LBTE=True, i.e. the
         full collision matrix is inverted and the coherence (off-diagonal)
         contribution is kept. Isotope scattering on, T = 300 K.
Output : k14_wigner_lbte.json  (written next to this script)
         TOT = total, P = population part, C = coherence part.
"""
import json
import time

import numpy as np
import phono3py
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5
from phonopy.interface.vasp import read_vasp

D = "/home/wangchao/repo_thermal/test/MnIn2Se4"
OUT = D
uc = read_vasp(f"{D}/fit/POSCAR")
ph3 = phono3py.Phono3py(uc, supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 1]],
                        phonon_supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 2]],
                        primitive_matrix=[[1, 1, 2], [-1, -2, -1], [1, 1, 1]],
                        log_level=0)
ph3.fc3 = read_fc3_from_hdf5(f"{D}/fit_cut3_6p45/fc3.hdf5")
ph3.fc2 = read_fc2_from_hdf5(f"{D}/fc2_442_phonon/fc2.hdf5")
ph3.mesh_numbers = [14, 14, 14]
ph3.init_phph_interaction()
print("Wigner-LBTE mesh 14 started", flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=True,
                             is_isotope=True, conductivity_type="wigner",
                             write_kappa=False, log_level=0)
dt = time.time() - t0
tc = ph3.thermal_conductivity


def ip(k):
    a = np.asarray(k)[0, 0]
    return float(0.5 * (a[0] + a[1])), float(a[2])


tot = ip(tc.kappa_TOT_exact)
pop = ip(tc.kappa_P_exact)
coh = (tot[0] - pop[0], tot[1] - pop[1])
print(f"Wigner-LBTE {dt:.0f}s", flush=True)
print(f"  TOT in-plane={tot[0]:.4f}  kappa_zz={tot[1]:.4f}", flush=True)
print(f"  P   in-plane={pop[0]:.4f}  kappa_zz={pop[1]:.4f}", flush=True)
print(f"  C   in-plane={coh[0]:.4f}  kappa_zz={coh[1]:.4f}", flush=True)
json.dump({"method": "Wigner-LBTE", "mesh": "14 14 14", "cutoff_A": 6.45, "T_K": [300.0],
           "TOT": {"inplane": tot[0], "kappa_zz": tot[1]},
           "population": {"inplane": pop[0], "kappa_zz": pop[1]},
           "coherence": {"inplane": coh[0], "kappa_zz": coh[1]},
           "seconds": round(dt, 1)},
          open(f"{OUT}/k14_wigner_lbte.json", "w"), indent=2)
print("DONE")
