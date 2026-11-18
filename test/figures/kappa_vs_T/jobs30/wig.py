#!/usr/bin/env python3
"""Wigner coherence term of MnIn2Se4, for comparing phono3py versions (report §7).

  python wig.py <mesh> [SMM19|NJC23|IBDB19]

phono3py 3.x : conductivity_type="wigner" (SMM19); the variant argument is ignored.
phono3py 4.x : transport_type=<variant>; prints kappa_intra (population) and
               kappa_inter (coherence).
Prints in-plane (xx) and zz at 300 K, isotopes on, 4x4x2 FC2 + 6.45 A FC3.
Reference result (11^3): 3.24.0 P 3.554/1.250, C 0.106/0.067;
4.4.0 SMM19 C 1.242/1.165, NJC23 C 1.242/1.165, IBDB19 C 0.114/0.067.
"""
import os
import sys

import numpy as np
import phono3py
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5
from phonopy.interface.vasp import read_vasp

R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "MnIn2Se4"))
N = int(sys.argv[1])
VARIANT = sys.argv[2] if len(sys.argv) > 2 else "SMM19"

ph3 = phono3py.Phono3py(read_vasp(f"{R}/fit/POSCAR"),
                        supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 1]],
                        phonon_supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 2]],
                        primitive_matrix=[[1, 1, 2], [-1, -2, -1], [1, 1, 1]], log_level=0)
ph3.fc3 = read_fc3_from_hdf5(f"{R}/fit_cut3_6p45/fc3.hdf5")
ph3.fc2 = read_fc2_from_hdf5(f"{R}/fc2_442_phonon/fc2.hdf5")
ph3.mesh_numbers = [N] * 3
ph3.init_phph_interaction()
kw = dict(temperatures=[300.0], is_isotope=True, log_level=0)
major = int(phono3py.__version__.split(".")[0])
if major < 4:
    ph3.run_thermal_conductivity(conductivity_type="wigner", **kw)
else:
    ph3.run_thermal_conductivity(transport_type=VARIANT, **kw)
tc = ph3.thermal_conductivity


def xz(v):
    a = np.asarray(v, dtype=float).reshape(-1, 6)[0]
    return f"{a[0]:.4f} / {a[2]:.4f}"


if major < 4:
    print(f"{phono3py.__version__} SMM19 N={N}: P {xz(tc.kappa_P_RTA)}  C {xz(tc.kappa_C)}  "
          f"TOT {xz(tc.kappa_TOT_RTA)}   (xx / zz)")
else:
    ex = tc.get_extra_kappa_output()
    print(f"{phono3py.__version__} {VARIANT} N={N}: P {xz(ex['kappa_intra'])}  "
          f"C {xz(ex['kappa_inter'])}  TOT {xz(tc.kappa)}   (xx / zz)")
