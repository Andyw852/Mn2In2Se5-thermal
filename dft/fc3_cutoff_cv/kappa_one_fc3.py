#!/usr/bin/env python3
"""kappa for one fitted FC3 of the CV study (protocol step C).

16^3, 300 K, isotopes on, primitive_matrix=P, production FC2 fixed and only FC3
varying -- exactly the settings of test/Mn2In2Se5/cutoff_ext. Writes the per-mode
hdf5 so kappa_zz_decompose.py can split the rest term.
"""
import argparse
import json
import os

import numpy as np
import phono3py
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

ap = argparse.ArgumentParser()
ap.add_argument('--fc3', required=True)
ap.add_argument('--prod-fc2', required=True)
ap.add_argument('--yaml', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--mesh', type=int, default=16)
a = ap.parse_args()

base = phono3py.load(a.yaml, produce_fc=False, log_level=0)
ph3 = phono3py.Phono3py(base.unitcell,
                        supercell_matrix=np.array(base.supercell_matrix, int),
                        primitive_matrix='P', log_level=0)
ph3.fc3 = read_fc3_from_hdf5(a.fc3)
ph3.fc2 = read_fc2_from_hdf5(a.prod_fc2)
ph3.mesh_numbers = [a.mesh] * 3
ph3.init_phph_interaction()
ph3.run_thermal_conductivity(temperatures=[300.0], is_isotope=True,
                             write_kappa=True, log_level=0)
k = np.asarray(ph3.thermal_conductivity.kappa).reshape(-1, 6)[0]
res = {'inplane': float(0.5 * (k[0] + k[1])), 'kappa_zz': float(k[2]),
       'mesh': a.mesh, 'fc3': a.fc3}
src = 'kappa-m' + str(a.mesh).zfill(2) * 3 + '.hdf5'
if os.path.exists(src):
    os.replace(src, a.out + '_kappa.hdf5')
with open(a.out + '.json', 'w') as fh:
    json.dump(res, fh, indent=1)
print(json.dumps(res))
