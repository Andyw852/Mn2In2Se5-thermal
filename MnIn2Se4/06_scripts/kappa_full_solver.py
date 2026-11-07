#!/usr/bin/env python3
"""phono3py RTA (MnIn2Se4, primitive = 7 atoms, supercell = 3x3x3 = 189 atoms).
Writes the combined kappa-m{mesh}.hdf5 (frequency, gamma, gamma_N, gamma_U,
gamma_isotope, group_velocity, heat_capacity, mode_kappa, kappa, ...)."""
import sys, numpy as np
from phonopy.interface.vasp import read_vasp
from phono3py import Phono3py
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5
mesh = int(sys.argv[1]) if len(sys.argv) > 1 else 11
uc = read_vasp('POSCAR')
M = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
ph3 = Phono3py(uc, supercell_matrix=M, primitive_matrix='P')
ph3.fc2 = read_fc2_from_hdf5('fc2.hdf5')
ph3.fc3 = read_fc3_from_hdf5('fc3.hdf5')
ph3.mesh_numbers = [mesh, mesh, mesh]
ph3.init_phph_interaction()
ph3.run_thermal_conductivity(temperatures=[300.0], is_isotope=True, write_kappa=True)
k = np.asarray(ph3.thermal_conductivity.kappa).ravel()
print('KAPPA %12.6f %12.6f %12.6f' % (k[0], k[1], k[2]))
