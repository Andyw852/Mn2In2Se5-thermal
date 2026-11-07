#!/usr/bin/env python3
"""Check that phonopy's 3x3x3 supercell atom order/direct coords == SPOSCAR."""
import numpy as np
from phonopy.interface.vasp import read_vasp
from phono3py import Phono3py
uc = read_vasp('POSCAR'); sc = read_vasp('SPOSCAR')
ph3 = Phono3py(uc, supercell_matrix=[[3,0,0],[0,3,0],[0,0,3]], primitive_matrix='P')
s = ph3.supercell
print('n atoms phonopy/sc', len(s), len(sc))
print('symbols equal:', list(s.symbols)==list(sc.symbols))
d = np.asarray(s.scaled_positions) - np.asarray(sc.scaled_positions)
d -= np.round(d)
print('max scaled-coord diff: %.3e' % np.abs(d).max())
print('species head phonopy:', list(s.symbols[:6]), 'sc:', list(sc.symbols[:6]))
