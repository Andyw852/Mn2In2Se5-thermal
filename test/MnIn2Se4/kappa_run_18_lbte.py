#!/usr/bin/env python3
"""MnIn2Se4: 18x18x18 LBTE (phono3py, fc2=4x4x2, fc3=6.45A)."""
import numpy as np, phono3py, time, json, os
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5
from phonopy.interface.vasp import read_vasp

D = "/home/wangchao/repo_thermal/test/MnIn2Se4"
OUT = "/home/wangchao/mn441_kappa"
uc = read_vasp(f"{D}/fit/POSCAR")
ph3 = phono3py.Phono3py(uc, supercell_matrix=[[4,0,0],[0,4,0],[0,0,1]],
                        phonon_supercell_matrix=[[4,0,0],[0,4,0],[0,0,2]],
                        primitive_matrix=[[1,1,2],[-1,-2,-1],[1,1,1]], log_level=0)
ph3.fc3 = read_fc3_from_hdf5(f"{D}/fit_cut3_6p45/fc3.hdf5")
ph3.fc2 = read_fc2_from_hdf5(f"{D}/fc2_442_phonon/fc2.hdf5")
ph3.mesh_numbers = [18,18,18]
ph3.init_phph_interaction()
print("18x18x18 LBTE 开始", flush=True)
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=True,
                             is_isotope=True, write_kappa=False, log_level=0)
dt = time.time()-t0
k = np.asarray(ph3.thermal_conductivity.kappa)[0,0]
ip = float(0.5*(k[0]+k[1])); kz = float(k[2])
print(f"18x18x18 LBTE {dt:.0f}s 面内={ip:.4f} kzz={kz:.4f}", flush=True)
json.dump({"method":"LBTE","mesh":"18 18 18","cutoff_A":6.45,"T_K":[300.0],
           "inplane":ip,"kappa_zz":kz,"seconds":round(dt,1)},
          open(f"{OUT}/k18_lbte.json","w"), indent=2)
print("DONE")
