import numpy as np, phono3py, time, json
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5
from phonopy.interface.vasp import read_vasp
d="/home/wangchao/repo_thermal/test/MnIn2Se4"
uc=read_vasp(f"{d}/fit/POSCAR")
FC2=f"{d}/fc2_442_phonon/fc2.hdf5"
for tag in ("5p60","6p23","6p45"):
    fc3=read_fc3_from_hdf5(f"/tmp/fc3scan/fc3_{tag}.hdf5")
    ph3=phono3py.Phono3py(uc,supercell_matrix=[[4,0,0],[0,4,0],[0,0,1]],
                          phonon_supercell_matrix=[[4,0,0],[0,4,0],[0,0,2]],
                          primitive_matrix=[[1,1,2],[-1,-2,-1],[1,1,1]],log_level=0)
    ph3.fc3=fc3; ph3.fc2=read_fc2_from_hdf5(FC2)
    ph3.mesh_numbers=[11,11,11]; ph3.init_phph_interaction()
    t0=time.time()
    ph3.run_thermal_conductivity(temperatures=[300.0],is_LBTE=True,
                                 is_isotope=True,write_kappa=False,log_level=0)
    dt=time.time()-t0
    k=np.asarray(ph3.thermal_conductivity.kappa)[0,0]
    inpl=0.5*(k[0]+k[1])
    print(f"LBTE {tag:6s} {dt:6.0f}s  in-plane={inpl:8.4f} kzz={k[2]:8.4f}",flush=True)
    json.dump({"cutoff_A":float(tag.replace('p','.')),"method":"LBTE","mesh":"11 11 11",
               "T_K":[300.0],"inplane":float(inpl),"kappa_zz":float(k[2]),
               "seconds":round(dt,1)},open(f"/tmp/klbte_{tag}.json","w"),indent=2)
print("DONE")
