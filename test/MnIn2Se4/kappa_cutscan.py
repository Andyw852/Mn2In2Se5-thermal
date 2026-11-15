import numpy as np, phono3py, time, json
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5
from phonopy.interface.vasp import read_vasp

d = "/home/wangchao/repo_thermal/test/MnIn2Se4"
uc = read_vasp(f"{d}/fit/POSCAR")
FC2_442 = f"{d}/fc2_442_phonon/fc2.hdf5"
PROD = f"{d}/fit_cut3_6p45/fc3.hdf5"

CASES = [("prod_6p45", PROD)] + [(f"scan_{t}", f"/tmp/fc3scan/fc3_{t}.hdf5")
                                 for t in ("5p60", "6p23", "6p45")]

for tag, f3 in CASES:
    fc3 = read_fc3_from_hdf5(f3)
    ph3 = phono3py.Phono3py(uc, supercell_matrix=[[4,0,0],[0,4,0],[0,0,1]],
                            phonon_supercell_matrix=[[4,0,0],[0,4,0],[0,0,2]],
                            primitive_matrix=[[1,1,2],[-1,-2,-1],[1,1,1]], log_level=0)
    ph3.fc3 = fc3
    ph3.fc2 = read_fc2_from_hdf5(FC2_442)
    ph3.mesh_numbers = [11,11,11]
    ph3.init_phph_interaction()
    t0=time.time()
    ph3.run_thermal_conductivity(temperatures=[300.0], is_LBTE=False,
                                 is_isotope=True, write_kappa=False, log_level=0)
    dt=time.time()-t0
    k=np.asarray(ph3.thermal_conductivity.kappa)[0,0]
    inpl=0.5*(k[0]+k[1])
    print(f"{tag:12s} fc3shape={fc3.shape} RTA {dt:5.0f}s  in-plane={inpl:8.4f} kzz={k[2]:8.4f}", flush=True)
    json.dump({"case":tag,"cutoff_A":float(tag.split('_')[-1].replace('p','.')) if tag!="prod_6p45" else 6.45,
               "method":"RTA","inplane":float(inpl),"kappa_zz":float(k[2]),"seconds":round(dt,1)},
              open(f"/tmp/kx_{tag}.json","w"), indent=2)
print("DONE")
