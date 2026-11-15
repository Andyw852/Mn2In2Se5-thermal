#!/usr/bin/env python3
"""
Standalone phono3py + symfc fitting script (fc2 + fc3)
Fitting only, submits nothing; writes fc2.hdf5, fc3.hdf5, phono3py_params.yaml
"""

import os
import sys
import numpy as np
from phonopy.interface.vasp import read_vasp
from phonopy.structure.cells import get_supercell
from phono3py import Phono3py
from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5

# ---------- configuration ----------
PRIM_FILE = "POSCAR"            # primitive cell
SPOSCAR_FILE = "SPOSCAR"        # supercell (the cell used by VASP)
DISPS_FILE = "dataset_disps.npy"
FORCES_FILE = "dataset_forces.npy"
FC3_CUTOFF = 5.5                # cutoff radius (Å); set to None for no cutoff
BORN_FILE = "BORN"              # optional; NAC is read in if the file exists

# ----- frame subsampling (can reduce the number of structures used) -----
NDATA = "all"                   # may be an integer, e.g. 20: first N frames or evenly sampled
MODE = "even"                   # "even" or "head"

# ---------------------------------------------

def main():
    # 1. load
    uc = read_vasp(PRIM_FILE)
    sc_vasp = read_vasp(SPOSCAR_FILE)
    disps_in = np.load(DISPS_FILE)      # (nframe+1, natom, 3) fractional coordinates
    forces_in = np.load(FORCES_FILE)    # (nframe+1, natom, 3) eV/Å
    natom = len(sc_vasp.numbers)
    assert disps_in.shape[1] == natom == forces_in.shape[1], \
        "atom count mismatch"

    # 2. infer the supercell matrix (pymatgen → phonopy convention)
    R = sc_vasp.cell @ np.linalg.inv(uc.cell)
    M_pmg = np.round(R).astype(int)
    if not np.allclose(R, M_pmg, atol=1e-4):
        sys.exit(f"the SPOSCAR lattice is not an integer combination of POSCAR\nR=\n{R}")
    M = M_pmg.T   # transpose to the phonopy column convention
    if np.linalg.det(M) < 0:
        M = M[:, [1, 0, 2]]   # fix handedness
    print(f"supercell matrix M (phonopy convention):\n{M}\ndet = {int(round(np.linalg.det(M)))}")

    # 3. displacements and forces: subtract the equilibrium frame
    eq_frac = disps_in[-1]
    d_frac = disps_in[:-1] - eq_frac[None]
    d_frac -= np.round(d_frac)          # minimum image
    u_vasp = d_frac @ sc_vasp.cell      # Cartesian displacements
    f_vasp = forces_in[:-1] - forces_in[-1][None]
    ntot = len(u_vasp)
    print(f"total training frames: {ntot}, max equilibrium residual force = {np.linalg.norm(forces_in[-1], axis=1).max():.5f} eV/Å")

    # 4. frame subsampling (optional)
    if NDATA == "all":
        idx = np.arange(ntot)
    else:
        n = int(NDATA)
        if MODE == "head":
            idx = np.arange(n)
        else:   # even
            idx = np.round(np.linspace(0, ntot-1, n)).astype(int)
            idx = np.unique(idx)
            while len(idx) < n:
                rest = np.setdiff1d(np.arange(ntot), idx)
                idx = np.sort(np.concatenate([idx, rest[:n-len(idx)]]))
    u, f = u_vasp[idx], f_vasp[idx]
    nuse = len(u)
    print(f"used in the fit: {nuse} frames → {nuse*natom*3} force equations")
    with open("ndata_used.txt", "w") as fo:
        fo.write(str(nuse))

    # 5. atom mapping: make the VASP and phonopy supercells use the same atom order
    sc_ph = get_supercell(uc, M)
    assert len(sc_ph.numbers) == natom
    # check that both cells describe the same superlattice
    T = sc_vasp.cell @ np.linalg.inv(sc_ph.cell)
    if not (np.allclose(T, np.round(T), atol=1e-5) and abs(abs(np.linalg.det(T))-1.0)<1e-6):
        sys.exit("the two supercells do not describe the same superlattice, check the matrix")
    # find the permutation from the coordinates
    cart_va = sc_vasp.scaled_positions @ sc_vasp.cell
    fp_va = cart_va @ np.linalg.inv(sc_ph.cell)
    fp_ph = sc_ph.scaled_positions
    perm = np.full(natom, -1, dtype=int)
    used = np.zeros(natom, dtype=bool)
    for i in range(natom):
        d = fp_va - fp_ph[i]
        d -= np.round(d)
        dist = np.linalg.norm(d @ sc_ph.cell, axis=1)
        dist[used] = np.inf
        j = int(dist.argmin())
        if dist[j] > 1e-3 or sc_vasp.numbers[j] != sc_ph.numbers[i]:
            sys.exit(f"atom mapping failed i={i}, min dist={dist[j]:.4f}")
        perm[i] = j; used[j] = True
    u = u[:, perm]
    f = f[:, perm]
    print("atom mapping done")

    # 6. fit
    ph3 = Phono3py(uc, supercell_matrix=M, primitive_matrix=np.eye(3))
    ph3.dataset = {"displacements": u, "forces": f}
    # fc2
    print("\nfitting fc2 ...")
    ph3.produce_fc2(fc_calculator="symfc")
    # fc3
    opts = None if FC3_CUTOFF is None else f"cutoff = {FC3_CUTOFF}"
    print(f"fitting fc3 (options: {opts}) ...")
    ph3.produce_fc3(fc_calculator="symfc", fc_calculator_options=opts)

    # 7. write out
    p2s = ph3.primitive.p2s_map
    write_fc2_to_hdf5(ph3.fc2, filename="fc2.hdf5", p2s_map=p2s)
    write_fc3_to_hdf5(ph3.fc3, filename="fc3.hdf5", p2s_map=p2s)

    # 8. NAC (if BORN is present)
    if os.path.isfile(BORN_FILE):
        from phonopy.file_IO import parse_BORN
        ph3.nac_params = parse_BORN(ph3.primitive, filename=BORN_FILE)
        print(f"read in {BORN_FILE}")
    ph3.save("phono3py_params.yaml")
    print("✅ wrote fc2.hdf5, fc3.hdf5, phono3py_params.yaml")

    # 9. check the force error
    fc2 = ph3.fc2 if ph3.fc2.shape[0] == natom else None
    fc3 = ph3.fc3 if ph3.fc3.shape[0] == natom else None
    if fc2 is not None and fc3 is not None:
        f_pred = -np.einsum("iajb,njb->nia", fc2.transpose(0,2,1,3), u)
        for n in range(len(u)):
            f_pred[n] -= 0.5 * np.einsum("ijkabc,jb,kc->ia", fc3, u[n], u[n])
        err = np.linalg.norm(f_pred - f) / np.linalg.norm(f)
        print(f"\nRelative error: {err:.6f}")
        with open("rel_err.txt", "w") as fo:
            fo.write(f"{err:.6f}")
    else:
        print("ℹ️  compact fc, skipping the force check")

if __name__ == "__main__":
    main()