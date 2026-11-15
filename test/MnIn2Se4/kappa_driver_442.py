#!/usr/bin/env python3
"""Compute κ with phono3py directly, bypassing AutoZT: fc2(4x4x2) + fc3(4x4x1).

Why the supercells can be mixed
-------------------------------
phono3py allows fc2 to use a larger supercell than fc3 (`phonon_supercell_matrix`).
Here fc3 is 4x4x1 (112 atoms) and fc2 is 4x4x2 (224 atoms).

⚠️ Why fc2 must be replaced
---------------------------
The old 4x4x1 fc2 gives a transverse sound velocity of only 322.9 m/s along the
layer normal (DFT elastic constants give 1338.3), so the phonon dispersion was
artificially softened → κzz was strongly underestimated. The same branch from the
4x4x2 fc2 gives 1348.3 m/s (+0.75%), i.e. a correct dispersion.

Usage
-----
    python kappa_driver_442.py            # run both RTA and LBTE
    python kappa_driver_442.py --rta      # run RTA only (about 2 min per case)
    python kappa_driver_442.py --lbte     # run LBTE only (about 15 min per case)

Output: kappa_results_442.json
"""
import argparse
import json
import time

import numpy as np
import phono3py
from phonopy.interface.vasp import read_vasp
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

# --- paths (following the repo test/MnIn2Se4 layout; pip install phonopy phono3py first) ---
ROOT = "."
POSCAR = f"{ROOT}/fit/POSCAR"
FC3 = f"{ROOT}/fit_cut3_6p45/fc3.hdf5"          # 4x4x1, cutoff 6.45 A
FC2_441 = f"{ROOT}/fit_cut3_6p45/fc2.hdf5"      # 4x4x1 (old, for comparison)
FC2_442 = f"{ROOT}/fc2_442_phonon/fc2.hdf5"     # 4x4x2 (new)

SUPERCELL_FC3 = [[4, 0, 0], [0, 4, 0], [0, 0, 1]]
PHONON_SUPERCELL = {                            # fc2 supercell
    "4x4x1": [[4, 0, 0], [0, 4, 0], [0, 0, 1]],
    "4x4x2": [[4, 0, 0], [0, 4, 0], [0, 0, 2]],
}
PRIMITIVE_MATRIX = [[1, 1, 2], [-1, -2, -1], [1, 1, 1]]
MESH = [11, 11, 11]
TEMPS = [300.0]


def run(fc2_path, phonon_supercell, lbte):
    ph3 = phono3py.Phono3py(
        read_vasp(POSCAR),
        supercell_matrix=SUPERCELL_FC3,
        phonon_supercell_matrix=phonon_supercell,
        primitive_matrix=PRIMITIVE_MATRIX,
        log_level=0,
    )
    ph3.fc3 = read_fc3_from_hdf5(FC3)
    ph3.fc2 = read_fc2_from_hdf5(fc2_path)
    ph3.mesh_numbers = MESH
    ph3.init_phph_interaction()
    t0 = time.time()
    ph3.run_thermal_conductivity(
        temperatures=TEMPS, is_LBTE=lbte, is_isotope=True,
        write_kappa=False, log_level=0,
    )
    dt = time.time() - t0
    k = np.asarray(ph3.thermal_conductivity.kappa)[0, 0]   # (6,) Voigt
    inplane = 0.5 * (k[0] + k[1])
    return {
        "method": "LBTE" if lbte else "RTA",
        "mesh": " ".join(map(str, MESH)),
        "T_K": TEMPS,
        "kappa_xx": float(k[0]), "kappa_yy": float(k[1]), "kappa_zz": float(k[2]),
        "inplane": float(inplane),
        "kzz_over_inplane": float(k[2] / inplane),
        "seconds": round(dt, 1),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rta", action="store_true")
    ap.add_argument("--lbte", action="store_true")
    a = ap.parse_args()
    do_rta = a.rta or not (a.rta or a.lbte)
    do_lbte = a.lbte or not (a.rta or a.lbte)

    cases = [("4x4x1_old", FC2_441, PHONON_SUPERCELL["4x4x1"]),
             ("4x4x2_new", FC2_442, PHONON_SUPERCELL["4x4x2"])]

    out = []
    for name, fc2, psc in cases:
        for lbte in ([False] if do_rta else []) + ([True] if do_lbte else []):
            r = run(fc2, psc, lbte)
            r["fc2_case"] = name
            r["fc2_path"] = fc2
            out.append(r)
            print(f"[{name}] {r['method']:<4} mesh={r['mesh']} 300K  "
                  f"in-plane={r['inplane']:.4f}  kzz={r['kappa_zz']:.4f}  "
                  f"ratio={r['kzz_over_inplane']:.3f}  ({r['seconds']}s)", flush=True)

    json.dump(out, open("kappa_results_442.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print("\nwrote kappa_results_442.json")


if __name__ == "__main__":
    main()
