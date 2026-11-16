#!/usr/bin/env python3
"""Compute kappa with phono3py directly, bypassing AutoZT: fc2(4x4x2) + fc3(4x4x1).

Why the supercells can be mixed
-------------------------------
phono3py allows fc2 to use a larger supercell than fc3 (`phonon_supercell_matrix`).
Here fc3 is 4x4x1 (112 atoms) and fc2 is 4x4x2 (224 atoms).

Why fc2 must be replaced
------------------------
The old 4x4x1 fc2 gives a transverse sound velocity of only 322.9 m/s along the
layer normal (DFT elastic constants give 1338.3), so the phonon dispersion was
artificially softened -> kappa_zz was strongly underestimated. The same branch
from the 4x4x2 fc2 gives 1348.3 m/s (+0.75%), i.e. a correct dispersion.

Usage
-----
    python kappa_driver_442.py                      # both fc2 cases, RTA + LBTE, 11^3
    python kappa_driver_442.py --rta                # RTA only
    python kappa_driver_442.py --lbte               # LBTE only
    python kappa_driver_442.py --mesh 18 18 18 --rta --write-kappa
    python kappa_driver_442.py --rta --boundary-mfp 1.0     # Casimir, micrometres

The kappa-*.hdf5 written by `--write-kappa` is what `kappa_zz_decompose.py`
needs; it lands in the current directory, so run each mesh in its own directory.

Output: kappa_results_442.json (override with --out).
"""
import argparse
import json
import os
import time

import numpy as np
import phono3py
from phonopy.interface.vasp import read_vasp
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

# 2026-11-16, explicit user instruction: no phono3py on this workstation (the 12
# cores get saturated and the Windows host plus the dsh web UI freeze). To
# re-enable local runs: rm ~/.no-local-compute. The file does not exist on a
# cluster, so cluster jobs are unaffected.
if os.path.exists(os.path.expanduser("~/.no-local-compute")):
    raise SystemExit("[BLOCKED] local computation is disabled (~/.no-local-compute exists); "
                     "run phono3py on a cluster")

# --- paths (following the repo test/MnIn2Se4 layout; pip install phonopy phono3py first) ---
ROOT = os.path.dirname(os.path.abspath(__file__))
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
DEFAULT_MESH = [11, 11, 11]
TEMPS = [300.0]


def run(fc2_path, phonon_supercell, lbte, mesh, boundary_mfp, write_kappa, fc3_path):
    ph3 = phono3py.Phono3py(
        read_vasp(POSCAR),
        supercell_matrix=SUPERCELL_FC3,
        phonon_supercell_matrix=phonon_supercell,
        primitive_matrix=PRIMITIVE_MATRIX,
        log_level=0,
    )
    ph3.fc3 = read_fc3_from_hdf5(fc3_path)
    ph3.fc2 = read_fc2_from_hdf5(fc2_path)
    ph3.mesh_numbers = mesh
    ph3.init_phph_interaction()
    t0 = time.time()
    ph3.run_thermal_conductivity(
        temperatures=TEMPS, is_LBTE=lbte, is_isotope=True,
        boundary_mfp=boundary_mfp, write_kappa=write_kappa, log_level=0,
    )
    dt = time.time() - t0
    k = np.asarray(ph3.thermal_conductivity.kappa)[0, 0]   # (6,) Voigt
    inplane = 0.5 * (k[0] + k[1])
    return {
        "method": "LBTE" if lbte else "RTA",
        "mesh": " ".join(map(str, mesh)),
        "boundary_mfp_um": boundary_mfp,
        "fc3": os.path.basename(fc3_path),
        "T_K": TEMPS,
        "kappa_xx": float(k[0]), "kappa_yy": float(k[1]), "kappa_zz": float(k[2]),
        "inplane": float(inplane),
        "kzz_over_inplane": float(k[2] / inplane),
        "seconds": round(dt, 1),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rta", action="store_true", help="run RTA only")
    ap.add_argument("--lbte", action="store_true", help="run LBTE only")
    ap.add_argument("--mesh", nargs=3, type=int, default=DEFAULT_MESH,
                    metavar=("N1", "N2", "N3"), help="q-mesh (default 11 11 11)")
    ap.add_argument("--fc2", choices=["4x4x1", "4x4x2", "both"], default="both",
                    help="which fc2 supercell to use (default both)")
    ap.add_argument("--boundary-mfp", type=float, default=None,
                    help="Casimir boundary mean free path in MICROMETRES")
    ap.add_argument("--write-kappa", action="store_true",
                    help="write kappa-*.hdf5 (needed by kappa_zz_decompose.py)")
    ap.add_argument("--fc3", default=FC3, help="fc3.hdf5 (default: production 6.45 A)")
    ap.add_argument("--out", default=f"{ROOT}/kappa_results_442.json")
    a = ap.parse_args()

    do_rta = a.rta or not (a.rta or a.lbte)
    do_lbte = a.lbte or not (a.rta or a.lbte)

    cases = []
    if a.fc2 in ("4x4x1", "both"):
        cases.append(("4x4x1_old", FC2_441, PHONON_SUPERCELL["4x4x1"]))
    if a.fc2 in ("4x4x2", "both"):
        cases.append(("4x4x2_new", FC2_442, PHONON_SUPERCELL["4x4x2"]))

    out = []
    for name, fc2, psc in cases:
        for lbte in ([False] if do_rta else []) + ([True] if do_lbte else []):
            r = run(fc2, psc, lbte, a.mesh, a.boundary_mfp, a.write_kappa, a.fc3)
            r["fc2_case"] = name
            r["fc2_path"] = fc2
            out.append(r)
            print(f"[{name}] {r['method']:<4} mesh={r['mesh']} "
                  f"boundary_mfp={r['boundary_mfp_um']} 300K  "
                  f"in-plane={r['inplane']:.4f}  kzz={r['kappa_zz']:.4f}  "
                  f"ratio={r['kzz_over_inplane']:.3f}  ({r['seconds']}s)", flush=True)

    json.dump(out, open(a.out, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
