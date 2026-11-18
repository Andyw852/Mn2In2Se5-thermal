#!/usr/bin/env python3
"""Extra random-displacement frames for the existing 4x4x1 FC3 datasets.

Purpose: allow an FC3 fit beyond the current 6.45 / 6.51 A cutoff *without a
larger supercell*. Only use it if the refit of the existing frames (README,
step 1) shows that the data/parameter ratio is too low at 7-8 A.

The cutoff limit of the 4x4x1 cells
-----------------------------------
Two criteria are printed for any supercell:
  * parallelepiped rule (AutoZT cut3_safe_cutoff): half the shortest interplanar
    distance - 0.1 A. It depends on the basis chosen for a non-orthogonal cell.
    It is only needed by codes that find periodic images by rounding fractional
    coordinates in that basis.
  * minimum-image rule: half the shortest supercell lattice vector - 0.1 A.
    It is basis-independent, and it is the real condition: below it, no atom
    pair has two periodic images inside the cutoff, so nothing aliases.
For the 4x4x1 cells: Mn2In2Se5 6.61 / 7.96 A, MnIn2Se4 6.51 / 8.01 A.
4x4x2 gives 6.61 / 7.96 and 6.55 / 8.05, i.e. no gain over 4x4x1 under either rule.

Protocol (identical to the existing datasets, so the new frames can be added
to them):
  * supercell and cell: the existing 4x4x1 ones (Mn2In2Se5/fc/phono3py_disp.yaml,
    test/MnIn2Se4/fit/POSCAR)
  * displacements: every atom, Gaussian, sigma = 0.017 A per component (measured
    on both existing datasets), no +/- pairs
  * INCAR / KPOINTS: the existing S4 inputs verbatim, only SYSTEM changed

  python make_frames.py --material Mn2In2Se5 --out new_frames --first 55 --count 30
--first must be past the last existing frame (Mn2In2Se5: 0..54, MnIn2Se4: check
its step4_disp). Frame k uses seed SEED + k. Existing directories are never
overwritten.
"""
import argparse
import itertools
import json
import os

import numpy as np
import phono3py
from phonopy.interface.calculator import write_crystal_structure
from phonopy.interface.vasp import read_vasp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MAT = {
    "Mn2In2Se5": dict(
        cell=("yaml", os.path.join(ROOT, "Mn2In2Se5", "fc", "phono3py_disp.yaml")),
        incar=os.path.join(ROOT, "Mn2In2Se5", "data", "dft_inputs", "step4_disp", "disp-00000", "INCAR"),
        kpoints=os.path.join(ROOT, "Mn2In2Se5", "data", "dft_inputs", "step4_disp", "disp-00000", "KPOINTS"),
        current_cutoff=6.51, pmat="P"),
    "MnIn2Se4": dict(
        cell=("poscar", os.path.join(ROOT, "test", "MnIn2Se4", "fit", "POSCAR")),
        incar=os.path.join(ROOT, "MnIn2Se4", "redo_441", "dft", "s4_disp_inputs", "INCAR_disp-00000"),
        kpoints=os.path.join(ROOT, "MnIn2Se4", "redo_441", "dft", "s4_disp_inputs", "KPOINTS_disp-00000"),
        current_cutoff=6.45, pmat=[[1, 1, 2], [-1, -2, -1], [1, 1, 1]]),
}
SIGMA = 0.017
SEED = 20261118


def ceilings(lat):
    v = abs(np.linalg.det(lat))
    d = [v / np.linalg.norm(np.cross(lat[(i + 1) % 3], lat[(i + 2) % 3])) for i in range(3)]
    t = min(np.linalg.norm(np.array(c) @ lat)
            for c in itertools.product(range(-3, 4), repeat=3) if c != (0, 0, 0))
    return min(d) / 2 - 0.1, t / 2 - 0.1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--material", required=True, choices=list(MAT))
    ap.add_argument("--out", default="new_frames")
    ap.add_argument("--first", type=int, required=True, help="first new frame index (past the existing ones)")
    ap.add_argument("--count", type=int, default=30)
    ap.add_argument("--dim", type=int, nargs=3, default=[4, 4, 1])
    a = ap.parse_args()
    m = MAT[a.material]

    if m["cell"][0] == "yaml":
        unit = phono3py.load(m["cell"][1], produce_fc=False, log_level=0).unitcell
    else:
        unit = read_vasp(m["cell"][1])
    ph3 = phono3py.Phono3py(unit, supercell_matrix=np.diag(a.dim), primitive_matrix=m["pmat"], log_level=0)
    sc = ph3.supercell
    par, mic = ceilings(np.array(sc.cell))
    print(f"{a.material} {a.dim}: {len(sc.symbols)} atoms; cutoff limit {par:.2f} A (parallelepiped rule), "
          f"{mic:.2f} A (minimum-image rule); current fit {m['current_cutoff']} A")

    os.makedirs(a.out, exist_ok=True)
    write_crystal_structure(os.path.join(a.out, "SPOSCAR"), sc, interface_mode="vasp")
    incar = open(m["incar"], encoding="utf-8").read()
    kpts = open(m["kpoints"], encoding="utf-8").read()
    man_path = os.path.join(a.out, "frames_manifest.json")
    man = json.load(open(man_path)) if os.path.exists(man_path) else {
        "material": a.material, "seed_base": SEED, "sigma_A": SIGMA, "dim": a.dim,
        "limit_parallelepiped_A": round(par, 3), "limit_minimum_image_A": round(mic, 3), "frames": {}}

    pos0 = sc.positions
    for k in range(a.first, a.first + a.count):
        tag = f"disp-{k:05d}"
        d = os.path.join(a.out, tag)
        if os.path.exists(os.path.join(d, "POSCAR")):
            raise SystemExit(f"{tag} exists -- refusing to overwrite")
        u = np.random.default_rng(SEED + k).normal(0.0, SIGMA, size=pos0.shape)
        cell = sc.copy()
        cell.positions = pos0 + u
        os.makedirs(d)
        write_crystal_structure(os.path.join(d, "POSCAR"), cell, interface_mode="vasp")
        lines = [f"SYSTEM  = {a.material}_441 fc3-extend {tag}" if ln.strip().upper().startswith("SYSTEM")
                 else ln for ln in incar.splitlines()]
        open(os.path.join(d, "INCAR"), "w").write("\n".join(lines) + "\n")
        open(os.path.join(d, "KPOINTS"), "w").write(kpts)
        n = np.linalg.norm(u, axis=1)
        man["frames"][tag] = {"seed": SEED + k, "mean_u": round(float(n.mean()), 4), "max_u": round(float(n.max()), 4)}
    json.dump(man, open(man_path, "w"), indent=1)
    print(f"wrote {a.count} frames {a.first}..{a.first + a.count - 1} in {a.out}/")


if __name__ == "__main__":
    main()
