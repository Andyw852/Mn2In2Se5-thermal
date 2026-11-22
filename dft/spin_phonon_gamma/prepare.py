#!/usr/bin/env python3
"""Set up the FM vs AFM Gamma-phonon test for Mn2In2Se5 (spin-lattice coupling screen).

Question: how much do the phonons depend on the Mn spin arrangement? The kappa calculation
assumes ferromagnetic (FM) order; the samples are paramagnetic at 300 K. If the force
constants hardly change between FM and antiferromagnetic (AFM) order, paramagnetic
force-constant disorder can be neglected; if they change by several per cent, it cannot.

Writes, under --out:

  FM/ref, AFM/ref                 9-atom primitive cell, undisplaced (energy, moments, residual force);
                                  LWAVE = .TRUE., its WAVECAR seeds the displaced runs
  FM/disp-001..054, AFM/disp-...  every atom displaced by +-0.01 A along a1, a2, a3 (no symmetry
                                  reduction, see below), ISTART = 1 from ../ref/WAVECAR
  energy_2x1/FM, stripeA, stripeB 18-atom 2x1x1 cell, in-plane stripe AFM (two interlayer variants)
                                  and FM in the same cell: total energies only
  FM/displacements.json, AFM/...  the displacement list read by analyse.py
  manifest.json                   configurations, MAGMOM, directory list

Spin configurations (atom order Mn1 Mn2 In In Se*5; Mn1 and Mn2 are the two Mn layers of one slab,
on either side of the central Se layer):
  FM      Mn1 +5, Mn2 +5                          (production setting)
  AFM     Mn1 +5, Mn2 -5                          (the two Mn layers antiparallel, each layer FM)
  stripeA 2x1x1: in each layer +- along a1, both layers in the same pattern (by cell index)
  stripeB 2x1x1: as stripeA with the Mn2 layer reversed

No symmetry reduction: phonopy's magnetic symmetry search fails for this cell (spglib
"spacegroup search failed", phonopy 4.4.0, tested), and treating FM and AFM identically lets the
splitting of the FM E doublets measure the numerical noise of the whole procedure.

DFT settings: those of the production S2 static run (Mn2In2Se5/data/dft_inputs/step2_static/INCAR:
PBE + D3(BJ), U(Mn) = 4 eV Dudarev, ENCUT 410, PREC Accurate, LREAL .FALSE., LASPH, ISMEAR 0 /
SIGMA 0.05, Gamma 8x8x2 for the 9-atom cell), plus the force-run settings of S4 (ISYM = 0,
ADDGRID), with EDIFF tightened to 1E-8 because 0.01 A displacements of soft modes give small
forces. The script checks the shared tags against the repository INCAR and stops on a mismatch.

  python prepare.py --repo <clone> --out <workdir> [--potcar <production POTCAR, Mn In Se>]
"""
import argparse
import json
import re
import shutil
from pathlib import Path

import numpy as np
import phonopy
from phonopy.interface.vasp import read_vasp, write_vasp
from phonopy.structure.atoms import PhonopyAtoms

DIST = 0.01
CONFIGS = {"FM": [5, 5], "AFM": [5, -5]}          # Mn1, Mn2

INCAR = """SYSTEM   = Mn2In2Se5 spin-phonon Gamma test: {tag}
ISTART   = {istart}
ICHARG   = {icharg}
GGA      = PE
IVDW     = 12
PREC     = Accurate
ENCUT    = 410
EDIFF    = 1E-8
LREAL    = .FALSE.
ADDGRID  = .TRUE.
LASPH    = .TRUE.
ALGO     = Normal
NELM     = 300
NELMIN   = 6
ISMEAR   = 0
SIGMA    = 0.05
IBRION   = -1
NSW      = 0
ISYM     = 0
LWAVE    = {lwave}
LCHARG   = .FALSE.
LORBIT   = 11
NCORE    = 4
KPAR     = 2
ISPIN    = 2
MAGMOM   = {magmom}
LMAXMIX  = 4
LDAU     = .TRUE.
LDAUTYPE = 2
LDAUL    = 2 -1 -1
LDAUU    = 4 0.0 0.0
LDAUJ    = 0.0 0.0 0.0
LDAUPRINT = 1
"""
# tags that must equal the production S2 static INCAR
SHARED = ["GGA", "IVDW", "PREC", "ENCUT", "LREAL", "LASPH", "ISMEAR", "SIGMA", "ISPIN",
          "LMAXMIX", "LDAU", "LDAUTYPE", "LDAUL", "LDAUU", "LDAUJ"]


def incar_tags(text):
    out = {}
    for line in text.splitlines():
        line = line.split("#")[0].split("!")[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip().upper()] = " ".join(v.split()).upper()
    return out


def magmom_str(m):
    """VASP MAGMOM: magnetic sites written out one by one, runs of zeros compressed (e.g. 5 -5 7*0)."""
    out, i = [], 0
    while i < len(m):
        if m[i] != 0:
            out.append(f"{m[i]:g}")
            i += 1
            continue
        j = i
        while j < len(m) and m[j] == 0:
            j += 1
        out.append(f"{j - i}*0" if j - i > 1 else "0")
        i = j
    return " ".join(out)


def kpoints(n):
    return "Gamma {0} {1} {2}\n0\nGamma\n   {0}   {1}   {2}\n0.0  0.0  0.0\n".format(*n)


def write_run(d, cell, magmom, tag, kmesh, istart, lwave, potcar):
    d.mkdir(parents=True, exist_ok=False)
    write_vasp(str(d / "POSCAR"), cell)
    (d / "INCAR").write_text(INCAR.format(tag=tag, istart=istart, icharg=0 if istart else 2,
                                          lwave=lwave, magmom=magmom_str(magmom)))
    (d / "KPOINTS").write_text(kpoints(kmesh))
    if potcar:
        shutil.copy(potcar, d / "POTCAR")


def supercell_2x1(prim, mn_pattern):
    """2x1x1 cell (2 a1, a2, a3), atoms grouped by species in the order of the primitive cell.
    mn_pattern[(mn_index, n1)] -> moment."""
    lat = np.array(prim.cell, float)
    lat2 = np.array([2 * lat[0], lat[1], lat[2]])
    syms, pos, mag = [], [], []
    for i, (s, f) in enumerate(zip(prim.symbols, prim.scaled_positions)):
        for n1 in (0, 1):
            syms.append(s)
            pos.append([(f[0] + n1) / 2.0, f[1], f[2]])
            mag.append(mn_pattern.get((i, n1), 0))
    return PhonopyAtoms(symbols=syms, cell=lat2, scaled_positions=pos), mag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="clone of Mn2In2Se5-thermal")
    ap.add_argument("--out", required=True, help="new work directory (must not exist)")
    ap.add_argument("--potcar", help="production POTCAR (species order Mn In Se); copied into every run")
    a = ap.parse_args()
    repo, out = Path(a.repo), Path(a.out)
    if out.exists():
        raise SystemExit(f"{out} exists; choose a new directory")

    # consistency with the production static INCAR
    prod = incar_tags((repo / "Mn2In2Se5/data/dft_inputs/step2_static/INCAR").read_text())
    mine = incar_tags(INCAR.format(tag="x", istart=0, icharg=2, lwave=".FALSE.", magmom="1"))
    bad = [k for k in SHARED if prod.get(k) != mine.get(k)]
    if bad:
        raise SystemExit("INCAR differs from production S2 in " + ", ".join(
            f"{k}: {prod.get(k)} vs {mine.get(k)}" for k in bad))
    if a.potcar:
        titles = re.findall(r"TITEL\s*=\s*\S+\s+(\S+)", Path(a.potcar).read_text())
        if [t.split("_")[0] for t in titles] != ["Mn", "In", "Se"]:
            raise SystemExit(f"POTCAR species {titles}, expected Mn, In, Se in that order")

    prim = read_vasp(str(repo / "Mn2In2Se5/POSCAR_relaxed"))
    assert list(prim.symbols) == ["Mn", "Mn", "In", "In", "Se", "Se", "Se", "Se", "Se"], prim.symbols
    manifest = {"distance_A": DIST, "kmesh_9": [8, 8, 2], "kmesh_18": [4, 8, 2], "configs": {}}

    for name, mn in CONFIGS.items():
        mag = mn + [0] * 7
        ph = phonopy.Phonopy(prim, np.eye(3, dtype=int), primitive_matrix="P", is_symmetry=False)
        ph.generate_displacements(distance=DIST, is_plusminus=True)
        disps = [{"number": int(d["number"]), "displacement": [float(x) for x in d["displacement"]]}
                 for d in ph.dataset["first_atoms"]]
        assert len(disps) == 54, len(disps)
        base = out / name
        write_run(base / "ref", prim, mag, f"{name} ref", (8, 8, 2), 0, ".TRUE.", a.potcar)
        dirs = []
        for k, sc in enumerate(ph.supercells_with_displacements, 1):
            d = base / f"disp-{k:03d}"
            write_run(d, sc, mag, f"{name} disp-{k:03d}", (8, 8, 2), 1, ".FALSE.", a.potcar)
            dirs.append(d.name)
        (base / "displacements.json").write_text(json.dumps(
            {"config": name, "magmom": mag, "dirs": dirs, "displacements": disps}, indent=1))
        manifest["configs"][name] = {"magmom": magmom_str(mag), "n_disp": len(dirs)}

    patterns = {
        "FM": {(0, 0): 5, (0, 1): 5, (1, 0): 5, (1, 1): 5},
        "stripeA": {(0, 0): 5, (0, 1): -5, (1, 0): 5, (1, 1): -5},
        "stripeB": {(0, 0): 5, (0, 1): -5, (1, 0): -5, (1, 1): 5},
    }
    for name, pat in patterns.items():
        cell, mag = supercell_2x1(prim, pat)
        write_run(out / "energy_2x1" / name, cell, mag, f"2x1 {name}", (4, 8, 2), 0, ".FALSE.", a.potcar)
        manifest["configs"][f"2x1_{name}"] = {"magmom": magmom_str(mag), "n_disp": 0}

    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    n = sum(v["n_disp"] + 1 for k, v in manifest["configs"].items() if not k.startswith("2x1")) + 3
    print(f"wrote {out}: {n} VASP runs (2 ref + 108 displaced, 9 atoms; 3 energies, 18 atoms)")
    for k, v in manifest["configs"].items():
        print(f"  {k:12s} MAGMOM = {v['magmom']:22s} displaced runs: {v['n_disp']}")
    if not a.potcar:
        print("no --potcar given: copy the production POTCAR (Mn In Se) into every run directory")


if __name__ == "__main__":
    main()
