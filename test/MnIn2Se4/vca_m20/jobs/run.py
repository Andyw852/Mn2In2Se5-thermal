#!/usr/bin/env python3
"""MnIn2Se4: VCA (virtual crystal) vs ordered baseline, at the largest FC3 cutoff.

The VCA recipe is the user's (E:\...\Mn2In2Se5\fix\run_phono3py_kappa_vca.sh):

  i = 0.44, i.e. octahedral site 0.56 Mn + 0.44 In, tetrahedral 0.22 Mn + 0.78 In
    Mn (octahedral) mass 81.29 amu, g = 0.1337
    In (tetrahedral) mass 101.65 amu, g = 0.0595
    Se natural isotope g = 4.63e-4 (phonopy isotope data; was 0.0003 in the first run)

The primitive cell is 7 atoms in the order  Mn In In Se Se Se Se  (checked against
the disp yaml), so the mass-variance list is
  [0.1337, 0.0595, 0.0595, 4.63e-4, 4.63e-4, 4.63e-4, 4.63e-4]

Two runs:
  baseline : original masses, is_isotope=True   (natural isotopes only)
  vca      : averaged masses, mass_variances    (disorder as mass variance)
Both at the largest available FC3 cutoff (6.45 A), FC2 4x4x2, T = 300 K, mesh 20^3
(the mesh the user's own script used). Per-mode output is written so the
scattering-rate figure can be drawn with analysis/plotting/plot_scattering_compare.py.
"""
import json
import os
import time

import numpy as np
import phono3py
from phonopy.interface.vasp import read_vasp
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

BASE = os.path.dirname(os.path.abspath(__file__))
MESH = [20, 20, 20]
MASS = {"Mn": 81.29, "In": 101.65}          # VCA averages at i = 0.44
G_SE = 4.63e-4                          # natural Se, same value is_isotope=True uses in the baseline
G = [0.1337, 0.0595, 0.0595, G_SE, G_SE, G_SE, G_SE]   # Mn In In Se Se Se Se


def build():
    uc = read_vasp(f"{BASE}/in/POSCAR")
    ph3 = phono3py.Phono3py(
        uc,
        supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 1]],
        phonon_supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 2]],
        primitive_matrix=[[1, 1, 2], [-1, -2, -1], [1, 1, 1]],
        log_level=0,
    )
    ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/in/fc3_645.hdf5")
    ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2_442.hdf5")
    ph3.mesh_numbers = MESH
    return ph3


def snap(k):
    a = np.asarray(k).reshape(-1, 6)[0]
    return {"inplane": float(0.5 * (a[0] + a[1])), "kappa_zz": float(a[2])}


meta = {"material": "MnIn2Se4", "mesh": MESH, "fc2": "4x4x2", "fc3_cutoff_A": 6.45,
        "T_K": 300.0, "vca": {"i": 0.44, "mass_Mn": 81.29, "mass_In": 101.65,
                              "mass_variances": G}}
out = dict(meta)

# ---- baseline: ordered masses + natural isotopes
print("start baseline (natural isotopes, ordered masses)", flush=True)
ph3 = build()
ph3.init_phph_interaction()
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_isotope=True,
                             write_kappa=True, log_level=0)
out["baseline"] = snap(ph3.thermal_conductivity.kappa)
out["baseline"]["seconds"] = round(time.time() - t0, 1)
print("baseline", out["baseline"], flush=True)
os.rename("kappa-m202020.hdf5", "kappa_baseline.hdf5") if os.path.exists("kappa-m202020.hdf5") else None
json.dump(out, open(f"{BASE}/result.json", "w"), indent=1)

# ---- VCA: averaged masses + mass-variance scattering
print("start VCA (averaged masses + mass variances)", flush=True)
uc = read_vasp(f"{BASE}/in/POSCAR")
mass_arr = np.array(uc.masses, dtype=float)
for i, s in enumerate(uc.symbols):
    if s in MASS:
        mass_arr[i] = MASS[s]
uc.masses = mass_arr
print("  VCA masses:", list(zip(uc.symbols, np.round(mass_arr, 3))), flush=True)
ph3v = phono3py.Phono3py(
    uc,
    supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 1]],
    phonon_supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 2]],
    primitive_matrix=[[1, 1, 2], [-1, -2, -1], [1, 1, 1]],
    log_level=0,
)
ph3v.fc3 = read_fc3_from_hdf5(f"{BASE}/in/fc3_645.hdf5")
ph3v.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2_442.hdf5")
ph3v.mesh_numbers = MESH
ph3v.init_phph_interaction()
t0 = time.time()
ph3v.run_thermal_conductivity(temperatures=[300.0], is_isotope=False,
                              mass_variances=G, write_kappa=True, log_level=0)
out["vca"] = snap(ph3v.thermal_conductivity.kappa)
out["vca"]["seconds"] = round(time.time() - t0, 1)
print("vca", out["vca"], flush=True)
if os.path.exists("kappa-m202020.hdf5"):
    os.replace("kappa-m202020.hdf5", "kappa_vca.hdf5")
json.dump(out, open(f"{BASE}/result.json", "w"), indent=1)
print("suppression in-plane %.2fx  kappa_zz %.2fx" % (
    out["baseline"]["inplane"] / out["vca"]["inplane"],
    out["baseline"]["kappa_zz"] / out["vca"]["kappa_zz"]), flush=True)
print("DONE", flush=True)
