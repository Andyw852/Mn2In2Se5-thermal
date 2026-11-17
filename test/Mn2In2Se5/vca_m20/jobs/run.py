#!/usr/bin/env python3
"""Mn2In2Se5: VCA (virtual crystal) vs ordered baseline, same settings as the MnIn2Se4 run.

VCA construction (same model as analysis/extrinsic/kappa_extrinsic.py, at the x the
kappa fit gives for this material, x = 0.13):

  Mn site : 0.87 Mn + 0.13 In   -> m = 62.73 amu, g = 0.1031
  In site : 0.87 In + 0.13 Mn   -> m = 107.03 amu, g = 0.0354   (f = x*n_Mn/n_In = x)
  Se      : natural             -> m = 78.96 amu, g = 4.63e-4

Primitive cell is the 9-atom cell order  Mn Mn In In Se Se Se Se Se  (with
primitive_matrix="P", as in every kappa run for this material), so the mass-variance
list is [g_Mn, g_Mn, g_In, g_In, g_Se x5].

Two runs at 20x20x20, 300 K, FC3 = 6.51 A, FC2 = 4x4x1:
  baseline : real masses, is_isotope=True
  vca      : averaged masses, mass_variances, is_isotope=False
"""
import json
import os
import time

import numpy as np
import phono3py
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

BASE = os.path.dirname(os.path.abspath(__file__))
MESH = [20, 20, 20]
G_SE = 4.63e-4
MASS = {"Mn": 62.73, "In": 107.03}      # VCA averages at x = 0.13
G = {"Mn": 0.1031, "In": 0.0354, "Se": G_SE}


def build(unitcell):
    p = phono3py.load(f"{BASE}/in/phono3py_disp.yaml", produce_fc=False, log_level=0)
    ph3 = phono3py.Phono3py(unitcell, supercell_matrix=np.array(p.supercell_matrix, int),
                            primitive_matrix="P", log_level=0)
    ph3.fc3 = read_fc3_from_hdf5(f"{BASE}/in/fc3.hdf5")
    ph3.fc2 = read_fc2_from_hdf5(f"{BASE}/in/fc2.hdf5")
    ph3.mesh_numbers = MESH
    return ph3


def snap(k):
    a = np.asarray(k).reshape(-1, 6)[0]
    return {"inplane": float(0.5 * (a[0] + a[1])), "kappa_zz": float(a[2])}


p = phono3py.load(f"{BASE}/in/phono3py_disp.yaml", produce_fc=False, log_level=0)
uc0 = p.unitcell
syms = list(uc0.symbols)
glist = [G[s] for s in syms]
print("primitive order:", syms, flush=True)
print("mass variances :", glist, flush=True)

out = {"material": "Mn2In2Se5", "mesh": MESH, "fc2": "4x4x1", "fc3_cutoff_A": 6.51,
       "T_K": 300.0, "vca": {"x": 0.13, "mass_Mn": 62.73, "mass_In": 107.03,
                             "mass_variances": glist}}

# ---- baseline
print("start baseline", flush=True)
ph3 = build(uc0)
ph3.init_phph_interaction()
t0 = time.time()
ph3.run_thermal_conductivity(temperatures=[300.0], is_isotope=True,
                             write_kappa=True, log_level=0)
out["baseline"] = snap(ph3.thermal_conductivity.kappa)
out["baseline"]["seconds"] = round(time.time() - t0, 1)
if os.path.exists("kappa-m202020.hdf5"):
    os.replace("kappa-m202020.hdf5", "kappa_baseline.hdf5")
print("baseline", out["baseline"], flush=True)
json.dump(out, open(f"{BASE}/result.json", "w"), indent=1)

# ---- VCA
print("start VCA", flush=True)
uc = phono3py.load(f"{BASE}/in/phono3py_disp.yaml", produce_fc=False, log_level=0).unitcell
masses = np.array(uc.masses, dtype=float)
for i, s in enumerate(uc.symbols):
    if s in MASS:
        masses[i] = MASS[s]
uc.masses = masses
ph3v = build(uc)
ph3v.init_phph_interaction()
t0 = time.time()
ph3v.run_thermal_conductivity(temperatures=[300.0], is_isotope=False,
                              mass_variances=glist, write_kappa=True, log_level=0)
out["vca"] = snap(ph3v.thermal_conductivity.kappa)
out["vca"]["seconds"] = round(time.time() - t0, 1)
if os.path.exists("kappa-m202020.hdf5"):
    os.replace("kappa-m202020.hdf5", "kappa_vca.hdf5")
print("vca", out["vca"], flush=True)
json.dump(out, open(f"{BASE}/result.json", "w"), indent=1)
print("suppression in-plane %.2fx kappa_zz %.2fx" % (
    out["baseline"]["inplane"] / out["vca"]["inplane"],
    out["baseline"]["kappa_zz"] / out["vca"]["kappa_zz"]), flush=True)
print("DONE", flush=True)
