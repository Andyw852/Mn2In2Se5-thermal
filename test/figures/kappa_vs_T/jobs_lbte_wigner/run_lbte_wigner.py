#!/usr/bin/env python3
"""kappa(T) of the ordered and the VCA crystal with RTA, LBTE and the Wigner coherence term.

Why: every kappa(T) curve so far (ordered and VCA) is RTA population only. The
300 K values to quote add the LBTE correction (in-plane) and the coherence term,
but only for the ordered crystal at 300 K. With strong disorder scattering the
LBTE correction is expected to shrink and the coherence term to grow, so the
ordered-crystal factors cannot be reused for VCA. This driver computes, on one
mesh and for each case:

  rta        RTA population                (reference, same mesh)
  lbte       direct solution of the linearised BTE (phono3py is_LBTE)
  coherence  Wigner coherence term in the RTA (3.x conductivity_type="wigner";
             4.x transport_type="IBDB19" — the 4.4.0 SMM19 value is wrong for
             MnIn2Se4, REPORT_MnIn2Se4_20261115.md §7)

and writes result_<MAT>_<case>_m<N>.json with kappa per T and the ratios
LBTE/RTA and coherence/RTA, which are then applied to the converged RTA curves.
kappa_zz here is the plain mesh total (Gamma-Z artifact included); use the
ratios, not the absolute kappa_zz.

  run_lbte_wigner.py --material MnIn2Se4  --cases ordered vca            # 14^3
  run_lbte_wigner.py --material Mn2In2Se5 --cases ordered vca --x 0.19   # 16^3
  run_lbte_wigner.py ... --mesh 6 --temps 100 300                        # quick test (~2 min)
  run_lbte_wigner.py --material Mn2In2Se5 --cases vca --mesh 20 --steps rta --write-kappa
                                             # 20^3 RTA VCA curve at the literature x

Checks at 300 K (ordered): MnIn2Se4 14^3 LBTE in-plane 4.5219 (job 412627),
coherence 0.0908 / 0.0661 (3.24.0); Mn2In2Se5 16^3 coherence 0.1466 / 0.1144
(4.4.0 SMM19, job 419208; IBDB19 may differ by a few per cent).
"""
import argparse
import json
import os
import time

import numpy as np
import phono3py
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5
from phonopy.interface.vasp import read_vasp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
TEMPS = [30.0, 60.0, 100.0, 150.0, 200.0, 300.0, 400.0, 500.0, 700.0]
M_MN, M_IN = 54.938045, 114.818
G_SE = 4.63e-4                                   # natural Se, as in all other runs

CFG = {
    "MnIn2Se4": {"poscar": "test/MnIn2Se4/fit/POSCAR",
                 "fc2": "test/MnIn2Se4/fc2_442_phonon/fc2.hdf5",
                 "fc3": "test/MnIn2Se4/fit_cut3_6p45/fc3.hdf5",
                 "sc": [[4, 0, 0], [0, 4, 0], [0, 0, 1]], "psc": [[4, 0, 0], [0, 4, 0], [0, 0, 2]],
                 "pmat": [[1, 1, 2], [-1, -2, -1], [1, 1, 1]], "mesh": 14, "x": 0.44,
                 "x_source": "Range et al., Z. Naturforsch. 46b, 1122 (1991)"},
    "Mn2In2Se5": {"yaml": "Mn2In2Se5/fc/phono3py_disp.yaml",
                  "fc2": "Mn2In2Se5/fc/fc2.hdf5", "fc3": "Mn2In2Se5/fc/fc3.hdf5",
                  "pmat": "P", "mesh": 16, "x": 0.19,
                  "x_source": "Williams et al., J. Mater. Chem. C 12, 1753 (2024), ESI Table S7"},
}


def site_mix(material, x):
    """VCA site occupancies -> {element: (mean mass, mass variance g)}.
    x = fraction of all Mn on In sites (analysis/extrinsic/kappa_extrinsic.py convention).
    MnIn2Se4: Mn (octahedral) site (1-x) Mn + x In; In (tetrahedral, 2 per f.u.) site x/2 Mn + (1-x/2) In
              -> x = 0.44 gives 81.29 / 101.65 amu and g = 0.1337 / 0.0595 (the 20^3 VCA runs).
    Mn2In2Se5: Mn sites (1-x) Mn + x In; In sites x Mn + (1-x) In."""
    c_in_on_mn = x
    c_mn_on_in = x / 2 if material == "MnIn2Se4" else x
    out = {}
    for el, cmn in (("Mn", 1 - c_in_on_mn), ("In", c_mn_on_in)):
        cin = 1 - cmn
        m = cmn * M_MN + cin * M_IN
        g = cmn * (1 - M_MN / m) ** 2 + cin * (1 - M_IN / m) ** 2
        out[el] = (m, g)
    return out


def build(cfg, case, material, x, mesh):
    if "poscar" in cfg:
        uc = read_vasp(os.path.join(ROOT, cfg["poscar"]))
        kw = dict(supercell_matrix=cfg["sc"], phonon_supercell_matrix=cfg["psc"])
    else:
        p = phono3py.load(os.path.join(ROOT, cfg["yaml"]), produce_fc=False, log_level=0)
        uc = p.unitcell
        kw = dict(supercell_matrix=p.supercell_matrix)
    gl = None
    if case == "vca":
        mix = site_mix(material, x)
        masses = np.array(uc.masses, float)
        for i, s in enumerate(uc.symbols):
            if s in mix:
                masses[i] = mix[s][0]
        uc.masses = masses
        gl = [mix[s][1] if s in mix else G_SE for s in uc.symbols]
    ph3 = phono3py.Phono3py(uc, primitive_matrix=cfg["pmat"], log_level=0, **kw)
    ph3.fc3 = read_fc3_from_hdf5(os.path.join(ROOT, cfg["fc3"]))
    ph3.fc2 = read_fc2_from_hdf5(os.path.join(ROOT, cfg["fc2"]))
    ph3.mesh_numbers = [mesh] * 3
    ph3.init_phph_interaction()
    return ph3, gl, [round(float(m), 3) for m in uc.masses]


def xz(k, nT):
    k = np.asarray(k, float).reshape(nT, 6)
    return [round(float(0.5 * (r[0] + r[1])), 5) for r in k], [round(float(r[2]), 5) for r in k]


def run(ph3, gl, temps, step, write_kappa=False):
    kw = dict(temperatures=temps, log_level=0)
    if write_kappa and step == "rta":
        kw["write_kappa"] = True                 # per-mode kappa-m<N>.hdf5 for kappa_zz_decompose.py
    kw.update({"is_isotope": True} if gl is None else {"mass_variances": gl})
    major = int(phono3py.__version__.split(".")[0])
    if step == "rta":
        ph3.run_thermal_conductivity(**kw)
        return {"population": ph3.thermal_conductivity.kappa}
    if step == "lbte":
        ph3.run_thermal_conductivity(is_LBTE=True, **kw)
        return {"population": ph3.thermal_conductivity.kappa}
    if major < 4:                                # coherence, Wigner-RTA
        ph3.run_thermal_conductivity(conductivity_type="wigner", **kw)
        tc = ph3.thermal_conductivity
        return {"population": tc.kappa_P_RTA, "coherence": tc.kappa_C}
    ph3.run_thermal_conductivity(transport_type="IBDB19", **kw)
    ex = ph3.thermal_conductivity.get_extra_kappa_output()
    return {"population": ex["kappa_intra"], "coherence": ex["kappa_inter"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--material", required=True, choices=sorted(CFG))
    ap.add_argument("--cases", nargs="+", default=["ordered", "vca"], choices=["ordered", "vca"])
    ap.add_argument("--steps", nargs="+", default=["rta", "lbte", "coherence"], choices=["rta", "lbte", "coherence"])
    ap.add_argument("--mesh", type=int, default=None)
    ap.add_argument("--x", type=float, default=None, help="antisite fraction for VCA (default: literature value)")
    ap.add_argument("--temps", nargs="+", type=float, default=TEMPS)
    ap.add_argument("--write-kappa", action="store_true",
                    help="RTA step also writes the per-mode hdf5 (renamed kappa_<case>_rta_m<N>.hdf5)")
    a = ap.parse_args()
    cfg = CFG[a.material]
    mesh = a.mesh or cfg["mesh"]
    x = cfg["x"] if a.x is None else a.x
    nT = len(a.temps)
    for case in a.cases:
        out = {"material": a.material, "case": case, "mesh": mesh, "T_K": a.temps,
               "phono3py": phono3py.__version__, "fc3": cfg["fc3"], "fc2": cfg["fc2"]}
        if case == "vca":
            out["vca"] = {"x": x, "x_source": cfg["x_source"] if a.x is None else "user",
                          "site": {k: {"mass": round(v[0], 3), "g": round(v[1], 5)}
                                   for k, v in site_mix(a.material, x).items()}}
        for step in a.steps:
            ph3, gl, masses = build(cfg, case, a.material, x, mesh)
            t0 = time.time()
            print(f"start {a.material} {case} {step} {mesh}^3 {nT} T  phono3py {phono3py.__version__}", flush=True)
            res = run(ph3, gl, a.temps, step, a.write_kappa)
            src = "kappa-m{0}{0}{0}.hdf5".format(mesh)
            if a.write_kappa and step == "rta" and os.path.exists(src):
                os.replace(src, f"kappa_{case}_rta_m{mesh}.hdf5")
            out[step] = {k: dict(zip(("inplane", "kzz_total"), xz(v, nT))) for k, v in res.items()}
            out[step]["seconds"] = round(time.time() - t0, 1)
            out.setdefault("masses", masses)
            if gl is not None:
                out["mass_variances"] = [round(g, 6) for g in gl]
            print(f"  done in {out[step]['seconds']} s", flush=True)
        r = out.get("rta", {}).get("population")
        if r:                                     # ratios to apply to the converged RTA curves
            ratios = {}
            if "lbte" in out:
                lb = out["lbte"]["population"]
                ratios["lbte_over_rta"] = {c: [round(b / a_, 4) for a_, b in zip(r[c], lb[c])] for c in r}
            if "coherence" in out:
                co = out["coherence"]["coherence"]
                ratios["coherence_over_rta"] = {c: [round(b / a_, 4) for a_, b in zip(r[c], co[c])] for c in r}
            out["ratios"] = ratios
        fn = f"result_{a.material}_{case}_m{mesh}.json"
        with open(fn, "w") as f:
            json.dump(out, f, indent=1)
        print("wrote", fn, flush=True)


if __name__ == "__main__":
    main()
