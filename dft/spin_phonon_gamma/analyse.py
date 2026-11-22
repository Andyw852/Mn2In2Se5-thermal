#!/usr/bin/env python3
"""Analyse the FM vs AFM Gamma-phonon test written by prepare.py.

  python analyse.py --repo <clone> --work <workdir>       # after all VASP runs
  python analyse.py --repo <clone> --selftest             # before submitting: checks this script

Reads only OUTCAR files (energy(sigma->0), TOTAL-FORCE, magnetization (x), iteration count).

Per configuration (FM, AFM): Gamma force constants of the 9-atom cell from the 54 displaced runs
(forces minus the ref forces, central differences, translational invariance imposed), Gamma
frequencies and eigenvectors. Then:

  - checks: every run finished, below NELM, Mn moments |m| > 4 muB with the intended signs;
  - noise:  splitting of the FM E doublets (exactly degenerate by symmetry; no symmetry is used,
            so the splitting is the numerical noise of the procedure);
  - FM against production: the same Gamma frequencies from the production 4x4x1 FC2
            (Mn2In2Se5/fc/fc2.hdf5), which should agree to about 1-2 %;
  - AFM against FM: each FM optical mode i gets the overlap-weighted AFM frequency
            w~_i = sum_j |<e_i(FM)|e_j(AFM)>|^2 w_j (robust against degenerate pairs and
            reordering), delta_i = w~_i / w_i - 1;
  - energies: E(AFM) - E(FM) per Mn (9-atom cell), E(stripe) - E(FM) per Mn (2x1x1 cell);
  - verdict (thresholds fixed before the runs):
        weak          max|delta| <= max(2 %, 3 x noise) and rms <= 1 %
        strong        max|delta| >= 5 % or rms >= 3 %
        intermediate  otherwise
    and a warning when the in-plane stripe energy exceeds the interlayer AFM energy by more than
    2x: the 9-atom AFM then probes the weaker exchange path, and a "weak" result is not conclusive.

Writes <workdir>/results.json.
"""
import argparse
import json
import re
from pathlib import Path

import h5py
import numpy as np
import phonopy
from phonopy.interface.vasp import read_vasp

NELM = 300
N_ACOUSTIC = 3


# ----------------------------------------------------------------------------- OUTCAR
def read_outcar(path, natom):
    txt = Path(path).read_text(errors="replace")
    out = {"finished": "General timing" in txt}
    e = re.findall(r"energy\(sigma->0\)\s*=\s*(-?\d+\.\d+)", txt)
    out["energy_eV"] = float(e[-1]) if e else None
    it = [int(m) for m in re.findall(r"Iteration\s+\d+\(\s*(\d+)\)", txt)]
    out["n_elec_steps"] = max(it) if it else None
    m = re.findall(r"number of electron\s+\S+\s+magnetization\s+(-?\d+\.\d+)", txt)
    out["total_magnetization"] = float(m[-1]) if m else None
    blocks = txt.split("magnetization (x)")
    if len(blocks) > 1:
        body = blocks[-1].split("tot", 1)[1] if "tot" in blocks[-1] else blocks[-1]
        rows = re.findall(r"^\s+\d+((?:\s+-?\d+\.\d+){4,5})\s*$", body, re.M)   # s p d (f) tot
        out["site_moments"] = [float(r.split()[-1]) for r in rows[:natom]]
    fb = txt.split("TOTAL-FORCE (eV/Angst)")
    if len(fb) > 1:
        lines = fb[-1].splitlines()[2:2 + natom]
        out["forces"] = np.array([[float(x) for x in ln.split()[3:6]] for ln in lines])
    return out


def check_run(r, magmom, label, problems):
    if not r["finished"]:
        problems.append(f"{label}: not finished")
    if r["n_elec_steps"] is not None and r["n_elec_steps"] >= NELM:
        problems.append(f"{label}: hit NELM ({r['n_elec_steps']})")
    sm = r.get("site_moments")
    if sm is None or "forces" not in r and "disp" in label:
        problems.append(f"{label}: no site moments / forces in OUTCAR")
        return
    for i, m0 in enumerate(magmom):
        if m0 != 0 and (abs(sm[i]) < 4.0 or np.sign(sm[i]) != np.sign(m0)):
            problems.append(f"{label}: atom {i + 1} moment {sm[i]:+.3f}, intended {m0:+g}")


# ----------------------------------------------------------------------------- phonons
def gamma_modes(cell, disps, forces):
    """Gamma frequencies (THz) and eigenvectors of the 9-atom cell from a displacement set."""
    ph = phonopy.Phonopy(cell, np.eye(3, dtype=int), primitive_matrix="P", is_symmetry=False)
    ph.dataset = {"natom": len(cell.symbols),
                  "first_atoms": [{"number": d["number"], "displacement": d["displacement"],
                                   "forces": f} for d, f in zip(disps, forces)]}
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    return modes_of(ph)


def modes_of(ph):
    ph.run_qpoints([[0, 0, 0]], with_eigenvectors=True)
    q = ph.qpoints
    return np.array(q.frequencies[0]), np.array(q.eigenvectors[0])


def production_modes(repo):
    cell = read_vasp(str(repo / "Mn2In2Se5/POSCAR_relaxed"))
    with h5py.File(repo / "Mn2In2Se5/fc/fc2.hdf5") as f:
        fc = f["fc2"][()] if "fc2" in f else f["force_constants"][()]
    ph = phonopy.Phonopy(cell, np.diag([4, 4, 1]), primitive_matrix="P")
    ph.force_constants = fc
    return ph, modes_of(ph)


def match(w_a, e_a, w_b, e_b):
    """Overlap-weighted frequency of set b for each optical mode of set a, and the largest overlap."""
    o = np.abs(e_a.conj().T @ e_b) ** 2                     # (mode a, mode b)
    opt = slice(N_ACOUSTIC, None)
    ob = o[opt, opt]
    ob = ob / ob.sum(axis=1, keepdims=True)
    return ob @ w_b[opt], ob.max(axis=1)


def doublet_noise(w_fm, e_fm, w_p, e_p):
    """Splitting of the FM modes that belong to the production E doublets."""
    opt = np.arange(N_ACOUSTIC, len(w_p))
    pairs = [(i, i + 1) for i in opt[:-1] if abs(w_p[i + 1] - w_p[i]) < 1e-3]
    o = np.abs(e_p.conj().T @ e_fm) ** 2                    # (production, FM)
    out = []
    for i, j in pairs:
        weight = o[i] + o[j]
        a, b = np.argsort(weight)[-2:]
        out.append({"production_THz": round(float(w_p[i]), 4),
                    "fm_THz": sorted([round(float(w_fm[a]), 4), round(float(w_fm[b]), 4)]),
                    "rel_split": float(abs(w_fm[a] - w_fm[b]) / (0.5 * (w_fm[a] + w_fm[b])))})
    return out


def compare(w_fm, e_fm, w_afm, e_afm, w_p, e_p):
    wt, ov = match(w_fm, e_fm, w_afm, e_afm)
    wf = w_fm[N_ACOUSTIC:]
    delta = wt / wf - 1
    wp_on_fm, _ = match(w_fm, e_fm, w_p, e_p)
    dprod = wf / wp_on_fm - 1
    soft = wf < 3.0
    noise = doublet_noise(w_fm, e_fm, w_p, e_p)
    return {
        "acoustic_FM_THz": [round(float(x), 4) for x in w_fm[:N_ACOUSTIC]],
        "acoustic_AFM_THz": [round(float(x), 4) for x in w_afm[:N_ACOUSTIC]],
        "modes": [{"fm_THz": round(float(a), 4), "afm_THz": round(float(b), 4), "delta": float(d),
                   "max_overlap": round(float(o), 3), "production_THz": round(float(p), 4),
                   "fm_vs_production": float(dp)}
                  for a, b, d, o, p, dp in zip(wf, wt, delta, ov, wp_on_fm, dprod)],
        "rms_delta": float(np.sqrt(np.mean(delta ** 2))),
        "max_abs_delta": float(np.max(np.abs(delta))),
        "rms_delta_below_3THz": float(np.sqrt(np.mean(delta[soft] ** 2))) if soft.any() else None,
        "max_abs_delta_below_3THz": float(np.max(np.abs(delta[soft]))) if soft.any() else None,
        "noise_doublets": noise,
        "noise_max_rel_split": max((n["rel_split"] for n in noise), default=0.0),
        "max_abs_fm_vs_production_above_1THz": float(np.max(np.abs(dprod[wf > 1.0]))),
    }


def verdict(c, energies):
    noise = c["noise_max_rel_split"]
    mx, rms = c["max_abs_delta"], c["rms_delta"]
    if mx <= max(0.02, 3 * noise) and rms <= 0.01:
        v = "weak"
    elif mx >= 0.05 or rms >= 0.03:
        v = "strong"
    else:
        v = "intermediate"
    notes = []
    if noise > 0.01:
        notes.append(f"numerical noise {100 * noise:.1f} % is large: tighten EDIFF or check the runs before using the result")
    if c["max_abs_fm_vs_production_above_1THz"] > 0.03:
        notes.append("FM differs from production by more than 3 % above 1 THz: check the setup before interpreting")
    if energies:
        ea = abs(energies.get("AFM_minus_FM_meV_per_Mn") or 0)
        es = [abs(energies[k]) for k in ("stripeA_minus_FM_meV_per_Mn", "stripeB_minus_FM_meV_per_Mn")
              if energies.get(k) is not None]
        if es and min(es) > 2 * ea:
            notes.append("in-plane stripe energy > 2x interlayer AFM energy: the 9-atom AFM probes the weaker "
                         "exchange path; a 'weak' result is then not conclusive (next: stripe phonons, 18 atoms)")
    return v, notes


# ----------------------------------------------------------------------------- main
def analyse(repo, work):
    prim = read_vasp(str(repo / "Mn2In2Se5/POSCAR_relaxed"))
    nat = len(prim.symbols)
    problems, runs, modes = [], {}, {}
    for name in ("FM", "AFM"):
        spec = json.loads((work / name / "displacements.json").read_text())
        mag = spec["magmom"]
        ref = read_outcar(work / name / "ref/OUTCAR", nat)
        check_run(ref, mag, f"{name}/ref", problems)
        forces = []
        for d in spec["dirs"]:
            r = read_outcar(work / name / d / "OUTCAR", nat)
            check_run(r, mag, f"{name}/{d}", problems)
            forces.append(r["forces"] - ref["forces"])
        runs[name] = {"energy_eV": ref["energy_eV"], "total_magnetization": ref["total_magnetization"],
                      "site_moments_Mn": ref["site_moments"][:2],
                      "max_residual_force_eV_A": float(np.abs(ref["forces"]).max())}
        modes[name] = gamma_modes(prim, spec["displacements"], forces)
    e2 = {}
    for name in ("FM", "stripeA", "stripeB"):
        p = work / "energy_2x1" / name / "OUTCAR"
        r = read_outcar(p, 2 * nat)
        if not r["finished"]:
            problems.append(f"energy_2x1/{name}: not finished")
        e2[name] = r
        runs[f"2x1_{name}"] = {"energy_eV": r["energy_eV"], "total_magnetization": r["total_magnetization"],
                               "site_moments_Mn": (r.get("site_moments") or [None] * 4)[:4]}
    energies = {
        "AFM_minus_FM_meV_per_Mn": 1e3 * (runs["AFM"]["energy_eV"] - runs["FM"]["energy_eV"]) / 2,
        "stripeA_minus_FM_meV_per_Mn": 1e3 * (e2["stripeA"]["energy_eV"] - e2["FM"]["energy_eV"]) / 4,
        "stripeB_minus_FM_meV_per_Mn": 1e3 * (e2["stripeB"]["energy_eV"] - e2["FM"]["energy_eV"]) / 4,
        "FM_2x1_minus_2xFM_9_meV": 1e3 * (e2["FM"]["energy_eV"] - 2 * runs["FM"]["energy_eV"]),
    }
    _, (w_p, e_p) = production_modes(repo)
    c = compare(*modes["FM"], *modes["AFM"], w_p, e_p)
    v, notes = verdict(c, energies)
    return {"problems": problems, "runs": runs, "energies": energies, "comparison": c,
            "verdict": v, "notes": notes}


def report(res):
    c = res["comparison"]
    if res.get("problems"):
        print("PROBLEMS (fix before interpreting):")
        for p in res["problems"]:
            print("  -", p)
    for k, r in res.get("runs", {}).items():
        print(f"  {k:12s} E = {r['energy_eV']}  M_tot = {r['total_magnetization']}  "
              f"m(Mn) = {r['site_moments_Mn']}" + (f"  max|F_res| = {r['max_residual_force_eV_A']:.4f}"
                                                  if "max_residual_force_eV_A" in r else ""))
    for k, v in res.get("energies", {}).items():
        print(f"  {k:30s} {v:+.2f}")
    print("\n  FM (THz)  AFM~ (THz)   delta   overlap   production   FM/prod-1")
    for m in c["modes"]:
        print(f"  {m['fm_THz']:8.4f}  {m['afm_THz']:9.4f}  {100 * m['delta']:+6.2f} %  {m['max_overlap']:6.3f}"
              f"   {m['production_THz']:9.4f}   {100 * m['fm_vs_production']:+6.2f} %")
    print(f"\n  acoustic FM {c['acoustic_FM_THz']}  AFM {c['acoustic_AFM_THz']}")
    print(f"  noise (FM doublet splitting, max): {100 * c['noise_max_rel_split']:.2f} %")
    print(f"  AFM vs FM: rms {100 * c['rms_delta']:.2f} %, max {100 * c['max_abs_delta']:.2f} %; below 3 THz: "
          f"rms {100 * (c['rms_delta_below_3THz'] or 0):.2f} %, max {100 * (c['max_abs_delta_below_3THz'] or 0):.2f} %")
    print(f"  FM vs production, max above 1 THz: {100 * c['max_abs_fm_vs_production_above_1THz']:.2f} %")
    print(f"\n  VERDICT: {res['verdict']}")
    for n in res["notes"]:
        print("  NOTE:", n)


def selftest(repo):
    """Synthetic forces from the production FC2: FM must reproduce production exactly; an AFM-like
    perturbation (a Mn2-Se spring, breaking the Mn1/Mn2 equivalence) must give the shifts of a
    direct diagonalisation."""
    prim = read_vasp(str(repo / "Mn2In2Se5/POSCAR_relaxed"))
    nat = len(prim.symbols)
    ph_p, (w_p, e_p) = production_modes(repo)
    dm = ph_p.dynamical_matrix
    dm.run([0, 0, 0])
    m = np.repeat(np.sqrt(prim.masses), 3)
    phi = np.real(dm.dynamical_matrix) * np.outer(m, m)        # Gamma force constants, eV/A^2
    ph = phonopy.Phonopy(prim, np.eye(3, dtype=int), primitive_matrix="P", is_symmetry=False)
    ph.generate_displacements(distance=0.01, is_plusminus=True)
    disps = [{"number": int(d["number"]), "displacement": list(d["displacement"])} for d in ph.dataset["first_atoms"]]

    def synth(p):
        return [-(p[:, 3 * d["number"]:3 * d["number"] + 3] @ np.array(d["displacement"])).reshape(nat, 3)
                for d in disps]

    w_fm, e_fm = gamma_modes(prim, disps, synth(phi))
    dev = np.max(np.abs(w_fm - w_p))
    print(f"selftest FM: max |w - w_production| = {dev:.1e} THz")
    assert dev < 1e-6
    # AFM-like perturbation: stiffen the Mn2 - nearest Se bond by k = 0.3 eV/A^2 (translation invariant)
    pos = prim.scaled_positions @ prim.cell
    d = [np.linalg.norm(pos[j] - pos[1] - np.round((pos[j] - pos[1]) @ np.linalg.inv(prim.cell)) @ prim.cell)
         for j in range(nat)]
    j = int(np.argmin([d[k] if prim.symbols[k] == "Se" else 1e9 for k in range(nat)]))
    p2 = phi.copy()
    k = 0.3 * np.eye(3)
    for a, b, s in ((1, 1, 1), (j, j, 1), (1, j, -1), (j, 1, -1)):
        p2[3 * a:3 * a + 3, 3 * b:3 * b + 3] += s * k
    w_afm, e_afm = gamma_modes(prim, disps, synth(p2))
    ev, vec = np.linalg.eigh(p2 / np.outer(m, m))
    w_dir = np.sign(ev) * np.sqrt(np.abs(ev)) * 15.633302      # VaspToTHz
    dev2 = np.max(np.abs(np.sort(w_afm) - np.sort(w_dir)))
    print(f"selftest AFM-like: max |w - direct diagonalisation| = {dev2:.1e} THz")
    assert dev2 < 1e-4
    c = compare(w_fm, e_fm, w_afm, e_afm, w_p, e_p)
    print(f"selftest AFM-like perturbation: rms delta {100 * c['rms_delta']:.2f} %, max {100 * c['max_abs_delta']:.2f} %, "
          f"noise {100 * c['noise_max_rel_split']:.1e} %")
    assert c["noise_max_rel_split"] < 1e-6 and c["max_abs_delta"] > 0.005
    print("selftest passed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--work")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    repo = Path(a.repo)
    if a.selftest:
        selftest(repo)
        return
    work = Path(a.work)
    res = analyse(repo, work)
    report(res)
    (work / "results.json").write_text(json.dumps(res, indent=1, default=float))
    print(f"\nwrote {work / 'results.json'}")


if __name__ == "__main__":
    main()
