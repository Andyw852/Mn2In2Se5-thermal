#!/usr/bin/env python3
"""Fit and compare the spin-configuration FC2 sets written by prepare.py.

  python analyse.py --repo <clone> --work <workdir> --refs-only   # after the ref runs (step 2)
  python analyse.py --repo <clone> --work <workdir>     # after all VASP runs
  python analyse.py --repo <clone> --selftest           # before submitting: synthetic data, no VASP

Reads OUTCAR only (energy(sigma->0), TOTAL-FORCE, magnetization (x), iteration count).

Per set: forces of each displaced frame minus the forces of the set's own undisplaced ref
(this removes the static force of the spin arrangement), then one symmetry-adapted
least-squares FC2 (phonopy + symfc, full R-3m symmetry of the cell, acoustic sum rule):

  FM_ctl  the 9 FM frames                 method control, compared with production fc2
  AFM     the 8 AFM frames                ordered bound
  PM      the 16 frames of PM01..PM08     paramagnetic FC2 = spin-space average (pooled fit)
  PM-k    leave-one-configuration-out     jackknife error of PM

Then the checks and comparisons, with gates fixed before the runs (see README):
  G1  FM_ctl vs production fc2 (Mn2In2Se5/fc/FORCE_CONSTANTS_2ND), 6x6x2 q mesh, modes > 0.5 THz:
      rms relative frequency difference <= 1 %
  G2  AFM vs FM_ctl at Gamma (overlap-weighted shifts) reproduces the 9-atom result of
      dft/spin_phonon_gamma within 1.5 percentage points for every optical mode below 3 THz
      (the fit alone gives up to ~0.9 pp on synthetic data with 5e-5 eV/A force noise)
  G3  jackknife: median relative error of the PM frequencies below 3 THz <= 1 %
Reported: PM and AFM vs FM_ctl frequency shifts on the mesh, sound velocities, Gamma shifts,
energies vs the Heisenberg prediction, static forces of each spin arrangement, fit residuals
(the excess PM residual measures the force-constant spread between spin configurations),
the FM_ctl even-order (anharmonic) check from frame 9 = -(frame 1).

Writes <workdir>/results/{results.json, fc2_FM_ctl.hdf5, fc2_AFM.hdf5, fc2_PM.hdf5}.
"""
import argparse
import json
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import phonopy
from phonopy.file_IO import parse_FORCE_CONSTANTS, write_force_constants_to_hdf5

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prepare as P  # noqa: E402

NELM = 300
NAT = 144
MESH = (6, 6, 2)
G1_RMS, G2_PP, G3_MED = 0.01, 0.015, 0.01


# ----------------------------------------------------------------------------- OUTCAR
def read_outcar(path, natom=NAT):
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
        rows = re.findall(r"^\s+\d+((?:\s+-?\d+\.\d+){4,5})\s*$", body, re.M)
        out["site_moments"] = [float(r.split()[-1]) for r in rows[:natom]]
    fb = txt.split("TOTAL-FORCE (eV/Angst)")
    if len(fb) > 1:
        lines = fb[-1].splitlines()[2:2 + natom]
        out["forces"] = np.array([[float(x) for x in ln.split()[3:6]] for ln in lines])
    return out


def check_run(r, spins, label, problems):
    if not r["finished"]:
        problems.append(f"{label}: not finished")
    if r["n_elec_steps"] is not None and r["n_elec_steps"] >= NELM:
        problems.append(f"{label}: hit NELM ({r['n_elec_steps']})")
    sm = r.get("site_moments")
    if sm is None or len(sm) < 32 or r.get("forces") is None or len(r["forces"]) != NAT:
        problems.append(f"{label}: no site moments / forces in OUTCAR")
        return
    bad = [i + 1 for i, s in enumerate(spins) if abs(sm[i]) < 4.0 or np.sign(sm[i]) != s]
    if bad:
        problems.append(f"{label}: Mn moments off the intended configuration at ions {bad[:8]}"
                        + (" ..." if len(bad) > 8 else ""))
    m_exp = 4.992 * int(np.sum(spins))        # production FM: 9.984 muB per 2 Mn
    if r["total_magnetization"] is not None and abs(r["total_magnetization"] - m_exp) > 0.5:
        problems.append(f"{label}: total magnetization {r['total_magnetization']:.3f}, expected {m_exp:.2f}")


# ----------------------------------------------------------------------------- force constants
def template(repo):
    return phonopy.load(str(Path(repo) / P.DISP_YAML), produce_fc=False, log_level=0)


def new_ph(tpl, fc=None):
    # primitive = the 9-atom unit cell (a, b, a3prime) as in kappa_extrinsic.py ("P"), so the q mesh is in-plane x out-of-plane
    ph = phonopy.Phonopy(tpl.unitcell, tpl.supercell_matrix, primitive_matrix="P", log_level=0)
    if fc is not None:
        ph.force_constants = fc
    return ph


def fit_fc2(tpl, U, F):
    ph = new_ph(tpl)
    ph.dataset = {"displacements": np.asarray(U), "forces": np.asarray(F)}
    ph.produce_force_constants(fc_calculator="symfc", show_drift=False)
    return ph.force_constants


def predict(fc, U):
    return -np.einsum("ijab,njb->nia", fc, np.asarray(U))


def rel_rms(a, b):
    return float(np.sqrt(np.mean((a - b) ** 2)) / np.sqrt(np.mean(b ** 2)))


def mesh_q():
    return np.array([[i / MESH[0], j / MESH[1], k / MESH[2]]
                     for i in range(MESH[0]) for j in range(MESH[1]) for k in range(MESH[2])])


def freqs(tpl, fc, q):
    ph = new_ph(tpl, fc)
    ph.run_qpoints(q)
    return np.array(ph.qpoints.frequencies)


def compare_freqs(w, w0, wmin=0.5):
    """Relative shifts of sorted frequencies (per q) with respect to w0, for w0 > wmin THz."""
    a, b = np.sort(w, axis=1), np.sort(w0, axis=1)
    m = b > wmin
    d = a[m] / b[m] - 1
    lo = m & (b < 3.0)
    dl = a[lo] / b[lo] - 1
    return {"rms": float(np.sqrt(np.mean(d ** 2))), "mean": float(np.mean(d)), "max_abs": float(np.max(np.abs(d))),
            "rms_below_3THz": float(np.sqrt(np.mean(dl ** 2))), "mean_below_3THz": float(np.mean(dl)),
            "max_abs_below_3THz": float(np.max(np.abs(dl)))}


def sound_velocities(tpl, fc):
    """Acoustic velocities (m/s) from omega/|q| at |q| = 0.004 1/A along x (in plane) and z."""
    ph = new_ph(tpl, fc)
    lat = np.array(ph.primitive.cell)
    out = {}
    for name, qc in (("x", [0.004, 0, 0]), ("z", [0, 0, 0.004])):
        qf = np.array(qc) @ lat.T                    # fractional q (1/A without 2 pi)
        ph.run_qpoints([qf])
        w = np.sort(ph.qpoints.frequencies[0])[:3]
        out[name] = [float(x / 0.004 * 100) for x in w]   # THz*A -> 100 m/s
    return out


def gamma_shifts(tpl, fc_ref, fc):
    """Overlap-weighted Gamma shifts of fc relative to fc_ref (optical modes of fc_ref)."""
    def modes(f):
        ph = new_ph(tpl, f)
        ph.run_qpoints([[0, 0, 0]], with_eigenvectors=True)
        return np.array(ph.qpoints.frequencies[0]), np.array(ph.qpoints.eigenvectors[0])
    w0, e0 = modes(fc_ref)
    w1, e1 = modes(fc)
    ov = np.abs(e0.conj().T @ e1) ** 2
    out = []
    for i in range(3, len(w0)):
        wt = float(ov[i] @ w1)
        out.append({"ref_THz": float(w0[i]), "THz": wt, "delta": wt / w0[i] - 1})
    return out


def jackknife(tpl, data, names, q):
    ws, fcs = [], []
    for leave in names:
        U = [u for n in names if n != leave for u in data[n]["U"]]
        F = [f for n in names if n != leave for f in data[n]["dF"]]
        fcs.append(fit_fc2(tpl, U, F))
        ws.append(np.sort(freqs(tpl, fcs[-1], q), axis=1))
    ws = np.array(ws)
    n = len(names)
    mean = ws.mean(axis=0)
    err = np.sqrt((n - 1) / n * np.sum((ws - mean) ** 2, axis=0))
    return mean, err, fcs


# ----------------------------------------------------------------------------- analysis
def check_refs(work):
    """Step-2 check: the ref runs only (spin state, energies vs Heisenberg, static forces)."""
    work = Path(work)
    man = json.loads((work / "manifest.json").read_text())
    problems, rows = [], {}
    for name, s in man["sets"].items():
        spins = np.array(json.loads((work / name / "displacements.json").read_text())["spins"])
        p = work / name / "ref" / "OUTCAR"
        if not p.exists():
            problems.append(f"{name}/ref: no OUTCAR")
            continue
        r = read_outcar(p)
        check_run(r, spins, f"{name}/ref", problems)
        rows[name] = r
    if "FM_ctl" in rows and rows["FM_ctl"]["energy_eV"] is not None:
        efm = rows["FM_ctl"]["energy_eV"]
        print("set      E-E_FM (meV/Mn)  Heisenberg   M_tot     min|m_Mn|  static F rms / max (eV/A)")
        for name, r in rows.items():
            f = r.get("forces")
            e = (r["energy_eV"] - efm) / 32 * 1000 if r["energy_eV"] is not None else float("nan")
            sm = r.get("site_moments") or [float("nan")] * 32
            print(f"  {name:7s} {e:+9.2f}      {man['sets'][name]['E_J_minus_FM_meV_per_Mn']:+9.2f}  "
                  f"{r['total_magnetization'] if r['total_magnetization'] is not None else float('nan'):+9.3f}  "
                  f"{np.min(np.abs(sm[:32])):8.3f}   "
                  + (f"{np.sqrt(np.mean(f ** 2)):.4f} / {np.max(np.abs(f)):.4f}" if f is not None else "-"))
    if problems:
        print("PROBLEMS")
        for p in problems:
            print("  " + p)
    else:
        print("ref checks passed")


def load_sets(work):
    man = json.loads((work / "manifest.json").read_text())
    data, problems, runs = {}, [], {}
    for name in man["sets"]:
        dj = json.loads((work / name / "displacements.json").read_text())
        spins = np.array(dj["spins"])
        ref = read_outcar(work / name / "ref" / "OUTCAR")
        check_run(ref, spins, f"{name}/ref", problems)
        U, dF, Fraw = [], [], []
        for d, u in zip(dj["dirs"], dj["displacements_A"]):
            r = read_outcar(work / name / d / "OUTCAR")
            check_run(r, spins, f"{name}/{d}", problems)
            if r.get("forces") is not None and ref.get("forces") is not None:
                U.append(np.array(u))
                dF.append(r["forces"] - ref["forces"])
                Fraw.append(r["forces"])
        data[name] = {"spins": spins, "U": U, "dF": dF, "ref": ref}
        f0 = ref.get("forces")
        runs[name] = {"energy_eV": ref["energy_eV"], "total_magnetization": ref["total_magnetization"],
                      "mn_moments_min_abs": float(np.min(np.abs(ref["site_moments"][:32]))) if ref.get("site_moments") else None,
                      "static_force_rms_eV_A": float(np.sqrt(np.mean(f0 ** 2))) if f0 is not None else None,
                      "static_force_max_eV_A": float(np.max(np.abs(f0))) if f0 is not None else None,
                      "n_frames": len(U)}
    return man, data, problems, runs


def analyse(repo, work, gamma9=None):
    repo, work = Path(repo), Path(work)
    man, data, problems, runs = load_sets(work)
    if problems:
        return {"problems": problems, "runs": runs}
    tpl = template(repo)
    pm_names = sorted(k for k in data if k.startswith("PM"))
    res = {"problems": [], "runs": runs, "notes": []}

    # energies
    efm = runs["FM_ctl"]["energy_eV"]
    en = {}
    for k in data:
        en[k] = {"dft_meV_per_Mn": (runs[k]["energy_eV"] - efm) / 32 * 1000,
                 "heisenberg_meV_per_Mn": man["sets"][k]["E_J_minus_FM_meV_per_Mn"]}
    pm_e = [en[k]["dft_meV_per_Mn"] for k in pm_names]
    res["energies"] = {"sets": en, "PM_mean_meV_per_Mn": float(np.mean(pm_e)), "PM_spread_meV_per_Mn": float(np.std(pm_e)),
                       "AFM_9atom_meV_per_Mn": -10.44}

    # fits
    fits, resid = {}, {}
    groups = {"FM_ctl": ["FM_ctl"], "AFM": ["AFM"], "PM": pm_names}
    for g, names in groups.items():
        U = [u for n in names for u in data[n]["U"]]
        F = [f for n in names for f in data[n]["dF"]]
        fits[g] = fit_fc2(tpl, U, F)
        resid[g] = rel_rms(predict(fits[g], U), np.array(F))
    res["fit_relative_residual"] = resid
    excess = np.sqrt(max(resid["PM"] ** 2 - resid["FM_ctl"] ** 2, 0.0))
    res["PM_excess_residual"] = float(excess)

    # even-order check of FM_ctl (frame 5 = -frame 1)
    fp, fm = data["FM_ctl"]["dF"][0], data["FM_ctl"]["dF"][-1]
    res["FM_ctl_even_over_odd"] = float(np.sqrt(np.mean(((fp + fm) / 2) ** 2) / np.mean(((fp - fm) / 2) ** 2)))

    # frequencies
    q = mesh_q()
    prod = parse_FORCE_CONSTANTS(str(repo / "Mn2In2Se5/fc/FORCE_CONSTANTS_2ND"))
    w = {g: freqs(tpl, fc, q) for g, fc in fits.items()}
    w["production"] = freqs(tpl, prod, q)
    cmp_ = {"FM_ctl_vs_production": compare_freqs(w["FM_ctl"], w["production"]),
            "AFM_vs_FM_ctl": compare_freqs(w["AFM"], w["FM_ctl"]),
            "PM_vs_FM_ctl": compare_freqs(w["PM"], w["FM_ctl"]),
            "PM_vs_production": compare_freqs(w["PM"], w["production"])}
    res["mesh_comparison"] = cmp_
    res["imaginary_modes_THz"] = {g: float(np.min(v)) for g, v in w.items()}
    res["sound_velocity_m_s"] = {g: sound_velocities(tpl, fc) for g, fc in list(fits.items()) + [("production", prod)]}

    # jackknife
    jm, je, jk_fcs = jackknife(tpl, data, pm_names, q)
    b = np.sort(w["FM_ctl"], axis=1)
    lo = (b > 0.5) & (b < 3.0)
    wpm = np.sort(w["PM"], axis=1)
    rel_err = je[lo] / wpm[lo]
    sig = np.abs(wpm[lo] - b[lo]) / np.maximum(je[lo], 1e-9)
    res["jackknife"] = {"median_rel_err_below_3THz": float(np.median(rel_err)),
                        "max_rel_err_below_3THz": float(np.max(rel_err)),
                        "fraction_PM_FM_shift_above_3sigma_below_3THz": float(np.mean(sig > 3))}

    # Gamma
    gs = {"AFM": gamma_shifts(tpl, fits["FM_ctl"], fits["AFM"]), "PM": gamma_shifts(tpl, fits["FM_ctl"], fits["PM"])}
    res["gamma_shifts"] = gs

    # gates
    g1 = cmp_["FM_ctl_vs_production"]["rms"] <= G1_RMS
    if gamma9 is None:
        gpath = repo / "dft/spin_phonon_gamma/results/results.json"
        gamma9 = [m["delta"] for m in json.loads(gpath.read_text())["comparison"]["modes"]]
    d4 = [m["delta"] for m in gs["AFM"]]
    lowmask = [m["ref_THz"] < 3.0 for m in gs["AFM"]]
    g2_dev = max(abs(a - b) for a, b, l in zip(d4, gamma9, lowmask) if l)
    g2 = g2_dev <= G2_PP
    g3 = res["jackknife"]["median_rel_err_below_3THz"] <= G3_MED
    res["gates"] = {"G1_FM_ctl_vs_production_rms": [cmp_["FM_ctl_vs_production"]["rms"], g1],
                    "G2_AFM_gamma_vs_9atom_max_dev": [float(g2_dev), g2],
                    "G3_jackknife_median_rel_err": [res["jackknife"]["median_rel_err_below_3THz"], g3]}
    res["gates_passed"] = bool(g1 and g2 and g3)

    # full 144 x 144 FC2 (independent of the primitive-cell choice of the reader)
    out = work / "results"
    out.mkdir(exist_ok=True)
    for g, fc in fits.items():
        write_force_constants_to_hdf5(fc, filename=str(out / f"fc2_{g}.hdf5"))
    jk = work / "jackknife"
    jk.mkdir(exist_ok=True)
    for n, fc in zip(pm_names, jk_fcs):
        write_force_constants_to_hdf5(fc, filename=str(jk / f"fc2_PM_without_{n}.hdf5"))
    return res


def report(res):
    if res.get("problems"):
        print("PROBLEMS")
        for p in res["problems"]:
            print("  " + p)
        return
    print("sets (ref runs):")
    for k, r in res["runs"].items():
        e = res["energies"]["sets"][k]
        print(f"  {k:7s} E-E_FM = {e['dft_meV_per_Mn']:+7.2f} meV/Mn (Heisenberg {e['heisenberg_meV_per_Mn']:+7.2f})  "
              f"M = {r['total_magnetization']:+7.3f}  min|m_Mn| = {r['mn_moments_min_abs']:.3f}  "
              f"static F rms/max = {r['static_force_rms_eV_A']:.4f}/{r['static_force_max_eV_A']:.4f} eV/A  frames {r['n_frames']}")
    e = res["energies"]
    print(f"  PM mean {e['PM_mean_meV_per_Mn']:+.2f} meV/Mn, spread {e['PM_spread_meV_per_Mn']:.2f}; AFM 9-atom {e['AFM_9atom_meV_per_Mn']:+.2f}")
    print("fit relative residual: " + ", ".join(f"{k} {100 * v:.2f} %" for k, v in res["fit_relative_residual"].items())
          + f";  PM excess {100 * res['PM_excess_residual']:.2f} %")
    print(f"FM_ctl even/odd force ratio (frame 9 = -frame 1): {100 * res['FM_ctl_even_over_odd']:.2f} %")
    print("mesh 6x6x2, modes > 0.5 THz (relative frequency shifts):")
    for k, c in res["mesh_comparison"].items():
        print(f"  {k:22s} mean {100 * c['mean']:+6.2f} %  rms {100 * c['rms']:5.2f} %  max {100 * c['max_abs']:5.2f} %   "
              f"< 3 THz: mean {100 * c['mean_below_3THz']:+6.2f} %  rms {100 * c['rms_below_3THz']:5.2f} %")
    print("lowest frequency on the mesh (THz): " + ", ".join(f"{k} {v:.3f}" for k, v in res["imaginary_modes_THz"].items()))
    print("sound velocities (m/s), x: TA TA LA / z: TA TA LA")
    for k, v in res["sound_velocity_m_s"].items():
        print(f"  {k:10s} x " + " ".join(f"{x:6.0f}" for x in v["x"]) + "   z " + " ".join(f"{x:6.0f}" for x in v["z"]))
    j = res["jackknife"]
    print(f"jackknife PM (< 3 THz): median rel. error {100 * j['median_rel_err_below_3THz']:.2f} %, max {100 * j['max_rel_err_below_3THz']:.2f} %, "
          f"PM-FM shift > 3 sigma for {100 * j['fraction_PM_FM_shift_above_3sigma_below_3THz']:.0f} % of modes")
    print("Gamma, overlap-weighted shifts vs FM_ctl (optical modes):")
    print("   FM_ctl     AFM      PM")
    for a, p in zip(res["gamma_shifts"]["AFM"], res["gamma_shifts"]["PM"]):
        print(f"  {a['ref_THz']:7.4f}  {100 * a['delta']:+6.2f} %  {100 * p['delta']:+6.2f} %")
    print("gates:")
    for k, (v, ok) in res["gates"].items():
        print(f"  {k:36s} {v:.4f}  {'pass' if ok else 'FAIL'}")
    print("ALL GATES PASSED" if res["gates_passed"] else "GATE FAILED: stop and report")


# ----------------------------------------------------------------------------- self-test
def spring_model(sc, layer, spins, lam):
    """Spin-dependent Mn-Mn springs: Phi += sum_bonds lam_shell s_i s_j (radial spring), translation invariant."""
    shells = P.mn_shells(sc, layer)
    d = np.zeros((NAT, NAT, 3, 3))
    for k, bonds in shells.items():
        if k not in lam:
            continue
        for i, j, v in bonds:
            n = v / np.linalg.norm(v)
            K = lam[k] * spins[i] * spins[j] * np.outer(n, n)
            d[i, j] -= K
            d[j, i] -= K
            d[i, i] += K
            d[j, j] += K
    return d


def fake_outcar(path, energy, moments, forces, mtot):
    lines = [" running on fake", f"   free  energy   TOTEN  =  {energy:.8f} eV",
             f"  energy  without entropy=  {energy:.8f}  energy(sigma->0) =  {energy:.8f}",
             "--- Iteration      1(  25)  ---",
             f" number of electron     1200.0000000 magnetization      {mtot:.7f}",
             " magnetization (x)", "", "# of ion       s       p       d       tot", "-" * 42]
    for i, m in enumerate(moments, 1):
        lines.append(f"   {i:3d}        0.010   0.020 {m:8.3f} {m + 0.03:8.3f}")
    lines += ["-" * 50, "tot          0.000   0.000   0.000   0.000", "",
              " POSITION                                       TOTAL-FORCE (eV/Angst)",
              " " + "-" * 83]
    for f in forces:
        lines.append(f"      0.00000      0.00000      0.00000     {f[0]:12.6f} {f[1]:12.6f} {f[2]:12.6f}")
    lines += [" " + "-" * 83, "", " General timing and accounting informations for this job:"]
    Path(path).write_text("\n".join(lines) + "\n")


def selftest(repo):
    """Synthetic data: production FC2 plus spin-dependent Mn-Mn springs, static spin forces and noise,
    written as fake OUTCARs into a prepare.py work directory, then the full analysis.
    FM_ctl and AFM must give the exact FC2 of their configuration; PM must give the spin-space
    average (springs weighted by the mean correlation, here 0); the energies must equal the Heisenberg values."""
    repo = Path(repo)
    tpl = template(repo)
    prod = parse_FORCE_CONSTANTS(str(repo / "Mn2In2Se5/fc/FORCE_CONSTANTS_2ND"))
    ph, sc, layer, shells, perms, sets = P.all_sets(repo)
    # spring sizes chosen so that the synthetic AFM-FM Gamma shifts (max ~11 %, rms < 3 THz ~3 %)
    # resemble the measured 9-atom ones (8 %, 4 %)
    lam = {"inter_3.89": 0.08, "intra_4.03": 0.03, "inter_5.60": 0.05}       # eV/A^2
    rng = np.random.default_rng(7)
    tmp = Path(tempfile.mkdtemp(prefix="spinfc2_selftest_"))
    work = tmp / "work"
    import subprocess
    subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "prepare.py"), "--repo", str(repo),
                    "--out", str(work)], check=True, stdout=subprocess.DEVNULL)
    man = json.loads((work / "manifest.json").read_text())
    J = P.J_GAMMA_FIT
    exact = {}
    for name, s in sets.items():
        spins = s["spins"]
        fc = prod + spring_model(sc, layer, spins, lam)
        exact[name] = fc
        e = -1500.0 + 32 * (P.heisenberg_per_mn(spins, shells, J)) / 1000
        f0 = np.zeros((NAT, 3)) if name == "FM_ctl" else rng.normal(scale=0.02, size=(NAT, 3))
        f0 -= f0.mean(axis=0)
        mom = 4.55 * spins
        mtot = float(4.992 * spins.sum())
        fake_outcar(work / name / "ref" / "OUTCAR", e, mom, f0, mtot)
        for d, u in zip([f"disp-{k:02d}" for k in range(1, len(s["frames"]) + 1)], s["frames"]):
            f = f0 + predict(fc, [u])[0] + rng.normal(scale=5e-5, size=(NAT, 3))
            fake_outcar(work / name / d / "OUTCAR", e, mom, f, mtot)

    # the Gamma shifts of the 9-atom test are replaced by those of the exact synthetic FCs
    g9 = [m["delta"] for m in gamma_shifts(tpl, exact["FM_ctl"], exact["AFM"])]
    res = analyse(repo, work, gamma9=g9)
    report(res)
    print("(G1 is expected to FAIL here: the synthetic FM differs from production by construction)")
    assert not res["problems"], res["problems"]

    q = mesh_q()
    out = work / "results"
    from phono3py.file_IO import read_fc2_from_hdf5  # noqa: F401  (the reader kappa_extrinsic.py uses)
    fit = {g: read_fc2_from_hdf5(str(out / f"fc2_{g}.hdf5")) for g in ("FM_ctl", "AFM", "PM")}
    pm_expected = prod + spring_model(sc, layer, np.ones(32), {k: v * 0.0 for k, v in lam.items()})
    checks = [("FM_ctl", exact["FM_ctl"]), ("AFM", exact["AFM"]), ("PM", pm_expected)]
    for g, fc_exact in checks:
        ph_fit = new_ph(tpl)
        ph_fit.force_constants = fit[g]
        ph_fit.run_qpoints(q)
        c = compare_freqs(np.array(ph_fit.qpoints.frequencies), freqs(tpl, fc_exact, q))
        print(f"selftest {g:6s}: fitted (read back from hdf5) vs exact synthetic FC2: rms |dw/w| "
              f"{100 * c['rms']:.2f} %, below 3 THz {100 * c['rms_below_3THz']:.2f} %")
        assert c["rms_below_3THz"] < (0.004 if g != "PM" else 0.008), g
    sig = res["mesh_comparison"]["PM_vs_FM_ctl"]["rms_below_3THz"]
    print(f"selftest: synthetic PM vs FM shift (rms < 3 THz) {100 * sig:.2f} % (must be resolved)")
    assert sig > 0.007
    for k, e in res["energies"]["sets"].items():
        assert abs(e["dft_meV_per_Mn"] - e["heisenberg_meV_per_Mn"]) < 0.01, k
    assert res["gates"]["G2_AFM_gamma_vs_9atom_max_dev"][1] and res["gates"]["G3_jackknife_median_rel_err"][1]
    assert res["FM_ctl_even_over_odd"] < 0.01
    shutil.rmtree(tmp)
    print("selftest passed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--work")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--refs-only", action="store_true", help="check the ref runs only (step 2)")
    a = ap.parse_args()
    if a.selftest:
        selftest(a.repo)
        return
    if a.refs_only:
        check_refs(a.work)
        return
    work = Path(a.work)
    res = analyse(a.repo, work)
    report(res)
    (work / "results").mkdir(exist_ok=True)
    (work / "results" / "results.json").write_text(json.dumps(res, indent=1, default=float))
    print(f"\nwrote {work / 'results'}")


if __name__ == "__main__":
    main()
