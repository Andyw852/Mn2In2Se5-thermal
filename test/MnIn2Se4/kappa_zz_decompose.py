#!/usr/bin/env python3
"""Separate the Gamma-Z quadrature error from kappa_zz in phono3py RTA output,
apply Casimir boundaries a posteriori, and estimate the mesh-converged kappa_zz.

Why this exists
---------------
On a Gamma-centred N^3 mesh, the first q-point on the Gamma-Z line (the layer
normal) carries the interlayer acoustic modes with the longest lifetimes. That
single point is given the weight of a whole mesh cell, so on coarse meshes it
dominates kappa_zz (50 % at 11^3) and the total falls with every refinement.
Everything else in the zone is already converged:

    kappa_zz(total) = kappa_zz(first Gamma-Z point) + kappa_zz(rest)

The first-point term falls as N^-s with s ~ 1.4-2.2 between successive meshes
(about 1.7 overall from 14^3 to 30^3; it is *not* 1/N). It is measured here, not
assumed. kappa_zz(rest) is the converged quantity, but it is biased slightly
low, because it drops the true contribution c(N) of the cell around that point.
This script estimates that bias from the drift of the rest term (see
`estimate_intrinsic`).

**Every derived number below uses the rest term, finite-grain values
included.** The finite-grain *totals* still contain the quadrature error and can
even come out above the intrinsic value (at 30^3 the 10 um total is 0.723, while
the intrinsic value is 0.646), which is physically impossible. The totals are
printed only so they can be compared with the old tables.

Gamma-Z line detection
----------------------
q-points are converted to Cartesian (minimum image) using the primitive cell,
and a point is on the Gamma-Z line when its in-plane Cartesian component is
zero. This works for any cell setting. For the rhombohedral MnIn2Se4 primitive
cell the line is (xi, xi, xi); for a hexagonal cell such as Mn2In2Se5 it is
(0, 0, xi). Pass --poscar/--pmat for another material.

Finite-grain kappa (RTA only)
-----------------------------
    gamma_b = |v| * Angstrom * 1e6 / (4 pi * L[um])   (phono3py convention)
    mode_kappa_L = mode_kappa * (gamma_anh + gamma_iso) / (gamma_anh + gamma_iso + gamma_b)

Validated against a real phono3py run with boundary_mfp=1.0 (kappa_zz_scan/m11_b1um,
kappa_zz 0.8976 both ways), and reproducing direct 18^3 runs (0.8849 / 0.6877).

Fixed-sphere check
------------------
An independent test of the quadrature diagnosis that does not depend on how the
first point is identified: kappa_zz is summed over the modes with |q| > q_c, for
a fixed physical radius q_c. For q_c a few mesh spacings wide, that sum must be
mesh-independent. As q_c -> 0 it approaches the intrinsic value from below.

Usage
-----
    python kappa_zz_decompose.py kappa_zz_scan/m*/kappa-m*.hdf5 --mfp 1 10
    python kappa_zz_decompose.py kappa_zz_scan/c*_m14/kappa-m*.hdf5 --no-estimate
    # another material:
    python kappa_zz_decompose.py run*/kappa-*.hdf5 --poscar POSCAR \
        --pmat 1 0 0 0 1 0 0 0 1

Writes kappa_zz_gamma_line.json (override with --json) and tables on stdout.
"""
import argparse
import glob
import json
import os
import re

import h5py
import numpy as np
from phonopy.interface.vasp import read_vasp
from phonopy.physical_units import get_physical_units

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_POSCAR = os.path.join(HERE, "fit", "POSCAR")
DEFAULT_PMAT = [1, 1, 2, -1, -2, -1, 1, 1, 1]       # same as kappa_driver_442.py


def reciprocal_rows(poscar, pmat):
    """Reciprocal vectors (rows, 1/Angstrom, no 2 pi) of the primitive cell."""
    cell = read_vasp(poscar)
    lat_p = (np.asarray(cell.cell).T @ np.reshape(pmat, (3, 3))).T   # rows = primitive vectors
    return np.linalg.inv(lat_p).T


def min_image_cart(q_red, rec):
    """Shortest Cartesian image of each reduced q (search the 27 neighbouring translations)."""
    q = np.asarray(q_red, float)
    q = q - np.round(q)
    shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    cand = (q[:, None, :] + shifts[None]) @ rec          # (nq, 27, 3)
    best = np.argmin(np.linalg.norm(cand, axis=2), axis=1)
    return cand[np.arange(len(q)), best]


def load(path):
    with h5py.File(path) as h:
        return {k: h[k][()] for k in h}


def decompose(path, mfps, rec, q_spheres, temp_index=0):
    d = load(path)
    mesh = [int(x) for x in d["mesh"]]
    w = np.asarray(d["weight"])
    n_bz = int(w.sum())
    if "boundary_mfp" in d:
        return None          # boundary already inside mode_kappa: use only as a validation file
    if "mode_kappa_RTA" in d:
        raise SystemExit(f"{path}: LBTE output. The per-mode split and the Casimir "
                         "rebuild are only valid for RTA runs.")
    mk = np.asarray(d["mode_kappa"])[temp_index] / n_bz      # (nir, nb, 6), already weighted
    kzz = mk[:, :, 2]
    kip = 0.5 * (mk[:, :, 0] + mk[:, :, 1])

    qc = min_image_cart(d["qpoint"], rec)
    qn = np.linalg.norm(qc, axis=1)
    q_inplane = np.linalg.norm(qc[:, :2], axis=1)
    on_line = (q_inplane < 1e-6 * max(qn.max(), 1)) & (qn > 1e-9)
    line = np.where(on_line)[0]
    line = line[np.argsort(qn[line])]
    if len(line) == 0:
        raise SystemExit(f"{path}: no q-point on the Gamma-Z line -- check --poscar/--pmat")
    first = int(line[0])

    gamma = np.asarray(d["gamma"])[temp_index]
    g_int = gamma + (np.asarray(d["gamma_isotope"]) if "gamma_isotope" in d else 0.0)
    vnorm = np.linalg.norm(np.asarray(d["group_velocity"]), axis=2)
    ang = get_physical_units().Angstrom

    def split(scale):
        z = kzz * scale
        tot = float(z.sum())
        f1 = float(z[first].sum())
        return {"kappa_zz_total": round(tot, 6),
                "kappa_zz_first_gamma_z_point": round(f1, 6),
                "kappa_zz_rest": round(tot - f1, 6),
                "inplane_total": round(float((kip * scale).sum()), 6)}

    out = {
        "file": os.path.relpath(path),
        "mesh": mesh,
        "n_bz_grid": n_bz,
        "T_K": float(np.asarray(d["temperature"])[temp_index]),
        "kappa_zz_stored": float(np.asarray(d["kappa"])[temp_index, 2]),
        "intrinsic": split(1.0),
        "gamma_z_line": [
            {"q_reduced": [round(float(x), 6) for x in d["qpoint"][i]],
             "q_cart_inv_A": round(float(qn[i]), 6),
             "weight": int(w[i]),
             "kappa_zz_star": round(float(kzz[i].sum()), 6)}
            for i in line[:6]],
        "finite_grain": {},
        "kappa_zz_outside_sphere": {},
    }
    for L in mfps:
        g_b = vnorm * ang * 1e6 / (4 * np.pi * L)
        scale = np.divide(g_int, g_int + g_b, out=np.zeros_like(g_int), where=(g_int + g_b) > 0)
        out["finite_grain"][f"{L:g}_um"] = split(scale)
    for qs in q_spheres:
        outside = float(kzz[qn > qs].sum())
        out["kappa_zz_outside_sphere"][f"{qs:g}"] = {
            "outside": round(outside, 6),
            "inside_without_first_point": round(float(kzz.sum()) - outside - float(kzz[first].sum()), 6)}
    return out


def estimate_intrinsic(rows, key=("intrinsic",), nmin=14, tol=0.01):
    """Mesh-converged kappa_zz from the rest term, with the cell-bias correction.

    Model: rest(N) = K - c0 N^-s, where s is the measured decay exponent of the
    first-point term over the meshes with N >= nmin. Then the bias of the finest
    mesh is c(Nmax) = [rest(Nmax) - rest(Nref)] / [(Nmax/Nref)^s - 1], with Nref
    the coarsest mesh >= nmin. The quoted uncertainty is the larger of |c| and
    the full range of rest over N >= nmin.
    """
    def pick(r):
        x = r
        for k in key:
            x = x[k]
        return x
    use = [r for r in rows if r["mesh"][0] >= nmin and len(set(r["mesh"])) == 1]
    if len(use) < 2:
        return None
    N = np.array([r["mesh"][0] for r in use], float)
    rest = np.array([pick(r)["kappa_zz_rest"] for r in use])
    f1 = np.array([pick(r)["kappa_zz_first_gamma_z_point"] for r in use])
    s = -np.polyfit(np.log(N), np.log(f1), 1)[0]
    step = [float(-np.log(f1[i + 1] / f1[i]) / np.log(N[i + 1] / N[i])) for i in range(len(N) - 1)]
    c = (rest[-1] - rest[0]) / ((N[-1] / N[0]) ** s - 1)
    out = {"meshes": N.astype(int).tolist(), "rest": rest.round(6).tolist(),
           "first_point_exponent_fit": round(float(s), 3),
           "first_point_exponent_steps": [round(x, 3) for x in step]}
    # The bias model only describes a rest term that RISES towards K (it misses
    # a positive cell contribution). A falling rest term means other near-Gamma
    # points are still overweighted; extrapolating it with the first-point
    # exponent has no basis. With two meshes the drift cannot be told from noise.
    last_step = (rest[-1] - rest[-2]) / rest[-2]
    if len(N) < 3 or (c < 0 and abs(last_step) > tol):
        why = "only two meshes" if len(N) < 3 else f"rest term still falling ({last_step:+.1%} in the last step)"
        out.update({"status": f"NOT CONVERGED: {why}",
                    "kappa_zz_rest_finest_mesh": round(float(rest[-1]), 4),
                    "kappa_zz_estimate": None, "uncertainty": None})
        return out
    unc = max(abs(c), float(rest.max() - rest.min()))
    out.update({"status": "converged", "bias_finest_mesh": round(float(c), 6),
                "kappa_zz_estimate": round(float(rest[-1] + c), 4), "uncertainty": round(unc, 4)})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("files", nargs="+", help="kappa-*.hdf5 from RTA runs (globs allowed)")
    ap.add_argument("--json", default="kappa_zz_gamma_line.json")
    ap.add_argument("--mfp", nargs="*", type=float, default=[1.0, 10.0],
                    help="Casimir boundary lengths in micrometres (default 1 10)")
    ap.add_argument("--poscar", default=DEFAULT_POSCAR)
    ap.add_argument("--pmat", nargs=9, type=float, default=DEFAULT_PMAT,
                    help="primitive matrix, row-major, as passed to phono3py")
    ap.add_argument("--sphere", nargs="*", type=float, default=[0.02, 0.03],
                    help="radii q_c (1/Angstrom, no 2 pi) for the fixed-sphere check")
    ap.add_argument("--nmin", type=int, default=14,
                    help="smallest mesh used in the intrinsic estimate (default 14)")
    ap.add_argument("--estimate-cases", default=r"^m\d+$",
                    help="regex on the case (directory) name selecting the mesh series for the "
                         "estimate; default picks m11, m14, ... and leaves out the cutoff cases")
    ap.add_argument("--no-estimate", action="store_true",
                    help="skip the mesh extrapolation (e.g. for a cutoff scan at one mesh)")
    a = ap.parse_args()

    rec = reciprocal_rows(a.poscar, a.pmat)
    files = []
    for f in a.files:
        files.extend(sorted(glob.glob(f)) or [f])
    rows = [decompose(f, a.mfp, rec, a.sphere) for f in files]
    skipped = [f for f, r in zip(files, rows) if r is None]
    rows = sorted((r for r in rows if r is not None), key=lambda r: (r["mesh"][0], r["file"]))
    for f in skipped:
        print(f"skipped {f}: run already contains a boundary term (validation file only)")

    def name(r):
        return os.path.basename(os.path.dirname(r["file"])) or "x".join(map(str, r["mesh"]))

    grains = [("intrinsic",)] + [("finite_grain", k) for k in rows[0]["finite_grain"]]
    for g in grains:
        label = "intrinsic" if len(g) == 1 else g[1].replace("_um", " um grain")
        print(f"\n== {label} ==")
        print(f"{'case':>10} {'in-plane':>9} {'kzz total':>10} {'first G-Z':>10} {'REST':>8} {'first/tot':>9}")
        for r in rows:
            x = r[g[0]] if len(g) == 1 else r[g[0]][g[1]]
            print(f"{name(r):>10} {x['inplane_total']:9.4f} {x['kappa_zz_total']:10.4f} "
                  f"{x['kappa_zz_first_gamma_z_point']:10.4f} {x['kappa_zz_rest']:8.4f} "
                  f"{x['kappa_zz_first_gamma_z_point'] / x['kappa_zz_total']:9.1%}")

    print("\n== Gamma-Z line, kappa_zz per star (first three points) ==")
    for r in rows:
        print(f"{name(r):>10} " + "  ".join(f"|q|={p['q_cart_inv_A']:.4f}: {p['kappa_zz_star']:.4f}"
                                             for p in r["gamma_z_line"][:3]))

    print("\n== fixed-sphere check (q_c in 1/A): outside |q|>q_c  /  inside minus first G-Z point ==")
    print("   both columns must be mesh-independent once q_c spans several mesh spacings;")
    print("   their sum is the rest term")
    print(f"{'case':>10} " + " ".join(f"{'q_c=' + k:>19}" for k in rows[0]["kappa_zz_outside_sphere"]))
    for r in rows:
        print(f"{name(r):>10} " + " ".join(f"{v['outside']:9.4f} {v['inside_without_first_point']:9.4f}"
                                           for v in r["kappa_zz_outside_sphere"].values()))

    est = {}
    if not a.no_estimate:
        print(f"\n== mesh-converged kappa_zz (rest term + cell-bias correction, N >= {a.nmin}) ==")
        for g in grains:
            series = [r for r in rows if re.search(a.estimate_cases, name(r))] or rows
            # the estimate only uses cubic meshes, so a material may legitimately
            # have several cases at the same N as long as only one of them is
            # cubic (e.g. Mn2In2Se5: 16x16x4 ... 16x16x16 plus an isotropic 24^3)
            cubic = [r for r in series if len(set(r["mesh"])) == 1]
            if len({r["mesh"][0] for r in cubic}) != len(cubic):
                raise SystemExit("estimate: more than one cubic case per mesh -- narrow --estimate-cases")
            e = estimate_intrinsic(series, g, a.nmin)
            if e is None:
                print("   need at least two cubic meshes with N >= nmin")
                break
            label = "intrinsic" if len(g) == 1 else g[1]
            est[label] = e
            if e["kappa_zz_estimate"] is None:
                print(f"{label:>10}: {e['status']} -- rest term {e['kappa_zz_rest_finest_mesh']:.4f} "
                      f"at N={e['meshes'][-1]}, still moving, not quotable (rest series {e['rest']})")
                continue
            print(f"{label:>10}: {e['kappa_zz_estimate']:.4f} +- {e['uncertainty']:.4f}   "
                  f"(rest at N={e['meshes'][-1]}: {e['rest'][-1]:.4f}, bias {e['bias_finest_mesh']:+.4f}, "
                  f"first-point exponent {e['first_point_exponent_fit']:.2f}, steps {e['first_point_exponent_steps']})")

    json.dump({"cases": rows, "estimate": est}, open(a.json, "w", encoding="utf-8"), indent=1)
    print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
