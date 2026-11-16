#!/usr/bin/env python3
"""Decompose kappa_zz of a phono3py kappa-*.hdf5 into the contribution of the
first Gamma-Z q-point and everything else, and rebuild the finite-grain values.

Why this exists
---------------
The q-points on the Gamma-Z line (the layer normal) are the only ones for which
kappa_zz keeps changing when the q-mesh is refined. Their contribution scales as
1/N with the mesh, so the total kappa_zz drifts downwards while the rest of the
Brillouin zone is already converged. Separating the two makes that visible:

    kappa_zz(total) = kappa_zz(first Gamma-Z point) + kappa_zz(rest)

and it is kappa_zz(rest) that is the mesh-converged quantity.

The Gamma-Z line is the (xi, xi, xi) direction in the *primitive* reciprocal
basis of this cell: the layer normal of the selected primitive cell is
(t, t, t), not (0, 0, t) -- see sound_velocity.py for the same pitfall in the
sound-velocity calculation.

Finite-grain kappa
------------------
RTA kappa is proportional to 1/gamma mode by mode, so a Casimir boundary can be
applied a posteriori without re-running phono3py:

    gamma_tot = gamma_anh + gamma_iso + gamma_b,
    gamma_b   = |v| * Angstrom * 1e6 / (4 pi * boundary_mfp[um])
    (phono3py/conductivity/base.py::_get_boundary_scattering)

    mode_kappa_bound = mode_kappa * gamma_anh+iso / gamma_tot

This reconstruction is validated against a real phono3py boundary run:
kappa_zz_scan/m11_b1um (boundary_mfp=1.0) gives kappa_zz = 0.8976, and the
a-posteriori value from kappa_zz_scan/m11 is also 0.8976 (12^3 apart at the 4th
decimal). So the finite-grain numbers follow from a single intrinsic run and no
extra phono3py runs are needed.

Usage
-----
    python kappa_zz_decompose.py kappa_zz_scan/m18/kappa-m181818.hdf5
    python kappa_zz_decompose.py kappa_zz_scan/*/kappa-m*.hdf5 \
        --json kappa_zz_gamma_line.json --mfp 1 10

Output: `kappa_zz_gamma_line.json` (override with --json), plus a table on stdout.
"""
import argparse
import glob
import json
import os

import h5py
import numpy as np
from phonopy.physical_units import get_physical_units


def load(path):
    with h5py.File(path) as h:
        d = {k: h[k][()] for k in h}
    return d


def reduce_q(q):
    """Map q to [0, 1) -- phono3py reports irreducible points with negative components."""
    a = np.asarray(q, dtype=float) % 1.0
    a[np.abs(a - 1.0) < 1e-6] = 0.0
    return a


def diag_mask(qn, tol=1e-6):
    """True for q-points on the (xi, xi, xi) line of the primitive reciprocal basis."""
    return (np.abs(qn[:, 0] - qn[:, 1]) < tol) & (np.abs(qn[:, 1] - qn[:, 2]) < tol)


def decompose(path, mfps, temp_index=0):
    d = load(path)
    mesh = [int(x) for x in d["mesh"]]
    n_bz = int(np.sum(d["weight"]))
    q = reduce_q(d["qpoint"])
    w = np.asarray(d["weight"])
    gp = np.asarray(d["grid_point"])

    # plain RTA runs write "mode_kappa"/"kappa"; LBTE runs write the *_RTA pair
    mk_key = "mode_kappa_RTA" if "mode_kappa_RTA" in d else "mode_kappa"
    k_key = "kappa_RTA" if "kappa_RTA" in d else "kappa"
    mk = np.asarray(d[mk_key])[temp_index]        # (n_irgp, n_band, 6)
    kzz_total = float(np.sum(mk[:, :, 2]) / n_bz)
    kzz_stored = float(np.asarray(d[k_key])[temp_index, 2])
    kip_total = float(0.5 * (np.sum(mk[:, :, 0]) + np.sum(mk[:, :, 1])) / n_bz)

    mask = diag_mask(q)
    line = np.where(mask & (np.abs(q[:, 0]) > 1e-9))[0]
    order = np.argsort(q[line, 0])
    line = line[order]

    per_point = np.array([np.sum(mk[i, :, 2]) / n_bz / w[i] for i in line])
    star = np.array([np.sum(mk[i, :, 2]) / n_bz for i in line])

    first = int(line[0])
    out = {
        "file": os.path.relpath(path),
        "mesh": mesh,
        "n_bz_grid": n_bz,
        "n_ir_grid": int(len(gp)),
        "T_K": float(np.asarray(d["temperature"])[temp_index]),
        "inplane_total": kip_total,
        "kappa_zz_total": kzz_total,
        "kappa_zz_stored": kzz_stored,
        "kappa_dataset": k_key,
        "gamma_z_line": [
            {"q_reduced": [round(float(x), 6) for x in q[i]],
             "weight": int(w[i]),
             "star_contribution": round(float(s), 6),
             "per_point_contribution": round(float(p), 6)}
            for i, s, p in zip(line, star, per_point)
        ],
        "kappa_zz_first_gamma_z_point": {
            "q_reduced": [round(float(x), 6) for x in q[first]],
            "weight": int(w[first]),
            "star_contribution": round(float(star[0]), 6),
            "per_point_contribution": round(float(per_point[0]), 6),
        },
        "kappa_zz_rest": round(float(kzz_total - star[0]), 6),
        "kappa_zz_gamma_z_line_total": round(float(np.sum(star)), 6),
        "kappa_zz_finite_grain": {},
    }

    # --- rebuild the finite-grain values from the same output file ---
    gamma = np.asarray(d["gamma"])[temp_index]                   # (n_irgp, n_band)
    g_iso = np.asarray(d["gamma_isotope"]) if "gamma_isotope" in d else np.zeros_like(gamma)
    gv = np.asarray(d["group_velocity"])                         # (n_irgp, n_band, 3)
    ang = get_physical_units().Angstrom
    g_int = gamma + g_iso
    for mfp in mfps:
        g_b = np.linalg.norm(gv, axis=2) * ang * 1e6 / (4 * np.pi * mfp)
        scale = np.divide(g_int, g_int + g_b, out=np.zeros_like(g_int), where=(g_int + g_b) > 0)
        mk_b = mk[:, :, 2] * scale
        tot = float(np.sum(mk_b) / n_bz)
        first_term = float(np.sum(mk_b[first]) / n_bz)
        ip_b = float(0.5 * (np.sum(mk[:, :, 0] * scale) + np.sum(mk[:, :, 1] * scale)) / n_bz)
        out["kappa_zz_finite_grain"][f"{mfp:g}_um"] = {
            "kappa_zz_total": round(tot, 6),
            "kappa_zz_first_gamma_z_point": round(first_term, 6),
            "kappa_zz_rest": round(tot - first_term, 6),
            "ratio_to_intrinsic": round(tot / kzz_total, 6),
            "inplane_total": round(ip_b, 6),
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+", help="kappa-*.hdf5 (globs allowed)")
    ap.add_argument("--json", default="kappa_zz_gamma_line.json")
    ap.add_argument("--mfp", nargs="*", type=float, default=[1.0, 10.0],
                    help="Casimir boundary MFPs in micrometres (default 1 10)")
    a = ap.parse_args()

    files = []
    for f in a.files:
        files.extend(sorted(glob.glob(f)) or [f])

    rows = [decompose(f, a.mfp) for f in files]
    rows.sort(key=lambda r: r["mesh"][0])

    keys = sorted({k for r in rows for k in r["kappa_zz_finite_grain"]})
    hdr = " ".join(f"{k.replace('_um', '') + 'um':>8}" for k in keys)
    print(f"{'case':>10} {'inplane':>9} {'kzz_tot':>9} {'first(g-Z)':>11} {'rest':>8} "
          f"{'rest/tot':>9} {hdr}")
    for r in rows:
        m = os.path.basename(os.path.dirname(r["file"])) or "x".join(map(str, r["mesh"]))
        f1 = r["kappa_zz_first_gamma_z_point"]["star_contribution"]
        gf = r["kappa_zz_finite_grain"]
        vals = " ".join(f"{gf[k]['kappa_zz_total']:8.4f}" for k in keys)
        print(f"{m:>10} {r['inplane_total']:9.4f} {r['kappa_zz_total']:9.4f} {f1:11.4f} {r['kappa_zz_rest']:8.4f} "
              f"{r['kappa_zz_rest'] / r['kappa_zz_total']:9.3f} {vals}")

    json.dump(rows, open(a.json, "w", encoding="utf-8"), indent=1)
    print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
