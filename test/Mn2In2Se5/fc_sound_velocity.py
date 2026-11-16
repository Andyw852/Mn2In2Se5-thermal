#!/usr/bin/env python3
"""Sound velocities and the implied elastic constants from a phono3py FC2.

Why this exists
---------------
Mn2In2Se5's FC2 and FC3 both live in a 4x4x1 supercell, and the cell is only
16.24 A along the layer normal, so that supercell is ONE LAYER THICK. That is the
defect that made MnIn2Se4's old 4x4x1 FC2 4.2x too soft between the layers and
suppressed its kappa_zz by ~3.5x. This script measures the same quantity -- the
interlayer transverse sound velocity -- from the FC2 itself, so that the question
can be answered without new DFT.

Conventions (the two pitfalls found for MnIn2Se4, see
../MnIn2Se4/fc2_442_phonon/sound_velocity.py):

1. phonopy q-points are in the reciprocal basis of the *selected primitive* cell,
   not of the POSCAR cell. Here the primitive matrix is the identity, so the two
   coincide and q along the layer normal is (0, 0, t) in reduced coordinates.
   The script nevertheless derives the reduced direction from the Cartesian layer
   normal, so it stays correct for any primitive matrix.
2. The Gamma point carries a numerical frequency offset (f0 ~ 1e-3 THz) which
   makes omega/q diverge at small q. Fit f^2 = f0^2 + (v q / 2pi / 100)^2 and take
   the square root of the slope.

For a hexagonal/rhombohedral layered crystal, with rho the density:
    q along c : TA1 = TA2 = sqrt(C44/rho)      (both polarised in the plane)
                LA        = sqrt(C33/rho)
    q in plane: TA (polarised along c) = sqrt(C44/rho)   <- same C44, a check
                TA (polarised in plane) = sqrt(C66/rho)
                LA                      = sqrt(C11/rho)

Usage
-----
    python fc_sound_velocity.py --fc2 F --poscar P --supercell 4 4 1 \
        [--label NAME] [--ref-vta M_S] [--ref-c44 GPa]
"""
import argparse

import numpy as np
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
from phono3py.file_IO import read_fc2_from_hdf5

# 1 THz.Angstrom = 100 m/s
THZ_A_TO_MS = 100.0


def fit_velocity(ph, direction, Bp, window, npts=25):
    """f^2 = f0^2 + (v q / 2pi / 100)^2 for the three acoustic branches.

    The branches are identified by ordering, not by eigenvectors: the eigenvector
    layout of phonopy's qpoint dictionaries differs between versions, and the
    acoustic branches are the three lowest frequencies at small q anyway. For
    q || c the two TAs are degenerate, which identifies them unambiguously.
    """
    qs, fs = [], []
    for t in np.linspace(0.001, window, npts):
        q = np.asarray(direction, float) * t
        ph.run_qpoints([q], with_eigenvectors=False)
        f = np.sort(ph.get_qpoints_dict()["frequencies"][0])[:3]
        qs.append(float(np.linalg.norm(q @ Bp)))
        fs.append(f)
    qs, fs = np.array(qs), np.array(fs)

    v, f0 = [], []
    for b in range(3):
        A = np.vstack([np.ones_like(qs), qs ** 2]).T
        (intercept, slope), *_ = np.linalg.lstsq(A, fs[:, b] ** 2, rcond=None)
        v.append(np.sqrt(max(slope, 0.0)) * 2 * np.pi * THZ_A_TO_MS)
        f0.append(np.sqrt(max(intercept, 0.0)))
    return np.array(v), np.array(f0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fc2", required=True)
    ap.add_argument("--poscar", required=True)
    ap.add_argument("--supercell", nargs=3, type=int, default=[4, 4, 1])
    ap.add_argument("--label", default="")
    ap.add_argument("--ref-vta", type=float, default=None,
                    help="reference interlayer TA velocity (m/s), e.g. from DFT C44")
    ap.add_argument("--ref-c44", type=float, default=None,
                    help="reference C44 in GPa")
    ap.add_argument("--window", type=float, default=0.02)
    ap.add_argument("--json", default=None,
                    help="append this case's result to a JSON list at this path")
    a = ap.parse_args()

    uc = read_vasp(a.poscar)
    ph = Phonopy(uc, supercell_matrix=np.diag(a.supercell),
                 primitive_matrix=np.eye(3), log_level=0)
    ph.force_constants = read_fc2_from_hdf5(a.fc2)
    pc = ph.primitive
    Bp = 2 * np.pi * np.linalg.inv(np.asarray(pc.cell)).T     # rows: b1,b2,b3
    rho = pc.masses.sum() / pc.volume * 1.66053906660e-27 / 1e-30

    a1, a2 = np.asarray(uc.cell)[0], np.asarray(uc.cell)[1]
    n = np.cross(a1, a2)
    n /= np.linalg.norm(n)                                   # layer normal
    # reduced coordinates of the layer normal in the primitive reciprocal basis
    u_c = np.linalg.solve(Bp.T, n)
    u_c /= np.linalg.norm(u_c)
    # an in-plane reduced direction: project b1 onto the plane
    b1 = Bp[0]
    ip = b1 - (b1 @ n) * n
    u_ip = np.linalg.solve(Bp.T, ip / np.linalg.norm(ip))
    u_ip /= np.linalg.norm(u_ip)

    tag = a.label or a.fc2
    print(f"== {tag}")
    print(f"   unitcell {len(uc)} atoms, supercell {np.prod(a.supercell)}x, "
          f"rho = {rho:.1f} kg/m^3")
    print(f"   layer normal (Cartesian) = {np.round(n, 4)}, "
          f"reduced direction = {np.round(u_c, 4)}")

    # --- q along the layer normal: interlayer TA (degenerate) and LA
    v_c, f0_c = fit_velocity(ph, u_c, Bp, a.window)
    order = np.argsort(v_c)
    v_ta = float(np.mean(v_c[order[:2]]))
    v_la = float(v_c[order[2]])
    print(f"   q || c       TA = {v_c[order[0]]:.1f} / {v_c[order[1]]:.1f} m/s "
          f"(degenerate pair), LA = {v_la:.1f} m/s   [f0 = {f0_c[order[0]]:.6f} THz]")
    c44 = rho * v_ta ** 2
    c33 = rho * v_la ** 2
    print(f"      -> interlayer TA = {v_ta:.1f} m/s")
    print(f"      -> C44 = {c44 / 1e9:.3f} GPa, C33 = {c33 / 1e9:.3f} GPa")
    if a.ref_vta:
        print(f"      -> vs reference TA {a.ref_vta:.1f} m/s: "
              f"{(v_ta / a.ref_vta - 1) * 100:+.1f} %")
    if a.ref_c44:
        print(f"      -> vs reference C44 {a.ref_c44:.3f} GPa: "
              f"{(c44 / 1e9 / a.ref_c44 - 1) * 100:+.1f} %")

    # --- q in the layer plane: one branch must be the c-polarised TA with the
    # same sqrt(C44/rho), which is an internal consistency check
    v_p, f0_p = fit_velocity(ph, u_ip, Bp, a.window)
    v_c44 = v_ta
    print(f"   q in plane   {', '.join(f'{x:.1f}' for x in np.sort(v_p))} m/s")
    closest = float(v_p[np.argmin(np.abs(v_p - v_c44))])
    c44_ip = rho * closest ** 2 / 1e9
    print(f"      -> one must equal sqrt(C44/rho) = {v_c44:.1f} m/s "
          f"(same C44 as above); closest branch {closest:.1f} m/s "
          f"({(closest / v_c44 - 1) * 100:+.1f} %) -> C44(route 2) = {c44_ip:.3f} GPa")
    print()

    if a.json:
        import json
        import os
        rec = {"label": tag, "fc2": a.fc2, "poscar": a.poscar,
               "supercell": a.supercell, "rho_kg_m3": round(float(rho), 1),
               "v_interlayer_TA_m_s": round(v_ta, 1),
               "v_interlayer_LA_m_s": round(v_la, 1),
               "C44_GPa": round(float(c44 / 1e9), 3),
               "C33_GPa": round(float(c33 / 1e9), 3),
               "v_inplane_m_s": sorted(round(float(x), 1) for x in v_p),
               "C44_GPa_from_inplane_route": round(float(c44_ip), 3),
               "gamma_point_offset_THz": round(float(f0_c[order[0]]), 6)}
        if a.ref_vta:
            rec["reference_TA_m_s"] = a.ref_vta
            rec["deviation_TA_pct"] = round((v_ta / a.ref_vta - 1) * 100, 2)
        if a.ref_c44:
            rec["reference_C44_GPa"] = a.ref_c44
            rec["deviation_C44_pct"] = round((c44 / 1e9 / a.ref_c44 - 1) * 100, 2)
        data = []
        if os.path.exists(a.json):
            data = json.load(open(a.json))
        data.append(rec)
        json.dump(data, open(a.json, "w", encoding="utf-8"), indent=1)
        print(f"   appended to {a.json}")


if __name__ == "__main__":
    main()
