#!/usr/bin/env python3
"""Extract the layer-normal sound velocity from the 4x4x2 fc2 and compare it with
the DFT elastic constants.

⚠️ Pitfalls hit before (user review 2026-11-13 + re-checked by this script):

1. **phonopy q points use the reciprocal-space coordinates of the *chosen
   primitive* cell, not of the POSCAR cell.**
   Here `primitive_matrix="auto"` switched the primitive cell to the rhombohedral
   one (7 atoms), whose b1/b2/b3 all make an angle of 84.91° with the layer normal
   (i.e. they lie almost in the layer plane).
   a1 and a2 of the POSCAR lie in the plane → layer normal = (0,0,1)_POSCAR.
   The same interlayer direction in the two sets of coordinates:
        (0,0,t)_POSCAR  ≡  (t,t,t)_phonopy_primitive
   Computing frequencies directly from (0,0,t) and then converting |q| with the
   POSCAR b3 gives bogus values of 4179/6256 m/s; the ratio of the two |b3|,
   3.760, is exactly the error factor.

2. **The Γ point carries a numerical offset** (f0 = 0.000719 THz) that pollutes
   ω/q at low q and makes q→0 diverge. One must fit
   f² = f0² + (v q / 2π/100)² and take the square root of the slope.

After the fix: TA = 1348.3 m/s (DFT 1338.3, +0.75%), LA = 3118.6 (DFT 3253.4,
-4.14%), C44 = 9.543 GPa (DFT 9.402, +1.50%). The result is stable for
t≤0.005~0.040.

Usage:
    python sound_velocity.py <directory containing fc2.hdf5 and phonopy_disp.yaml>
"""
import sys
import numpy as np
import phonopy
from phono3py.file_IO import read_fc2_from_hdf5

# reference sound velocities derived from the DFT elastic constants (same phonon
# structure, density 5249.3 kg/m^3)
REF_TA = 1338.319
REF_LA = 3253.408
REF_C44_GPa = 9.402


def main(d):
    ph = phonopy.load(f"{d}/phonopy_disp.yaml", produce_fc=False,
                      is_nac=False, log_level=0)
    ph.force_constants = read_fc2_from_hdf5(f"{d}/fc2.hdf5")

    uc = ph.unitcell
    pc = ph.primitive
    # ★ key point: use the primitive reciprocal vectors, not the unitcell ones
    Bp = 2 * np.pi * np.linalg.inv(np.asarray(pc.cell)).T
    rho = pc.masses.sum() / pc.volume * 1.66053906660e-27 / 1e-30

    a1, a2 = np.asarray(uc.cell)[0], np.asarray(uc.cell)[1]
    n = np.cross(a1, a2)
    n /= np.linalg.norm(n)
    print(f"phonopy primitive atom count = {len(pc)}")
    print(f"layer normal (Cartesian) = {np.round(n, 6)}")
    for i in range(3):
        ang = np.degrees(np.arccos(min(1.0, abs(Bp[i] @ n) / np.linalg.norm(Bp[i]))))
        print(f"  angle between primitive b{i+1} and the layer normal = {ang:.2f} deg")
    print(f"rho = {rho:.1f} kg/m^3\n")
    print("layer normal = the (t,t,t) direction of the phonopy primitive")
    print("fitting f^2 = f0^2 + (v·q/2π/100)^2 (removing the Γ-point offset)\n")

    ref = np.array([REF_TA, REF_TA, REF_LA])
    for win in (0.005, 0.010, 0.020, 0.040):
        qs, fs = [], []
        for t in np.linspace(0.001, win, 25):
            q = np.array([1.0, 1.0, 1.0])
            q = q / np.linalg.norm(q) * t
            ph.run_qpoints([q], with_eigenvectors=False)
            f = ph.get_qpoints_dict()["frequencies"][0]
            qs.append(np.linalg.norm(q @ Bp))
            fs.append(np.sort(f)[:3])
        qs, fs = np.array(qs), np.array(fs)

        v, f0 = [], []
        for b in range(3):
            A = np.vstack([np.ones_like(qs), qs ** 2]).T
            (intercept, slope), *_ = np.linalg.lstsq(A, fs[:, b] ** 2, rcond=None)
            v.append(np.sqrt(max(slope, 0.0)) * 2 * np.pi * 100.0)
            f0.append(np.sqrt(max(intercept, 0.0)))
        v = np.array(v)
        dev = (v / ref - 1) * 100
        print(f"window t<={win:.3f}: f0={f0[0]:.6f} THz  "
              f"TA1={v[0]:.1f} TA2={v[1]:.1f} LA={v[2]:.1f} m/s  "
              f"(deviation {dev[0]:+.2f}% {dev[1]:+.2f}% {dev[2]:+.2f}%)")
        print(f"                    C44={rho*v[0]**2/1e9:.3f} GPa "
              f"(DFT {REF_C44_GPa}, {(rho*v[0]**2/1e9/REF_C44_GPa-1)*100:+.2f}%)")
    print()
    print(f"DFT elastic-constant reference: TA={REF_TA:.1f} LA={REF_LA:.1f} m/s, C44={REF_C44_GPa} GPa")
    print("※ why (t,t,t): see note 1 in the header of this file.")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
