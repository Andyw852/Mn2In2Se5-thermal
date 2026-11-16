#!/usr/bin/env python3
"""
Compare baseline (ordered, real masses) vs VCA (averaged masses + mass-variance)
lattice thermal conductivity. Reads kappa_baseline.hdf5 / kappa_vca.hdf5 and plots
in-plane (kxx=kyy) and out-of-plane (kzz) vs temperature.

Usage (run inside test_vca/):
    python plot_kappa_compare.py
    python plot_kappa_compare.py kappa_baseline.hdf5 kappa_vca.hdf5
Output: kappa_compare.png and kappa_compare.pdf
"""
import sys
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

f_base = sys.argv[1] if len(sys.argv) > 1 else "kappa_baseline.hdf5"
f_vca  = sys.argv[2] if len(sys.argv) > 2 else "kappa_vca.hdf5"


def read_kappa(fn):
    with h5py.File(fn, "r") as f:
        T = f["temperature"][:]
        k = f["kappa"][:]        # (nT, 6): xx yy zz yz xz xy
    kin = 0.5 * (k[:, 0] + k[:, 1])          # in-plane (kxx=kyy)
    kout = k[:, 2]                            # out-of-plane (kzz)
    kiso = (k[:, 0] + k[:, 1] + k[:, 2]) / 3  # isotropic average
    return T, kin, kout, kiso


Tb, kin_b, kout_b, kiso_b = read_kappa(f_base)
Tv, kin_v, kout_v, kiso_v = read_kappa(f_vca)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))

C_BASE, C_VCA = "#1f77b4", "#d62728"

# Left: in-plane vs out-of-plane, baseline (solid) vs VCA (dashed)
ax1.plot(Tb, kin_b, "-o", color=C_BASE, ms=4, label="baseline in-plane")
ax1.plot(Tb, kout_b, "-s", color=C_BASE, ms=4, mfc="white", label="baseline out-of-plane")
ax1.plot(Tv, kin_v, "--o", color=C_VCA, ms=4, label="VCA in-plane")
ax1.plot(Tv, kout_v, "--s", color=C_VCA, ms=4, mfc="white", label="VCA out-of-plane")
ax1.set_xlabel("Temperature (K)")
ax1.set_ylabel(r"$\kappa_L$ (W m$^{-1}$ K$^{-1}$)")
ax1.set_title("In-plane / Out-of-plane")
ax1.set_yscale("log")
ax1.grid(True, which="both", ls=":", alpha=0.4)
ax1.legend(fontsize=8, frameon=False)

# Right: isotropic average
ax2.plot(Tb, kiso_b, "-o", color=C_BASE, ms=4, label="baseline (ordered)")
ax2.plot(Tv, kiso_v, "--o", color=C_VCA, ms=4, label="VCA (mass disorder)")
ax2.set_xlabel("Temperature (K)")
ax2.set_ylabel(r"$\kappa_L^{\mathrm{iso}}$ (W m$^{-1}$ K$^{-1}$)")
ax2.set_title("Isotropic average")
ax2.grid(True, ls=":", alpha=0.4)
ax2.legend(fontsize=9, frameon=False)

fig.tight_layout()
fig.savefig("kappa_compare.png", dpi=200)
fig.savefig("kappa_compare.pdf")
print("Saved kappa_compare.png / kappa_compare.pdf")

# 300 K table
def at300(T, *arrs):
    i = int(np.argmin(np.abs(T - 300)))
    return (T[i],) + tuple(a[i] for a in arrs)

tb, kin_b3, kout_b3, kiso_b3 = at300(Tb, kin_b, kout_b, kiso_b)
tv, kin_v3, kout_v3, kiso_v3 = at300(Tv, kin_v, kout_v, kiso_v)
print(f"\n{'T~300K':<10}{'in-plane':>10}{'out-plane':>11}{'iso':>10}")
print(f"{'baseline':<10}{kin_b3:>10.3f}{kout_b3:>11.3f}{kiso_b3:>10.3f}")
print(f"{'VCA':<10}{kin_v3:>10.3f}{kout_v3:>11.3f}{kiso_v3:>10.3f}")
print(f"{'ratio':<10}{kin_b3/kin_v3:>9.2f}x{kout_b3/kout_v3:>10.2f}x{kiso_b3/kiso_v3:>9.2f}x")
