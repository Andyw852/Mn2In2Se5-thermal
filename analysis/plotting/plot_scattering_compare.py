#!/usr/bin/env python3
"""
Compare baseline vs VCA scattering rates vs frequency.
Extracts from kappa_*.hdf5:
  gamma          = 3-phonon linewidth (THz)
  gamma_isotope  = isotope / mass-variance (inversion disorder) linewidth (THz)

phono3py convention: tau = 1/(2*Gamma), so tau^-1 = 2*Gamma.
Gamma is in THz, so tau^-1 = 2*Gamma is in ps^-1. NO factor of 2*pi.

Usage (run inside test_vca/):
    python plot_scattering_compare.py                 # default ~300K
    python plot_scattering_compare.py 300             # temperature (K)
    python plot_scattering_compare.py 300 kappa_baseline.hdf5 kappa_vca.hdf5
Output: scattering_compare.png / .pdf
"""
import sys
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Ttarget = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0
f_base = sys.argv[2] if len(sys.argv) > 2 else "kappa_baseline.hdf5"
f_vca  = sys.argv[3] if len(sys.argv) > 3 else "kappa_vca.hdf5"

FACTOR = 2.0   # tau^-1 = 2*Gamma ; Gamma(THz) -> ps^-1 ; no 2*pi


def load(fn, Ttarget):
    with h5py.File(fn, "r") as f:
        keys = list(f.keys())
        freq = f["frequency"][:]                 # (nq, nband) THz
        T = f["temperature"][:]
        it = int(np.argmin(np.abs(T - Ttarget)))
        gamma = f["gamma"][it]                    # (nq, nband) 3-phonon
        giso = f["gamma_isotope"][:] if "gamma_isotope" in keys else None
    fr = freq.ravel()
    g3 = gamma.ravel()
    order = np.argsort(fr)
    out = {"f": fr[order], "g3": g3[order] * FACTOR, "T": T[it], "keys": keys}
    if giso is not None:
        out["giso"] = giso.ravel()[order] * FACTOR
    return out


b = load(f_base, Ttarget)
v = load(f_vca, Ttarget)
print(f"baseline keys: {b['keys']}")
print(f"vca keys     : {v['keys']}")
print(f"Temperature point: baseline {b['T']:.0f}K, vca {v['T']:.0f}K")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))

# Left: 3-phonon scattering, baseline vs VCA (VCA shifted slightly by avg mass)
ax1.scatter(b["f"], b["g3"], s=6, alpha=0.4, color="#1f77b4", label="baseline 3-phonon")
ax1.scatter(v["f"], v["g3"], s=6, alpha=0.4, color="#d62728", label="VCA 3-phonon")
ax1.set_xlabel("Frequency (THz)")
ax1.set_ylabel(r"$\tau^{-1}$ (ps$^{-1}$)")
ax1.set_title(f"3-phonon scattering (~{b['T']:.0f}K)")
ax1.set_yscale("log")
ax1.grid(True, which="both", ls=":", alpha=0.4)
ax1.legend(fontsize=9, frameon=False)

# Right: disorder (mass-variance) scattering -- the core of VCA
ax2.scatter(v["f"], v.get("giso", np.zeros_like(v["f"])), s=6, alpha=0.5,
            color="#d62728", label="VCA disorder (mass-variance)")
if "giso" in b:
    ax2.scatter(b["f"], b["giso"], s=6, alpha=0.5, color="#1f77b4",
                label="baseline natural isotope")
# reference omega^2 line
ff = np.linspace(max(v["f"].min(), 0.1), v["f"].max(), 50)
if "giso" in v and np.any(v["giso"] > 0):
    mask = v["f"] > 1
    scale = np.median(v["giso"][mask] / v["f"][mask] ** 2)
    ax2.plot(ff, scale * ff ** 2, "k--", lw=1, alpha=0.6, label=r"$\propto\omega^2$")
ax2.set_xlabel("Frequency (THz)")
ax2.set_ylabel(r"$\tau^{-1}$ (ps$^{-1}$)")
ax2.set_title("Disorder elastic scattering")
ax2.set_yscale("log")
ax2.grid(True, which="both", ls=":", alpha=0.4)
ax2.legend(fontsize=9, frameon=False)

fig.tight_layout()
fig.savefig("scattering_compare.png", dpi=200)
fig.savefig("scattering_compare.pdf")
print("Saved scattering_compare.png / scattering_compare.pdf")
