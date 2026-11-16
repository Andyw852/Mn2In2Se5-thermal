#!/usr/bin/env python3
"""3-phonon vs disorder (VCA) scattering rates at the largest FC3 cutoff.

This is analysis/plotting/plot_scattering_compare.py (the project's own script)
with only the y-axis ranges added: the raw data contain near-Gamma acoustic modes
whose rates are ~1e-38, which would otherwise stretch both panels over 40 decades.
Floor set to 1e-21 as in the reference figure.

  left  : 3-phonon scattering, ordered baseline vs VCA
  right : disorder elastic scattering, VCA mass-variance vs natural isotope, omega^2 guide

  python plot_scattering_figure.py [T] [kappa_baseline.hdf5] [kappa_vca.hdf5]
"""
import sys

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Ttarget = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0
f_base = sys.argv[2] if len(sys.argv) > 2 else "kappa_baseline.hdf5"
f_vca = sys.argv[3] if len(sys.argv) > 3 else "kappa_vca.hdf5"

FACTOR = 2.0          # tau^-1 = 2*Gamma, Gamma in THz -> ps^-1 (no 2*pi)
YMIN, YMAX = 1e-21, 2.0


def load(fn, Ttarget):
    with h5py.File(fn, "r") as f:
        keys = list(f.keys())
        freq = f["frequency"][:]
        T = f["temperature"][:]
        it = int(np.argmin(np.abs(T - Ttarget)))
        gamma = f["gamma"][it]
        giso = f["gamma_isotope"][:] if "gamma_isotope" in keys else None
    fr = freq.ravel()
    order = np.argsort(fr)
    out = {"f": fr[order], "g3": gamma.ravel()[order] * FACTOR, "T": T[it], "keys": keys}
    if giso is not None:
        out["giso"] = giso.ravel()[order] * FACTOR
    return out


b = load(f_base, Ttarget)
v = load(f_vca, Ttarget)
print("baseline keys:", b["keys"])
print("vca keys     :", v["keys"])
print("Temperature point: baseline %.0f K, vca %.0f K" % (b["T"], v["T"]))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))

ax1.scatter(b["f"], b["g3"], s=6, alpha=0.4, color="#1f77b4", label="baseline 3-phonon")
ax1.scatter(v["f"], v["g3"], s=6, alpha=0.4, color="#d62728", label="VCA 3-phonon")
ax1.set_xlabel("Frequency (THz)")
ax1.set_ylabel(r"$\tau^{-1}$ (ps$^{-1}$)")
ax1.set_title("3-phonon scattering (~%.0f K)" % b["T"])
ax1.set_yscale("log")
ax1.set_ylim(YMIN, YMAX)
ax1.grid(True, which="both", ls=":", alpha=0.4)
ax1.legend(fontsize=9, frameon=False, loc="lower right")

ax2.scatter(v["f"], v.get("giso", np.zeros_like(v["f"])), s=6, alpha=0.5,
            color="#d62728", label="VCA disorder (mass-variance)")
if "giso" in b:
    ax2.scatter(b["f"], b["giso"], s=6, alpha=0.5, color="#1f77b4",
                label="baseline natural isotope")
ff = np.linspace(max(v["f"].min(), 0.1), v["f"].max(), 50)
if "giso" in v and np.any(v["giso"] > 0):
    mask = v["f"] > 1
    scale = np.median(v["giso"][mask] / v["f"][mask] ** 2)
    ax2.plot(ff, scale * ff ** 2, "k--", lw=1, alpha=0.6, label=r"$\propto\omega^2$")
ax2.set_xlabel("Frequency (THz)")
ax2.set_ylabel(r"$\tau^{-1}$ (ps$^{-1}$)")
ax2.set_title("Disorder elastic scattering")
ax2.set_yscale("log")
ax2.set_ylim(YMIN, YMAX)
ax2.grid(True, which="both", ls=":", alpha=0.4)
ax2.legend(fontsize=9, frameon=False, loc="lower right")

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig("scattering_compare." + ext, dpi=200)
plt.close(fig)
print("Saved scattering_compare.png / .pdf")
