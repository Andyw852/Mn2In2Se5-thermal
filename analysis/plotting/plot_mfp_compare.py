#!/usr/bin/env python3
"""
Compare mean free path (MFP) between baseline and VCA from phono3py kappa_*.hdf5.
MFP = |v_group| * tau, with tau = 1/(2*Gamma_total).
Gamma_total = gamma (3-phonon) + gamma_isotope (mass-variance disorder).

Two panels:
  (left)  MFP vs frequency, baseline vs VCA
  (right) cumulative kappa vs MFP (how far-travelling phonons carry the heat)

Usage (inside test_vca/):
    python plot_mfp_compare.py            # default ~300K
    python plot_mfp_compare.py 300
Output: mfp_compare.png / .pdf
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


def load(fn, Ttarget):
    with h5py.File(fn, "r") as f:
        keys = list(f.keys())
        T = f["temperature"][:]
        it = int(np.argmin(np.abs(T - Ttarget)))
        freq = f["frequency"][:]                       # (nq, nb) THz
        gamma = f["gamma"][it].copy()                  # (nq, nb) THz
        if "gamma_isotope" in keys:
            gamma = gamma + f["gamma_isotope"][:]      # add disorder channel
        gv = f["group_velocity"][:]                    # (nq, nb, 3) THz*Ang
        mk = f["mode_kappa"][it]                       # (nq, nb, 6) W/m/K per mode
        w = f["weight"][:] if "weight" in keys else np.ones(freq.shape[0])
    # tau (ps): Gamma in THz -> 2*Gamma in 1/ps -> tau = 1/(2*Gamma)
    with np.errstate(divide="ignore", invalid="ignore"):
        tau = 1.0 / (2.0 * gamma)                      # ps
        tau[~np.isfinite(tau)] = 0.0
    vmag = np.linalg.norm(gv, axis=2)                  # THz*Ang = 100 m/s
    # MFP = |v|*tau: (THz*Ang)*(ps) ; 1 THz*Ang = 1e12*1e-10 m/s = 100 m/s
    #   *1e-12 s = 100 m/s * 1e-12 s ... do it in nm:
    # |v| in THz*Ang = 100 m/s ; tau in ps=1e-12 s -> MFP(m)=100*|v|*tau*1e-12
    mfp_nm = 100.0 * vmag * tau * 1e-12 * 1e9          # nm
    kmode = mk[:, :, :3].sum(axis=2) / 3.0             # isotropic per-mode kappa
    return {"f": freq.ravel(), "mfp": mfp_nm.ravel(),
            "k": kmode.ravel(), "T": T[it], "keys": keys}


b = load(f_base, Ttarget)
v = load(f_vca, Ttarget)
print(f"baseline keys: {b['keys']}")
print(f"Temperature point: {b['T']:.0f}K")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))

# Left: MFP vs frequency
mb = b["mfp"] > 0
mv = v["mfp"] > 0
ax1.scatter(b["f"][mb], b["mfp"][mb], s=6, alpha=0.35, color="#1f77b4", label="baseline")
ax1.scatter(v["f"][mv], v["mfp"][mv], s=6, alpha=0.35, color="#d62728", label="VCA")
ax1.set_xlabel("Frequency (THz)")
ax1.set_ylabel("Mean free path (nm)")
ax1.set_title(f"MFP vs frequency (~{b['T']:.0f}K)")
ax1.set_yscale("log")
ax1.grid(True, which="both", ls=":", alpha=0.4)
ax1.legend(fontsize=9, frameon=False)


def cumulative(mfp, k):
    m = mfp > 0
    mfp, k = mfp[m], k[m]
    order = np.argsort(mfp)
    return mfp[order], np.cumsum(k[order])


xb, cb = cumulative(b["mfp"], b["k"])
xv, cv = cumulative(v["mfp"], v["k"])
ax2.plot(xb, cb, "-", color="#1f77b4", lw=1.8, label="baseline")
ax2.plot(xv, cv, "-", color="#d62728", lw=1.8, label="VCA")
ax2.set_xlabel("Mean free path (nm)")
ax2.set_ylabel(r"Cumulative $\kappa_L$ (W m$^{-1}$ K$^{-1}$)")
ax2.set_title("Cumulative kappa vs MFP")
ax2.set_xscale("log")
ax2.grid(True, which="both", ls=":", alpha=0.4)
ax2.legend(fontsize=9, frameon=False)

fig.tight_layout()
fig.savefig("mfp_compare.png", dpi=200)
fig.savefig("mfp_compare.pdf")
print("Saved mfp_compare.png / mfp_compare.pdf")
