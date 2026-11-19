#!/usr/bin/env python3
"""Export every dataset plotted in the HTML report as plain CSV / JSON.

Reads the repository files (no recomputation except the Mn2In2Se5 dispersion,
which is interpolated from its committed fc2.hdf5) and writes report/plot_data/.
plot_figures.py draws all figures from that folder alone.

Needs numpy, h5py, pyyaml, phonopy (+ seekpath).   python3 report/export_plot_data.py
"""
import csv
import json
import shutil
from pathlib import Path

import h5py
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "test" / "figures"
OUT = ROOT / "report" / "plot_data"
MAT = ("MnIn2Se4", "Mn2In2Se5")
SC_DIR = {"MnIn2Se4": FIG / "scattering", "Mn2In2Se5": FIG / "scattering" / "mn2in2se5"}
QUOTED = [  # 300 K values to quote (PROGRESS_REPORT_20261118.md §2, REPORT §9.1)
    ("MnIn2Se4", "inplane", 4.67, "LBTE + coherence; FC3 cutoff about -10 %"),
    ("MnIn2Se4", "kzz", 0.714, "mesh-converged + coherence; FC3 cutoff about -10 %"),
    ("Mn2In2Se5", "inplane", 3.37, "upper bound (FC3 cutoff not converged)"),
    ("Mn2In2Se5", "kzz", 1.089, "upper bound (FC3 cutoff not converged)"),
]


def load(rel):
    return json.loads((FIG / rel).read_text())


def write_csv(name, header, rows, comment=None):
    with open(OUT / name, "w", newline="", encoding="utf-8") as f:
        if comment:
            for line in comment.strip().splitlines():
                f.write(f"# {line}\n")
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def tex_label(s):
    return (s.replace("$", "").replace("\\mathrm{", "").replace("}", "")
             .replace("\\Gamma", "Γ").replace("_", ""))


# --------------------------------------------------------------------------- bands
def bands_mnin2se4():
    b = yaml.safe_load((ROOT / "test/MnIn2Se4/fc2_442_phonon/band-dft-cpu.yaml").read_text())
    segs, i = [], 0
    for n in b["segment_nqpoint"]:
        blk = b["phonon"][i:i + n]
        i += n
        segs.append((np.array([p["distance"] for p in blk]),
                     np.array([[bb["frequency"] for bb in p["band"]] for p in blk])))
    labels = []
    for k, (a, c) in enumerate(b["labels"]):
        a, c = tex_label(a), tex_label(c)
        if k == 0:
            labels.append(a)
        elif a != labels[-1]:
            labels[-1] += "|" + a
        labels.append(c)
    return segs, labels, "test/MnIn2Se4/fc2_442_phonon/band-dft-cpu.yaml (4x4x2 FC2)"


def bands_mn2in2se5():
    import phonopy
    d = ROOT / "Mn2In2Se5" / "fc"
    ph = phonopy.load(str(d / "phono3py_disp.yaml"), force_constants_filename=str(d / "fc2.hdf5"),
                      produce_fc=False, is_nac=False, log_level=0)
    ph.auto_band_structure(npoints=101, write_yaml=False)
    bs = ph.band_structure
    segs = [(np.array(x), np.array(f)) for x, f in zip(bs.distances, bs.frequencies)]
    lab, conn = [tex_label(x) for x in bs.labels], bs.path_connections
    labels, j = [lab[0]], 1
    for k, c in enumerate(conn):
        labels.append(lab[j])
        j += 1
        if not c and k < len(conn) - 1:
            labels[-1] += "|" + lab[j]
            j += 1
    return segs, labels, "Mn2In2Se5/fc/fc2.hdf5 (4x4x1 FC2), phonopy interpolation, seekpath path"


def export_bands():
    for m, fn in (("MnIn2Se4", bands_mnin2se4), ("Mn2In2Se5", bands_mn2in2se5)):
        segs, labels, src = fn()
        nb = segs[0][1].shape[1]
        rows = [[k, f"{x:.6f}"] + [f"{v:.5f}" for v in fr]
                for k, (xs, fs) in enumerate(segs) for x, fr in zip(xs, fs)]
        write_csv(f"bands_{m}.csv", ["segment", "distance"] + [f"band{i + 1}_THz" for i in range(nb)], rows,
                  f"{m} phonon dispersion, NAC off. Source: {src}.\n"
                  "Path labels and their positions: bands_" + m + "_ticks.csv")
        ticks = [s[0][0] for s in segs] + [segs[-1][0][-1]]
        write_csv(f"bands_{m}_ticks.csv", ["distance", "label"],
                  [[f"{t:.6f}", l] for t, l in zip(ticks, labels)])


# --------------------------------------------------------------------------- kappa
def export_kappa():
    for m, f in (("MnIn2Se4", "data_mnin2se4.json"), ("Mn2In2Se5", "data_mnin2se5.json")):
        d = load(f"kappa_vs_T/{f}")
        write_csv(f"kappa_intrinsic_{m}.csv", ["T_K", "kappa_inplane", "kappa_zz"],
                  zip(d["T_K"], d["inplane"], d["kzz_rest"]),
                  f"{m} intrinsic kappa_L (W/mK): ordered crystal, 3-phonon + natural isotopes, RTA, "
                  f"mesh {d['mesh']}.\nkappa_zz = mesh-converged term (Gamma-Z quadrature error removed). "
                  f"Source: test/figures/kappa_vs_T/{f}")
    write_csv("kappa_quoted_300K.csv", ["material", "component", "kappa_WmK", "note"], QUOTED,
              "300 K values to quote (stars in the intrinsic kappa(T) figure)")
    for m, f in (("MnIn2Se4", "data_mnin2se4_vca.json"), ("Mn2In2Se5", "data_mn2in2se5_vca.json")):
        d = load(f"kappa_vs_T/{f}")
        v = d["vca"]
        write_csv(f"kappa_ordered_vca_{m}.csv",
                  ["T_K", "ordered_inplane", "ordered_kzz", "vca_inplane", "vca_kzz"],
                  zip(d["T_K"], d["baseline_inplane"], d["baseline_kzz_rest"], d["vca_inplane"], d["vca_kzz_rest"]),
                  f"{m} single-crystal kappa (W/mK), RTA, mesh {d['mesh']}, natural isotopes, kappa_zz mesh-converged.\n"
                  f"ordered = real masses; vca = averaged site masses + mass-variance scattering, "
                  f"antisite fraction {v.get('x')}, masses {v.get('mass')}, g {v.get('g')}.\n"
                  f"Source: test/figures/kappa_vs_T/{f}")
        e = load(f"calc_vs_experiment/extrinsic_{m.lower()}.json")
        s = e["s"]
        rows = []
        for r in e["rows"]:
            kx, kz = r["kxx"], r["kzz"]
            rows.append([r["T"], kx, kz, round((1 - s / 2) * kx + s / 2 * kz, 5), round(s * kx + (1 - s) * kz, 5), s])
        write_csv(f"kappa_L3_{m}.csv", ["T_K", "kxx", "kzz", "Perp", "Para", "s"], rows,
                  f"{m} L3 extrinsic model (fitted to the 300-623 K experiment), texture s = {s}.\n"
                  "Perp = (1-s/2) kxx + (s/2) kzz ; Para = s kxx + (1-s) kzz.\n"
                  f"Source: test/figures/calc_vs_experiment/extrinsic_{m.lower()}.json")
    shutil.copy(ROOT / "experiment" / "kappa_tot_digitized.csv", OUT / "experiment_kappa_tot_digitized.csv")
    shutil.copy(ROOT / "experiment" / "sound_velocity_MnIn2Se4.csv", OUT / "experiment_sound_velocity_MnIn2Se4.csv")
    shutil.copy(ROOT / "report" / "independent_checks.json", OUT / "independent_checks.json")


# --------------------------------------------------------------------------- scattering
def export_scattering():
    for m in MAT:
        cols, T = {}, None
        for run in ("baseline", "vca"):
            with h5py.File(SC_DIR[m] / f"kappa_{run}.hdf5", "r") as f:
                cols[f"freq_{run}_THz"] = f["frequency"][:].ravel()
                cols[f"tau_inv_3ph_{run}_ps"] = f["gamma"][0].ravel() * 2.0
                cols[f"tau_inv_disorder_{run}_ps"] = f["gamma_isotope"][:].ravel() * 2.0
                T = float(f["temperature"][0])
        keys = list(cols)
        rows = [[f"{cols[k][i]:.5g}" for k in keys] for i in range(len(cols[keys[0]]))]
        write_csv(f"scattering_{m}_{T:.0f}K.csv", keys, rows,
                  f"{m} per-mode scattering rates at {T:.0f} K, 20^3 mesh, largest FC3 cutoff; one row per "
                  "(irreducible q-point, band).\ntau^-1 = 2 Gamma (ps^-1). baseline: real masses, disorder = natural "
                  "isotopes. vca: averaged masses, disorder = mass variance.\n"
                  f"Source: {SC_DIR[m].relative_to(ROOT)}/kappa_baseline.hdf5, kappa_vca.hdf5")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    export_bands()
    export_kappa()
    export_scattering()
    print(f"wrote {len(list(OUT.iterdir()))} files to {OUT.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
