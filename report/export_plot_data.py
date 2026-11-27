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
    # 2026-11-27: in-plane value carries the spin-disorder factor R_xx(300 K) = 0.857
    # (2.89 = 3.37 x 0.857, dft/spin_disorder_fc2/results/R_xx_vs_T.json); kappa_zz unchanged
    ("Mn2In2Se5", "inplane", 2.89, "LBTE + coherence, x R_xx = 0.857 for spin disorder;"
     " next FC3-cutoff step >= -2.1 % (hold-out test)"),
    ("Mn2In2Se5", "kzz", 1.089, "mesh-converged + coherence; next FC3-cutoff step >= -3.3 % (hold-out test)"),
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
    vca_x = {}
    for m in MAT:
        x = json.loads((SC_DIR[m] / "vca_x.json").read_text())["x"]   # antisite fraction of kappa_vca.hdf5
        vca_x[m] = x
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
                  f"isotopes. vca: averaged masses, disorder = mass variance, antisite fraction {x}.\n"
                  f"Source: {SC_DIR[m].relative_to(ROOT)}/kappa_baseline.hdf5, kappa_vca.hdf5, vca_x.json")
    (OUT / "scattering_vca_x.json").write_text(json.dumps(vca_x, indent=1) + "\n")


# --------------------------------------------------------------------------- convergence / solvers
# MnIn2Se4 q-mesh series: REPORT_MnIn2Se4_20261115.md §6.3 (RTA, 300 K, FC3 6.45 A, isotopes on;
# rest term from test/MnIn2Se4/kappa_zz_decompose.py on test/MnIn2Se4/kappa_zz_scan/).
QMESH_MNIN2SE4 = [  # N, in-plane, kzz total, kzz rest
    (11, 3.5375, 1.2495, 0.6251), (14, 3.6474, 1.0934, 0.6444), (18, 3.6945, 0.9435, 0.6461),
    (24, 3.7020, 0.8052, 0.6457), (30, 3.6946, 0.7645, 0.6473)]
SOLVERS = [  # material, quantity, value, source
    ("MnIn2Se4", "inplane_RTA_14", 3.6474, "REPORT §6.3 / §9, 14^3, RTA"),
    ("MnIn2Se4", "inplane_LBTE_14", 4.5219, "REPORT §9.1, 14^3, LBTE (HanHai job 412627)"),
    ("MnIn2Se4", "kzz_LBTE_minus_RTA", 4e-5, "REPORT §9: LBTE = RTA for kappa_zz to 4e-5"),
    ("MnIn2Se4", "mesh_correction_inplane", 1.013, "REPORT §9.1: RTA series 14^3 -> 30^3"),
    ("MnIn2Se4", "coherence_inplane", 0.0908, "REPORT §7, Wigner-RTA, phono3py 3.24.0, 14^3"),
    ("MnIn2Se4", "coherence_kzz", 0.0661, "REPORT §7, Wigner-RTA, phono3py 3.24.0, 14^3"),
    ("MnIn2Se4", "population_kzz_rest", 0.648, "REPORT §6.5, mesh-converged rest term"),
    ("Mn2In2Se5", "inplane_RTA_10", 2.8220, "analysis/shengbte_crosscheck/README.md, 10^3, job 415925"),
    ("Mn2In2Se5", "inplane_LBTE_10", 2.9879, "analysis/shengbte_crosscheck/README.md, 10^3, job 415925"),
    ("Mn2In2Se5", "inplane_RTA_32", 3.048, "test/Mn2In2Se5/iso_rerun/, 32^3, isotopes on"),
    ("Mn2In2Se5", "coherence_inplane", 0.1466, "Wigner-RTA (SMM19, phono3py 4.4.0), 16^3, job 419208"),
    ("Mn2In2Se5", "coherence_kzz", 0.112, "REPORT §9.1 (0.1115 at 16^3, soft-mode check job 440011)"),
    ("Mn2In2Se5", "population_kzz_rest", 0.977, "rest-term plateau 24^3-32^3 x isotope ratio 0.9961"),
]


def export_convergence():
    rows = []
    for m, f in (("MnIn2Se4", "data_mnin2se4.json"), ("Mn2In2Se5", "data_mnin2se5.json")):
        d = load(f"cutoff/{f}")
        for i, c in enumerate(d["cutoff_A"]):
            rows.append([m, c, d["mesh"], d["inplane_population"][i], d["kzz_rest"][i],
                         d["inplane_coherence"][i], d["kzz_coherence"][i]])
    write_csv("convergence_fc3_cutoff.csv",
              ["material", "cutoff_A", "mesh", "inplane_population", "kzz_rest", "inplane_coherence_4.4.0",
               "kzz_coherence_4.4.0"], rows,
              "FC3-cutoff scan, 300 K, RTA population, natural isotopes, kappa_zz = mesh-converged rest term.\n"
              "Coherence columns from phono3py 4.4.0 SMM19 (not valid for MnIn2Se4, see REPORT §7).\n"
              "Source: test/figures/cutoff/data_mnin2se4.json, data_mnin2se5.json")
    rows = [["MnIn2Se4", n, a, b, c, "on"] for n, a, b, c in QMESH_MNIN2SE4]
    dec = json.loads((ROOT / "test/Mn2In2Se5/kappa_zz_decomposition/kappa_zz_mesh_decomposition.json").read_text())
    for c in sorted(dec["cases"], key=lambda c: c["mesh"][0]):
        if len(set(c["mesh"])) == 1:
            i = c["intrinsic"]
            rows.append(["Mn2In2Se5", c["mesh"][0], round(i["inplane_total"], 4), round(i["kappa_zz_total"], 4),
                         round(i["kappa_zz_rest"], 4), "off"])
    write_csv("convergence_qmesh.csv", ["material", "N", "inplane", "kzz_total", "kzz_rest", "isotopes"], rows,
              "Isotropic N^3 q-mesh scans, 300 K, RTA. kzz_rest = kappa_zz without the first Gamma-Z q-point "
              "(quadrature error).\nMnIn2Se4: REPORT_MnIn2Se4_20261115.md §6.3 (FC3 6.45 A). Mn2In2Se5: "
              "test/Mn2In2Se5/kappa_zz_decomposition/kappa_zz_mesh_decomposition.json (FC3 6.51 A, isotopes off;\n"
              "isotopes on: x0.9884 in-plane, x0.9961 kzz_rest, test/Mn2In2Se5/iso_rerun/)")
    write_csv("solver_terms.csv", ["material", "quantity", "value", "source"], SOLVERS,
              "RTA / LBTE / Wigner terms behind the quoted 300 K values (W/mK unless a ratio)")


def export_disorder_sweep():
    """kappa at 300 K vs antisite fraction x from the L3 bases (24^3, ordered masses + mass disorder,
    no boundary). Shows that mass disorder lowers in-plane kappa strongly but kappa_zz only weakly."""
    import sys
    sys.path.insert(0, str(ROOT / "analysis" / "extrinsic"))
    from kappa_extrinsic import Model
    rows = []
    for m, d, xs in (("MnIn2Se4", "MnIn2Se4_442_m24", (0.0, 0.1, 0.2, 0.3, 0.44, 0.55, 0.6667)),
                     ("Mn2In2Se5", "Mn2In2Se5_m24", (0.0, 0.1, 0.13, 0.19, 0.3, 0.4, 0.5))):
        mod = Model(str(ROOT / "analysis" / "extrinsic" / "results" / d / "bases.npz"))
        it = int(np.argmin(np.abs(mod.T - 300.0)))
        for x in xs:
            a, z = mod.kappa(x, None)
            rows.append([m, x, round(float(a[it]), 4), round(float(z[it]), 4)])
    write_csv("disorder_sweep_300K.csv", ["material", "x", "inplane", "kzz_total"], rows,
              "300 K, 24^3, three-phonon + Mn/In antisite mass disorder (ordered masses, Tamura), no boundary.\n"
              "x = fraction of all Mn on In sites; random limit 2/3 (MnIn2Se4), 1/2 (Mn2In2Se5).\n"
              "kzz_total still contains the Gamma-Z quadrature error; use relative changes only.\n"
              "Source: analysis/extrinsic/results/*_m24/bases.npz via kappa_extrinsic.Model")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    export_convergence()
    export_disorder_sweep()
    export_bands()
    export_kappa()
    export_scattering()
    print(f"wrote {len(list(OUT.iterdir()))} files to {OUT.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
