#!/usr/bin/env python3
"""Write the VASP runs for the paramagnetic (spin-disordered) FC2 of Mn2In2Se5.

  python prepare.py --repo <clone> --out <new workdir> --potcar <production 4x4x1 S4 POTCAR>
  python prepare.py --repo <clone> --list-configs        # print the spin configurations only

Cell: the production 4x4x1 supercell (144 atoms, 32 Mn) of Mn2In2Se5/fc/phono3py_disp.yaml,
atom order identical to Mn2In2Se5/fc/SPOSCAR, at the production FM-relaxed geometry.

Sets (one directory each, every set = 1 undisplaced "ref" + displaced frames):
  FM_ctl    all Mn +5                      ref + 8 random frames + 1 frame = -(frame 1)   (method control)
  AFM       Mn1 layer +5, Mn2 layer -5     ref + 8 random frames                          (ordered bound)
  PM01..08  spin-disordered (SQS-like)     ref + 2 random frames each                     (paramagnet)
Random frames: every atom displaced by 0.01 A in a random direction (fixed seed).
Spin configurations: 8 up / 8 down in each Mn layer; Mn-Mn pair correlations <s_i s_j> close
to 0 in the five shells below 7.2 A (two intralayer, three interlayer); no two configurations
related by a space-group operation of the supercell or by a global spin flip.

DFT: the production S4 force-run INCAR (Mn2In2Se5/data/dft_inputs/step4_disp/disp-00000/INCAR),
except MAGMOM, EDIFF 1E-8, NELM 300, LORBIT 11 (site moments), ISTART/ICHARG/LWAVE (displaced
frames restart from the ref WAVECAR so that they stay in the same spin state). The script
stops if any other tag differs from production.
"""
import argparse
import json
import re
import shutil
from pathlib import Path

import numpy as np
import phonopy
from phonopy.interface.vasp import write_vasp
from phonopy.structure.atoms import PhonopyAtoms

DIST = 0.01
SEED = 20261123
N_PM = 8
N_ORDERED = 8              # random frames of FM_ctl and AFM (4 is the minimum for a fit; 8 makes single modes good to ~0.3 %)
N_PM_FRAMES = 2
SHELL_RMAX = 7.2           # A; supercell half-width in plane is 2a = 8.06 A
# exchange from the 9- and 18-atom energies (dft/spin_phonon_gamma/results), meV per bond,
# E = E0 + sum_bonds J s_i s_j, |s| = 1; used only to predict the ref energies (a check)
J_GAMMA_FIT = {"inter_3.89": 1.58, "intra_4.03": 0.175, "inter_5.60": 1.90}

INCAR = """SYSTEM   = Mn2In2Se5 spin-disorder FC2: {tag}
ISTART   = {istart}
ICHARG   = {icharg}
GGA      = PE
IVDW     = 12
PREC     = Accurate
ENCUT    = 410
EDIFF    = 1E-8
LREAL    = .FALSE.
ADDGRID  = .TRUE.
LASPH    = .TRUE.
ALGO     = Normal
NELM     = 300
NELMIN   = 6
ISMEAR   = 0
SIGMA    = 0.05
IBRION   = -1
NSW      = 0
ISYM     = 0
LWAVE    = {lwave}
LCHARG   = .FALSE.
LORBIT   = 11
NCORE    = 4
KPAR     = 8
ISPIN    = 2
MAGMOM   = {magmom}
LMAXMIX  = 4
LDAU     = .TRUE.
LDAUTYPE = 2
LDAUL    = 2 -1 -1
LDAUU    = 4 0.0 0.0
LDAUJ    = 0.0 0.0 0.0
LDAUPRINT = 1
"""
KPOINTS = "Gamma 2 2 2 (production S4, KSPACING 0.038)\n0\nGamma\n   2   2   2\n0.0  0.0  0.0\n"
# tags that must equal the production S4 force-run INCAR
SHARED = ["GGA", "IVDW", "PREC", "ENCUT", "LREAL", "ADDGRID", "LASPH", "ALGO", "ISMEAR", "SIGMA",
          "IBRION", "NSW", "ISYM", "NCORE", "KPAR", "ISPIN", "LMAXMIX", "LDAU", "LDAUTYPE", "LDAUL",
          "LDAUU", "LDAUJ"]
PROD_INCAR = "Mn2In2Se5/data/dft_inputs/step4_disp/disp-00000/INCAR"
PROD_KPOINTS = "Mn2In2Se5/data/dft_inputs/step4_disp/disp-00000/KPOINTS"
DISP_YAML = "Mn2In2Se5/fc/phono3py_disp.yaml"


# ----------------------------------------------------------------------------- helpers
def incar_tags(text):
    out = {}
    for line in text.splitlines():
        line = line.split("#")[0].split("!")[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip().upper()] = " ".join(v.split()).upper()
    return out


def magmom_str(m):
    """VASP MAGMOM, runs of equal values compressed (e.g. 5 -5 -5 112*0)."""
    out, i = [], 0
    while i < len(m):
        j = i
        while j < len(m) and m[j] == m[i]:
            j += 1
        out.append(f"{j - i}*{m[i]:g}" if j - i > 1 else f"{m[i]:g}")
        i = j
    return " ".join(out)


def load_supercell(repo):
    ph = phonopy.load(str(Path(repo) / DISP_YAML), produce_fc=False, log_level=0)
    sc = ph.supercell
    sym = list(sc.symbols)
    assert len(sym) == 144 and sym[:32] == ["Mn"] * 32, "unexpected supercell"
    layer = np.array([int(ph.supercell.s2u_map[i] != ph.supercell.s2u_map[0]) for i in range(32)])
    assert set(layer.tolist()) == {0, 1} and layer.sum() == 16, layer
    return ph, sc, layer


def mn_shells(sc, layer, rmax=SHELL_RMAX):
    """Mn-Mn pairs (i < j, minimum image) grouped into shells by distance and inter/intra layer."""
    lat = np.array(sc.cell)
    f = sc.scaled_positions[:32]
    images = np.array([[a, b, c] for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)])
    pairs = []
    for i in range(32):
        for j in range(i + 1, 32):
            d = f[j] - f[i]
            d -= np.floor(d)
            vs = (d + images) @ lat
            rs = np.linalg.norm(vs, axis=1)
            o = np.argsort(rs)
            r = rs[o[0]]
            if r < rmax:
                # a pair must have a unique shortest image, else the shell is not resolved by the supercell
                assert rs[o[1]] - r > 0.05, f"Mn {i}-{j}: two images at {r:.3f} A"
                pairs.append((i, j, r, vs[o[0]]))
    shells = {}
    for i, j, r, v in pairs:
        key = ("inter" if layer[i] != layer[j] else "intra", round(r, 2))
        shells.setdefault(key, []).append((i, j, v))
    names = {}
    for kind, r in sorted(shells, key=lambda k: k[1]):
        names[f"{kind}_{r:.2f}"] = shells[(kind, r)]
    return names


def correlations(spins, shells):
    return {k: float(np.mean([spins[i] * spins[j] for i, j, _ in v])) for k, v in shells.items()}


def mn_permutations(ph):
    """Permutations of the 32 Mn sites under the space-group operations of the supercell."""
    sc = ph.supercell
    from phonopy.structure.symmetry import Symmetry
    s = Symmetry(sc, symprec=1e-4)
    ops = s.symmetry_operations
    f = sc.scaled_positions[:32]
    perms = []
    for R, t in zip(ops["rotations"], ops["translations"]):
        g = f @ R.T + t
        d = g[:, None, :] - f[None, :, :]
        d -= np.rint(d)
        m = np.argmin(np.linalg.norm(d, axis=2), axis=1)
        if len(set(m.tolist())) == 32:
            perms.append(m)
    return np.array(perms)


def canonical(spins, perms):
    keys = []
    for sgn in (1, -1):
        s = sgn * np.asarray(spins)
        for p in perms:
            t = np.empty(32, int)
            t[p] = s
            keys.append(tuple(t.tolist()))
    return min(keys)


def make_pm_configs(ph, layer, shells, n=N_PM, seed=SEED, restarts=400):
    """Simulated annealing (spin swaps inside a layer) towards zero correlation in every shell."""
    rng = np.random.default_rng(seed)
    perms = mn_permutations(ph)
    idx = [np.where(layer == l)[0] for l in (0, 1)]
    bonds = [(np.array([b[0] for b in v]), np.array([b[1] for b in v])) for v in shells.values()]

    def cost(s):
        return sum(float(np.mean(s[a] * s[b])) ** 2 for a, b in bonds)

    found, keys = [], set()
    for _ in range(restarts):
        s = np.zeros(32, int)
        for ii in idx:
            s[ii] = rng.permutation([1] * 8 + [-1] * 8)
        c = cost(s)
        T = 0.05
        for step in range(4000):
            ii = idx[rng.integers(2)]
            up, dn = rng.choice(ii[s[ii] == 1]), rng.choice(ii[s[ii] == -1])
            s[up], s[dn] = -1, 1
            c2 = cost(s)
            if c2 <= c or rng.random() < np.exp(-(c2 - c) / T):
                c = c2
            else:
                s[up], s[dn] = 1, -1
            T *= 0.999
        k = canonical(s, perms)
        if k not in keys:
            keys.add(k)
            found.append((c, s.copy()))
    found.sort(key=lambda x: x[0])
    return [s for _, s in found[:n]], perms


def heisenberg_per_mn(spins, shells, J):
    """sum_bonds J s_i s_j / 32 (meV per Mn) for the shells that have a J."""
    e = 0.0
    for k, v in shells.items():
        if k in J:
            e += J[k] * sum(spins[i] * spins[j] for i, j, _ in v)
    return e / 32.0


def random_frames(n, natom, rng):
    out = []
    for _ in range(n):
        v = rng.normal(size=(natom, 3))
        out.append(DIST * v / np.linalg.norm(v, axis=1)[:, None])
    return out


HERE = Path(__file__).resolve().parent
CONFIG_FILE = HERE / "spin_configs.json"     # committed; written once by --regenerate
FRAMES_FILE = HERE / "frames.npz"            # committed; written once by --regenerate


def generate_sets(repo):
    """Spin configurations and displacements of every set, from the fixed seeds."""
    ph, sc, layer = load_supercell(repo)
    shells = mn_shells(sc, layer)
    pm, _ = make_pm_configs(ph, layer, shells)
    rng = np.random.default_rng(SEED + 1)
    natom = len(sc.symbols)
    spins, frames = {}, {}
    f = random_frames(N_ORDERED, natom, rng)
    spins["FM_ctl"], frames["FM_ctl"] = np.ones(32, int), f + [-f[0]]
    spins["AFM"], frames["AFM"] = np.where(layer == 0, 1, -1), random_frames(N_ORDERED, natom, rng)
    for k, s in enumerate(pm, 1):
        spins[f"PM{k:02d}"], frames[f"PM{k:02d}"] = s, random_frames(N_PM_FRAMES, natom, rng)
    return spins, frames


def all_sets(repo):
    """The committed spin configurations and displacements, checked against the cell."""
    ph, sc, layer = load_supercell(repo)
    shells = mn_shells(sc, layer)
    perms = mn_permutations(ph)
    cfg = json.loads(CONFIG_FILE.read_text())
    fr = np.load(FRAMES_FILE)
    sets = {}
    for name, s in cfg["spins"].items():
        u = fr[name]
        assert u.shape[1:] == (len(sc.symbols), 3) and np.allclose(np.linalg.norm(u, axis=2), DIST)
        s = np.array(s, int)
        assert all(s[layer == l].sum() == 0 for l in (0, 1)) or name in ("FM_ctl", "AFM")
        sets[name] = {"spins": s, "frames": list(u)}
    sets["FM_ctl"]["frame_note"] = f"frame {N_ORDERED + 1} = -(frame 1): even-order (anharmonic) check"
    assert np.allclose(sets["FM_ctl"]["frames"][N_ORDERED], -sets["FM_ctl"]["frames"][0])
    keys = {canonical(s["spins"], perms) for k, s in sets.items() if k.startswith("PM")}
    assert len(keys) == sum(k.startswith("PM") for k in sets), "two PM configurations are equivalent"
    return ph, sc, layer, shells, perms, sets


def write_run(d, cell, positions, spins, tag, istart, lwave, potcar):
    d.mkdir(parents=True, exist_ok=False)
    atoms = PhonopyAtoms(symbols=list(cell.symbols), cell=cell.cell, positions=positions)
    write_vasp(str(d / "POSCAR"), atoms)
    mag = [5 * int(x) for x in spins] + [0] * (len(cell.symbols) - 32)
    (d / "INCAR").write_text(INCAR.format(tag=tag, istart=istart, icharg=0 if istart else 2,
                                          lwave=lwave, magmom=magmom_str(mag)))
    (d / "KPOINTS").write_text(KPOINTS)
    if potcar:
        shutil.copy(potcar, d / "POTCAR")


def config_table(shells, sets):
    names = list(shells)
    print("shell (bonds in supercell):  " + "  ".join(f"{k}({len(shells[k])})" for k in names))
    for name, s in sets.items():
        c = correlations(s["spins"], shells)
        print(f"  {name:7s} <s s> = " + "  ".join(f"{c[k]:+.3f}" for k in names)
              + f"   E_J - E_FM = {heisenberg_per_mn(s['spins'], shells, J_GAMMA_FIT) - heisenberg_per_mn(np.ones(32), shells, J_GAMMA_FIT):+.2f} meV/Mn")
    pm = [s["spins"] for k, s in sets.items() if k.startswith("PM")]
    c = {k: np.mean([correlations(s, shells)[k] for s in pm]) for k in names}
    print("  PM mean <s s> = " + "  ".join(f"{c[k]:+.3f}" for k in names))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="clone of Mn2In2Se5-thermal")
    ap.add_argument("--out", help="new work directory (must not exist)")
    ap.add_argument("--potcar", help="production 4x4x1 S4 POTCAR (species Mn In Se); copied into every run")
    ap.add_argument("--list-configs", action="store_true")
    ap.add_argument("--regenerate", action="store_true",
                    help="rewrite spin_configs.json and frames.npz from the seeds (done once, committed)")
    a = ap.parse_args()
    repo = Path(a.repo)
    if a.regenerate:
        spins, frames = generate_sets(repo)
        CONFIG_FILE.write_text(json.dumps({"seed": SEED, "note": "8 up / 8 down per Mn layer for PM; "
                                           "Mn order = supercell atoms 1-32 of Mn2In2Se5/fc/SPOSCAR",
                                           "spins": {k: [int(x) for x in v] for k, v in spins.items()}}, indent=1))
        np.savez_compressed(FRAMES_FILE, **{k: np.array(v) for k, v in frames.items()})
        print(f"wrote {CONFIG_FILE.name}, {FRAMES_FILE.name}")
        return

    ph, sc, layer, shells, perms, sets = all_sets(repo)
    if a.list_configs:
        print(f"{len(perms)} space-group permutations of the Mn sites")
        config_table(shells, sets)
        return
    if not a.out:
        raise SystemExit("--out is required")
    out = Path(a.out)
    if out.exists():
        raise SystemExit(f"{out} exists; choose a new directory")

    prod = incar_tags((repo / PROD_INCAR).read_text())
    mine = incar_tags(INCAR.format(tag="x", istart=0, icharg=2, lwave=".FALSE.", magmom="1"))
    bad = [k for k in SHARED if prod.get(k) != mine.get(k)]
    if bad:
        raise SystemExit("INCAR differs from production S4 in " + ", ".join(
            f"{k}: {prod.get(k)} vs {mine.get(k)}" for k in bad))
    kp = (repo / PROD_KPOINTS).read_text().split("\n")
    if kp[2].strip().lower()[0] != "g" or kp[3].split() != ["2", "2", "2"]:
        raise SystemExit(f"production KPOINTS is not Gamma 2 2 2: {kp[2:4]}")
    if a.potcar:
        titles = re.findall(r"TITEL\s*=\s*\S+\s+(\S+)", Path(a.potcar).read_text())
        if [t.split("_")[0] for t in titles] != ["Mn", "In", "Se"]:
            raise SystemExit(f"POTCAR species {titles}, expected Mn, In, Se in that order")

    pos0 = sc.scaled_positions @ np.array(sc.cell)
    manifest = {"distance_A": DIST, "seed": SEED, "kmesh": [2, 2, 2], "supercell_matrix": [[4, 0, 0], [0, 4, 0], [0, 0, 1]],
                "shells": {k: len(v) for k, v in shells.items()}, "J_gamma_fit_meV": J_GAMMA_FIT, "sets": {}}
    for name, s in sets.items():
        base = out / name
        write_run(base / "ref", sc, pos0, s["spins"], f"{name} ref", 0, ".TRUE.", a.potcar)
        dirs = []
        for k, u in enumerate(s["frames"], 1):
            d = base / f"disp-{k:02d}"
            write_run(d, sc, pos0 + u, s["spins"], f"{name} disp-{k:02d}", 1, ".FALSE.", a.potcar)
            dirs.append(d.name)
        (base / "displacements.json").write_text(json.dumps(
            {"set": name, "spins": [int(x) for x in s["spins"]], "dirs": dirs,
             "displacements_A": [np.round(u, 10).tolist() for u in s["frames"]],
             "note": s.get("frame_note", "")}))
        manifest["sets"][name] = {"magmom": magmom_str([5 * int(x) for x in s["spins"]] + [0] * 112),
                                  "n_disp": len(dirs),
                                  "correlations": correlations(s["spins"], shells),
                                  "E_J_minus_FM_meV_per_Mn": heisenberg_per_mn(s["spins"], shells, J_GAMMA_FIT)
                                  - heisenberg_per_mn(np.ones(32), shells, J_GAMMA_FIT)}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    n = sum(v["n_disp"] + 1 for v in manifest["sets"].values())
    print(f"wrote {out}: {n} VASP runs of 144 atoms in {len(sets)} sets")
    config_table(shells, sets)
    if not a.potcar:
        print("no --potcar given: copy the production POTCAR (Mn In Se) into every run directory")


if __name__ == "__main__":
    main()
