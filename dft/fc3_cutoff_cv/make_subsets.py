#!/usr/bin/env python3
"""Build the subset datasets and per-run directories for the hold-out /
learning-curve test of the Mn2In2Se5 FC3 cutoff (protocol 2026-11-20).

Subset method: a directory of symlinks to the original vasprun frame folders
(disp-00000 ... disp-00054) plus the equilibrium frame, which is in every
training set. The driver re-reads that directory during `prep`, so no pheasy
internals are touched.
"""
import json
import os
from pathlib import Path

SRC = Path('/data/wangchaoyue852/taskflow/work/Mn2In2Se5_441/kl-dft-cpu/step4_disp')
X = Path('/data/wangchaoyue852/taskflow/work/Mn2In2Se5_441/kl-dft-cpu/step5_fc/cutoff_cv')
NF = 54                      # displaced frames (disp-00001 .. disp-00054)
SEED_FOLDS = 20261120
SEED_N36 = (20261122, 20261123)
CUTOFFS_HOLDOUT = [6.51, 6.90, 7.11, 7.43]
CUTOFFS_N36 = [6.51, 7.11, 7.43]
BASE_CFG = Path('/data/wangchaoyue852/taskflow/work/Mn2In2Se5_441/kl-dft-cpu/step5_fc/cutoff_scan_ext/fit_config.json')


def frame_dir(i):
    """dataset index i (0 = equilibrium) -> disp-XXXXX directory name"""
    return 'disp-%05d' % i


def link_subset(name, idxs):
    d = X / 'subsets' / name
    d.mkdir(parents=True, exist_ok=True)
    for f in d.iterdir():
        f.unlink()
    for i in [0] + sorted(idxs):
        os.symlink(SRC / frame_dir(i), d / frame_dir(i))
    return d


def make_run(name, subset_dir, cutoff):
    r = X / 'runs' / name
    r.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(BASE_CFG.read_text())
    cfg['dataset_dir'] = str(subset_dir)
    cfg['dataset_dir_abs'] = str(subset_dir)
    cfg['cut3_candidates'] = [cutoff]
    cfg['cut3_safe_cutoff'] = 8.0
    cfg['cut3_scan'] = True
    cfg['_cv_note'] = 'cutoff_cv run %s: subset %s, cutoff %.2f A' % (name, subset_dir.name, cutoff)
    (r / 'fit_config.json').write_text(json.dumps(cfg, indent=1))
    for f in ('POSCAR', 'SPOSCAR'):
        p = r / f
        if p.is_symlink() or p.exists():
            p.unlink()
        os.symlink(SRC / f, p)
    return r


def main():
    rng = __import__('numpy').random.default_rng(SEED_FOLDS)
    perm = rng.permutation(NF) + 1          # dataset indices 1..54
    folds = [perm[0:10], perm[10:20], perm[20:30]]
    manifest = {'seeds': {'folds': SEED_FOLDS, 'n36': list(SEED_N36)},
                'nframes_displaced': NF,
                'note': 'dataset index 0 is the equilibrium frame, present in every training set',
                'folds': []}
    for k, held in enumerate(folds):
        held = sorted(int(v) for v in held)
        train = [i for i in range(1, NF + 1) if i not in held]
        sub = link_subset('fold%d' % k, train)
        entry = {'fold': k, 'held_out': held, 'train': train, 'subset': str(sub)}
        manifest['folds'].append(entry)
        for c in CUTOFFS_HOLDOUT:
            make_run('c%04d_f%d' % (round(c * 100), k), sub, c)
    manifest['n36'] = []
    for si, seed in enumerate(SEED_N36):
        rng2 = __import__('numpy').random.default_rng(seed)
        pick = sorted(int(v) for v in (rng2.permutation(NF)[:36] + 1))
        sub = link_subset('n36_%d' % si, pick)
        manifest['n36'].append({'subset': si, 'seed': seed, 'frames': pick, 'path': str(sub)})
        for c in CUTOFFS_N36:
            make_run('c%04d_n36_%d' % (round(c * 100), si), sub, c)
    (X / 'subsets_manifest.json').write_text(json.dumps(manifest, indent=1))
    print('runs:', len(list((X / 'runs').iterdir())))
    print('subsets:', len(list((X / 'subsets').iterdir())))
    print('manifest ->', X / 'subsets_manifest.json')


if __name__ == '__main__':
    main()
