#!/usr/bin/env python3
"""Rebuild the cutoff_cv subsets on A800 from the authoritative 3090 split.

Dataset array convention (checked): indices 0..53 are the displaced frames, index
54 is the equilibrium frame (zero displacement, zero force). The 3090 manifest
names frames by their disp-XXXXX id, i.e. 1..54, so frame id f -> array index f-1.
"""
import json
from pathlib import Path
import numpy as np

B = Path('/fs0/home/wangcch/work/mn_cv')
CB = Path('/data/wangchaoyue852/taskflow/work/Mn2In2Se5_441/kl-dft-cpu/step5_fc/cutoff_cv')
MAN = json.loads((B / 'subsets_manifest_3090.json').read_text())
CUT_HOLDOUT = [6.51, 6.90, 7.11, 7.43]
CUT_N36 = [6.51, 7.11, 7.43]
D = np.load(B / 'full/dataset_disps.npy')
F = np.load(B / 'full/dataset_forces.npy')
NEQ = len(D) - 1
rms = np.sqrt((D ** 2).sum(axis=2).mean(axis=1))
assert rms[NEQ] == 0.0, 'last frame must be the equilibrium'


def make_subset(name, frame_ids):
    d = B / 'subsets' / name
    d.mkdir(parents=True, exist_ok=True)
    keep = sorted(int(f) - 1 for f in frame_ids) + [NEQ]
    np.save(d / 'dataset_disps.npy', D[keep])
    np.save(d / 'dataset_forces.npy', F[keep])
    return d, keep


BASE = json.loads((B / 'fit_config_template.json').read_text())


def make_run(name, sub, cut):
    r = B / 'runs' / name
    r.mkdir(parents=True, exist_ok=True)
    cfg = dict(BASE)
    cfg['dataset_dir'] = str(sub)
    cfg['dataset_dir_abs'] = str(sub)
    cfg['dataset_signature'] = 'arrays'
    cfg['cut3_candidates'] = [cut]
    cfg['cut3_safe_cutoff'] = 8.0
    cfg['cut3_scan'] = True
    (r / 'fit_config.json').write_text(json.dumps(cfg, indent=1))
    for f in ('POSCAR', 'SPOSCAR'):
        p = r / f
        if p.exists() or p.is_symlink():
            p.unlink()
        p.symlink_to(B / f)
    return r


man = {'source': 'subsets_manifest.json from the 3090 cutoff_cv (authoritative)',
       'seeds': MAN['seeds'], 'array_convention': '0..53 displaced, 54 equilibrium', 'folds': [], 'n36': []}
for f in MAN['folds']:
    held = f['held_out']
    train = [i for i in range(1, 55) if i not in held]
    sub, keep = make_subset('fold%d' % f['fold'], train)
    man['folds'].append({'fold': f['fold'], 'held_out': held, 'train': train, 'array_indices': keep})
    for c in CUT_HOLDOUT:
        make_run('c%04d_f%d' % (round(c * 100), f['fold']), sub, c)
for n in MAN['n36']:
    sub, keep = make_subset('n36_%d' % n['subset'], n['frames'])
    man['n36'].append({'subset': n['subset'], 'seed': n['seed'], 'frames': n['frames'], 'array_indices': keep})
    for c in CUT_N36:
        make_run('c%04d_n36_%d' % (round(c * 100), n['subset']), sub, c)
(B / 'subsets_manifest.json').write_text(json.dumps(man, indent=1))
print('subsets rebuilt; runs:', len(list((B / 'runs').iterdir())))
for n in ('fold0', 'n36_0'):
    d = np.load(B / 'subsets' / n / 'dataset_disps.npy')
    r = np.sqrt((d ** 2).sum(axis=2).mean(axis=1))
    print(n, 'frames', len(d), 'last-frame rms', r[-1], 'training frames', len(d) - 1)
