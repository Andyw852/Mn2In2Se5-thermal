#!/usr/bin/env python3
"""Hold-out errors for every finished CV fit (protocol step B).

Reads the split manifest and the fitted fc2/fc3 of each run under fits/, computes
the total-force relative error and eps_anh on the training and the held-out frames,
and writes holdout.json.
"""
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from predict_forces import predict_forces, rel_error
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5

B = sys.argv[1] if len(sys.argv) > 1 else '.'
OUT = os.path.join(B, 'holdout.json')

man = json.load(open(os.path.join(B, 'subsets_manifest_arrays.json')))
D = np.load(os.path.join(B, 'dataset_disps.npy'))
F = np.load(os.path.join(B, 'dataset_forces.npy'))
pf2 = read_fc2_from_hdf5(os.path.join(B, 'prod_fc2_b611f139.hdf5'))

folds = {f['fold']: f for f in man['folds']}
n36 = {n['subset']: n for n in man['n36']}


def harm(disps):
    return np.array([-np.einsum('ijab,jb->ia', pf2, u) for u in disps])


def errs(fitdir, train_ids, held_ids):
    fc2 = read_fc2_from_hdf5(os.path.join(fitdir, 'fc2.hdf5'))
    fc3 = read_fc3_from_hdf5(os.path.join(fitdir, 'fc3.hdf5'))
    res = {}
    for name, ids in (('train', train_ids), ('held', held_ids)):
        if not ids:
            continue
        sel = [i - 1 for i in ids]
        pred = predict_forces(D[sel], fc2, fc3)
        dft = F[sel]
        fh = harm(D[sel])
        res[name] = {
            'n_frames': len(sel),
            'frame_ids': ids,
            're_total': rel_error(pred, dft),
            'eps_anh': float(np.linalg.norm(pred - dft) / np.linalg.norm(dft - fh)),
        }
    return res


out = {'note': 'eps_anh = ||Fp-Fd|| / ||Fd-Fh,prod||, production FC2 b611f139 as the harmonic denominator', 'fits': {}}
fitroot = os.path.join(B, 'fits')
for name in sorted(os.listdir(fitroot)) if os.path.isdir(fitroot) else []:
    m = re.match(r'c(\d{4})_(f(\d)|n36_(\d))$', name)
    if not m:
        continue
    cutoff = int(m.group(1)) / 100.0
    fitdir = os.path.join(fitroot, name)
    if not (os.path.exists(os.path.join(fitdir, 'fc3.hdf5'))
            and os.path.exists(os.path.join(fitdir, 'fc2.hdf5'))):
        continue
    if m.group(2).startswith('f'):
        fold = folds[int(m.group(3))]
        train, held = fold['train'], fold['held_out']
        kind = 'fold%d' % fold['fold']
    else:
        sub = n36[int(m.group(4))]
        train = sub['frames']
        held = [i for i in range(1, 55) if i not in train]
        kind = 'n36_%d' % sub['subset']
    r = errs(fitdir, train, held)
    r.update({'cutoff_A': cutoff, 'split': kind})
    out['fits'][name] = r
    print(name, 'cut %.2f' % cutoff, kind,
          're_held %.5f eps_anh_held %.5f re_train %.5f' % (
              r['held']['re_total'], r['held']['eps_anh'], r['train']['re_total']), flush=True)
json.dump(out, open(OUT, 'w'), indent=1)
print('wrote', OUT, len(out['fits']), 'fits')
