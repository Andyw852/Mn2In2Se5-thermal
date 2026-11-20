#!/usr/bin/env python3
"""Hold-out errors for the FC3-cutoff CV runs (protocol step B).

Reports, on the training frames and on the held-out frames:
  * total-force relative error, pheasy definition ||Fp - Fd|| / ||Fd||;
  * eps_anh = ||Fp - Fd|| / ||Fd - Fh,prod||, where Fh,prod is the harmonic force
    from the PRODUCTION FC2 (md5 b611f139...), the same denominator for every fit,
    because the total error is dominated by FC2 and is insensitive to FC3.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from predict_forces import predict_forces, rel_error
from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5


def harm_forces(disps, fc2):
    return np.array([-np.einsum('ijab,jb->ia', fc2, u) for u in disps])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fitdir', required=True)
    ap.add_argument('--disps', required=True)
    ap.add_argument('--forces', required=True)
    ap.add_argument('--train', required=True)
    ap.add_argument('--held', required=True)
    ap.add_argument('--prod-fc2', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--tag', default='')
    a = ap.parse_args()
    D = np.load(a.disps)
    F = np.load(a.forces)
    fc2 = read_fc2_from_hdf5(os.path.join(a.fitdir, 'fc2.hdf5'))
    fc3 = read_fc3_from_hdf5(os.path.join(a.fitdir, 'fc3.hdf5'))
    pf2 = read_fc2_from_hdf5(a.prod_fc2)
    def ids(s):
        return [int(x) - 1 for x in s.split(',') if x.strip()]
    res = {'tag': a.tag, 'fitdir': a.fitdir, 'sets': {}}
    for name, sel in (('train', ids(a.train)), ('held', ids(a.held))):
        if not sel:
            continue
        pred = predict_forces(D[sel], fc2, fc3)
        dft = F[sel]
        fh = harm_forces(D[sel], pf2)
        res['sets'][name] = {
            'n_frames': len(sel),
            'frame_ids': [i + 1 for i in sel],
            're_total': rel_error(pred, dft),
            'eps_anh': float(np.linalg.norm(pred - dft) / np.linalg.norm(dft - fh)),
            'rmse_total_eV_per_A': float(np.sqrt(((pred - dft) ** 2).mean())),
        }
    json.dump(res, open(a.out, 'w'), indent=1)
    print(json.dumps(res['sets'], indent=1))


if __name__ == '__main__':
    main()
