#!/usr/bin/env python3
"""Build R_xx(T) = kappa_xx(PM) / kappa_xx(FM_ctl) with a jackknife standard error.

The error uses the delete-one jackknife estimator of the standard error,
    SE = sqrt((n-1)/n * sum_i (R_i - mean(R))^2),
which is the convention used in kappa_summary.json (a plain ddof=1 std would be
sqrt(n/(n-1)) times smaller and would have understated it by ~2.5x for n = 8).
"""
import json
import sys
from pathlib import Path

import numpy as np

B = Path(sys.argv[1] if len(sys.argv) > 1 else '.')
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else 'R_xx_vs_T.json')


def kxx(base):
    d = np.load(B / ('bases_' + base) / 'bases.npz')
    T = np.asarray(d['T'], float)
    k = np.asarray(d['kappa0'], float)          # (nT, 6)
    return T, 0.5 * (k[:, 0] + k[:, 1])


T, fm = kxx('FM_ctl')
Tp, pm = kxx('PM')
assert np.allclose(T, Tp)
R = pm / fm
jk = []
for i in range(1, 9):
    p = B / ('bases_PM_jk%02d' % i) / 'bases.npz'
    if not p.is_file():
        continue
    d = np.load(p)
    jk.append(0.5 * (np.asarray(d['kappa0'])[:, 0] + np.asarray(d['kappa0'])[:, 1]) / fm)
jk = np.array(jk)
n = len(jk)
err = np.sqrt((n - 1) / n * ((jk - jk.mean(axis=0)) ** 2).sum(axis=0)) if n > 1 else np.zeros_like(R)
out = {'what': 'kappa_xx(PM) / kappa_xx(FM_ctl) for Mn2In2Se5, 16x16x4, production fc3',
       'source': 'analysis/extrinsic/kappa_extrinsic.py bases on the committed FC2 files',
       'mesh': [16, 16, 4],
       'T_K': [float(t) for t in T],
       'R_xx': [float(v) for v in R],
       'R_xx_err': [float(v) for v in err],
       'kappa_xx_FM_ctl': [float(v) for v in fm],
       'kappa_xx_PM': [float(v) for v in pm],
       'jackknife_sets': n,
       'error_estimator': 'delete-one jackknife SE = sqrt((n-1)/n * sum (R_i - mean R)^2)',
       'note': 'R_xx is an RTA ratio. Applying it to the quoted LBTE + coherence value assumes the same ratio holds there; that assumption is not verified directly.'}
OUT.write_text(json.dumps(out, indent=1) + '\n')
print('wrote', OUT)
for t, r, e in zip(T, R, err):
    print('  %6.1f K  R = %.4f +- %.4f' % (t, r, e))
