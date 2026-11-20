#!/usr/bin/env python3
"""Predict supercell forces from fitted FC2 + FC3 and compute pheasy-style errors.

F_i = -sum_j Phi2[i,j] u_j - 1/2 sum_jk Phi3[i,j,k] u_j u_k

Used by the FC3-cutoff CV study (protocol/local_AI_next_steps_2026-11-20.md):
the hold-out metric eps_anh needs predicted forces on frames that were not in the
fit, and the mandatory sanity check needs them to reproduce pheasy's own training
relative error.

pheasy's definition (pheasy-gpu-gh/core/optimizer.py:9160):
    re = ||F_pred - F_dft|| / ||F_dft||
over all training frames and components, with the equilibrium frame not included.
"""
import numpy as np


def predict_forces(disps, fc2, fc3):
    """disps (n, nat, 3), fc2 (nat,nat,3,3), fc3 (nat,nat,nat,3,3,3) -> (n, nat, 3)."""
    n, nat = disps.shape[0], disps.shape[1]
    # (i, a) x (j, k, b, c) so that one matmul per frame gives the FC3 term
    m3 = fc3.transpose(0, 3, 1, 2, 4, 5).reshape(nat * 3, nat * nat * 9)
    out = np.empty_like(disps)
    for s in range(n):
        u = disps[s]
        f = -np.einsum('ijab,jb->ia', fc2, u)
        v = np.einsum('jb,kc->jkbc', u, u).ravel()
        out[s] = f - 0.5 * (m3 @ v).reshape(nat, 3)
    return out


def rel_error(pred, dft):
    return float(np.linalg.norm(pred - dft) / np.linalg.norm(dft))


if __name__ == '__main__':
    import sys, time
    from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5
    if len(sys.argv) < 4:
        sys.exit('usage: predict_forces.py <disps.npy> <forces.npy> <fc2.hdf5> <fc3.hdf5> [n_skip]')
    disps = np.load(sys.argv[1])
    forces = np.load(sys.argv[2])
    fc2 = read_fc2_from_hdf5(sys.argv[3])
    fc3 = read_fc3_from_hdf5(sys.argv[4])
    skip = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    t0 = time.time()
    pred = predict_forces(disps[skip:], fc2, fc3)
    print('predicted %d frames in %.1f s' % (len(pred), time.time() - t0))
    print('re = %.6f = %.4f %%' % (rel_error(pred, forces[skip:]), 100 * rel_error(pred, forces[skip:])))