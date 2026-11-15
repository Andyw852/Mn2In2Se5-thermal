#!/usr/bin/env python3
"""从 phono3py fc2.hdf5 写 ShengBTE FORCE_CONSTANTS_2ND（与修正后的 CONTROL 位置一致）。
用法: python write_fc2_sheng.py <fcdir> <outdir>
fcdir 需含 phono3py_disp.yaml + fc2.hdf5
"""
import sys, itertools, numpy as np, h5py, phono3py

fd, od = sys.argv[1], sys.argv[2]
p = phono3py.load(fd + "/phono3py_disp.yaml", produce_fc=False, log_level=0)
uc, sc, prim = p.unitcell, p.supercell, p.primitive
Lp, Ls = np.asarray(uc.cell), np.asarray(sc.cell)
fu = uc.scaled_positions % 1.0
# ★ 修正：把接近 1 的坐标写成 1（而不是 0），与修正后的 CONTROL 一致
fu = np.where(np.abs(fu - 1.0) < 1e-8, 1.0, fu)
xu = fu @ Lp
rs = sc.scaled_positions @ Ls
s2p = np.array([list(prim.p2s_map).index(x) for x in prim.s2p_map])
fc2 = h5py.File(fd + "/fc2.hdf5")["force_constants"] if "force_constants" in h5py.File(fd+"/fc2.hdf5") else h5py.File(fd+"/fc2.hdf5")["fc2"]
n = len(sc)
shifts = np.array(list(itertools.product(range(-2, 3), repeat=3))) @ Ls

lines = ["%d %d" % (n, n)]
maxd = 0.0
for ia, a in enumerate(prim.p2s_map):
    F = fc2[a]
    nz = np.argwhere(np.abs(F.reshape(n, 9).max(axis=1)) > 1e-10)
    blocks = []
    for s in nz[:, 0]:
        cand = rs[s] - rs[a] + shifts
        dn = np.linalg.norm(cand, axis=1)
        k = np.argsort(dn)
        if dn[k[1]] - dn[k[0]] < 1e-4:
            raise SystemExit("degenerate image for %d-%d at %.3f A" % (a, s, dn[k[0]]))
        d = cand[k[0]]
        maxd = max(maxd, np.linalg.norm(d))
        j = s2p[s]
        R = xu[ia] + d - xu[j]
        blocks.append((j + 1, R, np.asarray(F[s])))
    lines.append(str(len(blocks)))
    for j1, R, m in blocks:
        lines.append("%d %d" % (a + 1, j1))
        lines.append(" %18.10f %18.10f %18.10f" % tuple(R))
        for r in range(3):
            lines.append(" %18.10f %18.10f %18.10f" % tuple(m[r]))
open(od + "/FORCE_CONSTANTS_2ND", "w").write("\n".join(lines) + "\n")
print("blocks 元组数 %d, max pair distance %.3f A" % (len(lines), maxd))
