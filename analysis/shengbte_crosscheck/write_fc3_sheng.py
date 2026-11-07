"""从 phono3py 的 fc3.hdf5（全超胞）写 ShengBTE 的 FORCE_CONSTANTS_3RD：
每个三元组的 2、3 号原子都取相对 1 号原子的最短周期像（截断 < 超胞安全截断时唯一）；
R 向量 = 像的位置 − CONTROL 里该原子的位置（单胞分数坐标 wrap 到 [0,1)）。
同时写 CONTROL（与 AutoZT 同口径）。用法：python write_fc3_sheng.py <fcdir> <outdir> n1 n2 n3 T"""
import sys, itertools, numpy as np, h5py, phono3py
fd, od = sys.argv[1], sys.argv[2]
ngrid = [int(x) for x in sys.argv[3:6]]; T = float(sys.argv[6])
p = phono3py.load(fd + "/phono3py_disp.yaml", produce_fc=False, log_level=0)
uc, sc, prim = p.unitcell, p.supercell, p.primitive
assert len(prim) == len(uc)
Lp, Ls = np.asarray(uc.cell), np.asarray(sc.cell)
fu = uc.scaled_positions % 1.0
xu = fu @ Lp                                    # CONTROL 的笛卡尔位置
rs = sc.scaled_positions @ Ls
s2p = np.array([list(prim.p2s_map).index(x) for x in prim.s2p_map])   # 超胞原子 -> 单胞原子序号
fc3 = h5py.File(fd + "/fc3.hdf5")["fc3"]
shifts = np.array(list(itertools.product(range(-2, 3), repeat=3))) @ Ls
blocks = []
maxd = 0.0
for ia, a in enumerate(prim.p2s_map):
    F = fc3[a]                                  # (144,144,3,3,3)
    nz = np.argwhere(np.abs(F).reshape(len(sc), len(sc), 27).max(axis=2) > 1e-10)
    def img(s):
        cand = rs[s] - rs[a] + shifts
        dn = np.linalg.norm(cand, axis=1)
        k = np.argsort(dn)
        if dn[k[1]] - dn[k[0]] < 1e-4:
            raise SystemExit("degenerate image for %d-%d at %.3f A" % (a, s, dn[k[0]]))
        return cand[k[0]]
    cache = {}
    for s2, s3 in nz:
        for s in (s2, s3):
            if s not in cache:
                cache[s] = img(s)
        d2, d3 = cache[s2], cache[s3]
        maxd = max(maxd, np.linalg.norm(d2), np.linalg.norm(d3), np.linalg.norm(d3 - d2))
        j, k = s2p[s2], s2p[s3]
        R2 = xu[ia] + d2 - xu[j]; R3 = xu[ia] + d3 - xu[k]
        for R in (R2, R3):
            n = R @ np.linalg.inv(Lp)
            assert np.allclose(n, np.rint(n), atol=1e-5), n
        blocks.append((R2, R3, ia, j, k, F[s2, s3]))
with open(od + "/FORCE_CONSTANTS_3RD", "w") as f:
    f.write("%d\n" % len(blocks))
    for n, (R2, R3, i, j, k, v) in enumerate(blocks, 1):
        f.write("\n%d\n%.10f %.10f %.10f\n%.10f %.10f %.10f\n%d %d %d\n" % ((n,) + tuple(R2) + tuple(R3) + (i + 1, j + 1, k + 1)))
        for a, b, c in itertools.product(range(3), repeat=3):
            f.write("%d %d %d %.12e\n" % (a + 1, b + 1, c + 1, v[a, b, c]))
els = list(dict.fromkeys(uc.symbols))
L = ["&allocations", "        nelements=%d," % len(els), "        natoms=%d," % len(uc),
     "        ngrid(:)=%d %d %d" % tuple(ngrid), "&end", "&crystal", "        lfactor=0.1,"]
for i, v in enumerate(Lp):
    L.append("        lattvec(:,%d)=%.10f %.10f %.10f," % ((i + 1,) + tuple(v)))
L.append("        elements=%s" % " ".join('"%s"' % e for e in els))
L.append("        types=%s," % " ".join(str(els.index(s) + 1) for s in uc.symbols))
for i, x in enumerate(fu):
    L.append("        positions(:,%d)=%.10f %.10f %.10f," % ((i + 1,) + tuple(x)))
dims = np.rint(Ls @ np.linalg.inv(Lp)).astype(int).diagonal()
L.append("        scell(:)=%d %d %d" % tuple(dims))
L += ["&end", "&parameters", "        T=%g" % T, "        scalebroad=1.0", "&end", "&flags",
      "        convergence=.FALSE.,", "        isotopes=.TRUE.,", "        nonanalytic=.FALSE.,",
      "        nanowires=.FALSE.", "&end", ""]
open(od + "/CONTROL", "w").write("\n".join(L))
print("blocks %d, max pair distance %.3f A, scell %s" % (len(blocks), maxd, dims.tolist()))
