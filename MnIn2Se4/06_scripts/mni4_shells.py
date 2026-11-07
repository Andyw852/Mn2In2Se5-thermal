
import numpy as np, itertools
from collections import Counter

def read_poscar(p):
    L=open(p).read().splitlines()
    sc=float(L[1]); lat=np.array([[float(x) for x in L[i].split()] for i in range(2,5)])*sc
    sp=L[5].split(); cnt=[int(x) for x in L[6].split()]; n=sum(cnt)
    pos=np.array([[float(x) for x in L[8+i].split()[:3]] for i in range(n)])
    species=[s for s,c in zip(sp,cnt) for _ in range(c)]
    return lat,pos,species

CUTS=[2.85,3.70,3.95,4.30,4.70,5.00,5.20,5.50,5.80,6.10,6.50,7.00]
lat,pos,species=read_poscar('/home/wangchao/.dsh/.dsh-drops/SPOSCAR')
n=len(pos); cart=pos@lat
shifts=np.array(list(itertools.product([-1,0,1],repeat=3)))
dists=[]
for i in range(n):
    dr=cart[i]-cart; best=None
    for s in shifts:
        c=np.linalg.norm(dr+s@lat,axis=1); best=c if best is None else np.minimum(best,c)
    for j in range(i+1,n):
        dists.append((best[j],species[i],species[j]))
dists.sort()
tol=1e-3; shells=[]; cur=[dists[0]]
for d in dists[1:]:
    if d[0]-cur[-1][0]<tol: cur.append(d)
    else: shells.append(cur); cur=[d]
shells.append(cur)

lines=[]
lines.append("# MnIn2Se4 超胞 (3x3x3, 189 原子) 原子间距壳层与截断选择\n")
lines.append("超胞晶格: |a|=|b|=|c|=40.209 A, alpha=17.463 deg (菱方原胞; 六方描述 a=4.069 A, c=39.586 A)。")
lines.append("最小镜像上限 (supercell/2) = 20.10 A。\n")
lines.append("## 壳层表 (按距离聚类, tol=1e-3 A)\n")
lines.append("| shell | d (A) | 多重度 | 原子对组成 | 到下一壳层间隙 (A) |")
lines.append("|---|---|---|---|---|")
prev=None
for k,sh in enumerate(shells):
    d=sh[0][0]; comp=Counter(tuple(sorted((a,b))) for _,a,b in sh)
    gap = "" if prev is None else "%.3f"%(d-prev)
    comps=", ".join("%s-%s x%d"%(k2[0],k2[1],v) for k2,v in sorted(comp.items()))
    mark=""
    for c in CUTS:
        if prev is not None and prev < c < d:
            mark=" **<-c3=%.2f**"%c
    lines.append("| %d | %.4f | %d | %s | %s%s |"%(k,d,len(sh),comps,gap,mark))
    prev=d
lines.append("\n## 选定截断 (全部落在壳层间隙中)\n")
lines.append("| c3 (A) | 落在间隙 (A) | 间隙宽度 (A) | 说明 |")
lines.append("|---|---|---|---|")
prev=None
shell_d=[s[0][0] for s in shells]
for c in CUTS:
    lo=max([d for d in shell_d if d<c], default=0)
    hi=min([d for d in shell_d if d>c], default=20.10)
    lines.append("| %.2f | %.4f - %.4f | %.3f | |"%(c,lo,hi,hi-lo))

open('/home/wangchao/software/taskflow-v2.0/tmp/mn2inse4/壳层与截断.md','w').write("\n".join(lines)+"\n")
print("\n".join(lines[:40]))
print("...")
print("\n".join(lines[-20:]))
