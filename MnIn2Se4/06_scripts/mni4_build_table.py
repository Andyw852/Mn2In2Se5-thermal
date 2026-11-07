
import os, re, glob, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"results")
CUTS=["2.85","3.70","3.95","4.30","4.70","5.00","5.20","5.50","5.80","6.10","6.50","7.00"]
METHODS=["OLS","LASSO","ALASSO"]

def rd(p):
    try: return open(p,errors="ignore").read()
    except OSError: return ""

def fit_row(m,c3):
    t=rd(os.path.join(RES,"run_%s_c3_%s.log"%(m,c3)))
    if not t: return {}
    def num(pat):
        mm=re.search(pat,t); return float(mm.group(1)) if mm else None
    d=dict(p=num(r"Free IFC terms:\s*(\d+)"),
           rmse_in=num(r"RMSE:\s*([0-9.eE+-]+) eV/A"),
           rmse_cv=num(r"RMSE_CV:\s*([0-9.eE+-]+) eV/A"),
           rel=num(r"Relative error:\s*([0-9.eE+-]+)"),
           nnz=num(r"Non-zero IFC terms:\s*(\d+)"),
           alpha=num(r"alpha_opt:\s*([0-9.eE+-]+)"))
    mm=re.search(r"fc3\s+max\s*([0-9.eE+-]+)",t); d["fc3max"]=float(mm.group(1)) if mm else None
    return d

kmap={}
for ln in rd(os.path.join(RES,"kappa_all.txt")).splitlines():
    parts=ln.split()
    if len(parts)>=5 and parts[0].startswith(("ghsc_","symfc_")):
        kmap[parts[0]]=tuple(float(x) for x in parts[2:5])

def km(m,c3):
    return kmap.get("ghsc_%s_c3_%s"%(m,c3))
def ks(c3):
    return kmap.get("symfc_c3_%s"%c3)

def holdout(c3):
    t=rd(os.path.join(RES,"holdout_c3_%s.log"%c3)); out={}
    if not t: return out
    lines=t.splitlines()
    for i,ln in enumerate(lines):
        if ln.startswith("method ") and "rmse_mean" in ln:
            for l2 in lines[i+1:]:
                if not l2.strip(): break
                p=l2.split()
                if len(p)>=7 and p[0] in ("OLS","RIDGE","LASSO","ALASSO"):
                    out[p[0]]=(float(p[1]),float(p[2]),float(p[3]))
            break
    return out

rows=[]
for c3 in CUTS:
    ho=holdout(c3)
    for m in METHODS:
        rows.append(dict(c3=float(c3),method=m,fit=fit_row(m,c3),kappa=km(m,c3),ho=ho))
json.dump(rows,open(os.path.join(RES,"mni4_table.json"),"w"),indent=1,default=str)

print("### A. 拟合质量 + held-out RMSE + kappa (pheasy: OLS / LASSO / ALASSO, HARM_DENSE)")
print()
print("| c3 (A) | method | p (free IFC) | nnz | alpha_opt | RMSE_CV (eV/A) | RMSE in-sample | rel err | held-out OLS | held-out LASSO | held-out ALASSO | kxx (W/mK) | kzz (W/mK) |")
print("|"+"---|"*13)
for r in rows:
    d=r["fit"]; k=r["kappa"]; h=r["ho"]
    hh=lambda key: ("%.4e"%h[key][0]) if key in h else "-"
    print("| %.2f | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |"%(
        r["c3"],r["method"],(("%d"%d["p"]) if d.get("p") else "-"),(("%d"%d["nnz"]) if d.get("nnz") else "-"),
        ("%.3e"%d["alpha"]) if d.get("alpha") else "-",
        ("%.4e"%d["rmse_cv"]) if d.get("rmse_cv") else "-",
        ("%.4e"%d["rmse_in"]) if d.get("rmse_in") else "-",
        ("%.4e"%d["rel"]) if d.get("rel") else "-",
        hh("OLS"),hh("LASSO"),hh("ALASSO"),
        ("%.4f"%k[0]) if k else "-",("%.4f"%k[2]) if k else "-"))

print()
print("### B. 独立方法 symfc (对称适应最小二乘, fc2 不截断 + fc3 截断)")
print()
print("| c3 (A) | kxx (W/mK) | kyy | kzz |")
print("|---|---|---|---|")
for c3 in CUTS:
    k=ks(c3)
    if k: print("| %s | %.4f | %.4f | %.4f |"%(c3,k[0],k[1],k[2]))
    else: print("| %s | - | - | - |"%c3)

print()
print("### C. mesh 收敛 (c3=4.30, LASSO+HARM_DENSE)")
for ln in rd(os.path.join(RES,"mesh_kappa.txt")).splitlines():
    print(ln)
