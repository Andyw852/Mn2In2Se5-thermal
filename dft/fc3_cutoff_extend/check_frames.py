#!/usr/bin/env python3
"""Check the VASP force runs in step4_disp/: finished, electronically converged,
forces sane. Prints one line per frame and a summary.
  python check_frames.py [new_frames]
"""
import glob
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

root = sys.argv[1] if len(sys.argv) > 1 else "new_frames"
ok = bad = pending = 0
for d in sorted(glob.glob(os.path.join(root, "disp-*"))):
    tag = os.path.basename(d)
    vr = os.path.join(d, "vasprun.xml")
    if not os.path.exists(vr) or "</modeling>" not in open(vr, errors="ignore").read()[-200:]:
        print(f"{tag}  pending"); pending += 1; continue
    t = ET.parse(vr).getroot()
    nelm = int(t.find(".//parameters//i[@name='NELM']").text)
    nsc = len(t.findall(".//calculation/scstep"))
    f = np.array([[float(x) for x in v.text.split()] for v in t.findall(".//calculation/varray[@name='forces']/v")])
    fmax = np.linalg.norm(f, axis=1).max()
    drift = np.linalg.norm(f.sum(axis=0))
    flag = "OK"
    if nsc >= nelm:
        flag = "NOT CONVERGED (hit NELM)"
    elif fmax > 2.0:
        flag = "force suspiciously large"
    print(f"{tag}  {flag:30s} scf {nsc:3d}  max|F| {fmax:.4f} eV/A  |sum F| {drift:.2e}")
    ok += flag == "OK"; bad += flag != "OK"
print(f"\n{ok} OK, {bad} problem, {pending} pending")
