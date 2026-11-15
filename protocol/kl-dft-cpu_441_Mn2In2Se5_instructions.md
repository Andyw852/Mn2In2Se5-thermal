# Mn2In2Se5 lattice thermal conductivity recalculation instructions (AutoZT kl-dft-cpu, 4×4×1 rebased supercell)

> To the executing AI: first, per iron rule 0 of the AutoZT `AGENTS.md`, read `AGENTS.md` and `skill/kl-dft-cpu/README.md` (including the two sections added by [0105]),
> then execute this file. The `conf --set` and `register` commands listed here are already authorised by the user through this file; `rerun`, `clean`, `-f`, `-y`
> still require an individual request each time. Whenever a checkpoint fails, stop and report; do not work around it on your own.

## 0. Background (why it is done this way)

- Old 3×3×1: the supercell preserves only 4 of the 12 point-group operations, the force constants break symmetry, and κ jumps back and forth with the q-mesh (0.93→0.78→0.96→0.74).
- Old 3×3×3: equilibrium-frame residual force 0.103 eV/Å, the structure is not at equilibrium, and there is an imaginary frequency of −0.44 THz.
- New approach: use the same relaxed structure with a different set of basis vectors, (a, b, a3′). In this basis **4×4×1 = 144 atoms**, all 12 operations are preserved,
  and it is still a diagonal supercell, so both pheasy and ShengBTE can use it directly; the safe cutoff is 6.6 Å (the 3×3×3 rhombohedral cell gives only 5.1 Å).
  A cross-check with phono3py and ShengBTE 2.1.0 has been done: with the same force constants, the κ from 4×4×1 and from 3×3×3 agree to within the output precision.
- AutoZT needs the [0105] patch (branch `claude/jolly-curie-jsvosv`); it fixes the places on this chain that would damage such a primitive cell or make the results asymmetric:
  1. S2/S3/S4 previously carried no spin: S1 relaxed with ISPIN=2 + MAGMOM, but S4 evaluated forces non-magnetically. This is one direct cause of the large residual force and the imaginary frequencies.
  2. Symmetrization used to replace the primitive cell with the 27-atom conventional cell (written back in place).
  3. The KSPACING of S1 could not be changed before.
  4. The ShengBTE automatic mesh did not preserve symmetry.
  5. pheasy silently took only the diagonal elements when given a 9-number supercell matrix.

### 0.1 Magnetic order and U (decided by the user)

- The magnetic order is fixed to **FM** (`MAGMOM = 2*5 7*0`); no AFM comparison is made, and the results must note "ferromagnetic assumed".
- **U_eff = 4 eV** is applied to Mn 3d (Dudarev, LDAUTYPE=2); none on In or Se.
- Criterion (revised 2026-11-02): what is required is an **insulating state, i.e. integer occupations**; 0.3 eV is only an empirical threshold.
  The S2 result with U=4: band gap 0.254 eV (about 5× SIGMA), magnetic moment 9.9844.
  The moment is not an integer, which the Gaussian broadening with SIGMA = 0.05 eV explains: the Mn 3d minority-spin conduction band is very flat and the broadening tail moves a small number of electrons across.
  How to check: take the S2 input, change only `SIGMA` to 0.01, and run a separate single point outside the AutoZT directory (an independent diagnostic, a few minutes).
  - Moment = 10.00 ± 0.005: judged insulating, **U=4 is kept**, continue with S3/S4.
    All AutoZT steps keep SIGMA = 0.05 unchanged, because S1–S4 must be consistent.
  - Still about 9.98: there really is partial occupation; stop and report.
- Do not raise U just to reach the band-gap threshold. U is set only from the literature or physical arguments; changing U changes the structure and all the phonons with it.
- The expected magnetic moment per S4 frame becomes 16 × 9.9844 ≈ 159.75 (with SIGMA = 0.05), and frames should differ by < 0.2.

## 1. Pre-flight checks (stop if any fails)

1. AutoZT is switched to a version containing [0105] and [0106] (`3d13297` of branch `claude/jolly-curie-jsvosv` or newer):
   - `grep -c "total drift" relax_common.py` returns ≥ 1 (the [0106] force-parsing fix);
   - `git -C ~/software/AutoZT log --oneline -3` contains `[0105]`;
   - `grep -c inherit_spin_u skill/kl-dft-cpu/gen_step4_disp.py` returns ≥ 1;
   - `ls scripts/rebase_for_diag_supercell.py` exists.
2. Local structure self-check (read-only): `python3 scripts/rebase_for_diag_supercell.py <attached POSCAR> --dim 4 4 1 -o /tmp/x`.
   It should print `R-3m (#166), 9 atoms, 12 point-group operations` and `[OK] 4 4 1 on the input basis vectors preserves all point-group operations, no rebasing needed`.

## 2. Create the material (use a new name; do not reuse the old directory)

- Material name: **`Mn2In2Se5_441`** (a new name. The old material directory contains the 3×3×1 `step4_disp/disp-*` and possibly fit-fc-thermal's
  `kl_bundle/`; when S5 sees `kl_bundle` it automatically enters "import old force constants" mode, and then the primitive cell and the supercell would not match.)
- Material-root POSCAR = the attached `POSCAR` (9 atoms, the relaxed structure rebased; atom order Mn Mn In In Se×5).
  ```
  autozt -tt kl-dft-cpu -p Mn2In2Se5_441 register --poscar <path to the attached POSCAR>
  ```
- **Never** run `reuse_structure.py`: it overwrites the material-root POSCAR with the CONTCAR of another skill and the basis vectors revert.
- `nac` stays at its default of on (a polar semiconductor needs the LO-TO splitting).

## 3. step.conf for each step (all written with `conf --set`; after writing, check the final values step by step with `conf`)

| step | setting | reason |
|---|---|---|
| step1_std_opt | `params.KSPACING=0.038` | primitive cell Γ **8×8×2**, symmetry-preserving (0.03 would give 10×10×2/10×10×3, which does not preserve symmetry) |
| step1_std_opt | `params.CELL_POLICY=none` | explicitly leave the input cell untouched (the default primitive also leaves an already-primitive input unchanged; writing none is safer) |
| step1_std_opt | `params.AUTO_MAG=true`, `params.MAGMOM_OVERRIDE="Mn:5 In:0 Se:0"` | high-spin Mn²⁺ FM starting point (S2–S4 inherit it automatically) |
| step1_std_opt | `params.U_OVERRIDE="Mn:4 In:0 Se:0"` | **user decision: FM + U**. Without U, PBE-D3 gives an effectively metallic state (moment 9.31, "gap" 54 meV). Writing a non-zero U explicitly for Mn bypasses the anion gating of selenides; S2–S4 inherit LDAU* automatically via [0105] |
| step1_std_opt | `params.CELL_FORCE_TOL_EV_A=0.002` | tightens the final-state maximum force to 2 meV/Å (the default 0.01 does not meet the S4 requirement). **Requires [0106]**: earlier versions always reported 0.0000 when reading forces, so this gate never took effect |
| step1_std_opt | `params.CELL_STAGE_EDIFFG=1E-6` | tightens the energy criterion of the variable-cell stages from 1E-4 to 1E-6 eV. Without it VASP stops at a force of about 0.009 eV/Å and all three passes fail the 0.002 force gate |
| step1_std_opt | `FUNC` stays `auto` | layered multi-component systems select pbe-d3; after S1 is generated, check FUNC in `workflow_method.txt`; S2–S4 inherit it automatically, do not change it midway |
| step2_static | `params.KSPACING=0.038` | 8×8×2 |
| step3_nac | `params.KSPACING=0.038` | 8×8×2 (12×12×3 can be used for higher accuracy, see section 5) |
| step4_disp | `params.SUPERCELL="4 4 1"` | three diagonal numbers, **do not write 9 numbers** |
| step4_disp | `params.KSPACING=0.038` | supercell Γ **2×2×2**, the same physical k-point set as the 8×8×2 of S1, so that the reference-frame residual force can really approach zero |
| step4_disp | the rest stays at default | `SUPERCELL_SYMMETRY=strict`, `METHOD=alm`, `FORCE_PREC=Accurate`, `FORCE_EDIFF=1E-7`, `FORCE_LREAL=.FALSE.`, `CUT3_MAX=7.0`, `MAX_DISP=500` |
| step5_fc | `params.EQ_FORCE_MAX=0.01` | stops immediately if the equilibrium-frame residual force is > 0.01 eV/Å (the default 0.2 is too loose; it is what let the 0.103 of the old 3×3×3 through) |
| step5_fc | `params.FC_IMPORT_DIR=off` | force a fresh fit, never import old force constants |
| step6_kappa | `params.NAC=off` (decided 2026-11-02) | the S3 DFPT failed: ε∞ = identity matrix, Mn Born charge −0.9, so the result is invalid. Moreover, with U=4 the moment 9.9844 does not change with SIGMA, which shows a small genuine partial occupation — not a clean insulator in this calculation — so NAC itself does not apply. Turn it off for this material only; do not change the project-level `nac: false` (that would affect every kl material in the project) |
| step6_kappa | the rest stays at default | `SOLVER=shengbte`, `BTE_METHOD=conv`, `MESH=auto` (since [0105] the mesh is raised automatically until it is symmetry-preserving), `ISOTOPE=true` |

Command examples:
```
autozt -tt kl-dft-cpu -p Mn2In2Se5_441 -j step1_std_opt conf --set params.KSPACING=0.038
autozt -tt kl-dft-cpu -p Mn2In2Se5_441 -j step4_disp conf --set params.SUPERCELL="4 4 1"
autozt -tt kl-dft-cpu -p Mn2In2Se5_441 -j step4_disp conf
```

**Do not** hand-write ISPIN/MAGMOM in the S4 `[incar]`: [0105] inherits them from S1 and expands them to the supercell as `32*5 112*0`. Hand-writing them may actually be inconsistent with S1.

## 4. Reference KPOINTS (AutoZT generates them with VASPKIT from KSPACING; check each one after generation — they must match exactly)

S1 / S2 / S3 (9-atom primitive cell), file `KPOINTS_primitive_S1-S3_8x8x2`:
```
Gamma 8x8x2
0
Gamma
   8   8   2
0.0  0.0  0.0
```
For each S4 `disp-*` (144-atom supercell), file `KPOINTS_supercell_S4_2x2x2`:
```
Gamma 2x2x2
0
Gamma
   2   2   2
0.0  0.0  0.0
```
- Verified: with KSPACING=0.038, VASPKIT gives the two meshes above whether it rounds or rounds up.
- In this basis, the primitive-cell Γ mesh n n m preserves symmetry iff **the first two numbers are equal and m is a multiple of n/gcd(n,4)**.
  Example: 8 8 2, 12 12 3, 16 16 4 preserve symmetry; 10 10 2, 10 10 3, 13 13 3 do not.
  For the supercell (4×4×1 reciprocal basis) the condition for n n m to preserve symmetry is that m is a multiple of n: 1 1 1, 2 2 2, 3 3 3 are fine.
- If the generated KPOINTS differ from the above, stop and report; do not submit.

## 5. Step-by-step checkpoints (check after each step; read-only diagnosis via ssh tail/grep per iron rule 4)

**S1_opt**
- CONTCAR still has 9 atoms; |a| = |b| ≈ 4.03 Å and the third basis vector along c ≈ 18.7 Å (tilted). Space group R-3m: run `symmetry_audit.py --mode audit`; it must be R-3m at both 1e-5 and 1e-4.
- Convergence is judged by AutoZT's `relax_injob`, in three stages (relax_common `STAGE_SPEC`):
  - stage a (ISIF=2): force < 0.05 eV/Å;
  - stages b/c (ISIF=3): `EDIFFG = 1E-4` is **positive, i.e. an energy criterion** (energy difference between adjacent ionic steps < 1e-4 eV), not a force criterion of 0.1 meV/Å;
  - the job additionally checks: final-state maximum force < `CELL_FORCE_TOL_EV_A` (default 0.01 eV/Å), and the lattice change and stress must be stable too, otherwise another pass is run.
  - the `EDIFFG = -0.001` in the template is overridden by the three-stage logic and **has no effect**.
- After S1 finishes, look at the maximum force of the last frame:
  `awk '/TOTAL-FORCE/{f=1;n=0;x=0;next} f&&/^ -+/{if(n>0)f=0;next} f{m=sqrt($4^2+$5^2+$6^2);if(m>x)x=m;n++} END{print x}' OUTCAR`
  - **≤ 0.003 eV/Å**: pass;
  - 0.003–0.01 eV/Å: set `CELL_FORCE_TOL_EV_A=0.002` on step1_std_opt, then `retry` (continue from CONTCAR, not rerun), and then `start`.
    Reason: S4 requires the equilibrium frame to be < 0.005 eV/Å, and a margin of 0.01 is not enough.
  An initial residual force of ~1 eV/Å is normal: the attached structure was relaxed under the old calculation settings, whereas this S1 uses PBE-D3 + FM, so the lattice will change noticeably.
- The INCAR contains `ISPIN = 2` and `MAGMOM = 2*5 2*0 5*0` (this is S1's own automatic injection; [0105] governs S2–S4).
- Only the **final** total magnetic moment counts; intermediate relaxation values do not. An insulating FM Mn²⁺ should be close to the integer 10.0 μB/cell:
  - ≈ 10.0: pass;
  - collapsed to ≈ 0: stop;
  - clearly not an integer (e.g. 9.3): this suggests that the functional may give a metal or semimetal. Report first, wait for the S2 band gap, and only then decide whether to continue with S3.
- Record the final total magnetic moment M₀ of S1/S2; the later S4 checks need it.

**S2_static**
- The gen log contains `[OK] S2_static: spin/DFT+U inherited from S1: ISPIN=2  MAGMOM=2*5 7*0 …`.
- The total moment equals that of S1. If it is judged metallic: stop and report first; do not simply turn nac off (it may be a collapse of the magnetic state).
- If the gen log reports "numerical micro-distortion detected": since [0105] symmetrization keeps the primitive cell. Verify that `step2_static/POSCAR` still has 9 atoms and unchanged basis vectors.

**S3_nac**
- The log has the same inheritance line; the OUTCAR contains `MACROSCOPIC STATIC DIELECTRIC TENSOR`, and ε_xx = ε_yy.
- An optional denser mesh (more accurate dielectric): set `KSPACING=0.024` for S3 and check the KPOINTS after generation. If it is `12 12 3`, keep it; if it is `13 13 3`, change to `0.026` and regenerate — it should be `12 12 3`.

**S4_disp (gen stage, before submission)**
- The log contains `[OK] supercell 4 4 1 preserves all 12 point-group operations of R-3m`. If it is blocked by `SUPERCELL_SYMMETRY`, the primitive-cell basis vectors have been changed: stop and report.
- `SPOSCAR` has 144 atoms and its lattice matches the attached `SPOSCAR_441_reference` (component-wise difference < 1e-3 Å; there will be tiny changes because the structure was re-relaxed in S1).
- `disp-00000/INCAR` has `ISPIN = 2` and `MAGMOM = 32*5 112*0` (the log contains a line `step4 00000: spin/DFT+U inherited from S1`); all `disp-*/KPOINTS` are `2 2 2`.
- Number of frames ≤ `MAX_DISP`; if it exceeds that, lower `CUT3_MAX` to 6.0 and then retry.

**S4_disp (after the calculation) — the most critical gate**
- Maximum residual force of the equilibrium frame `disp-00000`:
  ```
  awk '/TOTAL-FORCE/{f=1;n=0;next} f&&/^ -+/{if(n>0)f=0;next} f{m=sqrt($4^2+$5^2+$6^2);if(m>x)x=m;n++} END{print x}' disp-00000/OUTCAR
  ```
  **< 0.005 eV/Å passes** (the old 3×3×1 gave 0.004). For 0.005–0.01, report and continue; above 0.01, stop: first check whether the magnetic moment, KPOINTS and FUNC/IVDW are consistent with S1.
- Total magnetic moment of each frame (take the last occurrence in each OUTCAR):
  `for f in disp-*/OUTCAR; do grep "number of electron" $f | tail -1 | awk '{printf "%.1f\n",$NF}'; done | sort | uniq -c`，
  There should be only one value, and it should be ≈ 16 × M₀ (M₀ = the final total magnetic moment of S2; when M₀ = 10 this is 160).
  If a frame ends up in another magnetic state, list that frame and report it.

**S5_fc（3090，pheasy）**
- `fc_dataset.json`：`supercell_matrix = [[4,0,0],[0,4,0],[0,0,1]]`，`natom_super = 144`，`equilibrium_max_force < 0.01`。
- `phonon_summary.json`：`stable = true`、`fc2_symmetry.broken = false`（`n_ops = 12`）、`fit_rmse_relative < 0.02`；
  `min_frequency_THz` not below −0.1. If imaginary frequencies away from Γ remain, stop and report, including the q-point and branch index.
- `cutoff_scan.json`: the candidate cutoff ceiling is about 6.6 Å (cutoffs above the safe cutoff are dropped in S4).

**S5.1_plot**: look only at `band_nonac*`; the acoustic branches between 0 and 1 THz must have no negative frequencies near Γ.
The `step3_nac` directory still exists, so S5 still generates BORN from that invalid vasprun; `band_nac*` is therefore wrong — ignore it.
For S6, confirm that the log shows `NAC=False`.

**S6_kappa（ShengBTE conv）**
- The log **must not** contain `rotations are incompatible with the q-point grid` ([0105] already raises the automatic mesh until it is symmetry-preserving; if it appears, stop and report).
- `kappa_summary.json`: κ_xx = κ_yy, off-diagonal elements ≈ 0; `mesh_converged = true`. If it stopped because of `MESH_CONV_MAX_LENGTH`, report and ask whether to raise it.
- κ_zz is very sensitive to the interlayer-direction mesh (in the control test the κ_zz at 16×16×4 and 16×16×16 differ by 60%). After convergence, do one more confirmation:
  on the converged mesh n n m, double the third number (e.g. from `MESH="16 16 8"` to `"16 16 16"`; n n 2m also satisfies the symmetry-preserving condition),
  and check whether κ_zz changes by < 5%. If not, refine further and report.

## 5.5 How to drive the workflow (AGENTS.md iron rule 8: do not write your own bash/ssh monitoring loop)

- For monitoring use only `autozt summary --diff`, or AutoZT's scoped monitor:
  `autozt -tt kl-dft-cpu -p Mn2In2Se5_441 monitor -d -i 300`. Do not start another bash background loop; if one is already running, stop it by PID.
- **S1→S2→S3 and the S4 generation stage: human-checked, with auto-advance off.** After each step is OK, verify the checkpoints above first, then
  `autozt -tt kl-dft-cpu -p Mn2In2Se5_441 start`.
- **Generate S4 first, do not submit**: `autozt -tt kl-dft-cpu -p Mn2In2Se5_441 -j step4_disp init`, then verify SPOSCAR,
  KPOINTS (2 2 2), `disp-00000/INCAR` (ISPIN/MAGMOM) and the number of frames, and only then `start`. S4 is the most expensive step; a bad generation wastes dozens of frames.
- Once all of S4 is done and the equilibrium-frame residual force and the per-frame moments have been verified, auto-advance may be switched on for this material:
  `autozt -tt kl-dft-cpu -p Mn2In2Se5_441 auto on`. S5 has the `EQ_FORCE_MAX=0.01` hard gate as a backstop.

## 6. "Primitive cell and supercell mismatch" risks and the corresponding measures

| where it can go wrong | consequence | measure in this plan |
|---|---|---|
| material-root POSCAR overwritten by `reuse_structure.py` | basis vectors revert to the old ones and 4 4 1 no longer preserves symmetry | forbidden; the S4 symmetry gate is the backstop |
| S1 standardizes the primitive cell | basis vectors are changed | `CELL_POLICY=none` |
| symmetrization of small distortions in S2/S3/S4 | used to replace the cell with the 27-atom conventional cell | fixed by [0105] (keeps the primitive cell); the checkpoint verifies the atom count |
| S4 given a 9-number SUPERCELL | pheasy used to silently take only the diagonal | write only `4 4 1`; since [0105] it raises an error instead |
| S5 importing an old kl_bundle | the force constants belong to another cell | new material name + `FC_IMPORT_DIR=off` |
| S1 magnetic, S4 non-magnetic | the reference structure is not at equilibrium and imaginary frequencies appear | automatic inheritance via [0105]; the checkpoint verifies MAGMOM and the moment |
| asymmetric k-mesh / inconsistent sampling between S1 and S4 | symmetry-breaking bias in the forces, larger residual force | KSPACING=0.038 (8×8×2 ↔ 2×2×2, the same k-point set) |
| asymmetric ShengBTE mesh | κ_xx ≠ κ_yy | raised automatically by [0105]; check the log |

## 7. Reporting content

Report after each completed step using the format of section 6 of AGENTS.md. The final report should contain:
- the S4 equilibrium-frame residual force, the number of frames, and the per-frame moments;
- the S5 minimum frequency, the fit error, and the selected third-order cutoff;
- the S6 κ(300 K) xx/yy/zz, the converged mesh, and the result of the κ_zz refinement check;
- the comparison with the old 3×3×1 results (jumping back and forth between 0.74 and 0.96 W/mK).
