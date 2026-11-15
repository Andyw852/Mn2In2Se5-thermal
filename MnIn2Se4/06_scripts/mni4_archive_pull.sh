#!/bin/bash
# mni4_archive_pull.sh -- assemble the MnIn2Se4 archive from 3090 + hanhai.
set -u
T=/home/wangchao/software/taskflow-v2.0/tmp/mn2inse4
A=$T/MnIn2Se4_data_archive_20260928
mkdir -p $A/{00_report,01_structure,02_dataset,03_fit_results,04_kappa,05_scripts,06_heldout_RMSE,07_kappa_raw_output_hdf5}
# structure + dataset (local staged copies are authoritative for the input)
cp $T/stage/POSCAR $T/stage/SPOSCAR $A/01_structure/
cp $T/stage/dataset_disps.npy $T/stage/dataset_forces.npy $A/02_dataset/
scp -q wangchao_3090:/data/wangchaoyue852/pheasy_mninse4_189/disp_matrix.pkl wangchao_3090:/data/wangchaoyue852/pheasy_mninse4_189/force_matrix.pkl wangchao_3090:/data/wangchaoyue852/pheasy_mninse4_189/dataset_alignment.json wangchao_3090:/data/wangchaoyue852/pheasy_mninse4_189/ndata_total.txt $A/02_dataset/ 2>/dev/null
# scripts
cp $T/mni4_*.sh $T/mni4_*.py $T/kappa_full_solver.py $T/check_supercell.py $T/mni4_shells.py $A/05_scripts/ 2>/dev/null
cp $T/stage/pheasy_fit.sh $A/05_scripts/pheasy_fit.sh_v6c
cp $T/prepare_dataset.py $T/symfc_fit_test.py $T/hanhai_symfc_run.sh $A/05_scripts/ 2>/dev/null
# fit results: POSCAR + fc2/fc3 + run.log for every dir
echo "[pull] fits from 3090"
ssh wangchao_3090 'cd /data/wangchaoyue852/pheasy_mninse4_189 && tar czf - $(for d in ghsc_*; do echo $d/POSCAR $d/fc2.hdf5 $d/fc3.hdf5 $d/run.log; done) 2>/dev/null' | tar xzf - -C $A/03_fit_results/
echo "[pull] symfc"
ssh wangchao_3090 'true' 2>/dev/null
SOCK="$HOME/.ssh/tf-wangchao217@114.214.255.25:22.sock"
/home/wangchao/bin/hanhai25-connect >/dev/null 2>&1 || true
ssh -o ControlPath="$SOCK" -o BatchMode=yes hanhai25 'cd ~/work/mninse4_189_kappa && tar czf - $(for d in symfc_c3_*; do echo $d/POSCAR $d/fc2.hdf5 $d/fc3.hdf5 $d/fit_test.py; done) 2>/dev/null' | tar xzf - -C $A/03_fit_results/
echo "[pull] kappa hdf5"
ssh -o ControlPath="$SOCK" -o BatchMode=yes hanhai25 'cd ~/work/mninse4_189_kappa && tar czf - $(for d in ghsc_* symfc_c3_*; do [ -f $d/kappa-m111111.hdf5 ] && echo $d/kappa-m111111.hdf5 $d/POSCAR $d/kappa_run.log; done) 2>/dev/null' | tar xzf - -C $A/07_kappa_raw_output_hdf5/
mkdir -p $A/04_kappa
cp $T/results/kappa_all.txt $A/04_kappa/kappa_all_cutoffs_summary.txt
cp $T/results/mesh_kappa.txt $A/04_kappa/mesh_convergence_11_15.txt 2>/dev/null
du -sh $A; find $A -type f | wc -l
