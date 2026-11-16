#!/bin/bash
# Run the layer-normal sound-velocity check on every FC2 that is relevant to the
# "is the 4x4x1 supercell too soft between the layers?" question, and write one
# JSON with all the cases.
#
#   MnIn2Se4 4x4x1  - the FC2 that WAS too soft (validation: should come out bad)
#   MnIn2Se4 4x4x2  - the corrected FC2      (validation: should match DFT)
#   Mn2In2Se5 fc/    - the FC2 the kappa runs used  (the actual question)
#   Mn2In2Se5 fit/   - the other FC2 revision, for completeness
#
# Pure post-processing: phonopy reads an FC2, no DFT, seconds per case.
set -eu
cd "$(dirname "$(readlink -f "$0")")"
S="python3 ../fc_sound_velocity.py"
D4=../../MnIn2Se4
OUT=fc_sound_velocity.json
rm -f "$OUT"

$S --fc2 $D4/fit_cut3_6p45/fc2.hdf5 --poscar $D4/fit/POSCAR --supercell 4 4 1 \
     --label "MnIn2Se4 4x4x1 (the FC2 that was too soft)" \
     --ref-vta 1338.319 --ref-c44 9.402 --json "$OUT"

$S --fc2 $D4/fc2_442_phonon/fc2.hdf5 --poscar $D4/fit/POSCAR --supercell 4 4 2 \
     --label "MnIn2Se4 4x4x2 (the corrected FC2)" \
     --ref-vta 1338.319 --ref-c44 9.402 --json "$OUT"

$S --fc2 ../../../Mn2In2Se5/fc/fc2.hdf5 --poscar ../fit/POSCAR --supercell 4 4 1 \
     --label "Mn2In2Se5 fc/ 4x4x1 (used for the kappa runs)" --json "$OUT"

$S --fc2 ../fit/fc2.hdf5 --poscar ../fit/POSCAR --supercell 4 4 1 \
     --label "Mn2In2Se5 fit/ 4x4x1 (other revision)" --json "$OUT"

echo
echo "=== summary (fc_sound_velocity.json)"
python3 - "$OUT" <<'EOF'
import json, sys
rows = json.load(open(sys.argv[1]))
print(f"{'case':52s} {'v_TA':>7s} {'C44':>7s} {'C44 route2':>10s} {'C33':>7s}")
for r in rows:
    print(f"{r['label']:52s} {r['v_interlayer_TA_m_s']:7.1f} {r['C44_GPa']:7.2f} "
          f"{r['C44_GPa_from_inplane_route']:10.2f} {r['C33_GPa']:7.2f}")
EOF
