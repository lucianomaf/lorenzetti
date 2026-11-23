#!/bin/bash
# One block of a minimum-bias production for the free-running library, inside the Lorenzetti image with this repository built. It does not
# know where it runs: run_local.sh calls it on one machine, condor/block_job.sh in an HTCondor job. See README.md.
#
# Usage: bash block.sh <mode> <production> <block> <events> <seed> [threads (4)] [watchdog, minutes (110)]
#   inclusive  generates the block with Pythia 8 (A3, the card next to this script, seed <seed>) and simulates it:
#              <name>.EVT.root, <name>.FRHITS.root, <name>.ok
#   mark       generates the block again with the same seed, marks the hard events (selection of the ATLAS sample 800831) and checks,
#              event by event, that the particles are those of the stored EVT ($LZT_STORED_EVT): <name>.hard.txt (EventNumber and bits:
#              1 jet, 2 photon, 4 B hadron and lepton; 0 not hard), <name>.TRUTH.root (full Pythia record); nothing if the check fails
#   hard       generates hard events only (the same selection) and simulates them: <name>.EVT.root, <name>.TRUTH.root, <name>.FRHITS.root,
#              <name>.ok (with the efficiency of the filter)
# <name> = <production>_<block>; the files are written in the current folder.
#
# Environment: LZT_SETUP (build/lzt_setup.sh of the build), LZT_RUN4_GEOMETRY (geometry-ATLAS-P2-RUN4-01-01-00-sim.db), LZT_FIELD_MAP
# (bfieldmap_7730_20400_14m.root), LZT_PYTHIA_CARD (default: minbias_inelastic_A3_14TeV.cmnd next to this script), LZT_STORED_EVT (mark).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
MODE=$1; PROD=$2; B=$3; NOV=$4; SEED=$5; NT=${6:-4}; TMIN=${7:-110}
NAME=${PROD}_${B}
CARD=${LZT_PYTHIA_CARD:-$HERE/minbias_inelastic_A3_14TeV.cmnd}
case $MODE in inclusive|mark|hard) ;; *) echo "unknown mode $MODE"; exit 2 ;; esac
source "$LZT_SETUP" > /dev/null 2>&1 || { echo "cannot source $LZT_SETUP"; exit 10; }
ulimit -v unlimited 2>/dev/null

# generation: one interaction per event, all particles, the vertex spread of the HL-LHC (sigma 50 mm and 175 ps)
GEN="--nov $NOV -s $SEED --pileup-file $CARD --pileup-per-bunch-crossing 1 --bc-id-start 0 --bc-id-end 0 --eta-max 10 --pt-min-charged 0
  --pt-min-neutral 0 --sigma-z 50 --sigma-t 175"
FILTER=
[ $MODE = mark ] && FILTER="--hard-filter mark --truth"
[ $MODE = hard ] && FILTER="--hard-filter keep --truth"
OUTEVT=$NAME.EVT.root
# mark: the regenerated EVT goes to a folder of its own, under the same name (gen_minbias.py names the truth file after it)
[ $MODE = mark ] && { mkdir -p regenerated; OUTEVT=regenerated/$NAME.EVT.root; }
T0=$(date +%s)
gen_minbias.py -o $OUTEVT $GEN $FILTER > gen.log 2>&1
echo "gen_exit=$? gen_s=$(( $(date +%s) - T0 ))" > times.txt
[ -s $OUTEVT ] || { echo "generation failed"; tail -20 gen.log; exit 21; }

if [ $MODE = mark ]; then
  [ -s "${LZT_STORED_EVT:-}" ] || { echo "LZT_STORED_EVT (the stored EVT of the block) is missing"; exit 12; }
  python3 - "$LZT_STORED_EVT" $OUTEVT $NAME.hard.txt.part > check.txt 2>&1 << 'PY'
import sys
import numpy as np
import ROOT
stored, regen, out = sys.argv[1:4]
fa, fb = ROOT.TFile.Open(stored), ROOT.TFile.Open(regen)
ta, tb = fa.Get("particles"), fb.Get("particles")
cols = ["p_pdg_id", "p_bc_id", "p_px", "p_py", "p_pz", "p_e", "p_prod_x", "p_prod_y", "p_prod_z", "p_prod_t"]
n, diff = ta.GetEntries(), 0
with open(out, "w") as o:
    if tb.GetEntries() != n:
        diff = -1
    else:
        for i in range(n):
            ta.GetEntry(i), tb.GetEntry(i)
            if any(not np.array_equal(np.array(getattr(ta, c)), np.array(getattr(tb, c))) for c in cols):
                diff += 1
            h = set(np.array(tb.p_hard)) - {0}
            o.write("%d %d\n" % (tb.EventNumber, max(h) if h else 0))
print("events", n, tb.GetEntries(), "different", diff)
PY
  cat times.txt check.txt
  mv regenerated/$NAME.TRUTH.root . 2>/dev/null
  grep -q "different 0$" check.txt && [ -s $NAME.TRUTH.root ] || { echo "block $NAME does not match the stored EVT"; rm -f $NAME.hard.txt.part; exit 31; }
  mv $NAME.hard.txt.part $NAME.hard.txt
  rm -rf regenerated
  echo "block $NAME marked: $(awk '$2 > 0' $NAME.hard.txt | wc -l) hard"
  exit 0
fi

[ $MODE = hard ] && { [ -s $NAME.TRUTH.root ] || exit 21; grep -o "HardFilter: [0-9]* events kept of [0-9]* generated" gen.log >> times.txt; }
# simulation: the ATLAS Run 4 calorimeters of the v2 detector and the free-running hits in the layout of the ATLAS HITS
ALL="--enable-solenoid-field --tile-atlas-geometry --tile-atlas-cells --atlas-material-in-front --tile-atlas-itc --atlas-emb
  --atlas-endcap-cryostat --atlas-barrel-cryostat --atlas-hec --atlas-emec --atlas-emb-cells --atlas-fcal --tile-atlas-gap --atlas-fcal-cells
  --atlas-fcal-frozen-showers --atlas-g4"
RUN4="--atlas-run4 --atlas-run4-geometry $LZT_RUN4_GEOMETRY --atlas-field-map $LZT_FIELD_MAP"
MORE="--atlas-emec-charge-collection --atlas-forward-shielding --atlas-full-field --atlas-barrel-gap-services --tile-atlas-run4
  --atlas-emec-argon-fraction"
T1=$(date +%s)
simu_trf.py -i $NAME.EVT.root -o $NAME.FRHITS.root --nov $NOV -t $TMIN -nt $NT --detector v2 $ALL $RUN4 --free-running-hits \
  --free-running-atlas-format $MORE > simu.log 2>&1
e=$?
echo "simu_exit=$e simu_s=$(( $(date +%s) - T1 ))" >> times.txt
echo "block=$NAME events=$NOV seed=$SEED exit=$e total_s=$(( $(date +%s) - T0 )) watchdog=$(grep -c WatchDog simu.log) node=$(hostname)" \
  "$(tr '\n' ' ' < times.txt)" > $NAME.ok
cat $NAME.ok
[ $e = 0 ] && [ -s $NAME.FRHITS.root ] || { tail -20 simu.log; exit 31; }
exit 0
