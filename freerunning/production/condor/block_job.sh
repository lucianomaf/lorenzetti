#!/bin/bash
# One block of a production in an HTCondor job (example for lxplus, the CERN batch system): brings the package of the build, runs block.sh
# in the Lorenzetti image (apptainer) and copies the files to EOS with xrdcp (temporary name, size check, rename). The FRHITS (inclusive,
# hard) or the hard.txt (mark) goes last and marks the block as done; a block already done is skipped. See ../README.md.
# Usage: block_job.sh <mode> <production> <block> <events> <seed> <threads> <watchdog, minutes>
# Environment (set in the submit file): PROD_IMAGE (the .sif of the image), PROD_PACKAGE (the .tgz of the build, on EOS), PROD_HOME (the
# home folder where the build was made; the package is mounted there), PROD_REPO (the folder of the repository in the package), PROD_OUT (the folder of the productions on EOS), PROD_CODE (this
# freerunning/production folder), PROD_FIELD_MAP (default: the map on CVMFS), XRD (default root://eosuser.cern.ch/).
set -u
MODE=$1; PROD=$2; B=$3; NOV=$4; SEED=$5; NT=$6; TMIN=$7
XRD=${XRD:-root://eosuser.cern.ch/}
HOST=${XRD#root://}; HOST=${HOST%/}
FIELD=${PROD_FIELD_MAP:-/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/MagneticFieldMaps/bfieldmap_7730_20400_14m.root}
DEST=$PROD_OUT/$PROD
NAME=${PROD}_$B
W=${_CONDOR_SCRATCH_DIR:-$PWD}
trap 'rm -rf "$W/home" "$W/run"' EXIT

marker=$DEST/frhits/$NAME.FRHITS.root
[ $MODE = mark ] && marker=$DEST/hard/$NAME.hard.txt
if xrdfs $HOST stat $marker > /dev/null 2>&1; then echo "block $NAME already done"; exit 0; fi
mkdir -p "$W/home" "$W/run"
xrdcp -s $XRD$PROD_PACKAGE "$W/package.tgz" && tar xzf "$W/package.tgz" -C "$W/home" && rm -f "$W/package.tgz" || { echo "cannot bring the package"; exit 11; }
cp "$PROD_CODE/block.sh" "$PROD_CODE/minbias_inelastic_A3_14TeV.cmnd" "$W/run/"
if [ $MODE = mark ]; then xrdcp -s $XRD$DEST/evt/$NAME.EVT.root "$W/run/stored.EVT.root" || { echo "cannot bring the EVT"; exit 12; }; fi
apptainer exec --cleanenv -B "$W/home:$PROD_HOME" -B "$W/run:/run_block" -B /cvmfs "$PROD_IMAGE" bash -c "cd /run_block && \
  LZT_SETUP=$PROD_HOME/$PROD_REPO/build/lzt_setup.sh LZT_RUN4_GEOMETRY=$PROD_HOME/atlas_geometry/geometry-ATLAS-P2-RUN4-01-01-00-sim.db \
  LZT_FIELD_MAP=$FIELD LZT_STORED_EVT=/run_block/stored.EVT.root bash block.sh $MODE $PROD $B $NOV $SEED $NT $TMIN"
e=$?
cd "$W/run"
[ $e = 0 ] || { tail -20 gen.log simu.log 2>/dev/null; exit 31; }

up() {  # local file, EOS folder
  local f=$1 d=$2 n
  xrdfs $HOST mkdir -p $d
  xrdcp -f -s $f $XRD$d/$f.part || return 1
  n=$(xrdfs $HOST stat $d/$f.part | awk '/Size:/{print $2}')
  [ "$n" = "$(stat -c %s $f)" ] || { echo "size differs for $f: $n on EOS"; return 1; }
  xrdfs $HOST mv $d/$f.part $d/$f
}
case $MODE in
  inclusive) up $NAME.EVT.root $DEST/evt && up $NAME.ok $DEST/ok && up $NAME.FRHITS.root $DEST/frhits || exit 41 ;;
  hard) up $NAME.EVT.root $DEST/evt && up $NAME.TRUTH.root $DEST/truth && up $NAME.ok $DEST/ok && up $NAME.FRHITS.root $DEST/frhits || exit 41 ;;
  mark) up $NAME.TRUTH.root $DEST/truth && up $NAME.hard.txt $DEST/hard || exit 41 ;;
esac
echo "block $NAME written"
