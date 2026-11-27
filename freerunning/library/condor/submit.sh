#!/bin/bash
# Writes and submits the HTCondor DAG that builds a free-running library: events -> blocks (one job per block) -> plan -> pieces (one job
# per piece) and cells -> check, with the geometry step before cells. Light: run it on the submit node; everything else runs in the
# jobs, inside the Lorenzetti image.
# The input list (see ../README.md) must use paths the jobs can read (/eos/...). The code is this checkout of the repository; the work
# folder (lists, logs, check.txt) must be readable and writable by the jobs (AFS with MY.SendCredential at CERN).
# The geometry step needs the package of the build (../../production/condor/make_package.sh): its path on EOS, the home folder where the
# build was made and the folder of the repository in it, as in ../../production/condor/production.sub.
# Usage: bash submit.sh <input list> <library folder on EOS> <work folder> <Lorenzetti .sif> <package .tgz> <home of the build>
#        <repository folder> [hits per piece] [extra submit lines file]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
LIST=$(readlink -f "$1"); OUT=${2%/}; WORK=$(mkdir -p "$3" && readlink -f "$3"); IMAGE=$4; PACKAGE=$5; BHOME=$6; REPO=$7
PER=${8:-160000000}; EXTRA=${9:-}
CODE=$(dirname "$HERE")
mkdir -p "$WORK/logs"
cp "$LIST" "$WORK/list.csv"
common() {
  printf '%s\n' "universe = vanilla" "executable = $HERE/lib_job.sh" "log = $WORK/logs/dag.log" \
    "environment = \"LIB_CODE=$CODE LIB_IMAGE=$IMAGE LIB_OUT=$OUT LIB_WORK=$WORK LIB_PACKAGE=$PACKAGE LIB_HOME=$BHOME LIB_REPO=$REPO\"" \
    "should_transfer_files = YES" "when_to_transfer_output = ON_EXIT" 'transfer_output_files = ""' "MY.SendCredential = True"
  [ -n "$EXTRA" ] && cat "$EXTRA"
  return 0
}
sub() {  # file, arguments, log name, memory MB, flavour, queue line
  { common; printf '%s\n' "arguments = $2" "output = $WORK/logs/$3.out" "error = $WORK/logs/$3.err" "request_memory = $4" \
      "+JobFlavour = \"$5\"" "$6"; } > "$WORK/$1"
}
sub geometry.sub "geometry" geometry 4000 microcentury "queue"
sub events.sub "events" events 4000 longlunch "queue"
sub blocks.sub "block \$(name) \$(frhits) \$(marking) \$(first) \$(n)" 'block_$(name)' 3000 microcentury \
  "queue name,frhits,marking,first,n from $WORK/blocks.csv"
sub plan.sub "plan $PER" plan 4000 longlunch "queue"
sub pieces.sub "piece \$(det) \$(j)" 'piece_$(det)_$(j)' 8000 workday "queue det,j from $WORK/pieces.csv"
sub cells.sub "cells" cells 4000 longlunch "queue"
sub check.sub "check" check 6000 workday "queue"
printf '%s\n' "JOB geometry $WORK/geometry.sub" "JOB events $WORK/events.sub" "JOB blocks $WORK/blocks.sub" "JOB plan $WORK/plan.sub" \
  "JOB pieces $WORK/pieces.sub" "JOB cells $WORK/cells.sub" "JOB check $WORK/check.sub" \
  "PARENT events CHILD blocks" "PARENT geometry CHILD cells" "PARENT blocks CHILD plan cells" "PARENT plan CHILD pieces" \
  "PARENT pieces cells CHILD check" \
  "RETRY geometry 2" "RETRY blocks 2" "RETRY pieces 2" > "$WORK/library.dag"
cd "$WORK" && condor_submit_dag -f -batch-name "library_$(basename "$OUT")" library.dag
echo "follow: condor_q -dag -nobatch; result in $WORK/check.txt"
