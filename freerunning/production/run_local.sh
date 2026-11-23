#!/bin/bash
# Runs the blocks of a production one after the other on one machine, inside the Lorenzetti image (see README.md), with block.sh, and
# keeps the files in the layout the library expects: <output>/<production>/{evt,frhits,truth,hard,ok}/. A block already done is skipped
# (inclusive and hard: its FRHITS exists; mark: its hard.txt exists); mark reads the EVT of the block from <output>/<production>/evt/.
# Usage: bash run_local.sh <mode> <production> <blocks file> <output folder> [threads (4)] [watchdog, minutes (110)]
# Blocks file: one block per line, "block events seed" (see README.md for the productions of the library).
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
MODE=$1; PROD=$2; LIST=$3; OUT=$(mkdir -p "$4/$2" && cd "$4/$2" && pwd); NT=${5:-4}; TMIN=${6:-110}
mkdir -p "$OUT"/{evt,frhits,truth,hard,ok}
bad=0
while read -r B NOV SEED; do
  [ -z "${B:-}" ] || [ "${B:0:1}" = "#" ] && continue
  NAME=${PROD}_$B
  done_file=$OUT/frhits/$NAME.FRHITS.root
  [ $MODE = mark ] && done_file=$OUT/hard/$NAME.hard.txt
  [ -s "$done_file" ] && { echo "$NAME: done"; continue; }
  W=$(mktemp -d "$OUT/.work_${NAME}_XXXX")
  ( cd "$W" && LZT_STORED_EVT=$OUT/evt/$NAME.EVT.root bash "$HERE/block.sh" $MODE $PROD $B $NOV $SEED $NT $TMIN ) > "$OUT/ok/$NAME.$MODE.log" 2>&1
  if [ $? = 0 ]; then
    for f in "$W"/$NAME.EVT.root "$W"/$NAME.TRUTH.root "$W"/$NAME.ok "$W"/$NAME.hard.txt; do
      [ -e "$f" ] || continue
      case $f in *.EVT.root) d=evt ;; *.TRUTH.root) d=truth ;; *.ok) d=ok ;; *.hard.txt) d=hard ;; esac
      mv "$f" "$OUT/$d/"
    done
    # the FRHITS last: it marks the block as done
    [ -e "$W/$NAME.FRHITS.root" ] && mv "$W/$NAME.FRHITS.root" "$OUT/frhits/"
    rm -rf "$W"
    echo "$NAME: ok"
  else
    bad=$((bad + 1))
    echo "$NAME: failed (log $OUT/ok/$NAME.$MODE.log, work folder $W)"
  fi
done < "$LIST"
[ $bad = 0 ]
