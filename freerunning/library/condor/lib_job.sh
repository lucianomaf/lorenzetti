#!/bin/bash
# One step of the free-running library as an HTCondor job, inside the Lorenzetti image (submit.sh writes the jobs).
# Environment: LIB_CODE (this freerunning/library folder), LIB_IMAGE (the Lorenzetti .sif), LIB_OUT (the library folder on EOS),
# LIB_WORK (the work folder with the lists, on AFS). The jobs read EOS by its mount (/eos), with the Kerberos ticket of the job passed
# into the image (the ROOT of the image has no xrootd), and write it with xrdcp.
# Usage: lib_job.sh events | block <name> <FRHITS> <marking> <first index> <events> | plan <hits per piece> | piece <subdetector> <n> |
#        cells | check
set -u
STEP=$1; shift
W=${_CONDOR_SCRATCH_DIR:-$PWD}/libjob
mkdir -p "$W/lib/support" "$W/lib/per_channel"
URL=root://eosuser.cern.ch/
py() { apptainer exec --cleanenv --env KRB5CCNAME="${KRB5CCNAME:-}" -B "${_CONDOR_SCRATCH_DIR:-$W}" -B /eos -B /afs -B /cvmfs -B "$W" "$LIB_IMAGE" python3 "$@"; }
put() {  # local file, EOS path: copy under a temporary name, compare the size, rename
  xrdfs eosuser.cern.ch mkdir -p "$(dirname "$2")"
  xrdcp -f -s "$1" "$URL$2.part" || return 1
  local n; n=$(xrdfs eosuser.cern.ch stat "$2.part" | awk '/Size:/{print $2}')
  [ "$n" = "$(stat -c %s "$1")" ] || { echo "size differs for $2: $n on EOS"; return 1; }
  xrdfs eosuser.cern.ch mv "$2.part" "$2"
}
case $STEP in
  events)
    py "$LIB_CODE/make_support.py" events "$LIB_WORK/list.csv" "$W/lib" "$W/blocks.csv" || exit 2
    for f in samples.json productions.json manifest.root; do put "$W/lib/support/$f" "$LIB_OUT/support/$f" || exit 3; done
    awk -F, -v d="$LIB_OUT/work/blocks" '{print d "/" $1}' "$W/blocks.csv" > "$LIB_WORK/block_files.txt"
    cp "$W/blocks.csv" "$LIB_WORK/blocks.csv" ;;
  block)
    NAME=$1; DEST=$LIB_OUT/work/blocks/$NAME
    if xrdfs eosuser.cern.ch stat "$DEST" > /dev/null 2>&1; then echo "$NAME already on EOS"; exit 0; fi
    py "$LIB_CODE/block_to_channels.py" "$2" "$3" "$4" "$5" "$W/$NAME" || exit 2
    put "$W/$NAME" "$DEST" || exit 3 ;;
  plan)
    py "$LIB_CODE/merge_pieces.py" plan "$LIB_WORK/block_files.txt" "$W/lib" "$W/pieces.csv" "$1" || exit 2
    put "$W/lib/support/per_channel_index.json" "$LIB_OUT/support/per_channel_index.json" || exit 3
    cp "$W/pieces.csv" "$LIB_WORK/pieces.csv" ;;
  piece)
    xrdcp -f -s "$URL$LIB_OUT/support/per_channel_index.json" "$W/lib/support/" || exit 1
    py "$LIB_CODE/merge_pieces.py" piece "$LIB_WORK/block_files.txt" "$W/lib" "$1" "$2" || exit 2
    F=$(printf '%s_%03d.root' "$1" "$2")
    put "$W/lib/per_channel/$F" "$LIB_OUT/per_channel/$F" || exit 3 ;;
  cells)
    py "$LIB_CODE/make_support.py" cells "$LIB_WORK/block_files.txt" "$W/lib" > "$W/cells.txt"; e=$?
    cat "$W/cells.txt"; cp "$W/cells.txt" "$LIB_WORK/cells.txt"
    for f in tile_cells.csv lar_cells.csv; do put "$W/lib/support/$f" "$LIB_OUT/support/$f" || exit 3; done
    exit $e ;;
  check)
    for f in samples.json per_channel_index.json; do xrdcp -f -s "$URL$LIB_OUT/support/$f" "$W/lib/support/" || exit 1; done
    py "$LIB_CODE/check_library.py" "$LIB_WORK/list.csv" "$W/lib" "$LIB_OUT/per_channel" > "$W/check.txt"; e=$?
    cat "$W/check.txt"; cp "$W/check.txt" "$LIB_WORK/check.txt"; exit $e ;;
  *) echo "unknown step $STEP"; exit 1 ;;
esac
