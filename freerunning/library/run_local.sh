#!/bin/bash
# Builds a free-running library on one machine, all steps in sequence (the same scripts the HTCondor jobs run), inside the Lorenzetti
# image. The geometry step needs the build of Lorenzetti and the external files of the v2 detector (LZT_SETUP, LZT_RUN4_GEOMETRY,
# LZT_FIELD_MAP, as ../production/block.sh); the other steps only ROOT and numpy. See README.md.
# Usage: bash run_local.sh <input list> <library folder> <work folder> [hits per piece]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
LIST=$1; LIB=$2; WORK=$3; PER=${4:-160000000}
mkdir -p "$LIB/support" "$WORK/blocks"
( source "$LZT_SETUP" > /dev/null 2>&1 || { echo "cannot source LZT_SETUP=$LZT_SETUP"; exit 10; }
  python3 "$HERE/geometry_channels.py" "$LIB/support/geometry_channels.npz" )
python3 "$HERE/make_support.py" events "$LIST" "$LIB" "$WORK/blocks.csv"
: > "$WORK/block_files.txt"
while IFS=, read -r name frhits marking first n; do
  python3 "$HERE/block_to_channels.py" "$frhits" "$marking" "$first" "$n" "$WORK/blocks/$name"
  echo "$WORK/blocks/$name" >> "$WORK/block_files.txt"
done < "$WORK/blocks.csv"
python3 "$HERE/merge_pieces.py" plan "$WORK/block_files.txt" "$LIB" "$WORK/pieces.csv" "$PER"
while IFS=, read -r det j; do
  python3 "$HERE/merge_pieces.py" piece "$WORK/block_files.txt" "$LIB" "$det" "$j"
done < "$WORK/pieces.csv"
python3 "$HERE/make_support.py" cells "$WORK/block_files.txt" "$LIB"
python3 "$HERE/check_library.py" "$LIST" "$LIB"
