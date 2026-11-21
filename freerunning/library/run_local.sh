#!/bin/bash
# Builds a free-running library on one machine, all steps in sequence (the same scripts the HTCondor jobs run), inside the Lorenzetti
# image (ROOT and numpy; no other package). See README.md.
# Usage: bash run_local.sh <input list> <library folder> <work folder> [hits per piece]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
LIST=$1; LIB=$2; WORK=$3; PER=${4:-160000000}
mkdir -p "$LIB" "$WORK/blocks"
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
