#!/bin/bash
# The package of the build for the HTCondor jobs (block_job.sh, PROD_PACKAGE): a .tgz of the home folder where the build was made, with the
# repository built (without .git), the GeoModel libraries and the Run 4 geometry dump, at the paths the build keeps (see ../README.md). The
# jobs mount it at that same home folder (PROD_HOME). Run on the machine of the build, from anywhere.
# Usage: bash make_package.sh <output .tgz> [home folder (default $HOME)] [repository folder in it (lorenzetti-geo)] [GeoModel folder
#        (geomodel)] [geometry dump folder (atlas_geometry)]
# It stops if a file the jobs need is missing: the setup of the build, the EMEC field map and the FCal shower libraries in geometry/data
# (../../../README.md, "Version 2"), the GeoModel libraries and the geometry dump. The ATLAS field map is not in the package: the jobs
# read it from CVMFS (PROD_FIELD_MAP). Then copy it to EOS and check it there (md5 of a copy brought back with xrdcp; just after xrdcp the
# /eos mount may not see the file yet), e.g.
#   xrdcp -f -s lzt_build.tgz root://eosuser.cern.ch//eos/.../lzt_build.tgz
#   xrdcp -f -s root://eosuser.cern.ch//eos/.../lzt_build.tgz /tmp/back.tgz && md5sum lzt_build.tgz /tmp/back.tgz
set -eu
OUT=$(realpath -m "$1"); H=${2:-$HOME}; REPO=${3:-lorenzetti-geo}; GM=${4:-geomodel}; GEO=${5:-atlas_geometry}
cd "$H"
need() { [ -e "$1" ] || { echo "missing: $H/$1"; exit 1; }; }
need $REPO/build/lzt_setup.sh
need $REPO/geometry/data/ec.fieldzv03
for f in FCAL1.11 FCAL1.22 FCAL1.2112 FCAL2.11 FCAL2.22 FCAL2.2112; do need $REPO/geometry/data/LArG4ShowerLibData/MC23_v1/LArG4ShowerLib.$f.root; done
need $GM/lib
need $GEO/geometry-ATLAS-P2-RUN4-01-01-00-sim.db
echo "commit $(git -C $REPO rev-parse --short HEAD)$(git -C $REPO diff --quiet HEAD -- . ':!geometry/data/ec.fieldzv03' || echo ', with local changes')"
tar czf "$OUT" --exclude=$REPO/.git $REPO $GM $GEO
echo "$OUT: $(du -h "$OUT" | cut -f1), md5 $(md5sum < "$OUT" | cut -c1-32)"
