"""The Tile sampling fraction of the ATLAS Run 4 simulation per ring (section, sampling, tower; both sides and both PMTs together), from
public sources only. Output: WORK/tables/tile_run4_sampling_fraction_by_ring.csv (median over the PMTs of the ring and how many are more than
5% away from it).

Sources: the simulation COOL folder /TILE/OFL02/CALIB/SFR (COOLOFL_TILE/OFLP200), tag TileOfl02CalibSfr-SIM-07 (the tag Athena 25.0.21
forces for the Run 4 geometry, TileConditions/python/TileSamplingFractionConfig.py), run 350000, one blob per drawer, through the public
Frontier; and the open Athena code of release 25.0.21 (gitlab.cern.ch/atlas/athena, branch release/25.0.21, readable without login):
blob format (TileCalibBlobObjs/TileCalibDrawerBase.h), drawer index and default drawers (TileCalibUtils.cxx: getDrawerIdx,
getDefaultDrawerIdx), channel policy (TileCalibDrawerDat.h), empty drawer -> default drawer (TileConditions/src/TileCondProxyCool.icc) and the
channel -> cell map (TileConditions/src/TileCablingService.cxx: hwid2section, hwid2sample, hwid2tower, with the E-cell channels of the Run 2
and Run 3 cabling: E1 12, E2 13, E3 0, E4 1). The special drawers (without C10 or D4) are not treated; the median of a ring does not change.
Usage: python3 tile_sfr.py
"""
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frontier import q  # noqa: E402
from sources import WORK  # noqa: E402

S = "ATLAS_COOLOFL_TILE.OFLP200"
TAG, RUN = "TileOfl02CalibSfr-SIM-07", 350000
MAX_DRAWR0, MAX_DRAWER, MAX_CHAN = 20, 64, 48
E1, E2, E3, E4 = 12, 13, 0, 1

TWR = [0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 4, 4, 5, 5, 5, 6, -1, -1, 6, 6, 6, 7, 8, 8, 7, 7, 7, 6, 8, -1, 6, 9, 9, 8,
       15, 13, 8, 8, 9, 9, 11, 11, 10, 10, 12, 12, 11, 10, 11, 11, 10, 10, -1, -1, 13, 13, 12, 12, -1, -1, -1, -1, -1, -1, 13, 14, 14, -1, -1, 13,
       14, 12, 12, 14, 15, 15, -1, -1, -1, -1, -1, -1]
SMP = [2, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 2, 2, 0, 1, 1, 0, 0, 0, 1, 1, 0, 2, 2, 0, 1, 1, 0, -1, -1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 2, 1, -1, 2, 0, 0, 1,
       3, 3, 2, 2, 1, 1, 0, 0, 1, 1, 0, 0, 3, 3, 1, 1, 2, 2, -1, -1, 0, 0, 1, 1, -1, -1, -1, -1, -1, -1, 1, 0, 0, -1, -1, 1, 1, 2, 2, 1, 0, 0, -1, -1, -1,
       -1, -1, 4]
TWR[48 + E1], TWR[48 + E2], TWR[48 + E3], TWR[48 + E4] = 10, 11, 13, 15


def cell_of(ros, ch):
    """(section, sampling, tower) of a channel: TileCablingService::hwid2section, hwid2sample, hwid2tower (1 barrel, 2 extended, 3 gap)"""
    k = ch + 48 if ros > 2 else ch
    smp, twr = SMP[k], TWR[k]
    if smp < 0 or twr < 0 or smp == 4:
        return None
    sec = (3 if (ch < 6 or ch in (12, 13)) else 2) if ros > 2 else 1
    return sec, smp, twr


def default_drawer(idx):
    """TileCalibUtils::getDefaultDrawerIdx"""
    if idx < MAX_DRAWR0:
        if idx < 4:
            return 0
        mod = (idx - 4) % 8
        return idx - mod if mod > 0 else 0
    i = idx - MAX_DRAWR0
    if i < 2 * MAX_DRAWER:
        return 4
    i -= 2 * MAX_DRAWER
    eba = [0, 0, 0, 0, 0, 0, 3, 2, 0, 0, 0, 0, 7, 6, 5, 7, 7, 6, 6, 7, 0, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    ebc = [0, 0, 0, 0, 0, 0, 3, 2, 0, 0, 0, 0, 7, 6, 6, 7, 7, 5, 6, 7, 0, 0, 0, 2, 3, 0, 0, 3, 4, 0, 3, 4,
           0, 4, 3, 0, 4, 3, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    if i < MAX_DRAWER:
        return 12 + eba[i]
    return 12 + ebc[i - MAX_DRAWER]


def read_blob(b):
    """(channels, gains, values) of a TileCalibDrawerFlt blob"""
    if len(b) < 20:
        return None
    _tipo, _versao = struct.unpack("<2H", b[:4])
    tam, nobj = struct.unpack("<2I", b[4:12])
    nch, ng = struct.unpack("<2H", b[12:16])
    v = np.frombuffer(b, "<f4", count=tam * nobj, offset=20).reshape(nobj, tam)
    return nch, ng, v


nid = q(f"select TO_CHAR(NODE_ID) from {S}_NODES where NODE_FULLPATH='/TILE/OFL02/CALIB/SFR'")[0][0]
t = f"{S}_F{int(nid):04d}"
tid = q(f"select TO_CHAR(TAG_ID) from {t}_TAGS where TAG_NAME='{TAG}'")[0][0]
k = RUN << 32
blobs = {int(c): read_blob(b.encode("latin-1")) for c, b in
         q(f'select TO_CHAR(CHANNEL_ID), "TileCalibBlob" from {t}_IOVS where USER_TAG_ID={tid} and IOV_SINCE<={k} and IOV_UNTIL>{k}')}
print(f"{TAG}, run {RUN}: {len(blobs)} COOL channels, {sum(v is not None for v in blobs.values())} with data")


def drawer_blob(idx):
    while blobs.get(idx) is None:
        idx = default_drawer(idx)
    return blobs[idx]


rings = {}
for ros in (1, 2, 3, 4):
    for dr in range(MAX_DRAWER):
        nch, ng, v = drawer_blob(MAX_DRAWR0 + (ros - 1) * MAX_DRAWER + dr)
        for ch in range(MAX_CHAN):
            c_ = cell_of(ros, ch)
            if c_ is None:
                continue
            c = ch if ch < nch else 0  # TileCalibDrawerDat: a channel beyond the stored ones -> 0
            rings.setdefault(c_, []).append(float(v[c * ng + 0][0]))  # high gain, first value
os.makedirs(os.path.join(WORK, "tables"), exist_ok=True)
with open(os.path.join(WORK, "tables", "tile_run4_sampling_fraction_by_ring.csv"), "w", newline="\n") as o:
    o.write("section,sampling,tower,n_pmts,sampling_fraction,sampling_fraction_n_off5pct\n")
    for (sec, smp, twr) in sorted(rings):
        x = np.array(rings[(sec, smp, twr)])
        m = np.median(x)
        o.write(f"{sec},{smp},{twr},{len(x)},{m:.6g},{int(np.sum(np.abs(x - m) > 0.05 * abs(m)))}\n")
print(f"{len(rings)} rings -> tables/tile_run4_sampling_fraction_by_ring.csv")
