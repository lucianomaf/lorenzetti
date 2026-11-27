"""Common pieces of the Tile scripts: the channel -> cell map, the default drawers and the blob format of the ATLAS Tile conditions, from
the open Athena code of release 25.0.21 (gitlab.cern.ch/atlas/athena, branch release/25.0.21, readable without login):
TileCalibBlobObjs/TileCalibDrawerBase.h (blob format), TileCalibUtils.cxx (getDrawerIdx, getDefaultDrawerIdx), TileCalibDrawerDat.h (a
channel beyond the stored ones -> 0), TileConditions/src/TileCondProxyCool.icc (empty drawer -> default drawer) and
TileConditions/src/TileCablingService.cxx (hwid2section, hwid2sample, hwid2tower, with the E-cell channels of the Run 2 and Run 3 cabling:
E1 12, E2 13, E3 0, E4 1). The special drawers (without C10 or D4) are not treated; the median of a ring does not change.
"""
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sources import payload  # noqa: E402

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



def drawers(folder, tag, run):
    """{COOL channel = drawer index: (channels, gains, values)} of a Tile COOL folder of the simulation (COOLOFL_TILE/OFLP200) at a run,
    through the public Frontier"""
    return {c: read_blob(b.encode("latin-1")) for c, b in payload("COOLOFL_TILE", "OFLP200", folder, tag, run, "TileCalibBlob").items()}


def channels(blobs):
    """(ros, drawer, channel, (section, sampling, tower), values of the channel: array gains x fields) for every Tile channel of a cell"""
    def drawer_blob(idx):
        while blobs.get(idx) is None:
            idx = default_drawer(idx)
        return blobs[idx]
    for ros in (1, 2, 3, 4):
        for dr in range(MAX_DRAWER):
            nch, ng, v = drawer_blob(MAX_DRAWR0 + (ros - 1) * MAX_DRAWER + dr)
            for ch in range(MAX_CHAN):
                cell = cell_of(ros, ch)
                if cell is None:
                    continue
                c = ch if ch < nch else 0
                yield ros, dr, ch, cell, v[c * ng:(c + 1) * ng]
