"""The Tile sampling fraction of the ATLAS Run 4 simulation per ring (section, sampling, tower; both sides and both PMTs together), from
public sources only. Output: WORK/tables/tile_run4_sampling_fraction_by_ring.csv (median over the PMTs of the ring and how many are more than
5% away from it).

Source: the simulation COOL folder /TILE/OFL02/CALIB/SFR (COOLOFL_TILE/OFLP200), tag TileOfl02CalibSfr-SIM-07 (the tag Athena 25.0.21
forces for the Run 4 geometry, TileConditions/python/TileSamplingFractionConfig.py), run 350000, one blob per drawer, through the public
Frontier. Blob format, drawers and channel -> cell map: tile.py.
Usage: python3 tile_sfr.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sources import WORK  # noqa: E402
from tile import channels, drawers  # noqa: E402

TAG, RUN = "TileOfl02CalibSfr-SIM-07", 350000
blobs = drawers("/TILE/OFL02/CALIB/SFR", TAG, RUN)
print(f"{TAG}, run {RUN}: {len(blobs)} COOL channels, {sum(v is not None for v in blobs.values())} with data")
rings = {}
for _ros, _dr, _ch, cell, v in channels(blobs):
    rings.setdefault(cell, []).append(float(v[0][0]))  # high gain, first value
os.makedirs(os.path.join(WORK, "tables"), exist_ok=True)
with open(os.path.join(WORK, "tables", "tile_run4_sampling_fraction_by_ring.csv"), "w", newline="\n") as o:
    o.write("section,sampling,tower,n_pmts,sampling_fraction,sampling_fraction_n_off5pct\n")
    for (sec, smp, twr) in sorted(rings):
        x = np.array(rings[(sec, smp, twr)])
        m = np.median(x)
        o.write(f"{sec},{smp},{twr},{len(x)},{m:.6g},{int(np.sum(np.abs(x - m) > 0.05 * abs(m)))}\n")
print(f"{len(rings)} rings -> tables/tile_run4_sampling_fraction_by_ring.csv")
