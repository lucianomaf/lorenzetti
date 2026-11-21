"""The Tile energy scale and photoelectron yield of the ATLAS Run 4 simulation per ring (section, sampling, tower; both sides and both PMTs
together), from public sources only. Output: WORK/tables/tile_run4_scale_by_ring.csv.

How Athena 25.0.21 turns the energy of a Tile hit into ADC counts in the simulation (TileSimAlgs/src/TileDigitsMaker.cxx): counts = energy
in the scintillator x sampling fraction / (MeV per count), with MeV per count = TileEMScale::calibrateChannel(1 count, ADCcounts ->
MegaElectronVolts) (TileConditions/src/TileEMScale.cxx): the charge calibration (CIS), the laser and the cesium corrections and the EM scale, each
a slope of the channel that divides the amplitude (TileCalibBlobObjs/src/TileCalibDrawerFlt.cxx, getCalib, version 1) from the simulation COOL folders that TileConditions/python/TileEMScaleConfig.py
reads for MC: /TILE/OFL02/CALIB/CIS/FIT/LIN, /TILE/OFL02/CALIB/LAS/LIN, /TILE/OFL02/CALIB/CES, /TILE/OFL02/CALIB/EMS (COOLOFL_TILE/OFLP200,
the tags of the global tag OFLCOND-MC21-SDR-RUN4-02, run 350000). The non-linear corrections (CIS/FIT/NLN, LAS/NLN) are read too and must be
the identity in the simulation (checked here). These are the values of the present electronics (10-bit ADC); the Phase-II scale is derived
from them in the generator (NOTICE.md). The photoelectrons per GeV are the second value of the sampling-fraction folder
/TILE/OFL02/CALIB/SFR (tag TileOfl02CalibSfr-SIM-07), which TileHitVecToCntTool.cxx uses for the photoelectron statistics
(TileSamplingFraction::getNumberOfPhotoElectrons). The gain index of the blobs is 0 for the low gain and 1 for the high gain.
Columns: median over the PMTs of the ring of each factor, of the MeV per count of each gain and of the photoelectrons per GeV, and how many
PMTs are more than 5% away from the median of the MeV per count of the high gain.
Usage: python3 tile_scale.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sources import MC_GLOBAL, MC_RUN, WORK, tag_in_global  # noqa: E402
from tile import channels, drawers  # noqa: E402

LINEAR = {"cis": "/TILE/OFL02/CALIB/CIS/FIT/LIN", "las": "/TILE/OFL02/CALIB/LAS/LIN", "ces": "/TILE/OFL02/CALIB/CES",
          "ems": "/TILE/OFL02/CALIB/EMS"}
NONLINEAR = {"cis_nln": "/TILE/OFL02/CALIB/CIS/FIT/NLN", "las_nln": "/TILE/OFL02/CALIB/LAS/NLN"}
SFR_TAG = "TileOfl02CalibSfr-SIM-07"


def per_channel(folder, tag):
    """{(ros, drawer, channel): values (gains x fields)} and the cell of each channel"""
    out, cell = {}, {}
    for ros, dr, ch, c, v in channels(drawers(folder, tag, MC_RUN)):
        out[(ros, dr, ch)], cell[(ros, dr, ch)] = np.asarray(v, float), c
    return out, cell


tags, vals, cell = {}, {}, None
for k, folder in {**LINEAR, **NONLINEAR}.items():
    tags[k] = tag_in_global("COOLOFL_TILE", "OFLP200", folder, MC_GLOBAL)
    vals[k], cell = per_channel(folder, tags[k])
    shapes = {v.shape for v in vals[k].values()}
    print(f"{folder}: tag {tags[k]}, {len(vals[k])} channels, values per channel (gains x fields) {sorted(shapes)}")
for k in NONLINEAR:  # TileCalibDrawerFlt::getCalib, look-up table (version 100): fewer than 2 (x, y) pairs is the identity
    n = max(v.shape[1] for v in vals[k].values())
    assert n < 4, f"{NONLINEAR[k]}: a look-up table of {n} values; it is not the identity"
    print(f"{NONLINEAR[k]}: {n} values per channel and gain, fewer than 2 (x, y) pairs: the identity")
sfr, _ = per_channel("/TILE/OFL02/CALIB/SFR", SFR_TAG)
rings = {}
for key, c in cell.items():
    lin = {k: vals[k][key] for k in LINEAR}
    # TileCalibDrawerFlt::getCalib, single slope (version 1): the slope is the first value of the channel and gain, and the amplitude
    # is divided by it; one gain stored is the same for both
    gain = lambda v, g: float(v[g if v.shape[0] > 1 else 0, 0])  # noqa: E731
    row = {"cis_hg": gain(lin["cis"], 1), "cis_lg": gain(lin["cis"], 0), "las": gain(lin["las"], 1), "ces": gain(lin["ces"], 1),
           "ems": gain(lin["ems"], 1), "npe": float(sfr[key][0, 1])}
    row["mev_hg"] = 1.0 / (row["cis_hg"] * row["las"] * row["ces"] * row["ems"])
    row["mev_lg"] = 1.0 / (row["cis_lg"] * gain(lin["las"], 0) * gain(lin["ces"], 0) * gain(lin["ems"], 0))
    rings.setdefault(c, []).append(row)
os.makedirs(os.path.join(WORK, "tables"), exist_ok=True)
cols = ["cis_hg", "cis_lg", "las", "ces", "ems", "mev_hg", "mev_lg", "npe"]
names = ["cis_slope_hg", "cis_slope_lg", "laser_slope", "cesium_slope", "ems_slope", "mev_per_count_legacy_hg",
         "mev_per_count_legacy_lg", "photoelectrons_per_gev"]
with open(os.path.join(WORK, "tables", "tile_run4_scale_by_ring.csv"), "w", newline="\n") as o:
    o.write("section,sampling,tower,n_pmts," + ",".join(names) + ",mev_per_count_legacy_hg_n_off5pct\n")
    for k in sorted(rings):
        med = {c_: float(np.median([r[c_] for r in rings[k]])) for c_ in cols}
        off = sum(abs(r["mev_hg"] - med["mev_hg"]) > 0.05 * abs(med["mev_hg"]) for r in rings[k])
        o.write(f"{k[0]},{k[1]},{k[2]},{len(rings[k])}," + ",".join(f"{med[c_]:.6g}" for c_ in cols) + f",{off}\n")
allm = np.array([r["mev_hg"] for v in rings.values() for r in v])
print(f"{len(rings)} rings; MeV per count, high gain: median {np.median(allm):.4g} (min {allm.min():.4g}, max {allm.max():.4g})")
print("tags: " + ", ".join(f"{k} {t}" for k, t in tags.items()) + f", sfr {SFR_TAG}")
