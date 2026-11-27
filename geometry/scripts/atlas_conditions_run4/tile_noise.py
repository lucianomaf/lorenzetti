"""The Tile noise per sample for the Phase-II electronics (FENICS), per ring (section, sampling, tower; both sides, every phi, both PMTs) and
gain, from public sources only. Output: WORK/tables/tile_run4_noise_by_ring.csv.

Source: the sample noise of the ATLAS Run 4 simulation, COOL folder /TILE/OFL02/NOISE/SAMPLE (COOLOFL_TILE/OFLP200), the tag of the global
tag OFLCOND-MC21-SDR-RUN4-02 (sources.json), run 350000, through the public Frontier. Fields of a channel and gain (open Athena 25.0.21,
TileConditions/TileConditions/TileSampleNoise.h): 0 pedestal, 1 HFN, 2 LFN, 3 HFN1, 4 HFN2, 5 HFNNORM. The ATLAS digitisation
(TileSimAlgs/src/TileDigitsMaker.cxx) draws the noise of a sample from a double Gaussian: sigma HFN1 with probability
w = HFN1 / (HFN1 + HFN2 * HFNNORM), sigma HFN2 otherwise; a channel without HFN1 and HFN2 uses HFN with w = 1. These are values of the present
electronics: 10-bit ADC, high gain up to 12.5 pC (800 pC / 64; ATLAS-TDR-028, pp. 95-96), low gain up to 800 pC.

Conversion to the FENICS counts, keeping the charge (decisions of 21/11/2026): the quantisation of the 10-bit ADC (1/sqrt(12) count) is
taken out of each sigma in quadrature (the rounding of the new ADC adds its own); then each sigma is multiplied by the ratio of the charge
per count, legacy / FENICS. FENICS: 12-bit ADC (ATLAS-TDR-028, p. 96), gains 0.4 and 16, ratio 40, low gain up to 1050 pC (PoS(ICHEP2020)860,
p. 4), so the high gain is 1050 / 40 = 26.25 pC full scale: high gain (12.5 / 1023) / (26.25 / 4095), low gain (800 / 1023) / (1050 / 4095).
Per ring, each parameter (sigma 1, sigma 2, w) is the median over the PMTs, leaving out the dead channels (HFN = 0).
Usage: python3 tile_noise.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sources import WORK, load  # noqa: E402
from tile import channels, drawers  # noqa: E402

FACTOR = {1: (12.5 / 1023) / (26.25 / 4095), 0: (800.0 / 1023) / (1050.0 / 4095)}   # TileID: 0 low gain, 1 high gain
QUANT = 1.0 / 12   # variance of the rounding of the 10-bit ADC (counts^2)
src = load()["tile_noise"]
blobs = drawers(src["folder"], src["tag"], src["run"])
print(f"{src['folder']}, tag {src['tag']}, run {src['run']}: {len(blobs)} COOL channels, {sum(v is not None for v in blobs.values())} with data")


def params(f):
    """(pedestal, sigma 1, sigma 2, w) of one channel and gain in 10-bit counts, as TileDigitsMaker; None for a dead channel"""
    f = list(f) + [0.0] * (6 - len(f))
    ped, hfn, _lfn, h1, h2, norm = f[:6]
    if hfn <= 0:
        return None
    if h1 > 0 or h2 > 0:
        return ped, h1, h2, h1 / (h1 + h2 * norm)
    return ped, hfn, 0.0, 1.0


rings, dead = {}, {}
for _ros, _dr, _ch, cell, v in channels(blobs):
    for g in range(v.shape[0]):
        p = params(v[g])
        if p is None:
            dead[(cell, g)] = dead.get((cell, g), 0) + 1
        else:
            rings.setdefault((cell, g), []).append(p)
os.makedirs(os.path.join(WORK, "tables"), exist_ok=True)
out = os.path.join(WORK, "tables", "tile_run4_noise_by_ring.csv")
with open(out, "w", newline="\n") as o:
    o.write("section,sampling,tower,gain,n_pmts,n_dead,sigma1_counts,sigma2_counts,w_sigma1,"
            "pedestal_legacy_counts,sigma1_legacy_counts,sigma2_legacy_counts\n")
    for (cell, g) in sorted(rings):
        x = np.array(rings[(cell, g)])
        ped, s1, s2, w = np.median(x, axis=0)
        f1, f2 = (FACTOR[g] * np.sqrt(max(s * s - QUANT, 0.0)) if s > 0 else 0.0 for s in (s1, s2))
        o.write(f"{cell[0]},{cell[1]},{cell[2]},{'high' if g == 1 else 'low'},{len(x)},{dead.get((cell, g), 0)},{f1:.5g},{f2:.5g},{w:.5g},"
                f"{ped:.5g},{s1:.5g},{s2:.5g}\n")
print(f"{len(rings)} (ring, gain) -> tables/tile_run4_noise_by_ring.csv")
