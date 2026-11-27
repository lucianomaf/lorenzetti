"""Every channel of the free-running hits of the v2 detector, from the configuration of the detector, without simulation.

The detector is built as simu_trf.py --detector v2 builds it with the options of ../production/block.sh, and each channel gets the
identifier, eta and phi that its hits carry, by the formulas of CaloHitMaker (reconstruction/calorimeter/CaloHitBuilder/src):
  tile  freeRunningTile and tileNPmt: one channel per PMT; module from the centre of the phi bin, phi the centre of the module; one PMT in
        the cells with NPmt = 1 or -1 (E1-E4, the special C10); the D4 of the special modules (NPmt = 0) is read as the D5 of its module,
        whose eta is taken (the larger |eta|, as block_to_channels.py)
  LAr   freeRunningLAr: EMB, EMEC and HEC, one channel per cell and phi bin, eta of the cell and centre of the phi bin
  fcal  fcalStep and fcalSpot: one channel per tile of each module, eta and phi of the tile
A change in those formulas must be made here too; make_support.py checks that every channel with hits is in this list, with the same
eta and phi, so a divergence stops the library.
Needs the Lorenzetti environment (lzt_setup.sh) and the external files of the v2 detector: LZT_RUN4_GEOMETRY and LZT_FIELD_MAP (see
../production/README.md).
Usage: python3 geometry_channels.py <output .npz>   (arrays <subdetector>_channel, _eta, _phi, as make_support.py uses them)
"""
import math
import os
import sys

import numpy as np

TWO_PI = 2 * math.pi
MOD_WIDTH = TWO_PI / 64


def detector():
    """the v2 detector with the options of ../production/block.sh"""
    from geometry import DetectorConstruction_v2
    return DetectorConstruction_v2(
        "ATLAS", UseMagneticField=False, UseSolenoidField=True, TileAtlasGeometry=True, TileAtlasCells=True, AtlasMaterialInFront=True,
        TileAtlasItc=True, AtlasEmb=True, AtlasEmbCells=True, AtlasHec=True, AtlasEmec=True, AtlasFcal=True, AtlasFcalCells=True,
        AtlasFcalElectrodes=os.path.join(os.environ["LORENZETTI_GEOMETRY_DATA_DIR"], "LArFCalElectrodes-06.csv"),
        AtlasFcalFrozenShowers=True, AtlasG4=True, TileAtlasGap=True, AtlasEndcapCryostat=True, AtlasBarrelCryostat=True,
        AtlasRun4=True, AtlasRun4Geometry=os.environ["LZT_RUN4_GEOMETRY"], AtlasForwardShielding=True, AtlasFullField=True,
        AtlasBarrelGapServices=True, TileAtlasRun4=True, AtlasEmecArgonFraction=True, AtlasFieldMap=os.environ["LZT_FIELD_MAP"])


def tile_channels(cells, phi_bins, out):
    """the PMTs of one tile sampling: {pmt_id: (eta, phi)}"""
    special = dict(zip(cells.get("SpecialModule", []), cells.get("SpecialNPmt", [])))
    for i, eta in enumerate(cells["Eta"]):
        for b in range(len(phi_bins) - 1):
            centre = 0.5 * (phi_bins[b] + phi_bins[b + 1])
            module = int(math.floor((centre + TWO_PI if centre < 0 else centre) / MOD_WIDTH)) % 64
            npmt = special.get(module, cells["NPmt"][i] if "NPmt" in cells else 2)
            if npmt == 0:
                sec, twr, smp, z = cells["ReplaceSection"], cells["ReplaceTower"], cells["ReplaceSampling"], cells["ReplaceZCentre"]
            else:
                sec, twr, smp, z = cells["AtlasSection"][i], cells["AtlasTower"][i], cells["AtlasSampling"][i], cells["ZCentre"][i]
            side = 1 if (z >= 0 or (sec == 1 and smp == 2 and twr == 0)) else -1
            phi = math.remainder((module + 0.5) * MOD_WIDTH, TWO_PI)
            for pmt in {2: (0, 1), 0: (0, 1), 1: (0,), -1: (1,)}[npmt]:
                x = ((2 << 61) | (int(sec) << 58) | ((2 if side > 0 else 0) << 54) | (module << 46) | (int(twr) << 40) | (int(smp) << 36)
                     | (pmt << 34))
                old = out.get(x)
                if old is None or abs(eta) > abs(old[0]):
                    out[x] = (float(eta), phi)


def lar_channels(cells, phi_bins, out):
    """the cells of one EMB, EMEC or HEC sampling: {subdetector: {cell_id: (eta, phi)}}"""
    n = len(phi_bins) - 1
    be, smp, rg = int(cells["AtlasBarrelEc"]), int(cells["AtlasLArSampling"]), int(cells["AtlasLArRegion"])
    if int(cells["AtlasLArPart"]) == 1:
        sub = "emb" if abs(be) == 1 else "emec"
        base = (1 << 29) | (3 << 26) | ([-3, -2, -1, 1, 2, 3].index(be) << 23) | (smp << 21) | (rg << 18)
        field = lambda e, p: (e << 9) | (p << 1)  # noqa: E731
    else:
        sub = "hec"
        base = (1 << 29) | (4 << 26) | ((1 if be > 0 else 0) << 25) | (smp << 23) | (rg << 22)
        field = lambda e, p: (e << 18) | (p << 12)  # noqa: E731
    for i, eta in enumerate(cells["Eta"]):
        e = int(cells["AtlasEtaIndex"][i])
        for b in range(n):
            phi = math.remainder(0.5 * (phi_bins[b] + phi_bins[b + 1]), TWO_PI)
            out[sub][base | field(e, (b + n // 2) % n)] = (float(eta), phi)


def fcal_channels(fm, out):
    """the tiles of one FCal module; returns the largest difference of eta or phi between the electrodes of a tile"""
    base = (1 << 29) | (5 << 26) | ((1 if fm.Side > 0 else 0) << 25) | ((int(fm.Module) - 1) << 23)
    dev = 0.0
    for e, p, eta, phi in zip(fm.CopyEtaIndex, fm.CopyPhiIndex, fm.CopyEta, fm.CopyPhi):
        x = base | (int(e) << 17) | (int(p) << 13)
        if x in out:
            dev = max(dev, abs(out[x][0] - eta), abs(out[x][1] - phi))
        else:
            out[x] = (float(eta), float(phi))
    return dev


def channels():
    """{subdetector: (channel uint64, eta, phi)}, sorted by channel"""
    d = detector()
    got = {k: {} for k in ("tile", "emb", "emec", "hec", "fcal")}
    for s in d.samplings:
        sv = s.sensitive()
        cells = getattr(sv, "Cells", None) or {}
        if "AtlasSection" in cells:
            tile_channels(cells, list(sv.PhiBins), got["tile"])
        elif "AtlasLArPart" in cells:
            lar_channels(cells, list(sv.PhiBins), got)
        else:
            raise SystemExit(f"ERROR: sampling {s.name()} has neither tile nor LAr identifier fields")
    dev = max(fcal_channels(fm, got["fcal"]) for fm in d.fcal_modules)
    if dev > 1e-6:
        raise SystemExit(f"ERROR: the electrodes of an FCal tile have different eta or phi ({dev})")
    out = {}
    for k, v in got.items():
        x = np.array(sorted(v), dtype=np.uint64)
        out[k] = (x, np.array([v[int(i)][0] for i in x]), np.array([v[int(i)][1] for i in x]))
    return out


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    ch = channels()
    np.savez(sys.argv[1], **{f"{k}_{n}": a for k, v in ch.items() for n, a in zip(("channel", "eta", "phi"), v)})
    print(", ".join(f"{k} {len(v[0])}" for k, v in ch.items()))
