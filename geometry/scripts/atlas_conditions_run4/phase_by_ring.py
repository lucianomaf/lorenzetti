"""The sampling phase per ring, from the free-running hits of Lorenzetti itself: the time, after the time of flight, at which the energy of
fast particles arrives in each ring of cells, so that the samples of the sequence fall on the peak of the pulse for high signals.

The time of a hit in the free-running hits (simu_trf.py --free-running-hits) is already that of the ATLAS simulation, the global time with the
time of flight from the interaction point taken out (CaloHitMaker.cxx: larHitTime, tileTime); a particle at the speed of light that deposits
at the centre of a cell is at t = 0. The energy of a shower is spread over the cell and a little late, so the energy of high signals arrives
at a ring-dependent time near 0. Per ring (LAr: the rings of the LAr tables; Tile: section, sampling and tower; all phi, both sides, both
PMTs), the deposits are put in a histogram of time weighted by energy (0.25 ns in the LAr, 0.5 ns in the Tile, from -12.5 to +12.5 ns); the
phase is the energy-weighted mean time of the bins within 3 ns of the highest bin. The short window keeps the late tail of the energy out of
the mean. Use the hard sample (events with a jet or a photon of high transverse momentum, whose deposits are dominated by fast particles as a
signal); the same on the inclusive sample shows how much the phase depends on the sample.
Output: WORK/tables/lar_phase_by_ring<suffix>.csv and tile_phase_by_ring<suffix>.csv.
Usage: python3 phase_by_ring.py <suffix, '' for the tables of the repository> <FRHITS file> [<FRHITS file> ...]
"""
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lar  # noqa: E402
from sources import WORK  # noqa: E402

SUFFIX, FILES = sys.argv[1], sys.argv[2:]
TMAX, WINDOW = 12.5, 3.0
STEP = {"lar": 0.25, "tile": 0.5}
EDGES = {k: np.arange(-TMAX, TMAX + s / 2, s) for k, s in STEP.items()}


class Acc:
    """per ring: energy and energy x time per time bin, and the total energy of the ring at any time"""
    def __init__(self, kind):
        self.kind, self.e, self.et, self.tot, self.n = kind, {}, {}, {}, {}

    def add(self, keys, e, t):
        nb = len(EDGES[self.kind]) - 1
        b = np.digitize(t, EDGES[self.kind]) - 1
        inw = (b >= 0) & (b < nb)
        uk, inv = np.unique(keys, axis=0, return_inverse=True)
        inv = inv.ravel()
        for j, k in enumerate(map(tuple, uk)):
            m = inv == j
            self.tot[k] = self.tot.get(k, 0.0) + float(e[m].sum())
            self.n[k] = self.n.get(k, 0) + int(m.sum())
            mw = m & inw
            if k not in self.e:
                self.e[k], self.et[k] = np.zeros(nb), np.zeros(nb)
            np.add.at(self.e[k], b[mw], e[mw])
            np.add.at(self.et[k], b[mw], e[mw] * t[mw])

    def phase(self, k):
        e, et = self.e[k], self.et[k]
        if e.sum() <= 0:
            return np.nan, 0.0
        c = 0.5 * (EDGES[self.kind][:-1] + EDGES[self.kind][1:])
        w = np.abs(c - c[np.argmax(e)]) <= WINDOW
        return et[w].sum() / e[w].sum(), e[w].sum() / self.tot[k]


lacc, tacc = Acc("lar"), Acc("tile")
for i, fn in enumerate(FILES):
    t = uproot.open(fn)["CollectionTree"]
    for d in ("EMB", "EMEC", "HEC", "FCAL"):
        z = {n.split("_", 1)[1]: np.concatenate(v) for n, v in
             t.arrays([f"LArHit{d}_{x}" for x in ("cellID", "energy", "time", "eta")], library="np").items()}
        if not len(z["cellID"]):
            continue
        ids = (z["cellID"].astype(np.uint64) >> np.uint64(32)).astype(np.int64)
        uid, inv = np.unique(ids, return_inverse=True)
        keys = np.empty((len(uid), 2), dtype=object)
        for j, x in enumerate(uid.tolist()):
            dd, _s, b, smp, rg, e, _p = lar.campos(int(x))
            keys[j] = (lar.layer(dd, b, smp), 0 if dd == "FCAL" else b * 10000 + rg * 1000 + e)
        keys = keys[inv.ravel()]
        if d == "FCAL":
            keys[:, 1] = 1000 + np.floor(np.abs(z["eta"]) / 0.1).astype(int)
        lacc.add(keys.astype(str), z["energy"].astype(float), z["time"].astype(float))
    z = {n.split("_", 1)[1]: np.concatenate(v) for n, v in
         t.arrays([f"TileCalHit_{x}" for x in ("cellID", "energy", "time")], library="np").items()}
    cid = z["cellID"].astype(np.uint64)
    sec = ((cid >> np.uint64(58)) & np.uint64(7)).astype(int)
    twr = ((cid >> np.uint64(40)) & np.uint64(63)).astype(int)
    smp = ((cid >> np.uint64(36)) & np.uint64(15)).astype(int)
    tacc.add(np.stack([sec, smp, twr], axis=1), z["energy"].astype(float), z["time"].astype(float))
    if i % 50 == 0:
        print(f"{i + 1}/{len(FILES)} files", flush=True)
os.makedirs(os.path.join(WORK, "tables"), exist_ok=True)
with open(os.path.join(WORK, "tables", f"lar_phase_by_ring{SUFFIX}.csv"), "w", newline="\n") as o:
    o.write("layer,ring,n_hits,energy_fraction_in_window,phase_ns\n")
    for k in sorted(lacc.e, key=lambda k: (k[0], int(k[1]))):
        ph, fr = lacc.phase(k)
        o.write(f"{k[0]},{k[1]},{lacc.n[k]},{fr:.4f},{ph:.4f}\n")
with open(os.path.join(WORK, "tables", f"tile_phase_by_ring{SUFFIX}.csv"), "w", newline="\n") as o:
    o.write("section,sampling,tower,n_hits,energy_fraction_in_window,phase_ns\n")
    for k in sorted(tacc.e, key=lambda k: tuple(int(x) for x in k)):
        ph, fr = tacc.phase(k)
        o.write(f"{k[0]},{k[1]},{k[2]},{tacc.n[k]},{fr:.4f},{ph:.4f}\n")
print(f"{len(lacc.e)} LAr rings, {len(tacc.e)} Tile rings -> tables/*_phase_by_ring{SUFFIX}.csv")
