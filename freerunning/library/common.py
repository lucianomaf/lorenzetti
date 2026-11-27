"""Common pieces of the scripts that build the free-running library: the calorimeter hits of a minimum-bias production of Lorenzetti
(simu_trf.py --free-running-hits, one FRHITS file per block of events) reorganised by channel, so that a sequence generator reads, for
one channel, every hit of every event of the library without reading the others.

Layout of the library (one folder):
  per_channel/<subdetector>_<NNN>.root  pieces: tree "channels" (channel, start, n) and tree "hits" (event, energy, time), sorted by
                                        channel, then event, then time; a piece holds a contiguous range of whole channels
  support/                              samples.json, productions.json, manifest.root (global index -> production, block, event in the
                                        block, seed, hard-event bits), per_channel_index.json (the channel range of each piece),
                                        tile_cells.csv and lar_cells.csv (geometry and ring of each channel)
Channel: the identifier written by Lorenzetti itself; Tile: the 64-bit cellID of the hit, which is the PMT; LAr: the upper 32 bits of
the cellID (the ATLAS cell identifier; the lower 32 bits are the same in every hit). Event: global index in the library (int32), the
inclusive sample first and the hard sample after it. Energy: visible MeV as in the FRHITS (float32). Time: ns as in the FRHITS, with the
time of flight already removed (float32).
"""
import json
import os
import sys
import types

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CONDITIONS_SCRIPTS = os.path.join(REPO, "geometry", "scripts", "atlas_conditions_run4")
CONDITIONS_DATA = os.path.join(REPO, "geometry", "data", "atlas_conditions_run4")

SUBDETECTORS = ("tile", "emb", "emec", "hec", "fcal")
FRHITS_PREFIX = {"tile": "TileCalHit", "emb": "LArHitEMB", "emec": "LArHitEMEC", "hec": "LArHitHEC", "fcal": "LArHitFCAL"}
EVENTS_PER_BLOCK = 50
# hard-event marking of an inclusive block ("<event> <bits>" per line): bit 1 jet, bit 2 photon, bit 4 B hadron with a lepton. Events with a
# jet or a photon leave the inclusive sample, as the hard events of the hard sample stand for them; the B-with-lepton-only events stay.
VETO_BITS = 1 | 2


def channel_of(subdetector, cell_id):
    """channel identifier (uint64) from the cellID of the FRHITS (int64)"""
    c = np.asarray(cell_id, dtype=np.int64).view(np.uint64)
    return c if subdetector == "tile" else c >> np.uint64(32)


def read_marking(path):
    """{EventNumber: bits} of a hard-event marking file"""
    m = np.loadtxt(path, dtype=np.int64, ndmin=2)
    return {int(e): int(b) for e, b in m}


def read_block(frhits, columns):
    """per-event arrays of a FRHITS file: {column: list of numpy arrays}, and the EventNumber of each entry"""
    import ROOT
    d = ROOT.RDataFrame("CollectionTree", frhits).AsNumpy(["EventNumber"] + list(columns))
    out = {c: [np.asarray(v) for v in d[c]] for c in columns}
    return np.asarray(d["EventNumber"], dtype=np.int64), out


def kept_events(event_numbers, marking):
    """entries of a block that enter the library, in EventNumber order; marking None for the hard sample (every event enters)"""
    ev = np.asarray(event_numbers)
    assert len(set(ev.tolist())) == len(ev) and ev.min() >= 0 and ev.max() < EVENTS_PER_BLOCK, "EventNumber repeated or out of range"
    if marking is not None:
        assert sorted(marking) == sorted(ev.tolist()), "the marking does not cover the events of the block"
    order = np.argsort(ev, kind="stable")
    return [int(i) for i in order if marking is None or (marking[int(ev[i])] & VETO_BITS) == 0]


def lar_module():
    """the lar.py of the conditions scripts, for the cell identifier and the ring; it needs uproot only for the POOL readers, which are
    not used here, so a placeholder stands in when uproot is missing (the Lorenzetti image has ROOT but no uproot)"""
    if "uproot" not in sys.modules:
        try:
            import uproot  # noqa: F401
        except ImportError:
            sys.modules["uproot"] = types.ModuleType("uproot")
    if CONDITIONS_SCRIPTS not in sys.path:
        sys.path.insert(0, CONDITIONS_SCRIPTS)
    import lar
    return lar


def lar_ring(lar, channel, abs_eta):
    """(layer, ring) of a LAr cell, the keys of the LAr tables of geometry/data/atlas_conditions_run4/ (make_lar_tables.py): EMB, EMEC and
    HEC by the barrel_ec, region and eta index of the identifier; FCal by bands of |eta| 0.1 wide of the cell"""
    d, _s, b, smp, rg, e, _p = lar.campos(int(channel))
    if d == "FCAL":
        return lar.layer(d, b, smp), 1000 + int(np.floor(abs_eta / 0.1))
    return lar.layer(d, b, smp), b * 10000 + rg * 1000 + e


def tile_fields(pmt_id):
    """section, side, module, tower, sampling and pmt of a Tile PMT identifier (bits as written by Lorenzetti)"""
    x = np.asarray(pmt_id, dtype=np.uint64)

    def f(shift, mask):
        return ((x >> np.uint64(shift)) & np.uint64(mask)).astype(np.int64)
    return f(58, 7), f(54, 15), f(46, 255), f(40, 63), f(36, 15), f(34, 3)


def snapshot(path, trees):
    """write {tree name: {column: numpy array}} to one ROOT file"""
    import ROOT
    opts = ROOT.RDF.RSnapshotOptions()
    for i, (name, cols) in enumerate(trees.items()):
        opts.fMode = "RECREATE" if i == 0 else "UPDATE"
        ROOT.RDF.FromNumpy({k: np.ascontiguousarray(v) for k, v in cols.items()}).Snapshot(name, path, list(cols), opts)


def read_tree(path, tree, columns, begin=None, end=None):
    import ROOT
    d = ROOT.RDataFrame(tree, path)
    if begin is not None:
        d = d.Range(begin, end)
    return {k: np.asarray(v) for k, v in d.AsNumpy(list(columns)).items()}


def read_lines(path):
    return [ln.split(",") for ln in open(path).read().splitlines() if ln and not ln.startswith("#")]


def dump_json(obj, path):
    with open(path, "w", newline="\n") as o:
        json.dump(obj, o, indent=1)
        o.write("\n")
