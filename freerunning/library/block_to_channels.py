"""First step of the free-running library: one block (one FRHITS file) -> its hits reorganised by channel, one tree per subdetector.

Keeps the events of the block that enter the library (common.kept_events: without the jet or photon events of an inclusive block), gives
them their global index (from make_support.py events), and writes, per subdetector, tree "<subdetector>" (channel, event, energy, time;
sorted by channel, event and time) and tree "<subdetector>_channels" (channel, n, eta, phi: number of hits of the channel in the block and
the nominal eta and phi of the cell, which must be the same in every hit). Energy and time are copied from the FRHITS as float32.
Usage: python3 block_to_channels.py <FRHITS> <marking or -> <first global index> <number of kept events> <output .root>
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as c  # noqa: E402

frhits, marking, first, expected, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
cols = [f"{c.FRHITS_PREFIX[d]}_{x}" for d in c.SUBDETECTORS for x in ("cellID", "energy", "time", "eta", "phi")]
evn, a = c.read_block(frhits, cols)
keep = c.kept_events(evn, None if marking == "-" else c.read_marking(marking))
assert len(keep) == expected, f"{len(keep)} kept events, expected {expected}"
trees, nhits = {}, 0
for det in c.SUBDETECTORS:
    p = c.FRHITS_PREFIX[det]

    def cat(x, dtype):
        v = [a[f"{p}_{x}"][i] for i in keep]
        return np.concatenate(v).astype(dtype) if v else np.zeros(0, dtype)
    ch = c.channel_of(det, cat("cellID", np.int64))
    ev = np.concatenate([np.full(len(a[f"{p}_cellID"][i]), first + n, np.int32) for n, i in enumerate(keep)]) if keep else np.zeros(0, np.int32)
    en, tm, eta, phi = cat("energy", np.float32), cat("time", np.float32), cat("eta", np.float64), cat("phi", np.float64)
    o = np.lexsort((tm, ev, ch))
    ch, ev, en, tm, eta, phi = ch[o], ev[o], en[o], tm[o], eta[o], phi[o]
    u, i0, n = np.unique(ch, return_index=True, return_counts=True)
    rep = np.repeat(np.arange(len(u)), n)
    eta_ch = eta[i0]
    if det == "tile" and len(u):
        # the D5 PMTs of the two special extended-barrel modules (EBA15, EBC18) also read the D4 region, as in ATLAS: their hits carry
        # two cell centres; the channel takes the D5 one, the larger |eta|
        eta_ch = np.sign(eta[i0]) * np.maximum.reduceat(np.abs(eta), i0)
    dev = max(float(np.max(np.abs(np.abs(eta) - np.abs(eta_ch[rep])) * (det != "tile"), initial=0)),
              float(np.max(np.abs(phi - phi[i0][rep]), initial=0)))
    assert dev < 1e-6, f"{det}: eta or phi changes between hits of a channel ({dev})"
    trees[det] = {"channel": ch, "event": ev, "energy": en, "time": tm}
    trees[f"{det}_channels"] = {"channel": u, "n": n.astype(np.int64), "eta": eta_ch, "phi": phi[i0]}
    nhits += len(ch)
c.snapshot(out, trees)
print(f"{os.path.basename(frhits)}: {len(keep)} events (global {first}..{first + len(keep) - 1}), {nhits} hits -> {out}")
