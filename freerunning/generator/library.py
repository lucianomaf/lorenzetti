"""Reading the free-running library (freerunning/library/): the samples, the index of the pieces and the hits of a set of channels."""
import json
import os

import numpy as np

import cxx


class Library:
    def __init__(self, folder, support=None):
        self.per_channel = os.path.join(folder, "per_channel")
        self.support = support or os.path.join(folder, "support")
        self.samples = json.load(open(os.path.join(self.support, "samples.json")))
        self.index = json.load(open(os.path.join(self.support, "per_channel_index.json")))
        self.total = sum(s["n_events"] for s in self.samples)
        self._tables = {}

    def sample(self, name):
        return next(s for s in self.samples if s["name"] == name)

    def _table(self, sub, j):
        if (sub, j) not in self._tables:
            import ROOT
            f = os.path.join(self.per_channel, self.index[sub][j]["file"])
            d = ROOT.RDataFrame("channels", f).AsNumpy(["channel", "start", "n"])
            self._tables[(sub, j)] = (f, d["channel"], d["start"], d["n"])
        return self._tables[(sub, j)]

    def hits(self, sub, channels):
        """{channel: (event int32, energy float32, time float32)} of the channels of a subdetector; only the entries of those channels are
        read (one C++ pass over the piece, entries in increasing order)"""
        ROOT = cxx.load()
        out = {int(c): (np.zeros(0, np.int32), np.zeros(0, np.float32), np.zeros(0, np.float32)) for c in channels}
        want = np.array(sorted(out), dtype=np.uint64)
        for j, info in enumerate(self.index[sub]):
            sel = want[(want >= np.uint64(info["first_channel"])) & (want <= np.uint64(info["last_channel"]))]
            if not len(sel):
                continue
            f, ch, st, n = self._table(sub, j)
            k = np.searchsorted(ch, sel)
            ok = (k < len(ch)) & (ch[np.minimum(k, len(ch) - 1)] == sel)
            k, sel = k[ok], sel[ok]
            k, sel = k[n[k] > 0], sel[n[k] > 0]
            if not len(k):
                continue
            o = np.argsort(st[k])
            k, sel = k[o], sel[o]
            lo, hi = st[k].astype(np.int64), (st[k] + n[k]).astype(np.int64)
            s = ROOT.lzt_read_ranges(f, ROOT.std.vector["Long64_t"](lo.tolist()), ROOT.std.vector["Long64_t"](hi.tolist()))
            ev, en, tm = np.array(s.ev, np.int32), np.array(s.e, np.float32), np.array(s.t, np.float32)
            off = np.r_[0, np.cumsum(hi - lo)]
            for i in range(len(k)):
                out[int(sel[i])] = (ev[off[i]:off[i + 1]], en[off[i]:off[i + 1]], tm[off[i]:off[i + 1]])
        return out

