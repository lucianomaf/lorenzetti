"""Reading the free-running library (freerunning/library/): the samples, the index of the pieces and the hits of a set of channels."""
import json
import os

import numpy as np


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

    def hits(self, sub, channels, gap=200_000):
        """{channel: (event int32, energy float32, time float32)} of the channels of a subdetector; neighbouring entry ranges of a piece are
        read together (gaps under <gap> entries)"""
        import ROOT
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
            if not len(k):
                continue
            a, b = st[k], st[k] + n[k]
            order = np.argsort(a)
            groups, cur = [], [order[0]]
            for i in order[1:]:
                if a[i] - b[cur[-1]] < gap:
                    cur.append(i)
                else:
                    groups.append(cur)
                    cur = [i]
            groups.append(cur)
            for g in groups:
                lo, hi = int(a[g].min()), int(b[g].max())
                d = ROOT.RDataFrame("hits", f).Range(lo, hi).AsNumpy(["event", "energy", "time"])
                for i in g:
                    s, e = int(a[i]) - lo, int(b[i]) - lo
                    out[int(sel[i])] = (d["event"][s:e], d["energy"][s:e], d["time"][s:e])
        return out
