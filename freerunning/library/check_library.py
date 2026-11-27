"""Check of the free-running library against the FRHITS it came from, read again here without the code of the other steps.

Per channel: number of hits and sum of the energy over the events that enter the library (inclusive blocks: without the events marked with
a jet or a photon), from the FRHITS, against the same in the pieces of per_channel/. The number of hits must be equal; the energy sums
must agree within the float32 rounding of each hit (relative 1e-6 of the sum of |energy|). Also checks the number of events of each
sample (support/samples.json), that the events of the pieces are inside the library and that every piece of the index exists.
Usage: python3 check_library.py <input list of make_support.py events> <library folder> [folder or URL of the pieces, default <library>/per_channel]
"""
import json
import os
import sys

import numpy as np
import ROOT

PREFIX = {"tile": "TileCalHit", "emb": "LArHitEMB", "emec": "LArHitEMEC", "hec": "LArHitHEC", "fcal": "LArHitFCAL"}
listing, lib = sys.argv[1], sys.argv[2]
pieces = sys.argv[3] if len(sys.argv) > 3 else os.path.join(lib, "per_channel")
rows = [ln.split(",") for ln in open(listing).read().splitlines() if ln and not ln.startswith("#")]

ref = {d: {} for d in PREFIX}  # channel -> [hits, energy, |energy|]
n_sample = {}
for sample, _prod, _block, _seed, frhits, marking in rows:
    cols = ["EventNumber"] + [f"{p}_{x}" for p in PREFIX.values() for x in ("cellID", "energy")]
    a = ROOT.RDataFrame("CollectionTree", frhits).AsNumpy(cols)
    veto = set()
    if marking != "-":
        for ln in open(marking).read().split("\n"):
            if ln.strip():
                e, bits = (int(x) for x in ln.split())
                if bits & 3:
                    veto.add(e)
    use = [i for i, e in enumerate(a["EventNumber"]) if int(e) not in veto]
    n_sample[sample] = n_sample.get(sample, 0) + len(use)
    for d, p in PREFIX.items():
        if not use:
            continue
        cid = np.concatenate([np.asarray(a[f"{p}_cellID"][i]) for i in use]).astype(np.int64).view(np.uint64)
        ch = cid if d == "tile" else cid >> np.uint64(32)
        en = np.concatenate([np.asarray(a[f"{p}_energy"][i]) for i in use]).astype(np.float64)
        u, inv = np.unique(ch, return_inverse=True)
        h, s, sa = np.bincount(inv), np.bincount(inv, en), np.bincount(inv, np.abs(en))
        r = ref[d]
        for k, x in enumerate(u.tolist()):
            v = r.setdefault(x, [0, 0.0, 0.0])
            v[0] += int(h[k]); v[1] += float(s[k]); v[2] += float(sa[k])

samples = json.load(open(os.path.join(lib, "support", "samples.json")))
total = sum(s["n_events"] for s in samples)
bad = [f"sample {s['name']}: {s['n_events']} events in the library, {n_sample.get(s['name'])} in the FRHITS"
       for s in samples if s["n_events"] != n_sample.get(s["name"])]
index = json.load(open(os.path.join(lib, "support", "per_channel_index.json")))
for d in PREFIX:
    got = {}
    for info in index[d]:
        f = f"{pieces}/{info['file']}"
        tf = ROOT.TFile.Open(f)
        if not tf or tf.IsZombie():
            bad.append(f"missing piece {info['file']}")
            continue
        chn = ROOT.RDataFrame("channels", f).AsNumpy(["channel", "start", "n"])
        hit = ROOT.RDataFrame("hits", f).AsNumpy(["event", "energy"])
        ev, en = hit["event"], hit["energy"].astype(np.float64)
        if len(ev) and (ev.min() < 0 or ev.max() >= total):
            bad.append(f"{info['file']}: event index outside 0..{total - 1}")
        st, n = chn["start"], chn["n"]
        if len(n) and (st[0] != 0 or np.any(st[1:] != st[:-1] + n[:-1]) or st[-1] + n[-1] != len(ev)):
            bad.append(f"{info['file']}: start and n do not tile the hits")
        sums = np.add.reduceat(en, st) if len(st) else np.zeros(0)
        for x, k, s in zip(chn["channel"].tolist(), n.tolist(), sums.tolist()):
            got[x] = (k, s)
    r = ref[d]
    if set(got) != set(r):
        bad.append(f"{d}: {len(set(r) - set(got))} channels only in the FRHITS, {len(set(got) - set(r))} only in the library")
    worst, nbad = 0.0, 0
    for x in set(got) & set(r):
        k, s = got[x]
        h, e, ea = r[x]
        dev = abs(s - e) / max(ea, 1e-30)
        worst = max(worst, dev)
        if k != h or dev > 1e-6:
            nbad += 1
    if nbad:
        bad.append(f"{d}: {nbad} channels with a different number of hits or energy")
    print(f"{d}: {len(got)} channels, {sum(v[0] for v in got.values())} hits; largest relative energy difference {worst:.2e}")
print(f"{total} events: " + ", ".join(f"{s['name']} {s['n_events']}" for s in samples))
print("CHECK OK" if not bad else "CHECK FAILED:\n  " + "\n  ".join(bad))
sys.exit(1 if bad else 0)
