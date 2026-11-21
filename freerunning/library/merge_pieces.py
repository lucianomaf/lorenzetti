"""Second step of the free-running library: the block files -> the pieces of per_channel/. Two modes:

  plan   sums the number of hits of every channel over the blocks and cuts each subdetector into pieces of contiguous whole channels with at
         most <hits per piece> hits (a channel with more hits than that is a piece alone); writes support/per_channel_index.json (per
         subdetector, for each piece: file, first and last channel, number of channels and of hits) and the list of piece jobs
  piece  builds one piece: for each block, in the order of the global index, reads only the entries of the channels of the piece (the
         block trees are sorted by channel), joins them with a stable sort by channel (inside a channel the hits stay in event and time
         order, since the blocks come in index order) and writes per_channel/<subdetector>_<NNN>.root with tree "channels" (channel, start,
         n) and tree "hits" (event, energy, time); checks the number of hits of every channel against the plan

Default hits per piece: 160 million, about 1 GB per file (about 6 bytes per hit after compression; 12 bytes per hit in memory).
Usage: python3 merge_pieces.py plan <block files list> <library folder> <pieces.csv> [hits per piece]
       python3 merge_pieces.py piece <block files list> <library folder> <subdetector> <piece number>
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as c  # noqa: E402


def counts(files, det):
    ch, n = [], []
    for f in files:
        t = c.read_tree(f, f"{det}_channels", ["channel", "n"])
        ch.append(t["channel"]), n.append(t["n"])
    ch, n = np.concatenate(ch), np.concatenate(n)
    u, inv = np.unique(ch, return_inverse=True)
    return u, np.bincount(inv, weights=n).astype(np.int64)


def plan(files, lib, jobs_out, per_piece=160_000_000):
    per_piece = int(per_piece)
    index, jobs = {}, []
    for det in c.SUBDETECTORS:
        u, n = counts(files, det)
        pieces, start, acc = [], 0, 0
        for k in range(len(u)):
            if acc and acc + n[k] > per_piece:
                pieces.append((start, k))
                start, acc = k, 0
            acc += n[k]
        if len(u):
            pieces.append((start, len(u)))
        index[det] = []
        for j, (a, b) in enumerate(pieces):
            index[det].append({"file": f"{det}_{j:03d}.root", "first_channel": int(u[a]), "last_channel": int(u[b - 1]),
                               "n_channels": int(b - a), "n_hits": int(n[a:b].sum())})
            jobs.append((det, j))
        print(f"{det}: {len(u)} channels, {int(n.sum())} hits, {len(pieces)} pieces")
    os.makedirs(os.path.join(lib, "support"), exist_ok=True)
    c.dump_json(index, os.path.join(lib, "support", "per_channel_index.json"))
    with open(jobs_out, "w", newline="\n") as o:
        for det, j in jobs:
            o.write(f"{det},{j}\n")


def piece(files, lib, det, j):
    j = int(j)
    info = json.load(open(os.path.join(lib, "support", "per_channel_index.json")))[det][j]
    lo, hi = np.uint64(info["first_channel"]), np.uint64(info["last_channel"])
    parts = {k: [] for k in ("channel", "event", "energy", "time")}
    for f in files:
        t = c.read_tree(f, f"{det}_channels", ["channel", "n"])
        cum = np.concatenate([[0], np.cumsum(t["n"])])
        a, b = np.searchsorted(t["channel"], lo, "left"), np.searchsorted(t["channel"], hi, "right")
        if b > a:
            h = c.read_tree(f, det, list(parts), int(cum[a]), int(cum[b]))
            assert h["channel"][0] == t["channel"][a] and h["channel"][-1] == t["channel"][b - 1], f"{f}: entries out of place"
            for k in parts:
                parts[k].append(h[k])
    ch = np.concatenate(parts["channel"])
    o = np.argsort(ch, kind="stable")
    ch = ch[o]
    u, start, n = np.unique(ch, return_index=True, return_counts=True)
    assert len(u) == info["n_channels"] and int(n.sum()) == info["n_hits"], f"{det} piece {j}: channels or hits differ from the plan"
    ev = np.concatenate(parts["event"])[o]
    assert np.all((np.diff(ev) >= 0) | (np.diff(ch) > 0)), "events out of order inside a channel"
    out = os.path.join(lib, "per_channel", info["file"])
    os.makedirs(os.path.dirname(out), exist_ok=True)
    c.snapshot(out, {"channels": {"channel": u, "start": start.astype(np.int64), "n": n.astype(np.int64)},
                     "hits": {"event": ev, "energy": np.concatenate(parts["energy"])[o], "time": np.concatenate(parts["time"])[o]}})
    print(f"{info['file']}: {len(u)} channels, {len(ch)} hits from {len(files)} blocks")


if __name__ == "__main__":
    files = [r[0] for r in c.read_lines(sys.argv[2])]
    if sys.argv[1] == "plan":
        plan(files, *sys.argv[3:])
    elif sys.argv[1] == "piece":
        piece(files, *sys.argv[3:])
    else:
        sys.exit(__doc__)
