"""The support/ folder of the free-running library. Two steps:

  events  before the blocks: reads the EventNumber of every block and its hard-event marking, gives each kept event its global index and
          writes support/samples.json (name, first index and number of events of each sample; in the hard sample also the fraction of
          hard events per interaction, from the marking), support/productions.json,
          support/manifest.root (tree "manifest": global_index, sample, production, block, event_in_block, seed, hard_bits; hard_bits is -1
          in the hard sample, whose events are not marked) and the list of block jobs (blocks.csv: output name, FRHITS, marking, first
          global index, number of kept events)
  cells   after the blocks: every channel of the detector, from support/geometry_channels.npz (geometry_channels.py, the geometry step), with
          its ring: support/tile_cells.csv (pmt_id, section, side, module, tower, sampling, pmt, eta, phi) and support/lar_cells.csv
          (cell_id, subdetector, layer, ring, eta, phi). Checks that every channel with hits, in the channel trees of the block files, is in
          that list with the same eta and phi (Tile: phi only, since the D5 of EBA15 and EBC18 also takes the D4 centre in its hits), and
          that every ring is in the conditions tables of geometry/data/atlas_conditions_run4/; reports the channels per ring against the
          numbers of those tables

Input list (events): one block per line, "sample,production,block,seed,FRHITS,marking" (marking: the hard-event marking file of an inclusive
block, or - for a block of the hard sample), the inclusive sample first; the order of the lines is the order of the global index.
Usage: python3 make_support.py events <input list> <library folder> <blocks.csv>
       python3 make_support.py cells <block files list> <library folder>
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as c  # noqa: E402


def events(listing, lib, blocks_out):
    rows = c.read_lines(listing)
    sup = os.path.join(lib, "support")
    os.makedirs(sup, exist_ok=True)
    samples, prods, man = [], [], {k: [] for k in ("global_index", "sample", "production", "block", "event_in_block", "seed", "hard_bits")}
    jobs, idx, marked, marked_hard = [], 0, 0, 0
    for sample, prod, block, seed, frhits, marking in rows:
        if not samples or samples[-1]["name"] != sample:
            assert sample not in [s["name"] for s in samples], f"sample {sample} not contiguous in the list"
            samples.append({"name": sample, "first_index": idx, "n_events": 0})
        if prod not in prods:
            prods.append(prod)
        evn, _ = c.read_block(frhits, [])
        mk = None if marking == "-" else c.read_marking(marking)
        if mk is not None:
            marked += len(mk)
            marked_hard += sum(1 for b in mk.values() if b)
        keep = c.kept_events(evn, mk)
        for i in keep:
            e = int(evn[i])
            for k, v in (("global_index", idx), ("sample", len(samples) - 1), ("production", prods.index(prod)), ("block", int(block)),
                         ("event_in_block", e), ("seed", int(seed)), ("hard_bits", -1 if mk is None else mk[e])):
                man[k].append(v)
            idx += 1
        jobs.append((f"{prod}_{block}.channels.root", frhits, marking, idx - len(keep), len(keep)))
        samples[-1]["n_events"] += len(keep)
    # the fraction of hard events per interaction, for the mixing in the generator: the events with any bit of the marking over the marked
    # inclusive events
    if marked:
        for smp in samples:
            if smp["name"] == "hard":
                smp.update(hard_fraction=marked_hard / marked, marked_events=marked, marked_hard_events=marked_hard)
    c.dump_json(samples, os.path.join(sup, "samples.json"))
    c.dump_json(prods, os.path.join(sup, "productions.json"))
    types = {"global_index": np.int32, "sample": np.int32, "production": np.int32, "block": np.int32, "event_in_block": np.int32,
             "seed": np.int64, "hard_bits": np.int32}
    c.snapshot(os.path.join(sup, "manifest.root"), {"manifest": {k: np.asarray(v, dtype=types[k]) for k, v in man.items()}})
    with open(blocks_out, "w", newline="\n") as o:
        for j in jobs:
            o.write(",".join(str(x) for x in j) + "\n")
    print(f"{idx} events in {len(rows)} blocks: " + ", ".join(f"{s['name']} {s['n_events']}" for s in samples))


def cells(block_files, lib):
    files = [r[0] for r in c.read_lines(block_files)]
    sup = os.path.join(lib, "support")
    g = np.load(os.path.join(sup, "geometry_channels.npz"))
    geo = {det: (g[f"{det}_channel"], g[f"{det}_eta"], g[f"{det}_phi"]) for det in c.SUBDETECTORS}
    bad = []
    for det in c.SUBDETECTORS:
        u, eta, phi = geo[det]
        for f in files:
            t = c.read_tree(f, f"{det}_channels", ["channel", "eta", "phi"])
            k = np.searchsorted(u, t["channel"])
            found = (k < len(u)) & (u[np.minimum(k, len(u) - 1)] == t["channel"])
            if not found.all():
                bad.append(f"{f}: {int(np.sum(~found))} {det} channels with hits are not in the geometry list")
                continue
            dphi = np.abs(np.remainder(t["phi"] - phi[k] + np.pi, 2 * np.pi) - np.pi)
            dev = max(float(np.max(dphi, initial=0)), 0.0 if det == "tile" else float(np.max(np.abs(t["eta"] - eta[k]), initial=0)))
            if dev > 1e-6:
                bad.append(f"{f}: {det} channels with an eta or phi other than the geometry ({dev:.2e})")
    if bad:
        print("\n".join(bad))
        return 1
    lar = c.lar_module()
    tables = {}
    for r in csv.DictReader(open(os.path.join(c.CONDITIONS_DATA, "lar_run4_constants_by_ring.csv"))):
        tables[(r["layer"], int(r["ring"]))] = int(r["n_cells"])
    for r in csv.DictReader(open(os.path.join(c.CONDITIONS_DATA, "tile_run4_sampling_fraction_by_ring.csv"))):
        tables[(int(r["section"]), int(r["sampling"]), int(r["tower"]))] = int(r["n_pmts"])
    count, missing = {}, set()
    with open(os.path.join(sup, "lar_cells.csv"), "w", newline="\n") as o:
        o.write("cell_id,subdetector,layer,ring,eta,phi\n")
        rows = []
        for det in ("emb", "emec", "hec", "fcal"):
            for x, e, p in zip(*geo[det]):
                ly, rg = c.lar_ring(lar, x, abs(e))
                count[(ly, rg)] = count.get((ly, rg), 0) + 1
                if (ly, rg) not in tables:
                    missing.add((det, ly, rg))
                rows.append((int(x), det, ly, rg, e, p))
        for r in sorted(rows):
            o.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]:.6f},{r[5]:.6f}\n")
    u, eta, phi = geo["tile"]
    sec, side, mod, twr, smp, pmt = c.tile_fields(u)
    with open(os.path.join(sup, "tile_cells.csv"), "w", newline="\n") as o:
        o.write("pmt_id,section,side,module,tower,sampling,pmt,eta,phi\n")
        for k in range(len(u)):
            key = (int(sec[k]), int(smp[k]), int(twr[k]))
            count[key] = count.get(key, 0) + 1
            if key not in tables:
                missing.add(("tile",) + key)
            o.write(f"{int(u[k])},{sec[k]},{side[k]},{mod[k]},{twr[k]},{smp[k]},{pmt[k]},{eta[k]:.6f},{phi[k]:.6f}\n")
    other = sorted(((k, n, tables[k]) for k, n in count.items() if k in tables and n != tables[k]), key=str)
    print(f"{len(u)} Tile PMTs, {len(rows)} LAr cells; rings missing in the conditions tables: {sorted(missing) if missing else 'none'}")
    print(f"rings whose number of channels differs from the conditions tables: {len(other)}"
          + "".join(f"\n  {k}: {n} here, {m} in the table" for k, n, m in other))
    return 1 if missing else 0


if __name__ == "__main__":
    if sys.argv[1] == "events":
        events(*sys.argv[2:5])
    elif sys.argv[1] == "cells":
        sys.exit(cells(*sys.argv[2:4]))
    else:
        sys.exit(__doc__)
