"""The support/ folder of the free-running library. Two steps:

  events  before the blocks: reads the EventNumber of every block and its hard-event marking, gives each kept event its global index and
          writes support/samples.json (name, first index and number of events of each sample), support/productions.json,
          support/manifest.root (tree "manifest": global_index, sample, production, block, event_in_block, seed, hard_bits; hard_bits is -1
          in the hard sample, whose events are not marked) and the list of block jobs (blocks.csv: output name, FRHITS, marking, first
          global index, number of kept events)
  cells   after the blocks: the geometry of every channel with hits, from the channel trees of the block files, with its ring:
          support/tile_cells.csv (pmt_id, section, side, module, tower, sampling, pmt, eta, phi) and support/lar_cells.csv (cell_id,
          subdetector, layer, ring, eta, phi); checks that every ring is in the conditions tables of geometry/data/atlas_conditions_run4/

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
    jobs, idx = [], 0
    for sample, prod, block, seed, frhits, marking in rows:
        if not samples or samples[-1]["name"] != sample:
            assert sample not in [s["name"] for s in samples], f"sample {sample} not contiguous in the list"
            samples.append({"name": sample, "first_index": idx, "n_events": 0})
        if prod not in prods:
            prods.append(prod)
        evn, _ = c.read_block(frhits, [])
        mk = None if marking == "-" else c.read_marking(marking)
        keep = c.kept_events(evn, mk)
        for i in keep:
            e = int(evn[i])
            for k, v in (("global_index", idx), ("sample", len(samples) - 1), ("production", prods.index(prod)), ("block", int(block)),
                         ("event_in_block", e), ("seed", int(seed)), ("hard_bits", -1 if mk is None else mk[e])):
                man[k].append(v)
            idx += 1
        jobs.append((f"{prod}_{block}.channels.root", frhits, marking, idx - len(keep), len(keep)))
        samples[-1]["n_events"] += len(keep)
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
    geo = {}
    for det in c.SUBDETECTORS:
        ch, eta, phi = [], [], []
        for f in files:
            t = c.read_tree(f, f"{det}_channels", ["channel", "eta", "phi"])
            ch.append(t["channel"]), eta.append(t["eta"]), phi.append(t["phi"])
        ch, eta, phi = np.concatenate(ch), np.concatenate(eta), np.concatenate(phi)
        u, i, inv = np.unique(ch, return_index=True, return_inverse=True)
        if det == "tile":  # the D5 of EBA15 and EBC18 (see block_to_channels.py): the larger |eta|
            big = np.zeros(len(u))
            np.maximum.at(big, inv, np.abs(eta))
            eta_u = np.sign(eta[i]) * big
        else:
            eta_u = eta[i]
        dev = max(float(np.max(np.abs(eta - eta_u[inv]) * (det != "tile"), initial=0)), float(np.max(np.abs(phi - phi[i][inv]), initial=0)))
        assert dev < 1e-6, f"{det}: eta or phi of a channel changes between blocks ({dev})"
        geo[det] = (u, eta_u, phi[i])
    lar = c.lar_module()
    known_lar = {(r["layer"], int(r["ring"])) for r in csv.DictReader(open(os.path.join(c.CONDITIONS_DATA, "lar_run4_constants_by_ring.csv")))}
    known_tile = {(int(r["section"]), int(r["sampling"]), int(r["tower"]))
                  for r in csv.DictReader(open(os.path.join(c.CONDITIONS_DATA, "tile_run4_sampling_fraction_by_ring.csv")))}
    missing = set()
    with open(os.path.join(sup, "lar_cells.csv"), "w", newline="\n") as o:
        o.write("cell_id,subdetector,layer,ring,eta,phi\n")
        rows = []
        for det in ("emb", "emec", "hec", "fcal"):
            for x, e, p in zip(*geo[det]):
                ly, rg = c.lar_ring(lar, x, abs(e))
                if (ly, rg) not in known_lar:
                    missing.add((det, ly, rg))
                rows.append((int(x), det, ly, rg, e, p))
        for r in sorted(rows):
            o.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]:.6f},{r[5]:.6f}\n")
    u, eta, phi = geo["tile"]
    sec, side, mod, twr, smp, pmt = c.tile_fields(u)
    with open(os.path.join(sup, "tile_cells.csv"), "w", newline="\n") as o:
        o.write("pmt_id,section,side,module,tower,sampling,pmt,eta,phi\n")
        for k in range(len(u)):
            if (int(sec[k]), int(smp[k]), int(twr[k])) not in known_tile:
                missing.add(("tile", int(sec[k]), int(smp[k]), int(twr[k])))
            o.write(f"{int(u[k])},{sec[k]},{side[k]},{mod[k]},{twr[k]},{smp[k]},{pmt[k]},{eta[k]:.6f},{phi[k]:.6f}\n")
    print(f"{len(u)} Tile PMTs, {len(rows)} LAr cells; rings missing in the conditions tables: {sorted(missing) if missing else 'none'}")
    return 1 if missing else 0


if __name__ == "__main__":
    if sys.argv[1] == "events":
        events(*sys.argv[2:5])
    elif sys.argv[1] == "cells":
        sys.exit(cells(*sys.argv[2:4]))
    else:
        sys.exit(__doc__)
