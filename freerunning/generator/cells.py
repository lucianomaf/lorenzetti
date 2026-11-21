"""The channels of the library: names, the selection of the command line (--cell, --cluster, --tower, --layers, with the meaning of the
ATLAS generator of the free-running sequences) and the symmetry copies (rotation by k x 2 pi / 64 and reflection between the sides).

Names: Tile, as ATLAS ("LBA14 A8": partition, module from 1, cell); LAr, "<layer> <A|C> eta=<eta> phi=<phi>" with the layers EMB0-EMB3,
EMEC0-EMEC3, HEC0-HEC3, FCAL1-FCAL3. A copy (k, reflected) of a target channel takes the hits of its source channel: the channel k steps of
2 pi / 64 back in phi (LAr: k x N / 64 cells of a ring of N cells in phi, rounded down; FCal: no rotation), on the other side when
reflected (Tile: the D0, across z = 0, is its own reflection).
"""
import csv
import math
import os

import numpy as np

PHI_STEPS = 64
CLUSTER_DETA, CLUSTER_DPHI = 0.075, 0.175
TOWER_DETA, TOWER_LAST = 0.1, 24
LAR_LAYER_NAME = {"EMB_PS": "EMB0", "EMB_1": "EMB1", "EMB_2": "EMB2", "EMB_3": "EMB3", "EMEC_PS": "EMEC0", "EMEC_1": "EMEC1",
                  "EMEC_2": "EMEC2", "EMEC_3": "EMEC3", "HEC_0": "HEC0", "HEC_1": "HEC1", "HEC_2": "HEC2", "HEC_3": "HEC3",
                  "FCal_1": "FCAL1", "FCal_2": "FCAL2", "FCal_3": "FCAL3"}
TILE_LAYER_NAME = {0: "TILEA", 1: "TILEBC", 2: "TILED", 3: "TILEE"}
GROUPS = {"em": ("EMB0", "EMB2", "EMB3", "EMEC0", "EMEC2", "EMEC3"),
          "had": ("HEC0", "HEC1", "HEC2", "HEC3", "TILEA", "TILEBC", "TILED"),
          "tile": ("TILEA", "TILEBC", "TILED"),
          "all": ("EMB0", "EMB1", "EMB2", "EMB3", "EMEC0", "EMEC1", "EMEC2", "EMEC3", "HEC0", "HEC1", "HEC2", "HEC3", "TILEA", "TILEBC",
                  "TILED")}


def tile_cell_name(sec, smp, twr):
    if smp == 0:
        return f"A{twr + 1}"
    if smp == 1:
        if sec == 1:
            return f"BC{twr + 1}" if twr < 8 else "B9"
        return f"B{twr + 1}" if sec == 2 else "C10"
    if smp == 2:
        return f"D{twr // 2}" if sec == 1 else ({10: "D5", 12: "D6"}[twr] if sec == 2 else "D4")
    return {10: "E1", 11: "E2", 13: "E3", 15: "E4"}[twr]


def lar_fields(x):
    """detector, side, |barrel_ec|, sampling, region, eta and phi indices of a LAr cell identifier (the layout of lar.py)"""
    part = (x >> 26) & 7
    if part == 3:
        bec = [-3, -2, -1, 1, 2, 3][(x >> 23) & 7]
        return "EM", (1 if bec > 0 else -1), abs(bec), (x >> 21) & 3, (x >> 18) & 7, (x >> 9) & 511, (x >> 1) & 255
    if part == 4:
        return "HEC", (1 if (x >> 25) & 1 else -1), 2, (x >> 23) & 3, (x >> 22) & 1, (x >> 18) & 15, (x >> 12) & 63
    return "FCAL", (1 if (x >> 25) & 1 else -1), 2, ((x >> 23) & 3) + 1, 0, (x >> 17) & 63, (x >> 13) & 15


class Cells:
    def __init__(self, support):
        self.lar, self.tile = {}, {}
        with open(os.path.join(support, "lar_cells.csv")) as f:
            for r in csv.DictReader(f):
                x = int(r["cell_id"])
                self.lar[x] = dict(id=x, det="lar", sub=r["subdetector"], layer=r["layer"], ring=int(r["ring"]), eta=float(r["eta"]),
                                   phi=float(r["phi"]), lname=LAR_LAYER_NAME[r["layer"]])
        with open(os.path.join(support, "tile_cells.csv")) as f:
            for r in csv.DictReader(f):
                x = int(r["pmt_id"])
                d = {k: int(r[k]) for k in ("section", "side", "module", "tower", "sampling", "pmt")}
                d.update(id=x, det="tile", sub="tile", eta=float(r["eta"]), phi=float(r["phi"]), lname=TILE_LAYER_NAME[d["sampling"]],
                         ring=(d["section"], d["sampling"], d["tower"]))
                self.tile[x] = d
        sign = {}
        for d in self.tile.values():
            if d["eta"] != 0:
                sign[d["side"]] = 1 if d["eta"] > 0 else -1
        self.tile_side_letter = {s: ("A" if v > 0 else "C") for s, v in sign.items()}
        self.tile_sides = sorted(sign)
        for d in self.tile.values():
            part = "LB" if d["section"] == 1 else "EB"
            d["cell"] = tile_cell_name(d["section"], d["sampling"], d["tower"])
            d["name"] = f"{part}{self.tile_side_letter.get(d['side'], 'A')}{d['module'] + 1:02d} {d['cell']}"
        self.nphi = {}
        for x, d in self.lar.items():
            det, side, b, smp, rg, e, p = lar_fields(x)
            key = (det, b, smp, rg)
            self.nphi[key] = max(self.nphi.get(key, 0), p + 1)
            d["name"] = f"{d['lname']} {'A' if side > 0 else 'C'} eta={d['eta']:.4f} phi={d['phi']:.4f}"

    def get(self, x):
        return self.lar.get(x) or self.tile.get(x)

    # ---------------------------------------------------------------- copies
    def source(self, x, copy):
        """the channel whose hits the copy (index 0..127: k = copy // 2 steps, reflected if copy is odd) of channel x takes"""
        k, refl = copy >> 1, copy & 1
        if x in self.tile:
            d = self.tile[x]
            mod = (d["module"] - k) % PHI_STEPS
            side = d["side"]
            if refl and not (d["section"] == 1 and d["sampling"] == 2 and d["tower"] == 0):
                side = {self.tile_sides[0]: self.tile_sides[-1], self.tile_sides[-1]: self.tile_sides[0]}[side]
            return (x & ~((255 << 46) | (15 << 54))) | (mod << 46) | (side << 54)
        det, side, b, smp, rg, e, p = lar_fields(x)
        if det == "EM":
            n = self.nphi[(det, b, smp, rg)]
            p2 = (p - (k * n) // PHI_STEPS) % n
            idx = (x >> 23) & 7
            if refl:
                idx = 5 - idx
            return (x & ~((7 << 23) | (255 << 1))) | (idx << 23) | (p2 << 1)
        if det == "HEC":
            n = self.nphi[(det, b, smp, rg)]
            p2 = (p - (k * n) // PHI_STEPS) % n
            y = (x & ~(63 << 12)) | (p2 << 12)
            return y ^ (1 << 25) if refl else y
        return x ^ (1 << 25) if refl else x  # FCal: reflection only

    # ---------------------------------------------------------------- selection
    def layers(self, text, default):
        names = set()
        for w in (text or default).split(","):
            w = w.strip().upper()
            if w.lower() in GROUPS:
                names.update(GROUPS[w.lower()])
            elif w in LAR_LAYER_NAME.values() or w in TILE_LAYER_NAME.values():
                names.add(w)
            else:
                raise SystemExit(f"ERROR: unknown layer {w}")
        return names

    def all_cells(self):
        return list(self.lar.values()) + list(self.tile.values())

    def cell(self, text, rng):
        """the channels of --cell: a list of identifiers"""
        t = text.strip()
        if t.isdigit():
            x = int(t)
            if x not in self.lar and x not in self.tile:
                raise SystemExit(f"ERROR: no channel {x} in the library")
            if x in self.tile:
                d = self.tile[x]
                return sorted(y for y, e in self.tile.items() if (e["side"], e["module"], e["ring"]) == (d["side"], d["module"], d["ring"]))
            return [x]
        w = t.split()
        if w[0].lower() == "tile":  # short form, a Tile cell name: drawn among the sides and modules
            keys = sorted({(d["side"], d["module"]) for d in self.tile.values() if d["cell"] == w[1]})
            if not keys:
                raise SystemExit(f"ERROR: no Tile cell {w[1]}")
            side, mod = keys[int(rng.integers(len(keys)))]
            return sorted(x for x, d in self.tile.items() if d["cell"] == w[1] and (d["side"], d["module"]) == (side, mod))
        if w[0][:3].upper() in ("LBA", "LBC", "EBA", "EBC"):
            got = sorted(x for x, d in self.tile.items() if d["name"] == f"{w[0].upper()} {w[1]}")
            if not got:
                raise SystemExit(f"ERROR: no Tile cell {t}")
            return got
        layer = w[0].upper()
        cand = [d for d in self.lar.values() if d["lname"] == layer]
        if not cand:
            raise SystemExit(f"ERROR: unknown LAr layer {layer}")
        if len(w) == 2:  # short form, layer and |eta|: drawn among the cells of the nearest ring, both sides and every phi
            target = abs(float(w[1]))
            best = min(cand, key=lambda d: abs(abs(d["eta"]) - target))
            ring = [d for d in cand if d["ring"] == best["ring"]]
            return [ring[int(rng.integers(len(ring)))]["id"]]
        side = 1 if w[1].upper() == "A" else -1
        kv = dict(s.split("=") for s in w[2:])
        eta, phi = float(kv["eta"]), float(kv["phi"])
        best = min((d for d in cand if np.sign(d["eta"]) == side),
                   key=lambda d: (d["eta"] - eta) ** 2 + (math.remainder(d["phi"] - phi, 2 * math.pi)) ** 2)
        return [best["id"]]

    def window(self, layers, eta_lo, eta_hi, side, phi_lo, dphi, tower=None, module=None):
        out = []
        for d in self.all_cells():
            if d["lname"] not in layers:
                continue
            if d["det"] == "tile" and tower is not None:
                if d["tower"] != tower or d["module"] != module or (d["side"] != side and d["eta"] != 0):
                    continue
                out.append(d["id"])
                continue
            if d["det"] == "tile" and self.tile_side_letter.get(d["side"]) != ("A" if side > 0 else "C") and d["eta"] != 0:
                continue
            if d["det"] == "lar" and np.sign(d["eta"]) != side:
                continue
            if not (eta_lo <= abs(d["eta"]) < eta_hi):
                continue
            if (d["phi"] - phi_lo) % (2 * math.pi) < dphi:
                out.append(d["id"])
        return sorted(out)

    def cluster(self, eta, deta, dphi, layers, rng):
        step = 2 * math.pi / PHI_STEPS
        dphi = max(1, int(round(dphi / step))) * step
        side = 1 if rng.integers(2) else -1
        phi_c = float(rng.uniform(0, 2 * math.pi))
        return self.window(layers, eta - deta / 2, eta + deta / 2, side, phi_c - dphi / 2, dphi), side, phi_c, dphi

    def tower(self, n, layers, rng):
        step = 2 * math.pi / PHI_STEPS
        side_l = 1 if rng.integers(2) else -1
        module = int(rng.integers(PHI_STEPS))
        tile_side = next(s for s, v in self.tile_side_letter.items() if (v == "A") == (side_l > 0))
        lar = self.window({x for x in layers if not x.startswith("TILE")}, TOWER_DETA * n, TOWER_DETA * (n + 1), side_l, module * step, step)
        tile = self.window({x for x in layers if x.startswith("TILE")}, 0, 99, tile_side, 0, 0, tower=n, module=module)
        return sorted(lar + tile), side_l, module
