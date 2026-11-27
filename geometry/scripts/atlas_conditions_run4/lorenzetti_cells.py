"""The LAr cells seen by Lorenzetti, from its own output: the free-running hits files (simu_trf.py --free-running-hits) of a minimum-bias
production. Each hit carries the ATLAS identifier of the cell (written by CaloFreeRunningHitWriter), the nominal eta and phi of the cell
and the dictionary fields. Writes, per cell: detector, sampling (FCal: module), region, side field, eta and phi, and checks that eta and phi
do not change between hits of the same cell. A few tens of 50-event files reach almost every cell; cells that no event reaches take the
|eta| of their ring in make_lar_tables.py.
Usage: python3 lorenzetti_cells.py WORK/lorenzetti_lar_cells.csv <FRHITS file> [<FRHITS file> ...]
"""
import sys

import numpy as np
import uproot

OUT, FILES = sys.argv[1], sys.argv[2:]
CAMPOS = {"EMB": ("sampling", "region", "barrel_ec"), "EMEC": ("sampling", "region", "barrel_ec"), "HEC": ("sampling", "region", "pos_neg"),
          "FCAL": ("module", None, None)}
cells = {}
conflicts = 0
for k, fn in enumerate(FILES):
    t = uproot.open(fn)["CollectionTree"]
    for d, (a, b, c) in CAMPOS.items():
        br = [f"LArHit{d}_cellID", f"LArHit{d}_eta", f"LArHit{d}_phi", f"LArHit{d}_{a}"] + [f"LArHit{d}_{x}" for x in (b, c) if x]
        z = {n: np.concatenate(v) for n, v in t.arrays(br, library="np").items()}
        cid = z[f"LArHit{d}_cellID"]
        u, i = np.unique(cid, return_index=True)
        eta, phi = z[f"LArHit{d}_eta"], z[f"LArHit{d}_phi"]
        # the same cell has the same eta and phi in all its hits
        _, inv = np.unique(cid, return_inverse=True)
        conflicts += int(np.sum(np.abs(eta - eta[i][inv]) > 1e-9) + np.sum(np.abs(phi - phi[i][inv]) > 1e-9))
        fa = z[f"LArHit{d}_{a}"][i]
        fb = z[f"LArHit{d}_{b}"][i] if b else np.full(len(u), -1)
        fc = z[f"LArHit{d}_{c}"][i] if c else np.sign(eta[i]).astype(int)
        for j, cc in enumerate(u.tolist()):
            v = (d, int(fa[j]), int(fb[j]), int(fc[j]), float(eta[i][j]), float(phi[i][j]))
            old = cells.setdefault(cc, v)
            if old != v:
                conflicts += 1
    if k % 100 == 0:
        print(f"{k + 1}/{len(FILES)} files, {len(cells)} cells, {conflicts} conflicts", flush=True)
with open(OUT, "w", newline="\n") as o:
    o.write("offline_id,detector,sampling_or_module,region,side_field,eta,phi\n")
    for cc in sorted(cells):
        d, fa, fb, fc, e, p = cells[cc]
        o.write(f"{cc & 0xFFFFFFFFFFFFFFFF},{d},{fa},{fb},{fc},{e:.6f},{p:.6f}\n")
print(f"{len(cells)} cells from {len(FILES)} files, {conflicts} conflicts -> {OUT}")
