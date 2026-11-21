#!/usr/bin/env python3
"""Free-running sequences of the ATLAS-like calorimeter of Lorenzetti for the HL-LHC (Run 4): the digitised samples, at 40 MHz, of
selected channels with the pile-up of minimum-bias interactions, from the free-running library (freerunning/library/) and the public
conditions of geometry/data/atlas_conditions_run4/. The options follow those of the ATLAS generator of free-running sequences.

Method (README.md of this folder):
- Bunch crossings: the HL-LHC filling scheme (hllhc_bunch_train_ip1.csv), from --first-slot; a Poisson number of interactions of mean --mu
  in each colliding bunch crossing; each interaction is a hard event with probability HARD_FRACTION, an inclusive one otherwise.
- Events: in the kept orbits each event of a sample is used once per pass, in an order shuffled from the seed (a new shuffle in each pass);
  pass p uses the symmetry copy copy_order[p] (the event as it is in the first pass, then rotated by k x 2 pi / 64 and reflected between the
  sides), the same copy in the whole detector; --reuse off stops when a sample runs out. The warm-up before the kept samples and the few
  bunch crossings after them draw events with replacement, each with a random copy, from a random stream of their own.
- Signal: each hit adds counts x pulse(25 (n - m) + phase - t) to sample n (m its bunch crossing, t its time, phase the sampling phase of the
  ring, the pulse with its peak at 0). Tile: photoelectron statistics as Athena (TileHitVecToCntTool, type 2), per channel and bunch crossing.
- Noise, pedestal, rounding (floor(x + 0.5)) and saturation per gain; truth a[n]: each hit in the sample of the peak of its pulse.
- The same seed gives the same samples for a channel whatever the selection: the interactions and events come from one stream for the whole
  detector, the noise and the photoelectron statistics of a channel from the seed and its identifier.
"""
import argparse
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conditions as cnd  # noqa: E402
from cells import CLUSTER_DETA, CLUSTER_DPHI, TOWER_LAST, Cells  # noqa: E402
from library import Library  # noqa: E402

VERSION = "0.1.0"
HARD_FRACTION = 0.00338  # hard events per interaction: the fraction the hard filter keeps in the inclusive production (F61)
N_COPIES = 128
PHOTOSTAT_WINDOW = 62.5  # ns, TileHitVecToCntTool PhotostatWindow
TS, NBC = cnd.TS, cnd.NBC


def parse(argv=None):
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument("--version", action="version", version=VERSION)
    a.add_argument("--cell", help='one cell: "Tile A13" or "EMB2 0.5" (short form: the side and phi drawn with the seed), "LBA14 A8", '
                   '"EMB2 A eta=0.5125 phi=0.0245", or the identifier')
    a.add_argument("--cluster", type=float, metavar="ETA", help="the cells of --layers in a window around |eta| ETA (side and phi drawn)")
    a.add_argument("--cluster-deta", type=float, default=CLUSTER_DETA, help=f"full width of the window in eta (default {CLUSTER_DETA})")
    a.add_argument("--cluster-dphi", type=float, default=CLUSTER_DPHI, help=f"full width of the window in phi (default {CLUSTER_DPHI})")
    a.add_argument("--tower", type=int, metavar="N", help=f"the cells of --layers in the tower N of 0.1 x 2 pi / 64 (0 to {TOWER_LAST})")
    a.add_argument("--layers", help="LAr layers (EMB0-3, EMEC0-3, HEC0-3, FCAL1-3), Tile layers (TileA, TileBC, TileD, TileE) and the groups em "
                   "(default of --cluster), had, tile and all (default of --tower), comma-separated")
    a.add_argument("--list-cells", metavar="LAYER", nargs="?", const="", help="list the cells of a layer or Tile partition and exit")
    a.add_argument("--mu", type=float, default=200.0, help="mean number of interactions per colliding bunch crossing (default 200)")
    a.add_argument("--orbits", type=int, default=9, help="kept orbits of 3564 bunch crossings (default 9)")
    a.add_argument("--seed", type=int, default=1, help="seed (default 1)")
    a.add_argument("--first-slot", type=int, default=0, help="slot of the filling scheme of the first kept sample (default 0)")
    a.add_argument("--truth-unit", choices=("adc", "visible_mev", "total_mev"), default="adc", help="unit of the truth a[n] (default adc: "
                   "counts of the high gain)")
    a.add_argument("--gains", default="hg,lg", help="the gains written, hg and lg (default both)")
    a.add_argument("--reuse", choices=("on", "off"), default="on", help="on (default): when a sample runs out, go on with rotated and "
                   "reflected copies; off: stop with an error")
    a.add_argument("--tile-tail", choices=("on", "off"), default="off", help="long tail of the Tile pulse (AC coupling); default off")
    a.add_argument("--tile-tail-tau", type=float, default=50.0, help="time constant of the Tile tail, in microseconds (default 50)")
    a.add_argument("--no-hits", action="store_true", help="do not write the hits tree")
    a.add_argument("--library", required=True, help="the library folder (per_channel/ and support/)")
    a.add_argument("--support", help="the support folder (default: <library>/support)")
    a.add_argument("--output", help="output .root file")
    A = a.parse_args(argv)
    if A.list_cells is None:
        chosen = [o for o, v in (("--cell", A.cell), ("--cluster", A.cluster), ("--tower", A.tower)) if v is not None]
        if len(chosen) != 1:
            a.error("one of --cell, --cluster or --tower is required")
        if not A.output:
            a.error("--output is required")
    A.gains = [g.strip() for g in A.gains.split(",")]
    if not set(A.gains) <= {"hg", "lg"} or not A.gains:
        a.error("--gains takes hg and lg")
    if A.tower is not None and not 0 <= A.tower <= TOWER_LAST:
        a.error(f"--tower takes 0 to {TOWER_LAST}")
    return A


def copy_order(seed):
    """the copy of each pass: 0 (the event as it is) first, then the 127 others in an order from the seed"""
    return np.r_[0, 1 + np.random.default_rng([seed, 11]).permutation(N_COPIES - 1)]


def interactions(A, lib, colliding, W, K, P):
    """(bunch-crossing index, global event index, copy) of every interaction of the sequence (warm-up W, kept K, after P samples)"""
    samples = [lib.sample("inclusive"), lib.sample("hard")]
    order = copy_order(A.seed)
    out_bc, out_ev, out_cp = [], [], []
    # kept: one stream for the whole detector
    rk = np.random.default_rng([A.seed, 0])
    kept_bc = np.arange(W, W + K)
    n = rk.poisson(A.mu, K) * colliding[W:W + K]
    bc = np.repeat(kept_bc, n)
    hard = rk.random(len(bc)) < HARD_FRACTION
    passes = {}
    for s_i, s in enumerate(samples):
        b = bc[hard == bool(s_i)]
        N = s["n_events"]
        q = np.arange(len(b))
        p = q // N
        if len(b) > N and A.reuse == "off":
            raise SystemExit(f"ERROR: the sequence needs {len(b)} {s['name']} events and the library has {N}; use --reuse on, fewer orbits "
                             "or a lower mu")
        if len(b) and p[-1] >= N_COPIES:
            raise SystemExit(f"ERROR: the sequence needs more than {N_COPIES} passes over the {s['name']} events")
        ev = np.empty(len(b), np.int64)
        for pp in np.unique(p):
            m = p == pp
            perm = np.random.default_rng([A.seed, 10, s_i, int(pp)]).permutation(N)
            ev[m] = s["first_index"] + perm[q[m] - pp * N]
        out_bc.append(b), out_ev.append(ev), out_cp.append(order[p])
        passes[s["name"]] = int(p[-1]) + 1 if len(b) else 0
    # warm-up and after: with replacement, random copy
    rw = np.random.default_rng([A.seed, 1])
    idx = np.r_[np.arange(0, W), np.arange(W + K, W + K + P)]
    n = rw.poisson(A.mu, len(idx)) * colliding[idx]
    b = np.repeat(idx, n)
    hard = rw.random(len(b)) < HARD_FRACTION
    ev = np.where(hard, samples[1]["first_index"] + rw.integers(samples[1]["n_events"], size=len(b)),
                  samples[0]["first_index"] + rw.integers(samples[0]["n_events"], size=len(b)))
    out_bc.append(b), out_ev.append(ev), out_cp.append(rw.integers(N_COPIES, size=len(b)))
    bc, ev, cp = np.concatenate(out_bc), np.concatenate(out_ev), np.concatenate(out_cp).astype(np.int64)
    o = np.lexsort((ev, cp))
    return bc[o], ev[o], cp[o], passes


def photostat_scale(rng, pe):
    """TileHitVecToCntTool::applyPhotoStatistics, type 2: Gaussian above 20 photoelectrons, Poisson with the single-photoelectron
    spectrum (Gaussian of mean 1 and width 1, cut at 0, divided by 1.08332) below"""
    out = np.ones_like(pe)
    big = pe > 20
    out[big] = np.maximum(0.0, rng.normal(pe[big], np.sqrt(pe[big]))) / pe[big]
    for i in np.flatnonzero((pe > 0) & ~big):
        k = rng.poisson(pe[i])
        out[i] = (np.maximum(0.0, rng.normal(1.0, 1.0, k)).sum() / 1.08332 / k) if k > 0 else 0.0
    return out


def channel_rng(seed, stream, x, gain=0):
    return np.random.default_rng([seed, stream, gain, int(x) >> 32, int(x) & 0xFFFFFFFF])


def main(argv=None):
    A = parse(argv)
    t0 = time.time()
    support = A.support or os.path.join(A.library, "support")
    cells = Cells(support)
    if A.list_cells is not None:
        want = A.list_cells.upper()
        for d in sorted(cells.all_cells(), key=lambda d: d["name"]):
            if not want or d["lname"] == want or d["name"].startswith(want):
                print(f"{d['name']}  {d['id']}")
        return 0
    lib = Library(A.library, support)
    cond = cnd.Conditions()
    rsel = np.random.default_rng([A.seed, 4])
    if A.cell is not None:
        chan, label = cells.cell(A.cell, rsel), f"cell {A.cell}"
    elif A.cluster is not None:
        chan, side, phi_c, dphi = cells.cluster(A.cluster, A.cluster_deta, A.cluster_dphi, cells.layers(A.layers, "em"), rsel)
        label = f"cluster |eta| {A.cluster}, side {'A' if side > 0 else 'C'}, phi {phi_c:.4f}, {A.cluster_deta} x {dphi:.4f}"
    else:
        chan, side, module = cells.tower(A.tower, cells.layers(A.layers, "all"), rsel)
        label = f"tower {A.tower}, side {'A' if side > 0 else 'C'}, module {module + 1}"
    if not chan:
        raise SystemExit("ERROR: the selection has no channel")
    info = [cells.get(x) for x in chan]
    tail = cnd.TileTail(A.tile_tail_tau * 1000.0) if A.tile_tail == "on" else None
    setup = []
    for d in info:
        g = {}
        for gain in ("hg", "lg"):
            g[gain] = cond.lar(d["layer"], d["ring"], gain) if d["det"] == "lar" else cond.tile(d["ring"], gain)
        setup.append(g)
    # lengths: warm-up W covers the longest pulse after its peak; P after the kept samples covers the longest part before the peak
    tile_end = cnd.tile_pulse_end()
    ends = [float(s["hg"]["pulse_t"][-1]) if d["det"] == "lar" else tile_end for d, s in zip(info, setup)]
    starts = [float(min(s["hg"]["pulse_t"][0], s["lg"]["pulse_t"][0])) if d["det"] == "lar" else -cnd.TILE_N * cnd.TILE_TAU
              for d, s in zip(info, setup)]
    W = int(math.ceil((max(ends) + 25.0) / TS)) + 1
    if tail is not None and any(d["det"] == "tile" for d in info):
        W = max(W, int(math.ceil(10 * tail.tau / TS)))
    P = int(math.ceil((-min(starts) + 25.0) / TS)) + 1
    K = A.orbits * NBC
    Ntot = W + K + P
    slots = (A.first_slot - W + np.arange(Ntot)) % NBC
    colliding = cond.train[slots]
    ibc, iev, icp, passes = interactions(A, lib, colliding, W, K, P)
    copies = np.unique(icp)
    cut = np.searchsorted(icp, np.r_[copies, N_COPIES])
    # the hits of every source channel
    need = {}
    for d in info:
        for c in copies:
            src = cells.source(d["id"], int(c))
            need.setdefault(d["sub"], set()).add(src)
    hits = {}
    for sub, chs in need.items():
        hits.update(lib.hits(sub, chs))
    t_read = time.time() - t0
    nch = len(chan)
    r = {g: np.zeros((K, nch), np.int16) for g in A.gains}
    sat = {g: np.zeros((K, nch), np.uint8) for g in A.gains}
    a_n = np.zeros((K, nch), np.float32)
    hit_rows = []
    for c, (d, s) in enumerate(zip(info, setup)):
        hb, hE, ht, hev, hcp = [], [], [], [], []
        for j, cp in enumerate(copies):
            src = cells.source(d["id"], int(cp))
            ev, en, tm = hits.get(src, (np.zeros(0, np.int32),) * 3)
            if not len(ev):
                continue
            e_ev, e_bc = iev[cut[j]:cut[j + 1]], ibc[cut[j]:cut[j + 1]]
            lo, hi = np.searchsorted(e_ev, ev, "left"), np.searchsorted(e_ev, ev, "right")
            cnt = hi - lo
            if not cnt.any():
                continue
            ih = np.repeat(np.arange(len(ev)), cnt)
            ie = lo[ih] + (np.arange(int(cnt.sum())) - np.repeat(np.cumsum(cnt) - cnt, cnt))
            hb.append(e_bc[ie]), hE.append(en[ih].astype(float)), ht.append(tm[ih].astype(float)), hev.append(ev[ih])
            hcp.append(np.full(len(ih), int(cp)))
        if hb:
            hb, hE, ht, hev, hcp = (np.concatenate(x) for x in (hb, hE, ht, hev, hcp))
        else:
            hb, hev, hcp = (np.zeros(0, np.int64),) * 3
            hE = ht = np.zeros(0)
        hA = hE.copy()
        if d["det"] == "tile" and len(hb):  # photoelectron statistics per bunch crossing, on the amplitude only
            ub, inv = np.unique(hb, return_inverse=True)
            ewin = np.bincount(inv, weights=np.where(np.abs(ht) < PHOTOSTAT_WINDOW, hE, 0.0), minlength=len(ub))
            scale = photostat_scale(channel_rng(A.seed, 3, d["id"]), ewin * s["hg"]["npe_per_visible_mev"])
            hA = hE * scale[inv]
        phase = s["hg"]["phase"]
        # truth
        n_h = hb + np.rint((ht - phase) / TS).astype(np.int64) - W
        ok = (n_h >= 0) & (n_h < K)
        unit = {"adc": hE * s["hg"]["counts_per_visible_mev"], "visible_mev": hE, "total_mev": hE * s["hg"]["visible_to_total"]}[A.truth_unit]
        a_n[:, c] = np.bincount(n_h[ok], weights=unit[ok], minlength=K)
        if not A.no_hits:
            hit_rows.append((np.full(int(ok.sum()), c), n_h[ok], hE[ok], ht[ok], hb[ok] - W, hev[ok], hcp[ok]))
        tau0 = phase - ht
        for gi, gain in enumerate(A.gains):
            p = s[gain]
            amp = hA * p["counts_per_visible_mev"]
            y = np.zeros(Ntot)
            if d["det"] == "lar":
                t_p, w_p = p["pulse_t"], p["pulse_w"]
                t_lo, t_hi, fn = float(t_p[0]), float(t_p[-1]), (lambda t: np.interp(t, t_p, w_p, left=0.0, right=0.0))
            else:
                t_lo = -cnd.TILE_N * cnd.TILE_TAU
                t_hi = tail.t1 if tail is not None else tile_end
                fn = tail.short if tail is not None else cnd.tile_pulse
            if len(hb):
                for jj in range(int(math.floor((t_lo - tau0.max()) / TS)), int(math.ceil((t_hi - tau0.min()) / TS)) + 1):
                    tt = TS * jj + tau0
                    m = (tt >= t_lo) & (tt <= t_hi)
                    idx = hb[m] + jj
                    k = (idx >= 0) & (idx < Ntot)
                    y += np.bincount(idx[k], weights=(amp[m] * fn(tt[m]))[k], minlength=Ntot)
                if tail is not None and d["det"] == "tile":  # beyond T1: -u(T1) exp(-(t - T1)/tau), summed recursively
                    j1 = np.floor((tail.t1 - tau0) / TS).astype(np.int64) + 1
                    idx = hb + j1
                    k = idx < Ntot
                    inj = np.bincount(idx[k], weights=(-amp * tail.u1 * np.exp(-(TS * j1 + tau0 - tail.t1) / tail.tau))[k],
                                      minlength=Ntot)
                    f = math.exp(-TS / tail.tau)
                    acc = 0.0
                    for i in range(Ntot):
                        acc = acc * f + inj[i]
                        y[i] += acc
            rn = channel_rng(A.seed, 2, d["id"], gi + 1 if gain == "lg" else 0)
            noise = cnd.lar_noise(rn, K, p["sigma"], p["rho"]) if d["det"] == "lar" else rn.normal(0.0, p["sigma"], K)
            q = np.floor(p["pedestal"] + y[W:W + K] + noise + 0.5)
            sat[gain][:, c] = (q < 0) | (q > p["top"])
            r[gain][:, c] = np.clip(q, 0, p["top"])
    t_gen = time.time() - t0 - t_read
    write(A, label, info, setup, chan, r, sat, a_n, hit_rows, slots[W:W + K], colliding[W:W + K], ibc, W, K, passes, lib, tail)
    print(f"{label}: {nch} channels, {A.orbits} orbits, {len(ibc)} interactions (warm-up included); passes {passes}; "
          f"read {t_read:.1f} s, generated {t_gen:.1f} s -> {A.output}")
    return 0


def git_commit():
    try:
        return subprocess.run(["git", "-C", os.path.dirname(os.path.abspath(__file__)), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, timeout=10).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def write(A, label, info, setup, chan, r, sat, a_n, hit_rows, slots, colliding, ibc, W, K, passes, lib, tail):
    import ROOT
    from array import array
    f = ROOT.TFile(A.output, "RECREATE")
    par = dict(version=VERSION, commit=git_commit(), command=" ".join(sys.argv), selection=label, seed=A.seed, mu=A.mu, orbits=A.orbits,
               first_slot=A.first_slot, filling_scheme="25ns_2760b_2748_2492_2574_288bpi_13inj_800ns_bs200ns (hllhc_bunch_train_ip1.csv)",
               library=os.path.abspath(A.library), samples=lib.samples, hard_fraction=HARD_FRACTION, reuse=A.reuse, passes=passes,
               truth_unit=A.truth_unit, gains=A.gains, tile_tail=A.tile_tail, tile_tail_tau_us=A.tile_tail_tau if tail else None,
               conditions="geometry/data/atlas_conditions_run4", warmup_samples=W,
               tile_pulse=f"CR-RC^{cnd.TILE_N}, tau {cnd.TILE_TAU} ns", tile_noise="one Gaussian, sigma keeping the variance of the double Gaussian",
               lar_noise="AR(4) from the correlation of the ring at lags 1-4", photostatistics="TileHitVecToCntTool type 2, per channel and bunch crossing")
    js = json.dumps(par, indent=1)
    ROOT.TNamed("parameters_json", js).Write()
    tp = ROOT.TTree("parameters", "parameters of the run")
    k_, v_ = ROOT.std.string(), ROOT.std.string()
    tp.Branch("key", k_), tp.Branch("value", v_)
    for k, v in par.items():
        k_.replace(0, k_.size(), str(k)), v_.replace(0, v_.size(), json.dumps(v))
        tp.Fill()
    tp.Write()
    # channels
    tc = ROOT.TTree("channels", "the channels of the sequence")
    cols = dict(index=array("i", [0]), channel=array("L", [0]), eta=array("d", [0]), phi=array("d", [0]), phase_ns=array("d", [0]),
                pedestal=array("d", [0]), counts_per_visible_mev_hg=array("d", [0]), counts_per_visible_mev_lg=array("d", [0]),
                sigma_hg=array("d", [0]), sigma_lg=array("d", [0]), visible_to_total=array("d", [0]))
    tps = {"i": "I", "L": "l", "d": "D"}
    for k, v in cols.items():
        tc.Branch(k, v, f"{k}/{tps[v.typecode]}")
    name, det, layer, ring = ROOT.std.string(), ROOT.std.string(), ROOT.std.string(), ROOT.std.string()
    for k, v in (("name", name), ("subdetector", det), ("layer", layer), ("ring", ring)):
        tc.Branch(k, v)
    for c, (d, s) in enumerate(zip(info, setup)):
        cols["index"][0], cols["channel"][0], cols["eta"][0], cols["phi"][0] = c, int(d["id"]), d["eta"], d["phi"]
        cols["phase_ns"][0], cols["pedestal"][0] = s["hg"]["phase"], s["hg"]["pedestal"]
        cols["counts_per_visible_mev_hg"][0], cols["counts_per_visible_mev_lg"][0] = s["hg"]["counts_per_visible_mev"], s["lg"]["counts_per_visible_mev"]
        cols["sigma_hg"][0], cols["sigma_lg"][0], cols["visible_to_total"][0] = s["hg"]["sigma"], s["lg"]["sigma"], s["hg"]["visible_to_total"]
        for v, x in ((name, d["name"]), (det, d["sub"]), (layer, d["lname"]), (ring, str(d["ring"]))):
            v.replace(0, v.size(), x)
        tc.Fill()
    tc.Write()
    # samples
    nch = len(chan)
    ts = ROOT.TTree("samples", "one entry per 25 ns sample")
    n_, bcid_, orbit_, filled_, ni_ = array("i", [0]), array("i", [0]), array("i", [0]), array("B", [0]), array("i", [0])
    ts.Branch("n", n_, "n/I"), ts.Branch("bcid", bcid_, "bcid/I"), ts.Branch("orbit", orbit_, "orbit/I"), ts.Branch("filled", filled_, "filled/b")
    ts.Branch("n_interactions", ni_, "n_interactions/I")
    bufs = {"nch": array("i", [nch])}  # the buffers must live as long as the tree
    ts.Branch("nch", bufs["nch"], "nch/I")
    for g in A.gains:
        bufs[f"r_{g}"] = np.zeros(nch, np.int16)
        ts.Branch(f"r_{g}", bufs[f"r_{g}"], f"r_{g}[nch]/S")
        bufs[f"saturated_{g}"] = np.zeros(nch, np.uint8)
        ts.Branch(f"saturated_{g}", bufs[f"saturated_{g}"], f"saturated_{g}[nch]/b")
    bufs["a"] = np.zeros(nch, np.float32)
    ts.Branch("a", bufs["a"], "a[nch]/F")
    nint = np.bincount(np.clip(ibc - W, -1, K)[((ibc - W) >= 0) & ((ibc - W) < K)], minlength=K)
    for n in range(K):
        n_[0], bcid_[0], orbit_[0], filled_[0], ni_[0] = n, int(slots[n]), (A.first_slot + n) // NBC, int(colliding[n]), int(nint[n])
        for g in A.gains:
            bufs[f"r_{g}"][:] = r[g][n]
            bufs[f"saturated_{g}"][:] = sat[g][n]
        bufs["a"][:] = a_n[n]
        ts.Fill()
    ts.Write()
    if not A.no_hits:
        th = ROOT.TTree("hits", "the hits of the kept samples")
        hv = dict(channel=array("i", [0]), sample=array("i", [0]), visible_energy_mev=array("f", [0]), time_ns=array("f", [0]),
                  origin_bc=array("i", [0]), event=array("i", [0]), copy=array("i", [0]))
        for k, v in hv.items():
            th.Branch(k, v, f"{k}/{'F' if v.typecode == 'f' else 'I'}")
        for rows in hit_rows:
            for vals in zip(*rows):
                for (k, v), x in zip(hv.items(), vals):
                    v[0] = x
                th.Fill()
        th.Write()
    f.Close()


if __name__ == "__main__":
    sys.exit(main())
