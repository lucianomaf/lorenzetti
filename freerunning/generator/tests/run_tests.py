"""Tests of the free-running generator (no pytest in the Lorenzetti image: plain functions, run in sequence; exit 1 if one fails).
They need a small library built by freerunning/library/ (a few blocks of each sample are enough) and run with a low mu.
Usage: python3 run_tests.py <library folder> <work folder>
"""
import math
import os
import sys
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import conditions as cnd  # noqa: E402
import generate as g  # noqa: E402
from cells import Cells  # noqa: E402

LIB, WORK = sys.argv[1], sys.argv[2]
os.makedirs(WORK, exist_ok=True)
MU = "0.05"


def run(name, *opts):
    out = os.path.join(WORK, name + ".root")
    g.main(["--library", LIB, "--orbits", "1", "--mu", MU, "--output", out, *opts])
    return out


def read(fn, tree="samples", cols=("r_hg", "r_lg", "a")):
    import ROOT
    d = ROOT.RDataFrame(tree, fn).AsNumpy(list(cols))
    return {k: (np.array([np.asarray(x) for x in v]) if v.dtype == object else v) for k, v in d.items()}


def channels(fn):
    import ROOT
    return [int(x) for x in ROOT.RDataFrame("channels", fn).AsNumpy(["channel"])["channel"]]


# ---------------------------------------------------------------- tests
def test_same_seed_same_sequence():
    a, b, c = run("s1a", "--cell", "LBA14 A8", "--seed", "3"), run("s1b", "--cell", "LBA14 A8", "--seed", "3"), run("s1c", "--cell", "LBA14 A8", "--seed", "4")
    ra, rb, rc = read(a), read(b), read(c)
    assert all(np.array_equal(ra[k], rb[k]) for k in ra), "the same seed gave different samples"
    assert not np.array_equal(ra["r_hg"], rc["r_hg"]), "another seed gave the same samples"


def test_channel_alone_equals_channel_in_a_tower():
    t = run("tower", "--tower", "5", "--seed", "7")
    chs = channels(t)
    x = next(c for c in chs if c < (1 << 32))  # a LAr cell of the tower
    one = run("one", "--cell", str(x), "--seed", "7")
    rt, ro = read(t), read(one)
    k = chs.index(x)
    for col in ("r_hg", "r_lg", "a"):
        assert np.array_equal(rt[col][:, k], ro[col][:, 0]), f"{col} of the cell differs inside the tower"


def test_each_event_once_per_pass_and_a_new_order():
    class A:
        seed, mu, reuse = 5, 0.2, "on"
    lib = g.Library(LIB)
    cond = cnd.Conditions()
    K = 3 * cnd.NBC
    col = np.r_[np.zeros(10, bool), np.tile(cond.train, 3), np.zeros(10, bool)]
    bc, ev, cp, passes = g.interactions(A, lib, col, 10, K, 10)
    kept = (bc >= 10) & (bc < 10 + K)
    s = lib.sample("inclusive")
    m = kept & (ev >= s["first_index"]) & (ev < s["first_index"] + s["n_events"])
    assert passes["inclusive"] >= 2, "the test needs two passes"
    order = g.copy_order(A.seed)
    p0, p1 = m & (cp == order[0]), m & (cp == order[1])
    assert len(np.unique(ev[p0])) == p0.sum() == s["n_events"], "the first pass does not use each event once"
    assert len(np.unique(ev[p1])) == p1.sum(), "an event repeats inside the second pass"
    e0 = ev[p0][np.argsort(bc[p0], kind="stable")]
    e1 = ev[p1][np.argsort(bc[p1], kind="stable")]
    assert not np.array_equal(e0[:len(e1)], e1), "the second pass has the order of the first"
    A.reuse = "off"
    try:
        g.interactions(A, lib, col, 10, K, 10)
    except SystemExit:
        pass
    else:
        raise AssertionError("--reuse off did not stop when the events ran out")


def test_copies():
    cells = Cells(os.path.join(LIB, "support"))
    for x, d in list(cells.tile.items())[:200]:
        for c in range(g.N_COPIES):
            y = cells.source(x, c)
            k, refl = c >> 1, c & 1
            assert ((y >> 46) & 255) == (d["module"] - k) % 64, "Tile module of the source"
            if not refl:
                assert ((y >> 54) & 15) == d["side"]
    by = {}
    for x, d in cells.lar.items():
        by.setdefault(d["layer"], x)
    for layer, x in by.items():
        f = g.Cells.__dict__  # noqa: F841
        det, side, b, smp, rg, e, p = __import__("cells").lar_fields(x)
        n = cells.nphi[(det, b, smp, rg)]
        for c in range(g.N_COPIES):
            y = cells.source(x, c)
            det2, side2, b2, smp2, rg2, e2, p2 = __import__("cells").lar_fields(y)
            k, refl = c >> 1, c & 1
            assert (det2, b2, smp2, rg2, e2) == (det, b, smp, rg, e), f"{layer}: the copy left the ring"
            assert side2 == (-side if refl else side), f"{layer}: side of the copy"
            if det == "FCAL":
                assert p2 == p, "the FCal rotated"
            else:
                assert p2 == (p - (k * n) // 64) % n, f"{layer}: phi of the copy"
            if cells.get(y) is not None:
                assert abs(abs(cells.get(y)["eta"]) - abs(cells.get(x)["eta"])) < 1e-3, f"{layer}: |eta| of the copy"


def test_truth_is_the_sum_of_the_hits():
    fn = run("truth", "--cell", "LBA14 A8", "--seed", "2")
    a = read(fn)["a"]
    h = read(fn, "hits", ("channel", "sample", "visible_energy_mev"))
    cond = cnd.Conditions()
    cells = Cells(os.path.join(LIB, "support"))
    for k, x in enumerate(channels(fn)):
        cpm = cond.tile(cells.get(x)["ring"], "hg")["counts_per_visible_mev"]
        m = h["channel"] == k
        ref = np.bincount(h["sample"][m], weights=h["visible_energy_mev"][m] * cpm, minlength=len(a))
        assert np.allclose(ref, a[:, k], rtol=1e-5, atol=1e-4), "a[n] is not the sum of the hits"


def test_noiseless_samples_are_the_pulse_of_the_hits():
    sig = {}
    for name in ("lar", "tile"):
        f = getattr(cnd.Conditions, name)
        sig[name] = f

        def nonoise(self, *args, _f=f):
            p = _f(self, *args)
            p["sigma"] = 0.0
            return p
        setattr(cnd.Conditions, name, nonoise)
    keep = g.photostat_scale
    g.photostat_scale = lambda rng, pe: np.ones_like(pe)
    try:
        for cell in ("LBA14 A8", "EMB2 0.5"):
            fn = run("nonoise", "--cell", cell, "--seed", "9", "--mu", "0.3")
            r = read(fn)["r_hg"]
            h = read(fn, "hits", ("channel", "visible_energy_mev", "time_ns", "origin_bc"))
            cond = cnd.Conditions()
            cells = Cells(os.path.join(LIB, "support"))
            for k, x in enumerate(channels(fn)):
                d = cells.get(x)
                p = cond.lar(d["layer"], d["ring"], "hg") if d["det"] == "lar" else cond.tile(d["ring"], "hg")
                fn_p = (lambda t: np.interp(t, p["pulse_t"], p["pulse_w"], left=0, right=0)) if d["det"] == "lar" else cnd.tile_pulse
                m = h["channel"] == k
                y = np.zeros(len(r))
                for e, t, b in zip(h["visible_energy_mev"][m], h["time_ns"][m], h["origin_bc"][m]):
                    n = np.arange(len(r))
                    y += e * p["counts_per_visible_mev"] * fn_p(25.0 * (n - b) + p["phase"] - t)
                # the hits tree has the hits whose truth falls in the kept samples; hits of the warm-up are not in it, so compare from
                # the samples where the warm-up no longer reaches
                s = slice(40, len(r) - 40)
                dev = float(np.max(np.abs(r[s, k] - p["pedestal"] - y[s])))
                assert dev <= 0.5 + 1e-6, f"{cell}: noiseless samples differ from the pulse of the hits by {dev:.3f} (rounding: 0.5)"
                assert np.max(y[s]) > 1, f"{cell}: no signal in the test"
    finally:
        for name, f in sig.items():
            setattr(cnd.Conditions, name, f)
        g.photostat_scale = keep


def test_noise():
    fn = run("noise", "--cell", "EMB2 0.5", "--seed", "11", "--mu", "0", "--orbits", "3")
    r = read(fn)["r_hg"][:, 0].astype(float)
    cond = cnd.Conditions()
    cells = Cells(os.path.join(LIB, "support"))
    d = cells.get(channels(fn)[0])
    p = cond.lar(d["layer"], d["ring"], "hg")
    x = r - r.mean()
    rho = [float(np.mean(x[:-k] * x[k:]) / np.mean(x * x)) for k in (1, 2, 3, 4)]
    assert abs(r.mean() - p["pedestal"]) < 0.05, "LAr pedestal"
    assert abs(math.sqrt(np.var(r) - 1 / 12) / p["sigma"] - 1) < 0.05, "LAr sigma"
    assert np.allclose(rho, p["rho"], atol=0.05), f"LAr correlation {np.round(rho, 3)} against {p['rho']}"
    fn = run("noise_tile", "--cell", "LBA14 A8", "--seed", "11", "--mu", "0", "--orbits", "3")
    rr = read(fn)
    for gain in ("hg", "lg"):
        rt = rr[f"r_{gain}"][:, 0].astype(float)
        pt = cond.tile(Cells(os.path.join(LIB, "support")).get(channels(fn)[0])["ring"], gain)
        assert abs(rt.mean() - pt["pedestal"]) < 0.05, f"Tile pedestal {gain}"
        assert abs(math.sqrt(np.var(rt) - 1 / 12) / pt["sigma"] - 1) < 0.05, f"Tile sigma {gain}"


def test_tail_has_zero_area():
    t = cnd.TileTail(50_000.0)
    grid = np.arange(-40, t.t1, 0.5)
    area = np.trapezoid(t.short(grid), grid) if hasattr(np, "trapezoid") else np.trapz(t.short(grid), grid)
    beyond = -t.u1 * t.tau  # the integral of -u(T1) exp(-(t - T1)/tau) beyond T1
    peak = np.trapezoid(cnd.tile_pulse(grid), grid) if hasattr(np, "trapezoid") else np.trapz(cnd.tile_pulse(grid), grid)
    assert abs(area + beyond) < 1e-3 * peak, f"tail area {area + beyond} against the pulse area {peak}"


def test_saturation():
    f = cnd.Conditions.tile

    def big(self, *args):
        p = f(self, *args)
        p["counts_per_visible_mev"] *= 1e6
        return p
    cnd.Conditions.tile = big
    try:
        fn = run("sat", "--cell", "LBA14 A8", "--seed", "2", "--mu", "0.3")
        d = read(fn, "samples", ("r_hg", "saturated_hg"))
        assert d["r_hg"].max() == cnd.TILE_TOP and d["saturated_hg"].max() == 1, "no saturation at 4095"
        assert np.all(d["r_hg"][d["saturated_hg"] == 1] == cnd.TILE_TOP), "saturated samples not at the top"
    finally:
        cnd.Conditions.tile = f


if __name__ == "__main__":
    tests = [v for k, v in list(globals().items()) if k.startswith("test_")]
    bad = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
        except Exception:
            bad += 1
            print(f"FAIL {t.__name__}")
            traceback.print_exc()
    print(f"{len(tests) - bad} of {len(tests)} tests passed")
    sys.exit(1 if bad else 0)
