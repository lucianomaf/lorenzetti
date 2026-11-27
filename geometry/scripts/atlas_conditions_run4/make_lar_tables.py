"""The LAr tables of geometry/data/atlas_conditions_run4/, per layer and ring of |eta| (both sides together), from public sources only.

Inputs: WORK/sources.json (sources.py), WORK/lorenzetti_lar_cells.csv (lorenzetti_cells.py), the public CVMFS repository, the public Frontier
and CREST. Output: WORK/tables/lar_run4_constants_by_ring.csv, lar_run4_noise_autocorr_by_ring.csv, lar_run4_pulses_by_ring.csv.

Rings: EMB, EMEC and HEC, one per (detector, wheel, sampling, region, eta index) of the cell identifier, with the nominal |eta| of the
Lorenzetti cell (unique in each ring); FCal, bands of |eta| 0.1 wide per module. Constants: median over the cells of the ring, with the number
of cells more than 5% away from it. Pulse: mean of the waves of the ring, each normalised to its peak with the peak at t = 0, on a 0.5 ns
grid (STEP), leaving out waves still above 5% of their peak at their last point and cells whose peak time is more than 10 ns from the median of the
ring; a constant is subtracted over the window so that the area is zero, and the peak is brought back to 1. Noise autocorrelation: the vector
of the representative channel of the ring (FCal: the mean of the representatives of the band).
Usage: python3 make_lar_tables.py
"""
import base64
import json
import os
import sys
import urllib.request

import numpy as np

STEP = 0.5   # ns, time step of the pulse tables
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lar  # noqa: E402
from sources import WORK, load  # noqa: E402

src = load()
OUT = os.path.join(WORK, "tables")
os.makedirs(OUT, exist_ok=True)
off_mc = lar.cabling("COOLOFL_LAR/OFLP200", src["lar_cabling_mc"]["tag"], src["mc_run"])
off_data = lar.cabling("COOLONL_LAR/CONDBR2", src["lar_cabling_data"]["tag"], src["data_run"])
rank = lar.febs(lar.file_of(src["lar_febs"]["refs"]))
inv = {i: f for f, i in rank.items()}
cells = [int(x) for x in off_mc if x != 0xFFFFFFFF]
hash_of = {int(x): i for i, x in enumerate(off_mc) if x != 0xFFFFFFFF}

# rings
eta_lzt = {}
for ln in open(os.path.join(WORK, "lorenzetti_lar_cells.csv")).read().splitlines()[1:]:
    p = ln.split(",")
    eta_lzt[int(p[0]) >> 32] = abs(float(p[5]))
eta_ring = {}
for x in cells:
    d, s, b, smp, rg, e, p = lar.campos(x)
    if d != "FCAL" and x in eta_lzt:
        eta_ring.setdefault((d, b, smp, rg, e), set()).add(round(eta_lzt[x], 6))
assert all(len(v) == 1 for v in eta_ring.values()), "a ring with more than one |eta|"


def ring(x):
    d, s, b, smp, rg, e, p = lar.campos(x)
    if d == "FCAL":
        k = int(np.floor(eta_lzt[x] / 0.1))
        return (lar.layer(d, b, smp), 1000 + k), (round(k * 0.1, 5), round((k + 1) * 0.1, 5))
    ea = next(iter(eta_ring[(d, b, smp, rg, e)]))
    return (lar.layer(d, b, smp), b * 10000 + rg * 1000 + e), (ea, ea)


ring_of, lim = {}, {}
for x in cells:
    k, l_ = ring(x)
    ring_of[x], lim[k] = k, l_
order = sorted(lim, key=lambda k: (k[0], lim[k][0], k[1]))
print(f"{len(cells)} cells, {len(order)} rings", flush=True)

# simulation constants
M = src["lar_mc"]
one = lambda f: lar.file_of(M[f]["refs"])  # noqa: E731
ramp = lar.read_subset(one("Ramp"), "ConditionsContainerLArRampSubset_p1", [("m_vRamp", "f4", "m_vRampSize"), ("m_vRampSize", "size", None)], entry=0)[0][2]
d2u = lar.read_subset(one("DAC2uA"), "ConditionsContainerLArSingleFloatSubset_p1", [("m_values", "f4", 1)])[0][2]
u2m = lar.read_subset(one("uA2MeV"), "CollectionTreeLArSingleFloatSubset_p1", [("m_values", "f4", 1)])[0][2]
mpm = lar.read_subset(one("MphysOverMcal"), "ConditionsContainerLArMphysOverMcalSubset_p1", [("m_values", "f4", 1)], entry=0)[0][2]
hv = {}
for _c, _g, d in lar.read_subset(one("HVScaleCorr"), "ConditionsContainerLArSingleFloatSubset_p1", [("m_values", "f4", 1)]):
    hv.update(d)
fs = lar.read_fsampl(one("fSampl"))
rep = lar.representatives(ramp, rank, off_mc)
fc_of = lambda h: (inv[h // 128], h % 128)  # noqa: E731


def crest(folder):
    tag = src["crest"]["tags"][folder]
    get = lambda u: json.loads(urllib.request.urlopen(f"{src['crest']['server']}/{u}", timeout=300).read())  # noqa: E731
    iov = get(f"iovs?tagname={tag}")
    iov = iov.get("resources", iov)
    assert len(iov) == 1, len(iov)
    return get(f"payloads/data?hash={iov[0]['payloadHash']}")


fb = lambda b: np.frombuffer(base64.b64decode(b), "<f4")  # noqa: E731
N, R4, PED = crest("Noise"), crest("Ramp"), crest("Pedestal")
cr = {"ramp_r4_hg": fb(R4["0"][0]).reshape(-1, 2)[:, 1], "ramp_r4_lg": fb(R4["1"][0]).reshape(-1, 2)[:, 1],
      "noise_adc_r4_hg": fb(N["0"][0]), "noise_adc_r4_lg": fb(N["1"][0])}
print(f"CREST pedestal (value, second field): high {[float(fb(b)[0]) for b in PED['0'][:2]]}, "
      f"low {[float(fb(b)[0]) for b in PED['1'][:2]]}", flush=True)
COLS = ["fsampl_mc", "adc2mev_hg_mc", "ramp_hg_mc", "ramp_r4_hg", "ramp_r4_lg", "noise_adc_r4_hg", "noise_adc_r4_lg"]
val = {k: {c: [] for c in COLS} for k in lim}
for x in cells:
    h, r = hash_of[x], lar.rep_of(rep, x)
    v = {"fsampl_mc": fs[fc_of(h)], "ramp_hg_mc": float(ramp[r]["m_vRamp"][1]),
         "adc2mev_hg_mc": float(ramp[r]["m_vRamp"][1]) * float(d2u[r]["m_values"][0]) * float(u2m[r]["m_values"][0]) / float(mpm[r]["m_values"][0]) * float(hv[fc_of(h)]["m_values"][0])}
    v.update({c: float(cr[c][h]) for c in cr})
    for c in COLS:
        val[ring_of[x]][c].append(v[c])
with open(os.path.join(OUT, "lar_run4_constants_by_ring.csv"), "w", newline="\n") as o:
    o.write("layer,ring,abs_eta_min,abs_eta_max,n_cells," + ",".join(f"{c},{c}_n_off5pct" for c in COLS) + "\n")
    for k in order:
        parts = [k[0], str(k[1]), f"{lim[k][0]:.5f}", f"{lim[k][1]:.5f}", str(len(val[k][COLS[0]]))]
        for c in COLS:
            a = np.array(val[k][c])
            m = float(np.median(a))
            parts += [f"{m:.8g}", str(int(np.sum(np.abs(a - m) > 0.05 * abs(m))))]
        o.write(",".join(parts) + "\n")
print("constants: done", flush=True)

# noise autocorrelation
ac = {g: d for _c, g, d in lar.read_subset(one("AutoCorr"), "CollectionTreeLArAutoCorrSubset_p1", [("m_vAutoCorr", "f4", "m_vAutoCorrSize"), ("m_vAutoCorrSize", "size", None)])}
reps = {}
for x in cells:
    reps.setdefault(ring_of[x], set()).add(lar.rep_of(rep, x))
with open(os.path.join(OUT, "lar_run4_noise_autocorr_by_ring.csv"), "w", newline="\n") as o:
    o.write("layer,ring,abs_eta_min,abs_eta_max,gain,rho1,rho2,rho3,rho4\n")
    for g, name in ((0, "high"), (1, "medium")):
        for k in order:
            v = np.mean([np.asarray(ac[g][r]["m_vAutoCorr"])[:4] for r in reps[k]], axis=0)
            o.write(f"{k[0]},{k[1]},{lim[k][0]:.5f},{lim[k][1]:.5f},{name}," + ",".join(f"{y:.6g}" for y in v) + "\n")
print("autocorrelation: done", flush=True)

# pulses: RTM waves (gains high and medium), FCal completed with the LArPileupShape
WAVE = [("m_dt", "f4", 1), ("m_vAmplitudes", "f4", "m_samples"), ("m_vErrors", "f4", "m_samples"), ("m_vTriggers", "i4", "m_samples"),
        ("m_flag", "u4", 1), ("m_timeOffset", "i4", 1), ("m_samples", "size", None)]
shape = {g: d for _c, g, d in lar.read_subset(lar.file_of(src["pileup_shape"]["refs"]), "ConditionsContainerLArShapeSubset_p1",
                                               [("m_vShape", "f4", "m_vShapeSize"), ("m_vShapeDer", "f4", "m_vShapeDerSize"),
                                                ("m_vShapeSize", "size", None), ("m_vShapeDerSize", "size", None)])}
waves = {0: {}, 1: {}}
by_file = {}
for ch, (path, cnt, entry) in src["rtm"]["refs"].items():
    by_file.setdefault((path, cnt), {})[entry] = int(ch)
for (path, cnt), entries in by_file.items():
    br = uproot.open(path)[cnt][cnt.replace("ConditionsContainer", "").replace("CollectionTree", "")]
    if len(br.keys()) > 0:  # written split in branches (the FCal file)
        rows = [r for i, r in enumerate(lar.read_split(path, cnt)) if i in entries]
        rows = [(c, g, {fc: (dt, w) for fc, (dt, w) in d.items()}) for c, g, d in rows]
    else:
        rows = [(c, g, {fc: (float(v["m_dt"][0]), v["m_vAmplitudes"]) for fc, v in d.items()})
                for i, (c, g, d) in enumerate(lar.read_subset(path, cnt, WAVE)) if i in entries]
    for c, g, d in rows:
        if g not in (0, 1):
            continue
        for (f, ch), (dt, w) in d.items():
            x = int(off_data[128 * rank[int(f)] + ch])
            if x == 0xFFFFFFFF:
                continue
            if lar.campos(x)[0] == "FCAL":
                w = lar.fcal_tail(w, dt, np.asarray(shape[g][lar.rep_of(rep, x)]["m_vShape"], float))
            waves[g][x] = (dt, np.asarray(w, float))
    print(f"{os.path.basename(path)}: {len(entries)} COOL channels", flush=True)
with open(os.path.join(OUT, "lar_run4_pulses_by_ring.csv"), "w", newline="\n") as o:
    o.write("layer,ring,abs_eta_min,abs_eta_max,gain,n_cells,n_used,t_peak_mean_ns,t_peak_std_ns,t_start_ns,dt_ns,n_points,wave\n")
    for g, name in ((0, "high"), (1, "medium")):
        per = {}
        for x in cells:
            if x in waves[g]:
                per.setdefault(ring_of[x], []).append(x)
        for k in order:
            xs, info = per.get(k, []), []
            for x in xs:
                dt, w = waves[g][x]
                w = w[:int(np.isfinite(w).sum())]
                ip = int(np.argmax(w))
                pk = w[ip]
                if pk <= 0 or abs(w[-1] / pk) > 0.05:
                    continue
                info.append((dt, w / pk, ip * dt))
            if info:
                med = np.median([i[2] for i in info])
                info = [i for i in info if abs(i[2] - med) <= 10.0]
            if not info:
                o.write(f"{k[0]},{k[1]},{lim[k][0]:.5f},{lim[k][1]:.5f},{name},{len(xs)},0,nan,nan,0,{STEP},0,\n")
                continue
            st = max(-i[2] for i in info)
            en = min(-i[2] + (len(i[1]) - 1) * i[0] for i in info)
            grid = np.arange(np.ceil(st / STEP) * STEP, np.floor(en / STEP) * STEP + STEP / 2, STEP)
            m = np.mean([np.interp(grid, -i[2] + np.arange(len(i[1])) * i[0], i[1]) for i in info], axis=0)
            m = m - m.sum() / len(m)
            m = m / m[np.argmin(np.abs(grid))]
            tp = np.array([i[2] for i in info])
            o.write(f"{k[0]},{k[1]},{lim[k][0]:.5f},{lim[k][1]:.5f},{name},{len(xs)},{len(info)},{tp.mean():.4f},{tp.std():.4f},"
                    f"{grid[0]:.1f},{STEP},{len(m)}," + " ".join(f"{y:.6g}" for y in m) + "\n")
print(f"pulses: done -> {OUT}")
