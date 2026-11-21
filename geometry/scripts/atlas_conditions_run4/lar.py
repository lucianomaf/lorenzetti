"""Readers of the ATLAS LAr conditions POOL files of the public CVMFS repository, without Athena, and the links channel -> cell -> ring.

All formats were checked on 21 November 2026 against the values Athena 25.0.21 extracts from the same files (182 468 cells): equal.
- LArConditionsSubset_p1 (Ramp, DAC2uA, uA2MeV, MphysOverMcal, AutoCorr, Shape, HVScaleCorr, PhysWave): front-end boards (FEB, 32-bit
  identifier), boards with sparse data (identifier + 4 words of a 128-channel mask) and the values in board and channel order. Some files
  are read with uproot, the others by hand from the streamer (the uproot reader misplaces some unsplit objects).
- LArConditionsSubset<LArfSamplP> (fSampl): written member-wise, read by hand.
- FEB identifier: 001 | 110 | barrel_ec | pos_neg | feedthrough (5 bits) | slot - 1 (4) | channel (7) | 8 zeros; online hash = 128 x the
  position of the board in the order (barrel_ec, side, feedthrough, slot) + channel.
- cell identifier (32 bits): LArEM, LArHEC and LArFCAL layouts (see campos()).
- symmetrisation of the simulation constants: one representative channel per ring on side A, phi 0 (EM, HEC) or phi 0 to 7 (FCal, period
  8 in phi); a cell takes the one of its ring (FCal: phi mod 8).
"""
import os
import struct

import numpy as np
import uproot

from frontier import q


def cabling(schema_db, tag, run):
    """online hash -> 32-bit cell identifier (0xFFFFFFFF: disconnected) from /LAR/Identifier/OnOffIdMap (blob in the IOV table)"""
    schema, db = schema_db.split("/")
    s = f"ATLAS_{schema}.{db}"
    n = int(q(f"select TO_CHAR(NODE_ID) from {s}_NODES where NODE_FULLPATH='/LAR/Identifier/OnOffIdMap'")[0][0])
    t = f"{s}_F{n:04d}"
    tid = q(f"select TO_CHAR(TAG_ID) from {t}_TAGS where TAG_NAME='{tag}'")[0][0]
    k = run << 32
    r = q(f'select "OnlineHashToOfflineId" from {t}_IOVS where USER_TAG_ID={tid} and IOV_SINCE<={k} and IOV_UNTIL>{k}')
    assert len(r) == 1, len(r)
    return np.frombuffer(r[0][0].encode("latin-1"), "<u4").astype(np.int64)


def febs(pool):
    """{FEB identifier: position in the online hash order}, from the simulation cabling POOL file (LArOnOffIdMapAtlas)"""
    t = uproot.open(pool)["ConditionsContainerLArOnOffId_P_LArOnOffIdMapAtlas"]
    a = {k.split(".")[-1]: np.asarray(v[0]).astype(np.int64) for k, v in t.arrays(library="np").items() if "m_v." in k}
    p = a["pn"]
    be, side = np.where(np.abs(p) == 1, 0, 1), (p > 0).astype(np.int64)
    f = np.unique(np.c_[be, side, a["ft_num"], a["feb_slot"]], axis=0)
    hw = (1 << 29) | (6 << 26) | (f[:, 0] << 25) | (f[:, 1] << 24) | (f[:, 2] << 19) | ((f[:, 3] - 1) << 15)
    return {int(h): i for i, h in enumerate(hw)}


def campos(id32):
    """detector, side (+1/-1), |barrel_ec| (EM: 1 barrel, 2 outer wheel, 3 inner wheel; HEC and FCal: 2), sampling (FCal: module),
    region, eta and phi indices of a 32-bit cell identifier"""
    part = (id32 >> 26) & 7
    if part == 3:
        bec = [-3, -2, -1, 1, 2, 3][(id32 >> 23) & 7]
        return "EM", int(np.sign(bec)), abs(bec), (id32 >> 21) & 3, (id32 >> 18) & 7, (id32 >> 9) & 511, (id32 >> 1) & 255
    if part == 4:
        return "HEC", 1 if (id32 >> 25) & 1 else -1, 2, (id32 >> 23) & 3, (id32 >> 22) & 1, (id32 >> 18) & 15, (id32 >> 12) & 63
    if part == 5:
        return "FCAL", 1 if (id32 >> 25) & 1 else -1, 2, ((id32 >> 23) & 3) + 1, 0, (id32 >> 17) & 63, (id32 >> 13) & 15
    return None


def layer(d, b, smp):
    """layer name: EMB_PS..EMB_3, EMEC_PS..EMEC_3 (inner wheel: dictionary samplings 1 and 2 are layers 2 and 3), HEC_0..3, FCal_1..3"""
    if d == "EM":
        if b == 1:
            return ["EMB_PS", "EMB_1", "EMB_2", "EMB_3"][smp]
        if b == 2:
            return ["EMEC_PS", "EMEC_1", "EMEC_2", "EMEC_3"][smp]
        return ["", "EMEC_2", "EMEC_3"][smp]
    return f"HEC_{smp}" if d == "HEC" else f"FCal_{smp}"


def _channels(feb, sp):
    masc = {int(r[0]): r[1:] for r in np.asarray(sp).reshape(-1, 5)} if len(sp) else {}
    for f in np.asarray(feb).tolist():
        chs = [32 * w + b for w in range(4) for b in range(32) if (int(masc[f][w]) >> b) & 1] if f in masc else range(128)
        for c in chs:
            yield f, c


class _Reader:
    def __init__(self, b):
        self.b, self.p = b, 0

    def head(self):
        n, v = struct.unpack(">IH", self.b[self.p:self.p + 6])
        self.p += 6
        return v

    def u32(self):
        x = struct.unpack(">I", self.b[self.p:self.p + 4])[0]
        self.p += 4
        return x

    def vec(self, fmt):
        self.head()
        n = self.u32()
        x = np.frombuffer(self.b, dtype=">" + fmt, count=n, offset=self.p)
        self.p += 4 * n
        return x


def _entries(pool, tree, branch):
    br = uproot.open(pool)[tree][branch]
    for i in range(br.num_baskets):
        bk = br.basket(i)
        data, offs = bytes(bk.data), bk.byte_offsets
        for k in range(bk.num_entries):
            yield data[offs[k]:offs[k + 1]]


def read_subset(pool, tree, members, entry=None):
    """[(COOL channel, gain, {(FEB, channel): {member: values}})] of a LArConditionsSubset_p1 object.
    members: [(member name, numpy format or 'size', values per channel: 1 or the name of the size member)], after the subset, in streamer
    order. Objects written split in branches are read with uproot; unsplit ones by hand from the streamer (the uproot reader misplaces
    some of them)."""
    t = uproot.open(pool)[tree]
    top = tree.replace("ConditionsContainer", "").replace("CollectionTree", "")
    out = []
    if len(t[top].keys()) > 0:  # split in branches (the data are in the sub-branches)
        a = t.arrays(library="np")
        g = lambda k, i: a[[x for x in a if x.split(".")[-1].split("/")[-1] == k][0]][i]  # noqa: E731
        rng = range(t.num_entries) if entry is None else [entry]
        for i in rng:
            feb, sp = np.asarray(g("m_febIds", i)), np.asarray(g("m_febsWithSparseData", i))
            vals = {n: (int(g(n, i)) if f == "size" else np.asarray(g(n, i))) for n, f, _ in members}
            out.append((int(g("m_channel", i)), int(g("m_gain", i)), _rows(feb, sp, members, vals)))
        return out
    for idx, x in enumerate(_entries(pool, tree, top)):
        if entry is not None and idx != entry:
            continue
        L = _Reader(x)
        if L.head() == 0:
            L.u32()
        feb, _corr, sp = L.vec("u4"), L.vec("u4"), L.vec("u4")
        gain, chan, _group = L.u32(), L.u32(), L.u32()
        vals = {n: (L.u32() if f == "size" else L.vec(f)) for n, f, _ in members}
        out.append((chan, gain, _rows(feb, sp, members, vals)))
    return out


def _rows(feb, sp, members, vals):
    d = {}
    for j, fc in enumerate(_channels(feb, sp)):
        row = {}
        for name, fmt, per in members:
            if fmt == "size":
                continue
            n = vals[per] if isinstance(per, str) else 1
            row[name] = vals[name][j * n:(j + 1) * n]
        d[fc] = row
    return d


def read_split(pool, tree):
    """the same for subsets written split in branches (the FCal RTM file), with uproot: [(COOL channel, gain, {(FEB, ch): (dt, wave)})]"""
    a = uproot.open(pool)[tree].arrays(library="np")
    g = lambda k, i: a[[x for x in a if x.endswith(k)][0]][i]  # noqa: E731
    out = []
    for i in range(len(a[list(a)[0]])):
        ns = int(g("m_samples", i))
        dt, amp = np.asarray(g("m_dt", i)), np.asarray(g("m_vAmplitudes", i))
        d = {fc: (float(dt[j]), amp[j * ns:(j + 1) * ns]) for j, fc in enumerate(_channels(g("m_febIds", i), g("m_febsWithSparseData", i)))}
        out.append((int(g("m_channel", i)), int(g("m_gain", i)), d))
    return out


def read_fsampl(pool):
    """{(FEB, channel): fSampl} of the member-wise LArConditionsSubset<LArfSamplP>: a 12-byte header, the number of boards, their
    identifiers, another 12-byte header and, per board, the count and the float32 values (big-endian)"""
    b = uproot.open(pool)["ConditionsContainerLArConditionsSubset<LArfSamplP>"]["LArConditionsSubset_LArfSamplP_"]
    x = bytes(b.basket(0).data)
    n = struct.unpack(">I", x[12:16])[0]
    boards = np.frombuffer(x, ">u4", count=n, offset=16)
    p = 16 + 4 * n + 12
    out = {}
    for f in boards.tolist():
        c = struct.unpack(">I", x[p:p + 4])[0]
        v = np.frombuffer(x, ">f4", count=c, offset=p + 4)
        for ch in range(c):
            out[(f, ch)] = float(v[ch])
        p += 4 + 4 * c
    return out


def representatives(ramp, feb_rank, offline):
    """{(detector, |barrel_ec|, sampling, region, eta): {phi: (FEB, channel)}}: the representative channels of the symmetrised
    simulation constants, from the channels of the ramp"""
    rep = {}
    for (f, c) in ramp:
        x = int(offline[128 * feb_rank[f] + c])
        if x == 0xFFFFFFFF:
            continue
        d, s, b, smp, rg, e, p = campos(x)
        rep.setdefault((d, b, smp, rg, e), {})[p] = (f, c)
    return rep


def rep_of(rep, x):
    d, s, b, smp, rg, e, p = campos(x)
    m = rep[(d, b, smp, rg, e)]
    return m[0] if len(m) == 1 else m[p % 8]


def fcal_tail(w, dt, s, npts=768):
    """the complete FCal wave: the 156-point RTM wave (1 ns), then straight lines through the 32 samples (25 ns) of the LArPileupShape of
    the channel, placed on the time axis of the wave and scaled to it by the least squares of the samples inside the wave (position: from
    the alignment of the peaks, -5 to +5 ns in 0.05 ns steps); from the last point of the wave to the next sample a straight line; zero one
    sample after the last one"""
    w = np.asarray(w, float)
    n = int(np.isfinite(w).sum())
    w = w[:n]
    tw = np.arange(n) * dt
    ts0 = tw[np.argmax(w)] - 25.0 * np.argmax(s)
    best = None
    for t0 in ts0 + np.arange(-5.0, 5.0001, 0.05):
        tj = t0 + 25.0 * np.arange(len(s))
        inside = (tj >= tw[0]) & (tj <= tw[-1])
        if inside.sum() < 2:
            continue
        y, x = np.interp(tj[inside], tw, w), s[inside]
        a = float(np.dot(x, y) / np.dot(x, x))
        r = float(np.sum((a * x - y) ** 2))
        if best is None or r < best[0]:
            best = (r, t0, a)
    _, t0, a = best
    tj = t0 + 25.0 * np.arange(len(s))
    after = tj > tw[-1]
    tx, yx = np.r_[tw, tj[after], tj[-1] + 25.0], np.r_[w, a * s[after], 0.0]
    t = np.arange(npts) * dt
    return np.where(t <= tx[-1], np.interp(t, tx, yx), 0.0)


def file_of(refs, ch=None):
    """the POOL file of a folder (one file for all its COOL channels, or the file of channel ch)"""
    files = {v[0] for v in refs.values()} if ch is None else {refs[str(ch)][0]}
    assert len(files) == 1, files
    return files.pop()


def exists(path):
    return os.path.exists(path)
