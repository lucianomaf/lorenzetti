"""Where every input comes from, resolved without any ATLAS account, and the common paths of the chain.

- WORK (environment LZT_COND_WORK, default ./work): results of the chain.
- CVMFS (environment CVMFS_ROOT, default /cvmfs): the public ATLAS conditions repository atlas-condb.cern.ch, mounted with the CernVM-FS
  client (free, no account).
Run `python3 sources.py` first: it resolves, through the public Frontier, the COOL tag and the POOL file of every folder and writes
WORK/sources.json, which the other scripts read.
"""
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frontier import q  # noqa: E402

WORK = os.path.abspath(os.environ.get("LZT_COND_WORK", "work"))
CVMFS = os.environ.get("CVMFS_ROOT", "/cvmfs")
CATALOG = os.path.join(CVMFS, "atlas-condb.cern.ch", "repo", "conditions", "poolcond")
CREST = "https://crest.cern.ch/api-v6.0"
MC_GLOBAL, MC_RUN = "OFLCOND-MC21-SDR-RUN4-02", 350000  # the ATLAS Run 4 simulation with the COOL global tag
DATA_GLOBAL, DATA_RUN = "CONDBR2-BLKPA-2026-01", 520451  # the run of the RTM physics waves


def node(schema, db, path):
    r = q(f"select TO_CHAR(NODE_ID) from ATLAS_{schema}.{db}_NODES where NODE_FULLPATH='{path}'")
    assert r, f"{schema}/{db} {path}: no such folder"
    return int(r[0][0])


def tag_in_global(schema, db, path, glob):
    """the tag of a folder in a global tag, down the COOL tag hierarchy (TAG2TAG)"""
    s = f"ATLAS_{schema}.{db}"
    cur = [q(f"select TO_CHAR(TAG_ID) from {s}_TAGS where TAG_NAME='{glob}'")[0][0]]
    parts = path.strip("/").split("/")
    for i in range(1, len(parts) + 1):
        n = node(schema, db, "/" + "/".join(parts[:i]))
        cur = sum(([x[0] for x in q(f"select TO_CHAR(CHILD_TAGID) from {s}_TAG2TAG where PARENT_TAGID={t} and CHILD_NODEID={n}")]
                   for t in cur), [])
    assert len(cur) == 1, (path, cur)
    n = node(schema, db, path)
    return q(f"select TAG_NAME from {s}_F{n:04d}_TAGS where TAG_ID={cur[0]}")[0][0]


def payload(schema, db, path, tag, run, column):
    """{COOL channel: payload column} valid at run in a tag: per channel the interval of validity that starts last before run (the
    current version of a multi-version folder); the payload column is in the IOV table or in a separate payload table"""
    s = f"ATLAS_{schema}.{db}"
    n = node(schema, db, path)
    t = f"{s}_F{n:04d}"
    tid = q(f"select TO_CHAR(TAG_ID) from {t}_TAGS where TAG_NAME='{tag}'")[0][0]
    k = run << 32
    sep = bool(q(f"select COLUMN_NAME from ALL_TAB_COLUMNS where OWNER='ATLAS_{schema}' and TABLE_NAME='{db}_F{n:04d}_IOVS' "
                 f"and COLUMN_NAME='PAYLOAD_ID'"))
    col, extra, join = (f'p."{column}"', f", {t}_PAYLOAD p", " and p.PAYLOAD_ID = i.PAYLOAD_ID") if sep else (f'i."{column}"', "", "")
    rows = q(f"select TO_CHAR(i.CHANNEL_ID), TO_CHAR(i.IOV_SINCE), TO_CHAR(i.OBJECT_ID), {col} from {t}_IOVS i{extra} "
             f"where i.USER_TAG_ID={tid} and i.IOV_SINCE<={k} and i.IOV_UNTIL>{k}{join}")
    best = {}
    for ch, since, obj, val in rows:
        key = (int(since), int(obj))
        if ch not in best or key > best[ch][0]:
            best[ch] = (key, val)
    return {int(ch): v for ch, (key, v) in best.items()}


def pool_file(guid):
    """the path of a POOL file in the public CVMFS repository, from its PoolFileCatalog"""
    for name in sorted(os.listdir(CATALOG)):
        if not name.startswith("PoolFileCatalog.xml"):
            continue
        txt = open(os.path.join(CATALOG, name), errors="replace").read()
        m = re.search(rf'<File ID="{guid}">\s*<physical>\s*<pfn filetype="[^"]*" name="([^"]+)"', txt)
        if m:
            return m.group(1).replace("/cvmfs", CVMFS, 1)
    raise FileNotFoundError(f"POOL file {guid} is not in the CVMFS catalogs")


def pool_refs(refs):
    """{COOL channel: (file, container, entry)} from PoolRef strings"""
    out = {}
    for ch, ref in refs.items():
        guid = re.search(r"\[DB=([0-9A-F-]+)\]", ref).group(1)
        cnt = re.search(r"\[CNT=([^\]]+)\]", ref).group(1).split("(")[0]
        entry = int(re.search(r"\[OID=[0-9A-F]+-([0-9A-F]+)\]", ref).group(1), 16) & 0xFFFFFFFF
        out[ch] = (pool_file(guid), cnt, entry)
    return out


def resolve():
    src = {"mc_global": MC_GLOBAL, "mc_run": MC_RUN, "data_global": DATA_GLOBAL, "data_run": DATA_RUN, "lar_mc": {}}
    # LAr simulation folders (COOLOFL_LAR / OFLP200), tags from the global tag
    for f in ("fSampl", "DAC2uA", "uA2MeV", "MphysOverMcal", "Ramp", "HVScaleCorr", "AutoCorr"):
        path = f"/LAR/ElecCalibMC/{f}"
        tag = tag_in_global("COOLOFL_LAR", "OFLP200", path, MC_GLOBAL)
        refs = pool_refs(payload("COOLOFL_LAR", "OFLP200", path, tag, MC_RUN, "PoolRef"))
        src["lar_mc"][f] = {"folder": path, "tag": tag, "refs": refs}
    # LAr cabling of the simulation: Athena 25.0.21 LArCabling/python/LArCablingConfig.py sets this tag for MC ("MC folder-tag hack")
    src["lar_cabling_mc"] = {"folder": "/LAR/Identifier/OnOffIdMap", "schema": "COOLOFL_LAR/OFLP200", "tag": "LARIdentifierOnOffIdMap-012"}
    # the list of front-end boards, from the simulation cabling POOL file of the online folder
    tag = "LARIdentifierOnOffIdAtlas-013-MC"
    src["lar_febs"] = {"folder": "/LAR/Identifier/OnOffIdAtlas", "tag": tag,
                       "refs": pool_refs(payload("COOLONL_LAR", "OFLP200", "/LAR/Identifier/OnOffIdAtlas", tag, MC_RUN, "PoolRef"))}
    # LAr cabling of the data, for the RTM waves (from the data global tag)
    tag = tag_in_global("COOLONL_LAR", "CONDBR2", "/LAR/Identifier/OnOffIdMap", DATA_GLOBAL)
    src["lar_cabling_data"] = {"folder": "/LAR/Identifier/OnOffIdMap", "schema": "COOLONL_LAR/CONDBR2", "tag": tag}
    # RTM physics waves and the 25 ns shape that completes the FCal ones
    for key, path, tag in (("rtm", "/LAR/ElecCalibOfl/PhysWaves/RTM", "LARElecCalibOflPhysWavesRTM-RUN2-UPD3-00"),
                           ("pileup_shape", "/LAR/ElecCalibOfl/LArPileupShape", "LARElecCalibOflLArPileupShape-RUN2-2018")):
        src[key] = {"folder": path, "tag": tag, "refs": pool_refs(payload("COOLOFL_LAR", "CONDBR2", path, tag, DATA_RUN, "PoolRef"))}
    # LAr Phase-II front-end board (FEB 2) constants of the Run 4 simulation in CREST
    src["crest"] = {"server": CREST, "global": "COND-MC21-SDR-RUN4-06",
                    "tags": {p: f"LARElecCalibMC{p}-R4-00" for p in ("Noise", "Pedestal", "Ramp")}}
    # Tile sampling fraction: Athena 25.0.21 TileConditions/python/TileSamplingFractionConfig.py forces this tag for Run 4
    src["tile_sfr"] = {"folder": "/TILE/OFL02/CALIB/SFR", "schema": "COOLOFL_TILE/OFLP200", "tag": "TileOfl02CalibSfr-SIM-07", "run": MC_RUN}
    # Tile sample noise of the simulation (pedestal and the double Gaussian of TileDigitsMaker), tag from the global tag
    path = "/TILE/OFL02/NOISE/SAMPLE"
    src["tile_noise"] = {"folder": path, "schema": "COOLOFL_TILE/OFLP200", "tag": tag_in_global("COOLOFL_TILE", "OFLP200", path, MC_GLOBAL),
                         "run": MC_RUN}
    os.makedirs(WORK, exist_ok=True)
    json.dump(src, open(os.path.join(WORK, "sources.json"), "w"), indent=1)
    return src


def load():
    return json.load(open(os.path.join(WORK, "sources.json")))


if __name__ == "__main__":
    s = resolve()
    for k, v in s.items():
        if isinstance(v, dict) and "refs" in v:
            files = sorted({x[0] for x in v["refs"].values()})
            print(f"{k}: {v['folder']} {v['tag']}, {len(v['refs'])} COOL channels in {len(files)} files")
        elif k == "lar_mc":
            for f, w in v.items():
                print(f"lar_mc {f}: {w['tag']}, {len(w['refs'])} COOL channels")
    print(f"-> {os.path.join(WORK, 'sources.json')}")
