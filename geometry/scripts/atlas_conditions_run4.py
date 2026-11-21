#!/usr/bin/env python3
"""Rebuild and check the CREST columns of geometry/data/atlas_conditions_run4/lar_run4_constants.csv.

The LAr Phase-II (FEB 2) constants of the ATLAS Run 4 simulation are read from CREST, the public ATLAS
conditions service (https://crest.cern.ch/api-v6.0, no authentication), global tag COND-MC21-SDR-RUN4-06,
tags LARElecCalibMC{Noise,Pedestal,Ramp}-R4-00. Each payload holds one float32 per online channel, in the
order of the LAr online-channel hash; the CSV keeps that hash in the column online_hash.

Usage:
  python3 atlas_conditions_run4.py check  [csv]   compare the CSV with CREST (default csv: the one of the repository)
  python3 atlas_conditions_run4.py rebuild [csv]  rewrite the CREST columns of the CSV from CREST

Only the standard library and numpy are needed. See NOTICE.md in the data folder for the origin and the terms.
"""
import base64
import json
import os
import sys
import urllib.request

import numpy as np

API = "https://crest.cern.ch/api-v6.0"
TAGS = {p: f"LARElecCalibMC{p}-R4-00" for p in ("Noise", "Pedestal", "Ramp")}
HIGH, LOW = "0", "1"  # the payload keys of the two gains (high gain has the larger noise in ADC counts)
CREST_COLUMNS = ("ramp_r4_hg", "ramp_r4_lg", "noise_adc_r4_hg", "noise_adc_r4_lg")
HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(HERE, "..", "data", "atlas_conditions_run4", "lar_run4_constants.csv")


def get_json(url):
    with urllib.request.urlopen(url, timeout=300) as r:
        return json.loads(r.read())


def payload(folder):
    """the payload of one tag (each tag has a single interval of validity)"""
    iovs = get_json(f"{API}/iovs?tagname={TAGS[folder]}")
    iovs = iovs.get("resources", iovs)
    if len(iovs) != 1:
        sys.exit(f"{TAGS[folder]}: expected one interval of validity, found {len(iovs)}")
    return get_json(f"{API}/payloads/data?hash={iovs[0]['payloadHash']}")


def floats(blob):
    return np.frombuffer(base64.b64decode(blob), dtype="<f4")


def crest_values():
    noise, ramp, ped = payload("Noise"), payload("Ramp"), payload("Pedestal")
    values = {
        "noise_adc_r4_hg": floats(noise[HIGH][0]),
        "noise_adc_r4_lg": floats(noise[LOW][0]),
        "ramp_r4_hg": floats(ramp[HIGH][0]).reshape(-1, 2)[:, 1],  # c1 of the ramp (c0 is zero)
        "ramp_r4_lg": floats(ramp[LOW][0]).reshape(-1, 2)[:, 1],
    }
    pedestal = {g: [float(floats(b)[0]) for b in ped[g][:2]] for g in (HIGH, LOW)}
    return values, pedestal


def read_csv(path):
    with open(path) as f:
        header = f.readline().strip().split(",")
        rows = [line.strip().split(",") for line in f if line.strip()]
    return header, rows


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "check"
    path = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_CSV
    header, rows = read_csv(path)
    h = np.array([int(r[header.index("online_hash")]) for r in rows])
    values, pedestal = crest_values()
    print(f"CREST pedestal (value, second field) per gain: high {pedestal[HIGH]}, low {pedestal[LOW]}")
    if mode == "check":
        bad = 0
        for col in CREST_COLUMNS:
            have = np.array([float(r[header.index(col)]) for r in rows])
            want = values[col][h].astype(float)
            n = int(np.sum(~np.isclose(have, want, rtol=1e-6, atol=0)))
            bad += n
            print(f"{col}: {len(rows) - n} of {len(rows)} cells equal to CREST")
        sys.exit(1 if bad else 0)
    if mode == "rebuild":
        for col in CREST_COLUMNS:
            k = header.index(col)
            for r, v in zip(rows, values[col][h]):
                r[k] = f"{v:.8g}"
        with open(path, "w", newline="\n") as f:
            f.write(",".join(header) + "\n")
            f.writelines(",".join(r) + "\n" for r in rows)
        print(f"rewrote the CREST columns of {path}")
        return
    sys.exit(__doc__)


if __name__ == "__main__":
    main()
