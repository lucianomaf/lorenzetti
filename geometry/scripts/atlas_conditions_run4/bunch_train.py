"""The HL-LHC bunch pattern at the ATLAS interaction point (IP1), from a public filling scheme. Output: WORK/tables/hllhc_bunch_train_ip1.csv,
the 3564 bunch-crossing slots of one LHC orbit as runs: first slot, number of slots, 1 if the slots have colliding bunches at IP1 else 0.

Source: the HL-LHC baseline 25 ns filling scheme of the standard production mode, 25ns_2760b_2748_2492_2574_288bpi_13inj_800ns_bs200ns
(2760 bunches per beam, 2748 colliding pairs at IP1 and IP5; HL-LHC Technical Design Report, CERN-2020-010, Table 2-1), as the JSON file of
the LHC Programme Coordination filling scheme tool published in the examples of the FillingPatterns package
(https://github.com/giadarol/FillingPatterns, readable without login). Each beam is a list of 3564 slots, 1 for a filled bucket; at IP1
the bunch of beam 1 in slot i meets the bunch of beam 2 in slot i. The repository declares no licence, so the JSON file is not copied: this
script downloads it, writes the table and, when the table of the Lorenzetti repository is given, checks that they are equal.
Usage: python3 bunch_train.py [geometry/data/atlas_conditions_run4/hllhc_bunch_train_ip1.csv]
"""
import json
import os
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sources import WORK  # noqa: E402

NAME = "25ns_2760b_2748_2492_2574_288bpi_13inj_800ns_bs200ns"
URL = f"https://raw.githubusercontent.com/giadarol/FillingPatterns/master/examples/{NAME}.json"

d = json.load(urllib.request.urlopen(URL, timeout=60))
b1, b2 = d["beam1"], d["beam2"]
assert len(b1) == len(b2) == 3564, (len(b1), len(b2))
coll = [int(x and y) for x, y in zip(b1, b2)]
print(f"{NAME}: {sum(b1)} and {sum(b2)} bunches, {sum(coll)} colliding at IP1")
runs, i = [], 0
while i < len(coll):
    j = i
    while j < len(coll) and coll[j] == coll[i]:
        j += 1
    runs.append((i, j - i, coll[i]))
    i = j
lines = ["first_slot,n_slots,colliding_ip1\n"] + [f"{a},{n},{c}\n" for a, n, c in runs]
os.makedirs(os.path.join(WORK, "tables"), exist_ok=True)
out = os.path.join(WORK, "tables", "hllhc_bunch_train_ip1.csv")
open(out, "w", newline="\n").writelines(lines)
print(f"{len(runs)} runs -> tables/hllhc_bunch_train_ip1.csv")
if len(sys.argv) > 1:
    same = open(sys.argv[1]).read() == "".join(lines)
    print(f"table of the repository {'equal' if same else 'DIFFERENT'}")
    sys.exit(0 if same else 1)
