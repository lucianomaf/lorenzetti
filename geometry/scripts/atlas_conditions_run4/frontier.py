"""Read-only SQL queries to the ATLAS COOL conditions database through the public ATLAS Frontier server (no account needed).

Every number in the query must be turned into text (TO_CHAR); the reader parses the text fields of the Frontier BLOB answer.
Usage as a module: from frontier import q; q("select ... from ATLAS_COOLOFL_LAR.CONDBR2_NODES")
"""
import base64
import struct
import sys
import urllib.request
import zlib

URL = "http://atlasfrontier-ai.cern.ch:8000/atlr/Frontier"
_ENC, _DEC = str.maketrans("+/=", ".-_"), str.maketrans(".-_", "+/=")


def _raw(sql):
    enc = base64.b64encode(zlib.compress(sql.encode())).decode().translate(_ENC)
    url = f"{URL}?type=frontier_request:1:DEFAULT&encoding=BLOBzip5&p1={enc}"
    req = urllib.request.Request(url, headers={"X-Frontier-Id": "lorenzetti-atlas-conditions"})
    r = urllib.request.urlopen(req, timeout=120).read().decode("ascii", "replace")
    if "global_error" in r or 'error="0"' not in r:
        raise RuntimeError(r[r.find("<frontier"):][:900])
    d = r.split("<data>")[1].split("</data>")[0]
    return zlib.decompress(base64.b64decode(d.strip().translate(_DEC)))


def q(sql):
    """rows of the query; Frontier BLOB format: type (1 byte), size (4 bytes) and value; type 7 ends a row"""
    b, p, rows, row = _raw(sql), 0, [], []
    while p < len(b):
        t = b[p]
        p += 1
        if t & 0x7F == 7:
            rows.append(row)
            row = []
            continue
        if t & 0x80:
            row.append(None)
            continue
        n = struct.unpack(">I", b[p:p + 4])[0]
        p += 4
        row.append(b[p:p + n].decode("latin-1"))
        p += n
    return rows[1:]  # the first row has the column names and types


if __name__ == "__main__":
    for x in q(sys.argv[1]):
        print(x)
