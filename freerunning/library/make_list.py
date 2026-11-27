"""The input list of the library from the folders of the productions (../production/README.md): one line per block with a FRHITS,
"sample,production,block,seed,FRHITS,marking", the blocks of each production in order, the productions in the order given (the inclusive ones
first: the order of the lines is the order of the global index). A block of an inclusive production must have its marking
(<folder>/hard/<production>_<block>.hard.txt); the hard sample has none (-).
Usage: python3 make_list.py <output list> <sample>:<production folder>:<seed base> [...] [--max-blocks N]
  e.g. python3 make_list.py list.csv inclusive:<storage>/minbias/a3_inclusivo_v1:100000 inclusive:<storage>/minbias/a3_inclusivo_v2:102000 \
         hard:<storage>/minbias/a3_hard_v1:200000
"""
import os
import sys


def main(argv):
    max_blocks = 0
    if "--max-blocks" in argv:
        i = argv.index("--max-blocks")
        max_blocks = int(argv[i + 1])
        argv = argv[:i] + argv[i + 2:]
    out, specs = argv[0], argv[1:]
    lines, seen = [], []
    for spec in specs:
        sample, folder, seed0 = spec.rsplit(":", 2)
        if sample not in ("inclusive", "hard"):
            raise SystemExit(f"unknown sample {sample}")
        if seen and sample != seen[-1] and sample in seen:
            raise SystemExit(f"the productions of the sample {sample} must be contiguous")
        seen.append(sample)
        folder = folder.rstrip("/")
        prod = os.path.basename(folder)
        pre, suf = prod + "_", ".FRHITS.root"
        blocks = sorted(f[len(pre):-len(suf)] for f in os.listdir(os.path.join(folder, "frhits")) if f.startswith(pre) and f.endswith(suf))
        if max_blocks:
            blocks = blocks[:max_blocks]
        for b in blocks:
            mark = "-"
            if sample == "inclusive":
                mark = os.path.join(folder, "hard", f"{prod}_{b}.hard.txt")
                if not (os.path.exists(mark) and os.path.getsize(mark)):
                    raise SystemExit(f"block without marking: {mark} (mark it first: ../production/, mode mark)")
            lines.append(f"{sample},{prod},{b},{int(seed0) + int(b)},{folder}/frhits/{pre}{b}{suf},{mark}")
        print(f"{prod}: {len(blocks)} blocks")
    with open(out, "w", newline="\n") as o:
        o.write("\n".join(lines) + "\n")
    print(f"{len(lines)} blocks in {out}")


if __name__ == "__main__":
    main(sys.argv[1:])
