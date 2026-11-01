"""Cross-check primecount-library values (default Gourdon) in bigargs.txt against the
primecount CLI with a different algorithm (Deleglise-Rivat by default; --lmo optional).
Usage: python3 crosscheck_bigargs.py bigargs.txt [--lmo] [--max N]
"""
import subprocess, sys

path = sys.argv[1]
alg = "--lmo" if "--lmo" in sys.argv else "--deleglise-rivat"
mx = None
if "--max" in sys.argv:
    mx = int(sys.argv[sys.argv.index("--max") + 1])
vals = {}
for line in open(path):
    y, p = line.split()
    y, p = int(y), int(p)
    if y in vals and vals[y] != p:
        print("INTERNAL INCONSISTENCY", y, vals[y], p)
    vals[y] = p
ys = sorted(vals)
if mx is not None and len(ys) > mx:
    # keep the largest ones plus an even spread
    step = max(1, len(ys) // mx)
    ys = sorted(set(ys[::step] + ys[-20:]))
bad = 0
for i, y in enumerate(ys):
    out = subprocess.run(["primecount", str(y), alg, "--threads=3"], capture_output=True, text=True, check=True).stdout.strip()
    if int(out) != vals[y]:
        bad += 1
        print("MISMATCH", y, vals[y], out, flush=True)
print(f"{alg}: checked {len(ys)} distinct arguments (of {len(vals)}), range {ys[0]}..{ys[-1]}, mismatches {bad}")
