"""Compare primecount-free Lucy results with both grid CSVs.
primecount grid: L(x_i) and pi(x_{i+1}) must equal the CSV; 2L(x_i) - pi(x_{i+1}) >= 0.
Dusart grid (if Lucy values exist): exact L(x_i) >= LD(x_i) and exact pi(x_{i+1}) <= hi(x_{i+1})
(consistency of the Dusart bounds at the points actually used).
Usage: python3 -I compare_lucy.py lucy_out_*.txt
"""
import csv, re, sys

vals = {}
for f in sys.argv[1:]:
    for line in open(f):
        m = re.match(r"x=(\d+) pi=(\d+) L=(\d+)", line)
        if m:
            x, p, L = map(int, m.groups())
            vals[x] = (p, L)
rows = list(csv.DictReader(open("../runs/grid_1e9_1e16.csv")))
done = bad = 0
for r in rows:
    x, L, xn, pn = int(r["x_i"]), int(r["L_x_i"]), int(r["x_next"]), int(r["pi_x_next"])
    if x in vals and xn in vals:
        done += 1
        Lr, pr = vals[x][1], vals[xn][0]
        ok = (Lr == L) and (pr == pn) and (2 * Lr - pr >= 0)
        if not ok:
            bad += 1
            print("MISMATCH row", r["i"], x, L, Lr, xn, pn, pr)
print(f"primecount grid: rows checked by Lucy {done}/{len(rows)}, mismatches {bad}")
drows = list(csv.DictReader(open("../runs/dusart_grid_1e10_2e13.csv")))
dd = db = 0
for r in drows:
    x, xn = int(r["x_i"]), int(r["x_next"])
    LD, hi = float(r["LD_x_i"]), float(r["hi_pi_x_next"])
    if x in vals and xn in vals:
        dd += 1
        Lr, pr = vals[x][1], vals[xn][0]
        if not (Lr >= LD and pr <= hi):
            db += 1
            print("DUSART INCONSISTENCY row", r["i"], x, Lr, LD, xn, pr, hi)
        else:
            print(f"dusart row {r['i']}: L/LD = {Lr/LD:.6f}, pi/hi = {pr/hi:.6f}, true 2L-pi = {2*Lr-pr}")
print(f"dusart grid: rows checked {dd}/{len(drows)}, inconsistencies {db}")
for x in sorted(vals):
    if x in (10**9, 10**11, 10**12, 10**13, 10**14, 10**15, 10**16):
        p, L = vals[x]
        print(f"x=1e{len(str(x))-1}: L={L} pi={p} 2L/pi={2*L/p:.4f}")
