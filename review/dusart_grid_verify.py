"""Referee: independent verification of the primecount-free 'Dusart grid' certificate
(runs/dusart_grid_1e10_2e13.csv), written from the coordinator's description only.

Claim per row i:  hi(x_{i+1}) <= 2 LD(x_i), where
  LD(x) = sum_{k in G} max(0, lo(x/k) - pi(g(k))),
  lo(y) = y/(ln y - 1) for real y >= 5393 (Dusart), exact pi(floor y) for y < 5393,
  hi(y) = y/(ln y - 1.1) for y >= 60184 (Dusart).
Then for n in [x_i, x_{i+1}]: pi(n) <= pi(x_{i+1}) <= hi(x_{i+1}) <= 2 LD(x_i) <= 2 L(x_i) < 2 T(x_i) <= 2 T(n).

Rigour: Decimal at 60 digits, ln correctly rounded (<= 0.5 ulp); every lo term is
multiplied by (1 - 1e-45), hi by (1 + 1e-45), and the final sum is lowered by 1e-20
absolute (accumulated rounding is < 1e-35 absolute).  Generators: the binary file,
which gens_check.c verified record-by-record against two independent enumerations.
Usage: python3 -I dusart_grid_verify.py gens.bin grid.csv
"""
import sys, struct, csv
from decimal import Decimal as D, getcontext

getcontext().prec = 60
gens = open(sys.argv[1], "rb").read()
n = len(gens) // 16
rec = struct.unpack("<%dQ" % (2 * n), gens)
ks, gs = rec[0::2], rec[1::2]

# pi table up to max g (sieve)
M = max(max(gs), 6000)
sieve = bytearray([1]) * (M + 1)
sieve[0] = sieve[1] = 0
for i in range(2, int(M ** 0.5) + 1):
    if sieve[i]:
        sieve[i * i :: i] = bytearray(len(range(i * i, M + 1, i)))
pi = [0] * (M + 1)
c = 0
for i in range(M + 1):
    c += sieve[i]
    pi[i] = c
pig = [pi[g] for g in gs]

ONE_M = D(1) - D("1e-45")
ONE_P = D(1) + D("1e-45")


def LD(x):
    x = int(x)
    s = D(0)
    for k, pg in zip(ks, pig):
        q = x // k
        if q < 5393:
            # y = x/k < 5393 iff floor(x/k) < 5393 (5393 integer)
            t = pi[q] - pg
            if t > 0:
                s += t
        else:
            y = D(x) / D(k)
            lo = y / (y.ln() - 1) * ONE_M
            t = lo - pg
            if t > 0:
                s += t
    return s - D("1e-20")


def hi(y):
    y = D(int(y))
    assert y >= 60184
    return y / (y.ln() - D("1.1")) * ONE_P + D("1e-20")


rows = list(csv.DictReader(open(sys.argv[2])))
bad = 0
minslack = None
for i, r in enumerate(rows):
    x, xn = int(r["x_i"]), int(r["x_next"])
    ld = LD(x)
    h = hi(xn)
    h1 = hi(xn + 1) - D("2e-20")  # lower-ish estimate for x_next+1, only to report maximality
    slack = 2 * ld - h
    ok = slack >= 0
    chain = (i + 1 == len(rows)) or (int(rows[i + 1]["x_i"]) == xn)
    if not (ok and chain):
        bad += 1
    if minslack is None or slack < minslack[0]:
        minslack = (slack, i)
    print(f"{i} x={x} LD_ref={ld:.6f} LD_csv={r['LD_x_i']} x_next={xn} hi_ref={h:.6f} hi_csv={r['hi_pi_x_next']} "
          f"2LD-hi={slack:.6f} next_int_would_fail={2*ld < h1} chain={chain} {'ok' if ok and chain else 'FAIL'}", flush=True)
print("first x:", rows[0]["x_i"], " last x_next:", rows[-1]["x_next"], " >= 2e13:", int(rows[-1]["x_next"]) >= 2 * 10 ** 13)
print("rows:", len(rows), " failures:", bad, " min slack 2LD-hi:", f"{minslack[0]:.6f}", "at row", minslack[1])
