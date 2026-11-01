"""dusart_grid.py: cover [X0, X1] without any prime-counting software.

For the generator set G (k <= 28,653,696, 2 tau(k) | k) and real x:
    T(x) > L(x) >= LD(x) = sum_{k in G} max(0, lo_pi(x/k) - pi(g(k))),
where lo_pi(y) = y/(ln y - 1) for y >= 5393 (Dusart: pi(y) >= y/(ln y - 1) there),
and lo_pi(y) = pi(floor(y)) from our own sieve for y < 5393.
Upper side: pi(y) <= hi_pi(y) = y/(ln y - 1.1) for y >= 60184 (Dusart).
Grid: x_{i+1} = the largest integer with hi_pi(x_{i+1}) <= 2 LD(x_i); then for n in
[x_i, x_{i+1}]: pi(n) <= hi_pi(x_{i+1}) <= 2 LD(x_i) < 2 T(x_i) <= 2 T(n).

Arithmetic: Decimal at 50 digits; ln k is computed once per generator.  Each LD term is
computed to relative error far below 1e-40; the check subtracts a margin of 1e-30 * LD
and adds 1e-30 * hi_pi, which dominates all rounding.

Usage: python dusart_grid.py gens_28653696.bin X0 X1 > dusart_grid.csv
"""
import sys
from decimal import Decimal, getcontext

import numpy as np

getcontext().prec = 50
ONE = Decimal(1)
C_LO = Decimal(1)
C_HI = Decimal("1.1")
MARGIN = Decimal(10) ** -30


def main(path, X0, X1):
    rec = np.fromfile(path, dtype="<u8").reshape(-1, 2)
    ks = [int(v) for v in rec[:, 0]]
    gs = rec[:, 1].astype(np.int64)
    gmax = int(gs.max())
    lim = max(gmax, 5393) + 1
    is_p = np.ones(lim + 1, dtype=bool)
    is_p[:2] = False
    for p in range(2, int(lim ** 0.5) + 1):
        if is_p[p]:
            is_p[p * p::p] = False
    pi_tab = np.cumsum(is_p)
    pig = [int(v) for v in pi_tab[gs]]
    lnk = [Decimal(k).ln() for k in ks]
    Dk = [Decimal(k) for k in ks]

    def LD(x):
        X = Decimal(x)
        lx = X.ln()
        total = Decimal(0)
        for k, dk, lk, pg in zip(ks, Dk, lnk, pig):
            if x < 5393 * k:
                y = x // k
                v = int(pi_tab[y]) - pg if y <= lim else None
                if v is None:
                    raise RuntimeError("table too small")
                if v > 0:
                    total += v
                continue
            ly = lx - lk  # ln(x/k)
            v = (X / dk) / (ly - C_LO) - pg
            if v > 0:
                total += v
        return total * (ONE - MARGIN)

    def hi_pi(y):
        Y = Decimal(y)
        return Y / (Y.ln() - C_HI) * (ONE + MARGIN)

    print("i,x_i,LD_x_i,x_next,hi_pi_x_next,ratio")
    x = X0
    i = 0
    while True:
        assert x >= 60184
        L2 = 2 * LD(x)
        if hi_pi(x) > L2:
            print(f"FAIL at x={x}: hi_pi(x) > 2 LD(x)", file=sys.stderr)
            sys.exit(2)
        lo, hi = x, 2 * x
        if hi_pi(hi) <= L2:
            lo = hi
        else:
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if hi_pi(mid) <= L2:
                    lo = mid
                else:
                    hi = mid
        y = lo
        print(f"{i},{x},{L2 / 2:.6f},{y},{hi_pi(y):.6f},{y / x:.6f}", flush=True)
        if y <= x:
            print("no progress", file=sys.stderr)
            sys.exit(3)
        if y >= X1:
            break
        x = y
        i += 1


if __name__ == "__main__":
    main(sys.argv[1], int(float(sys.argv[2])), int(float(sys.argv[3])))
