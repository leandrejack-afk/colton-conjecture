"""tail_check.py: exact check of the tail lemma for Colton's conjecture.

Claim.  Let G be the set of generators k <= K = 28,653,696 (2 tau(k) | k), with
g(k) the largest prime factor of k, A = sum_{k in G} 1/k, B = sum_{k in G} pi(g(k)).
For every real x >= X1:   2 L(x) > pi(x),  where  L(x) = sum_{k in G} max(0, pi(x/k) - pi(g(k))) <= T(x).

Proof used (Dusart: pi(y) >= y/(ln y - 1) for y >= 5393, pi(y) <= y/(ln y - 1.1) for y >= 60184):
for x >= X1, every x/k >= X1/K >= 5393, so
    2 L(x) - pi(x) >= x [2A/(l-1) - 1/(l-1.1)] - 2B = F(x) - 2B,   l = ln x,
    F(x) = x (A(2l - 2.2) - l + 1) / ((l - 1)(l - 1.1)).
F is increasing once N = A(2l-2.2) - l + 1 > 0, i.e. for l > 7.35 (then
d ln F/dl = 1 + (2A-1)/N - 1/(l-1) - 1/(l-1.1) > 0), so it suffices that F(X1) > 2B.

This script (independent of gens.c and lbound.c):
  1. recomputes tau(k) for every k <= K with a divisor-count sieve and checks that
     the generator file holds exactly the generators <= K, with the right g(k);
  2. computes A as an exact rational lower bound and B exactly (own prime sieve);
  3. checks F(X1) > 2B for several X1: ln X1 is bracketed by +-1e-70, the other
     Decimal operations round to nearest at 80 digits (errors of order 1e-78, far
     below the margin F/2B - 1 >= 0.27 at X1 = 2e13).

Usage: python tail_check.py gens_28653696.bin
"""
import sys
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np

K = 28653696


def main(path):
    rec = np.fromfile(path, dtype="<u8").reshape(-1, 2)
    ks = rec[:, 0].astype(np.int64)
    gs = rec[:, 1].astype(np.int64)
    print(f"records: {len(ks)}, max k {ks.max()}, max g {gs.max()}")
    assert np.all(np.diff(ks) > 0), "k not strictly increasing"
    assert ks.max() <= K

    # 1. tau(n) for all n <= K by the divisor-count sieve (no factorization)
    tau = np.zeros(K + 1, dtype=np.int32)
    for d in range(1, K + 1):
        tau[d::d] += 1
    n = np.arange(K + 1, dtype=np.int64)
    gen_mask = np.zeros(K + 1, dtype=bool)
    gen_mask[1:] = (n[1:] % (2 * tau[1:].astype(np.int64))) == 0
    all_gens = np.nonzero(gen_mask)[0]
    print(f"generators <= K by divisor-count sieve: {len(all_gens)}")
    assert len(all_gens) == len(ks) and np.array_equal(all_gens, ks), "generator list differs"

    # largest prime factor by a largest-prime-factor sieve
    lpf = np.zeros(K + 1, dtype=np.int64)
    for p in range(2, K + 1):
        if lpf[p] == 0:  # p is prime: no smaller prime has marked it
            lpf[p::p] = p  # primes visited in increasing order: last write is the largest
    assert np.array_equal(lpf[ks], gs), "g(k) differs"
    print("generator list and g(k) confirmed")

    # 2. A (exact lower bound) and B (exact)
    SCALE = 10 ** 30
    S = sum(SCALE // int(k) for k in ks)
    A_low = Fraction(S, SCALE)
    print(f"A >= {float(A_low):.12f}  (A >= 0.508: {A_low >= Fraction(508, 1000)})")
    gmax = int(gs.max())
    is_p = np.ones(gmax + 1, dtype=bool)
    is_p[:2] = False
    for p in range(2, int(gmax ** 0.5) + 1):
        if is_p[p]:
            is_p[p * p::p] = False
    pi_tab = np.cumsum(is_p)
    B = int(pi_tab[gs].sum())
    print(f"B = sum pi(g(k)) = {B}")

    # 3. tail inequality
    getcontext().prec = 80
    A = Decimal(A_low.numerator) / Decimal(A_low.denominator)  # rounding error ~1e-80, absorbed below
    A = A - Decimal(10) ** -70
    for X1 in [10 ** 13, 2 * 10 ** 13, 5 * 10 ** 13, 10 ** 14, 10 ** 15, 10 ** 16]:
        assert Decimal(X1) / Decimal(K) >= 5393 and X1 >= 60184
        l = Decimal(X1).ln()
        eps = Decimal(10) ** -70
        l_lo, l_hi = l - eps, l + eps
        num = A * (2 * l_lo - Decimal("2.2")) - l_lo + 1  # increasing in l (2A - 1 > 0)
        den = (l_hi - 1) * (l_hi - Decimal("1.1"))       # increasing in l
        assert num > 0
        F_low = Decimal(X1) * num / den
        ok = F_low > 2 * B
        print(f"X1={X1:.0e}: ln X1={float(l):.6f}, F(X1) >= {F_low:.6e}, 2B = {2 * B:.6e}, "
              f"ratio {float(F_low / (2 * B)):.3f}, holds: {ok}")


if __name__ == "__main__":
    main(sys.argv[1])
