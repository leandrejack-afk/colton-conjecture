"""Referee check of Part 3 (analytic tail) and the corollary, independent of src/.

Uses Python's decimal module (ln is correctly rounded) at 60 digits, and interval
style margins.  A is enclosed by the integer sums printed by gens_check.c:
    10^30 * A in [508000047238841733452073056315, 508000047238841733452073470276].
B = 3348506095 (gens_check.c).
"""
from decimal import Decimal as D, getcontext, ROUND_FLOOR, ROUND_CEILING
import math

getcontext().prec = 60
A_LO = D(508000047238841733452073056315) / D(10) ** 30
A_HI = D(508000047238841733452073470276) / D(10) ** 30
B = D(3348506095)
K = 28653696


def F(x, A):
    """F(x) = 2A x/(l-1) - x/(l-1.1), l = ln x."""
    x = D(x)
    l = x.ln()
    return 2 * A * x / (l - 1) - x / (l - D("1.1"))


def F_readme(x, A):
    x = D(x)
    l = x.ln()
    return x * (A * (2 * l - D("2.2")) - l + 1) / ((l - 1) * (l - D("1.1")))


print("5393*K =", 5393 * K, " (x/k >= 5393 for all k<=K needs x >= this)")
for x in [D("1.5e13"), D("1.55e13"), D("1.6e13"), D("2e13"), D("1e14"), D("1e16")]:
    f1, f2 = F(x, A_LO), F_readme(x, A_LO)
    print(f"x={x:.3e} F={f1:.6e} (README form {f2:.6e}) 2B={2*B:.6e} F/2B={f1/(2*B):.6f}")

# threshold: smallest real x with F(x) >= 2B, by bisection (F increasing there)
for A, name in [(A_LO, "A_lo"), (D("0.508"), "A=0.508")]:
    lo, hi = D("1e13"), D("2e13")
    for _ in range(200):
        mid = (lo + hi) / 2
        if F(mid, A) > 2 * B:
            hi = mid
        else:
            lo = mid
    print(f"threshold with {name}: F(x) > 2B iff x > {lo:.12e}")

# monotonicity: N(l) = (2A-1) l - 2.2A + 1 > 0 iff l > (2.2A-1)/(2A-1)
for A, name in [(A_LO, "A_lo"), (D("0.508"), "0.508")]:
    lstar = (D("2.2") * A - 1) / (2 * A - 1)
    print(f"N(l)>0 iff l > {lstar:.10f} ({name}); e^l = {lstar.exp():.6e}")
    # derivative of ln F wrt l: 1 + (2A-1)/N - 1/(l-1) - 1/(l-1.1); check at a grid of l >= 7.4
    worst = None
    l = D("7.4")
    while l < 60:
        N = (2 * A - 1) * l - D("2.2") * A + 1
        d = 1 + (2 * A - 1) / N - 1 / (l - 1) - 1 / (l - D("1.1"))
        if worst is None or d < worst[0]:
            worst = (d, l)
        l += D("0.01")
    print(f"  min of dlnF/dl on l in [7.4,60] (step .01): {worst[0]:.6f} at l={worst[1]}")

# Dusart thresholds vs use
x0 = D("2e13")
print("x/k at x=2e13, k=K:", x0 / K, ">= 5393:", x0 / K >= 5393, "; x >= 60184:", x0 >= 60184)

# Corollary: T(x) >= x/(2 ln x) for real x >= 2.
# For x >= 599: Dusart (6.5) lower bound pi(x) >= x/ln x (1 + 1/ln x) > x/ln x, so T >= pi/2 > x/(2 ln x).
# For 2 <= x < 599: direct check on each [n, n+1): need T(n) >= sup f, f(x) = x/(2 ln x).
def tau(n):
    c, d = 0, 1
    while d * d <= n:
        if n % d == 0:
            c += 1 if d * d == n else 2
        d += 1
    return c

T = 0
P = 0
bad = []
badRS = []
for n in range(1, 700):
    t = tau(n)
    if n % t == 0:
        T += 1
    if t == 2:
        P += 1
    if n >= 2:
        f = lambda x: x / (2 * math.log(x))
        sup = max(f(n), f(n + 1))  # f decreasing on [2,e], increasing after; n=2 handled by f(2)
        if not (T >= sup + 1e-12):
            bad.append((n, T, sup))
    if n >= 17:
        g = lambda x: x / math.log(x)
        # pi(x) > x/ln x on [n, n+1): need P > sup g = g(n+1) (left limit)
        if not (P > g(n + 1) + 1e-12):
            badRS.append((n, P, g(n + 1)))
print("corollary small range failures (2<=x<700):", bad[:5], "count", len(bad))
print("pi(x) > x/ln x failures on 17<=x<700:", badRS[:5], "count", len(badRS))
# Colton's literal 'less than x' version at x = 2: #{tau numbers < 2} = 1 vs 2/(2 ln 2)
print("count of tau numbers < 2 = 1, 2/(2 ln 2) =", 2 / (2 * math.log(2)))
