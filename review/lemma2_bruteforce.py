"""Brute-force check of Lemmas 1-2 at x = 2*10^6: every kp (k generator, p prime,
g(k) < p <= x/k) is a tau number, all are distinct, and their number equals L(x);
also T(x) >= 1 + L(x) and the tau numbers NOT of this form are listed by count."""
import struct, sys
X = 2 * 10**6
d = open(sys.argv[1], 'rb').read(); n = len(d) // 16
rec = struct.unpack('<%dQ' % (2 * n), d); ks, gs = rec[0::2], rec[1::2]
tau = [0] * (X + 1)
for a in range(1, X + 1):
    for m in range(a, X + 1, a):
        tau[m] += 1
isp = [t == 2 for t in tau]
pi = [0] * (X + 1); c = 0
for i in range(X + 1):
    c += isp[i]; pi[i] = c
prods = []
for k, g in zip(ks, gs):
    if k * 2 > X: break
    for p in range(g + 1, X // k + 1):
        if isp[p]: prods.append(k * p)
assert all(tau[m] and m % tau[m] == 0 for m in prods), "a product is not a tau number"
L = sum(max(0, pi[X // k] - pi[g]) for k, g in zip(ks, gs) if k <= X)
T = sum(1 for m in range(1, X + 1) if m % tau[m] == 0)
print("x", X, "pairs", len(prods), "distinct", len(set(prods)), "L", L, "T", T, "T-1-L", T - 1 - L)
