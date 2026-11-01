# Colton's conjecture on refactorable numbers holds for every n

By Leandre Jack, with programs written and run with the assistance of a large language model.

A positive integer n is a *tau number* (refactorable number) if its number of divisors, tau(n), divides n: 1, 2, 8, 9, 12, 18, 24, ... ([OEIS A033950](https://oeis.org/A033950)). Write T(x) for the number of tau numbers up to x and pi(x) for the number of primes up to x.

In 1999 Simon Colton observed from his tables that "the number of refactorables is always at least half the number of primes" ([J. Integer Seq. 2 (1999), Article 99.1.2](https://cs.uwaterloo.ca/journals/JIS/colton/joisol.html)). In 2002 Joshua Zelinsky proved this for every n > 7.42 x 10^13 ([J. Integer Seq. 5 (2002), Article 02.2.8](https://cs.uwaterloo.ca/journals/JIS/VOL5/Zelinsky/zelinsky9.html), Theorem 8); the range below was left open, and the statement is listed as open in [google-deepmind/formal-conjectures](https://github.com/google-deepmind/formal-conjectures/blob/main/FormalConjectures/OEIS/111291.lean) (`colton_conjecture`).

## Result

**Theorem.** For every integer n >= 1, pi(n) <= 2 T(n), with equality only at n = 7 (pi(7) = 4, T(7) = 2).

Consequently T(x) >= x / (2 log x) for every real x >= 2: for x >= 599 by Dusart's bound pi(x) >= (x/ln x)(1 + 1/ln x) ([arXiv:1002.0442](https://arxiv.org/abs/1002.0442), (6.5)), and below 599 by direct check (`review/part3_check.py`). Colton phrased this form with "less than x"; read that way it holds for every real x > 2 and fails at x = 2.

## Proof

**Generators.** Call k a generator if 2 tau(k) divides k (Zelinsky, Lemma 6): 8, 12, 36, 72, 80, 96, 128, 180, ... For a generator k and any prime p not dividing k, tau(kp) = 2 tau(k) divides kp, so kp is a tau number.

**Lower bound.** Let G be a finite set of generators and g(k) the largest prime factor of k. The pairs (k, p) with k in G, p prime and g(k) < p <= x/k give distinct tau numbers kp <= x (p is the largest prime factor of kp and occurs once, so kp determines p and k), and none of them is 1. Hence

    T(x) >= 1 + L(x),   L(x) = sum over k in G of max(0, pi(x/k) - pi(g(k))).

We use Zelinsky's own set: G = all 413,980 generators k <= 28,653,696 (28,653,696 = 2^7 x 3^3 x 8291 is itself one), with A = sum 1/k = 0.508000047... and B = sum pi(g(k)) = 3,348,506,095.

**Part 1, n <= 10^10 (exact).** A segmented sieve computes tau(n) for every n <= 10^10, hence T(n) and pi(n) exactly; 2T(n) - pi(n) is evaluated at every prime, the only places it can drop. Its minimum is 0, attained only at n = 7 (`results/pointwise_1e10.txt`). The same program run to 10^12 gives the same minimum and T(10^12) = 26,568,322,474, a value not previously published (`results/pointwise_1e12.txt`, produced before the wording of its last output line was changed; a second, independent count to 10^12 with the referee's divisor-pair program is in progress).

**Part 2, 10^10 <= n <= 2.07 x 10^13 (grid, published prime bounds only).** Dusart's bounds pi(y) >= y/(ln y - 1) for y >= 5393 and pi(y) <= y/(ln y - 1.1) for y >= 60184 ([arXiv:1002.0442](https://arxiv.org/abs/1002.0442), Theorem 6.9, (6.6); Ramanujan J. 45 (2018) 227-251; independently, both were checked directly for every real y up to 2 x 10^10 (`review/dusart_numeric.c`) and follow from Axler, [arXiv:1703.08032](https://arxiv.org/abs/1703.08032), Theorems 1.2 and 1.3, for y >= 1.91 x 10^10) give LD(x) = sum over k in G of max(0, lo(x/k) - pi(g(k))) <= L(x), with lo(y) = y/(ln y - 1) for y >= 5393 and lo(y) = pi(y) from a small sieve below that. 26 points x_0 = 10^10 < ... < x_25 = 20,723,943,136,807 satisfy x_{i+1}/(ln x_{i+1} - 1.1) <= 2 LD(x_i) (50-digit decimal arithmetic with a 10^-30 relative margin on each side). For n in [x_i, x_{i+1}]: pi(n) <= x_{i+1}/(ln x_{i+1} - 1.1) <= 2 LD(x_i) < 2 T(x_i) <= 2 T(n). `src/dusart_grid.py`, `results/dusart_grid_1e10_2e13.csv`.

**Part 3, x >= 2 x 10^13 (analytic).** The same bounds apply to every x/k once x >= 2 x 10^13, and with l = ln x and 0 < ln(x/k) - 1 <= l - 1 they give

    2 L(x) - pi(x) >= F(x) - 2B,   F(x) = x (A(2l - 2.2) - l + 1) / ((l - 1)(l - 1.1)).

F is increasing for l > 7.4, and F(2 x 10^13) > 8.5 x 10^9 > 2B = 6,697,012,190. `src/tail_check.py` re-derives the generator list with a divisor-count sieve, computes A and B exactly, and checks the inequality with directed bounds on the logarithm (`results/tail_check.txt`).

The three ranges overlap and cover every n >= 1. The proof uses no prime-counting software.

**Cross-check with exact prime counts.** Evaluating L with exact prime counts ([primecount](https://github.com/kimwalisch/primecount) 8.8) shows the real margin: 57 points from 10^9 to 11,032,993,982,389,318 with pi(x_{i+1}) <= 2 L(x_i) (each x_{i+1} is the largest integer allowed, so the slack column is 0 by construction), and 2L(x)/pi(x) between 1.26 and 1.38 throughout (`results/grid_1e9_1e16.csv`, columns i, x_i, L(x_i), x_{i+1}, pi(x_{i+1}), 2L(x_i) - pi(x_{i+1}), x_{i+1}/x_i).

## Checks

- T(10^k) for k = 0..11 agrees with [OEIS A111291](https://oeis.org/A111291); T(10^12) = 26,568,322,474 is new.
- The generator count 413,980 and A > 0.508 reproduce Zelinsky's Theorem 8 numbers exactly; his bound B < 8,694,520,815 holds (the true sum of pi(g(k)) + 1 is 3,348,920,075).
- The lower bound L is tight: L(10^11) = 2,785,736,202 against T(10^11) = 2,823,898,245.
- An independent referee re-implemented the computations without seeing these programs (a factorization-free divisor-pair sieve, two enumerations of the generators, its own prime sieve and a library-free Lucy_Hedgehog prime count, Dusart's bounds checked numerically to 2 x 10^10): verdict sound, with minor edits that have been made. Report and programs: `review/REVIEW.md`.

## Reproduce

Needs a C compiler, Python 3 with NumPy (for `dusart_grid.py` and `tail_check.py`) and, for the generator list and the cross-check only, [primecount](https://github.com/kimwalisch/primecount).

    cc -O3 -o gens src/gens.c -lprimecount -lm && ./gens 28653696 gens_28653696.bin
    cc -O3 -o grid src/grid.c -lprimecount -lm && ./grid gens_28653696.bin 1e9 1e16 > grid.csv
    cc -O3 -o pointwise src/pointwise.c -lpthread -lm && ./pointwise 1000000000000 8
    python3 src/dusart_grid.py gens_28653696.bin 1e10 2e13 > dusart_grid.csv
    python3 src/tail_check.py gens_28653696.bin

`gens_28653696.bin` (413,980 records of two little-endian uint64: k, g(k)) has the SHA-256 in `results/gens_28653696.sha256`. Times on a 10-core laptop: generators 1 s, Dusart grid 1 min, primecount grid 11 min, pointwise to 10^10 34 s and to 10^12 about an hour.
