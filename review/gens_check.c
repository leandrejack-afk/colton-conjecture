/* Independent generator check (referee code, not derived from src/).
 *
 * Generators: k >= 1 with 2*tau(k) | k.  For each k <= K we need tau(k) and
 * g(k) = largest prime factor.  Two unrelated methods:
 *   M1: linear smallest-prime-factor sieve, then factor k by repeated spf.
 *   M2: depth-first enumeration of every k <= K by its factorization
 *       (primes in increasing order), carrying tau and the largest prime.
 * Both lists are compared with each other and with a binary file of
 * (uint64 k, uint64 g) little-endian records.  Also prints
 *   count, B = sum pi(g(k)), Zelinsky's sum (pi(g)+1),
 *   floor/ceil sums of 10^30/k (rigorous enclosure of A = sum 1/k).
 * Usage: gens_check K file.bin
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static uint32_t K;
static uint32_t *primes; static size_t np;
static uint8_t *isgen2;      /* M2 marks */
static uint32_t *g2;         /* M2 largest prime factor */

static void dfs(size_t pi_start, uint64_t k, uint64_t tau, uint32_t gpf) {
    /* k is the current number; record it */
    if (k <= K) {
        if (k % (2 * tau) == 0) { isgen2[k] = 1; g2[k] = gpf; }
    }
    for (size_t i = pi_start; i < np; i++) {
        uint64_t p = primes[i];
        if (k * p > K) break;
        uint64_t m = k * p; uint64_t e = 1;
        while (m <= K) {
            dfs(i + 1, m, tau * (e + 1), (uint32_t)p);
            if (m > K / p) break;
            m *= p; e++;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s K file.bin\n", argv[0]); return 2; }
    K = (uint32_t)strtoul(argv[1], 0, 10);
    uint32_t *spf = calloc((size_t)K + 1, sizeof(uint32_t));
    primes = malloc(sizeof(uint32_t) * ((size_t)K / 5 + 1000));
    if (!spf || !primes) { fprintf(stderr, "alloc\n"); return 1; }
    for (uint32_t i = 2; i <= K; i++) {
        if (!spf[i]) { spf[i] = i; primes[np++] = i; }
        for (size_t j = 0; j < np && primes[j] <= spf[i] && (uint64_t)i * primes[j] <= K; j++)
            spf[(size_t)i * primes[j]] = primes[j];
    }
    fprintf(stderr, "pi(K) = %zu\n", np);
    /* M1 */
    uint8_t *isgen1 = calloc((size_t)K + 1, 1);
    uint32_t *g1 = calloc((size_t)K + 1, sizeof(uint32_t));
    size_t c1 = 0;
    for (uint32_t k = 1; k <= K; k++) {
        uint32_t m = k; uint64_t tau = 1; uint32_t gp = 1;
        while (m > 1) {
            uint32_t p = spf[m]; uint32_t e = 0;
            while (m % p == 0) { m /= p; e++; }
            tau *= (e + 1); if (p > gp) gp = p;
        }
        if (k % (2 * tau) == 0) { isgen1[k] = 1; g1[k] = gp; c1++; }
    }
    /* M2 */
    isgen2 = calloc((size_t)K + 1, 1);
    g2 = calloc((size_t)K + 1, sizeof(uint32_t));
    dfs(0, 1, 1, 1);
    size_t c2 = 0, mism12 = 0;
    for (uint32_t k = 1; k <= K; k++) {
        c2 += isgen2[k];
        if (isgen1[k] != isgen2[k] || (isgen1[k] && g1[k] != g2[k])) {
            if (mism12 < 10) fprintf(stderr, "M1/M2 mismatch at %u\n", k);
            mism12++;
        }
    }
    printf("M1 count %zu, M2 count %zu, M1-vs-M2 mismatches %zu\n", c1, c2, mism12);
    /* compare with file */
    FILE *f = fopen(argv[2], "rb");
    if (!f) { perror("file"); return 1; }
    uint64_t rec[2]; size_t nf = 0, mismf = 0; uint64_t prevk = 0;
    uint8_t *seen = calloc((size_t)K + 1, 1);
    while (fread(rec, 8, 2, f) == 2) {
        nf++;
        if (rec[0] <= prevk) { mismf++; fprintf(stderr, "not increasing at record %zu\n", nf); }
        prevk = rec[0];
        if (rec[0] > K || !isgen1[rec[0]] || g1[rec[0]] != rec[1]) {
            if (mismf < 10) fprintf(stderr, "file record %zu (k=%llu,g=%llu) wrong\n", nf,
                                    (unsigned long long)rec[0], (unsigned long long)rec[1]);
            mismf++;
        } else seen[rec[0]] = 1;
    }
    fclose(f);
    size_t missing = 0;
    for (uint32_t k = 1; k <= K; k++) if (isgen1[k] && !seen[k]) { if (missing < 10) fprintf(stderr, "missing from file: %u\n", k); missing++; }
    printf("file records %zu, bad records %zu, generators missing from file %zu\n", nf, mismf, missing);
    /* B and A enclosure */
    uint64_t B = 0, Bz = 0; uint32_t maxg = 0;
    unsigned __int128 S = 1; for (int i = 0; i < 30; i++) S *= 10;
    unsigned __int128 lo = 0, hi = 0;
    for (uint32_t k = 1; k <= K; k++) if (isgen1[k]) {
        uint32_t g = g1[k]; if (g > maxg) maxg = g;
        /* pi(g) by binary search in primes[] */
        size_t a = 0, b = np; while (a < b) { size_t mid = (a + b) / 2; if (primes[mid] <= g) a = mid + 1; else b = mid; }
        B += a; Bz += a + 1;
        unsigned __int128 q = S / k, r = S % k;
        lo += q; hi += q + (r != 0);
    }
    char buf[64], buf2[64];
    /* print 128-bit */
    #define P128(x, out) do { unsigned __int128 _v = (x); char _t[64]; int _n = 0; if (!_v) _t[_n++] = '0'; while (_v) { _t[_n++] = '0' + (int)(_v % 10); _v /= 10; } for (int _i = 0; _i < _n; _i++) out[_i] = _t[_n - 1 - _i]; out[_n] = 0; } while (0)
    P128(lo, buf); P128(hi, buf2);
    printf("max g %u\nB = sum pi(g(k)) = %llu\nsum (pi(g)+1) = %llu\n", maxg, (unsigned long long)B, (unsigned long long)Bz);
    printf("10^30 * A in [%s, %s]\n", buf, buf2);
    /* counts of generators below some thresholds, useful for planning */
    size_t cnt = 0; uint32_t th[] = {1000, 10000, 100000, 1000000, 10000000};
    int ti = 0;
    for (uint32_t k = 1; k <= K && ti < 5; k++) { if (k > th[ti]) { printf("generators <= %u: %zu\n", th[ti], cnt); ti++; } cnt += isgen1[k]; }
    return (mism12 || mismf || missing || c1 != nf) ? 1 : 0;
}
