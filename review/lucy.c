/* Referee: primecount-free L(x) and pi(x) via the Lucy_Hedgehog / Legendre-type
 * dynamic programme, which yields pi(floor(x/i)) for EVERY i at once.
 * (Independent of src/ and of primecount/primesieve.)
 *   S(v) initialised to v-1; for each prime p <= sqrt(x), for v >= p^2:
 *   S(v) -= S(v/p) - S(p-1).
 * Prints, for each x given: pi(x), L(x) = sum_{k in G} max(0, pi(x/k) - pi(g(k))),
 * with G = generators k <= 28653696 (2 tau(k) | k) found by an spf sieve.
 * Usage: lucy x1 [x2 ...]   (single-threaded, memory 12*sqrt(x) bytes; small table int32, valid for sqrt(x) < 2^31)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define K 28653696u

static uint64_t isqrt64(uint64_t x) { uint64_t r = (uint64_t)sqrtl((long double)x); while (r * r > x) r--; while ((r + 1) * (r + 1) <= x) r++; return r; }

int main(int argc, char **argv) {
    /* generators */
    uint32_t *spf = calloc((size_t)K + 1, 4); uint32_t *pr = malloc(4 * ((size_t)K / 5 + 1000)); size_t np = 0;
    for (uint32_t i = 2; i <= K; i++) {
        if (!spf[i]) { spf[i] = i; pr[np++] = i; }
        for (size_t j = 0; j < np && pr[j] <= spf[i] && (uint64_t)i * pr[j] <= K; j++) spf[(size_t)i * pr[j]] = pr[j];
    }
    uint32_t *gk = malloc(4 * 500000), *gpig = malloc(4 * 500000); size_t ng = 0;
    for (uint32_t k = 2; k <= K; k++) {
        uint32_t m = k; uint64_t tau = 1; uint32_t gp = 1;
        while (m > 1) { uint32_t p = spf[m]; uint32_t e = 0; while (m % p == 0) { m /= p; e++; } tau *= e + 1; if (p > gp) gp = p; }
        if (k % (2 * tau) == 0) {
            size_t a = 0, b = np; while (a < b) { size_t mid = (a + b) / 2; if (pr[mid] <= gp) a = mid + 1; else b = mid; }
            gk[ng] = k; gpig[ng] = (uint32_t)a; ng++;
        }
    }
    free(spf); free(pr);
    for (int ai = 1; ai < argc; ai++) {
        clock_t t0 = clock();
        uint64_t x = strtoull(argv[ai], 0, 10);
        uint64_t r = isqrt64(x);
        int32_t *sm = malloc(sizeof(int32_t) * (r + 2)); int64_t *lg = malloc(sizeof(int64_t) * (r + 2));
        if (!sm || !lg) { fprintf(stderr, "alloc\n"); return 1; }
        for (uint64_t v = 0; v <= r; v++) sm[v] = (int32_t)v - 1;
        sm[0] = 0;
        for (uint64_t i = 1; i <= r; i++) lg[i] = (int64_t)(x / i) - 1;
        for (uint64_t p = 2; p <= r; p++) {
            if (sm[p] == sm[p - 1]) continue;          /* not prime */
            int64_t sp = sm[p - 1]; uint64_t p2 = p * p;
            uint64_t imax = x / p2; if (imax > r) imax = r;
            uint64_t ilim = r / p; if (ilim > imax) ilim = imax;
            for (uint64_t i = 1; i <= ilim; i++) lg[i] -= lg[i * p] - sp;
            for (uint64_t i = ilim + 1; i <= imax; i++) lg[i] -= sm[x / (i * p)] - sp;
            if (p2 <= r) for (uint64_t v = r; v >= p2; v--) sm[v] -= (int32_t)(sm[v / p] - sp);
        }
        /* pi(floor(x/k)) = lg[k] if k <= r else sm[x/k] */
        int64_t L = 0;
        for (size_t j = 0; j < ng; j++) {
            uint64_t k = gk[j];
            int64_t pk = (k <= r) ? lg[k] : sm[x / k];
            int64_t t = pk - (int64_t)gpig[j]; if (t > 0) L += t;
        }
        printf("x=%llu pi=%lld L=%lld secs=%.1f\n", (unsigned long long)x, (long long)lg[1], (long long)L, (double)(clock() - t0) / CLOCKS_PER_SEC);
        fflush(stdout);
        free(sm); free(lg);
    }
    return 0;
}
