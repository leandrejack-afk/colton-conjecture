/* Referee Part 1 check (independent of src/): exact T(n), pi(n) for all n <= N.
 *
 * Method: divisor-PAIR counting, no factorization at all.
 *   tau(n) = 2 * #{d | n : d*d < n} + [n is a square].
 * For each block [lo, hi) and each d <= sqrt(hi-1), add 2 at every multiple
 * m of d with m > d*d, and 1 at m = d*d.  n is a tau number iff n % tau(n) == 0,
 * prime iff tau(n) == 2.
 *
 * Reports: T(q), pi(q), D(q) = sum tau for query points q (10^k and any extra
 * points given on the command line); the global minimum of 2T(n) - pi(n), the
 * number of n attaining it and the first such n; the minimum over n >= 8 and
 * where; every n <= 2^20 with 2T(n) - pi(n) <= 3.
 * Usage: part1_divpairs N threads [extra query points...]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>
#include <math.h>

#define S (1u << 20)

typedef struct {
    uint64_t ntau, nprime, sumtau;
    int64_t minpref; uint64_t mincount, minpos;   /* over all n in block */
    int64_t minpref8; uint64_t minpos8;           /* over n >= 8 in block */
} blockres;

static uint64_t N; static uint64_t nb;
static blockres *R;
static atomic_ullong next_block;
static uint64_t *Q; static int nq;                 /* sorted query points */
static uint64_t *Qt, *Qp, *Qd;                     /* in-block partial counts */
static int64_t *Qblock;

static uint64_t isqrt64(uint64_t x) {
    uint64_t r = (uint64_t)sqrtl((long double)x);
    while (r * r > x) r--;
    while ((r + 1) * (r + 1) <= x) r++;
    return r;
}

static void *worker(void *arg) {
    (void)arg;
    uint16_t *cnt = malloc(sizeof(uint16_t) * S);
    for (;;) {
        uint64_t b = atomic_fetch_add(&next_block, 1);
        if (b >= nb) break;
        uint64_t lo = b * (uint64_t)S + 1, hi = lo + S; if (hi > N + 1) hi = N + 1;
        uint64_t len = hi - lo;
        memset(cnt, 0, sizeof(uint16_t) * len);
        uint64_t sq = isqrt64(hi - 1);
        for (uint64_t d = 1; d <= sq; d++) {
            uint64_t d2 = d * d;
            uint64_t m = ((lo + d - 1) / d) * d;
            if (m < d2) m = d2;
            if (m == d2) { cnt[m - lo] += 1; m += d; }
            uint16_t *p = cnt + (m - lo), *end = cnt + len;
            for (; p < end; p += d) *p += 2;
        }
        blockres r; memset(&r, 0, sizeof r);
        int64_t s = 0; r.minpref = INT64_MAX; r.minpref8 = INT64_MAX;
        for (uint64_t i = 0; i < len; i++) {
            uint64_t n = lo + i; uint32_t t = cnt[i];
            int tn = (n % t) == 0, pr = (t == 2);
            r.ntau += tn; r.nprime += pr; r.sumtau += t;
            s += 2 * tn - pr;
            if (s < r.minpref) { r.minpref = s; r.mincount = 1; r.minpos = n; }
            else if (s == r.minpref) r.mincount++;
            if (n >= 8 && s < r.minpref8) { r.minpref8 = s; r.minpos8 = n; }
            if (b == 0 && s <= 3) fprintf(stderr, "low: n=%llu 2T-pi=%lld\n", (unsigned long long)n, (long long)s);
        }
        /* query partials */
        for (int j = 0; j < nq; j++) if (Q[j] >= lo && Q[j] < hi) {
            uint64_t a = 0, c = 0, dd = 0;
            for (uint64_t i = 0; i <= Q[j] - lo; i++) { uint64_t n = lo + i; uint32_t t = cnt[i]; a += (n % t) == 0; c += (t == 2); dd += t; }
            Qt[j] = a; Qp[j] = c; Qd[j] = dd; Qblock[j] = (int64_t)b;
        }
        R[b] = r;
    }
    free(cnt);
    return NULL;
}

static int cmpu64(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return x < y ? -1 : x > y; }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage\n"); return 2; }
    N = strtoull(argv[1], 0, 10); int th = atoi(argv[2]);
    nb = (N + S - 1) / S;
    R = calloc(nb, sizeof(blockres));
    Q = malloc(sizeof(uint64_t) * (32 + argc)); nq = 0;
    for (uint64_t p = 1; p <= N; p *= 10) Q[nq++] = p;
    for (int i = 3; i < argc; i++) { uint64_t q = strtoull(argv[i], 0, 10); if (q <= N) Q[nq++] = q; }
    qsort(Q, nq, sizeof(uint64_t), cmpu64);
    Qt = calloc(nq, 8); Qp = calloc(nq, 8); Qd = calloc(nq, 8); Qblock = calloc(nq, 8);
    pthread_t tid[64];
    for (int i = 0; i < th; i++) pthread_create(&tid[i], 0, worker, 0);
    for (int i = 0; i < th; i++) pthread_join(tid[i], 0);
    /* combine */
    uint64_t T = 0, P = 0, D = 0;
    int64_t gmin = INT64_MAX, gmin8 = INT64_MAX; uint64_t gcount = 0, gpos = 0, gpos8 = 0;
    uint64_t *Tb = malloc(8 * nb), *Pb = malloc(8 * nb), *Db = malloc(8 * nb);
    for (uint64_t b = 0; b < nb; b++) {
        Tb[b] = T; Pb[b] = P; Db[b] = D;
        int64_t off = 2 * (int64_t)T - (int64_t)P;
        int64_t v = off + R[b].minpref;
        if (v < gmin) { gmin = v; gcount = R[b].mincount; gpos = R[b].minpos; }
        else if (v == gmin) gcount += R[b].mincount;
        int64_t v8 = off + R[b].minpref8;
        if (v8 < gmin8) { gmin8 = v8; gpos8 = R[b].minpos8; }
        T += R[b].ntau; P += R[b].nprime; D += R[b].sumtau;
    }
    printf("N=%llu T(N)=%llu pi(N)=%llu D(N)=%llu\n", (unsigned long long)N, (unsigned long long)T, (unsigned long long)P, (unsigned long long)D);
    printf("min_{1<=n<=N} 2T(n)-pi(n) = %lld, attained at %llu n, first at n=%llu\n", (long long)gmin, (unsigned long long)gcount, (unsigned long long)gpos);
    printf("min_{8<=n<=N} 2T(n)-pi(n) = %lld, first at n=%llu\n", (long long)gmin8, (unsigned long long)gpos8);
    for (int j = 0; j < nq; j++) {
        int64_t b = Qblock[j];
        printf("q=%llu T=%llu pi=%llu D=%llu 2T-pi=%lld\n", (unsigned long long)Q[j],
               (unsigned long long)(Tb[b] + Qt[j]), (unsigned long long)(Pb[b] + Qp[j]), (unsigned long long)(Db[b] + Qd[j]),
               (long long)(2 * (int64_t)(Tb[b] + Qt[j]) - (int64_t)(Pb[b] + Qp[j])));
    }
    return 0;
}
