/* Referee Part 2 check (independent of src/).
 *
 * L(x) = sum_{k in G} max(0, pi(floor(x/k)) - pi(g(k))), G = generators <= K.
 * Generators and g(k) from a linear smallest-prime-factor sieve (as in gens_check.c).
 * pi(y) for y <= YS: the referee's own segmented sieve of Eratosthenes (odd-only
 *   bytes, pre-sieved pattern for 3..13), streaming over [0, YS] in parallel chunks
 *   and answering all queries after sorting them by y.  No primesieve/primecount.
 * pi(y) for y > YS: primecount library (default algorithm), each argument also
 *   written to bigargs.txt for an external cross-check with other algorithms.
 * Input: grid CSV (i,x_i,L,x_next,pi_next,...).  Output, per row: recomputed L(x_i),
 *   recomputed pi(x_next), CSV values, and the verdict pi(x_next) <= 2 L(x_i).
 * Extra L points may follow on the command line.
 * Usage: lcheck grid.csv YS threads [extra x ...]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>
#include <math.h>
#include <primecount.h>

#define K 28653696u

typedef struct { uint64_t y; uint32_t row; uint32_t pig; } query;   /* pig==UINT32_MAX: pi-point */

static uint32_t *gk, *gpig; static size_t ng;
static uint32_t *bprimes; static size_t nbp;           /* base primes up to sqrt(YS) */
static uint64_t YS;
static query *Qs; static size_t nqs; static uint64_t *Qans;  /* answers (local, then global) */
#define CH (1ull << 32)
#define SEGB (1u << 21)          /* bytes per segment = odd numbers per segment */
#define PER 15015u
static uint8_t *pattern;
static uint64_t nch; static uint64_t *chcount; static size_t *chq0, *chq1;
static atomic_ullong nextch;

static uint64_t isqrt64(uint64_t x) { uint64_t r = (uint64_t)sqrtl((long double)x); while (r * r > x) r--; while ((r + 1) * (r + 1) <= x) r++; return r; }

static void *chunk_worker(void *arg) {
    (void)arg;
    uint8_t *seg = malloc(SEGB);
    uint64_t *nextm = malloc(sizeof(uint64_t) * nbp);
    for (;;) {
        uint64_t c = atomic_fetch_add(&nextch, 1);
        if (c >= nch) break;
        uint64_t A = c * CH, Bend = A + CH; if (Bend > YS + 1) Bend = YS + 1;   /* numbers in [A, Bend) */
        /* init next-multiple index for each base prime >= 17 (as odd multiples) */
        for (size_t j = 0; j < nbp; j++) {
            uint64_t p = bprimes[j];
            uint64_t m = p * p;
            if (m < A + 1) { m = ((A + 1 + p - 1) / p) * p; if (!(m & 1)) m += p; }
            nextm[j] = m;     /* odd multiple of p, >= max(p^2, A+1) */
        }
        uint64_t cnt = 0; size_t qi = chq0[c];
        for (uint64_t lo = A; lo < Bend; lo += 2ull * SEGB) {
            uint64_t hi = lo + 2ull * SEGB; if (hi > Bend) hi = Bend;     /* [lo, hi) */
            uint64_t nodd = (hi - lo) / 2;     /* lo even, hi even except possibly last */
            if ((hi - lo) & 1) nodd = (hi - lo + 1) / 2;
            /* index i <-> n = lo + 2i + 1 */
            memcpy(seg, pattern + ((lo / 2) % PER), nodd);
            if (lo == 0) { seg[0] = 0; seg[1] = 1; seg[2] = 1; seg[3] = 1; seg[5] = 1; seg[6] = 1; } /* 1 not prime; 3,5,7,11,13 prime */
            uint64_t lim = hi - 1;
            for (size_t j = 0; j < nbp; j++) {
                uint64_t p = bprimes[j];
                if (p * p > lim) break;
                uint64_t m = nextm[j];
                if (m >= hi) continue;
                uint64_t i = (m - lo - 1) / 2;
                for (; i < nodd; i += p) seg[i] = 0;
                nextm[j] = lo + 2 * i + 1;
            }
            /* answer queries inside [lo, hi) */
            uint64_t run = 0, pos = 0;
            while (qi < chq1[c] && Qs[qi].y < hi) {
                uint64_t y = Qs[qi].y;
                /* count odd primes in [lo, y]: indices 0..(y-lo-1)/2 when y > lo */
                uint64_t upto = (y > lo) ? (y - lo - 1) / 2 + 1 : 0;   /* exclusive bound */
                if (upto > nodd) upto = nodd;
                for (; pos < upto; pos++) run += seg[pos];
                Qans[qi] = cnt + run;
                qi++;
            }
            uint64_t s = 0;
            for (uint64_t i = 0; i < nodd; i++) s += seg[i];
            cnt += s;
        }
        chcount[c] = cnt;
    }
    free(seg); free(nextm);
    return NULL;
}

static int cmpq(const void *a, const void *b) { const query *x = a, *y = b; return x->y < y->y ? -1 : x->y > y->y; }

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage\n"); return 2; }
    YS = strtoull(argv[2], 0, 10); int th = atoi(argv[3]);
    primecount_set_num_threads(th);
    /* generators */
    uint32_t *spf = calloc((size_t)K + 1, 4); uint32_t *pr = malloc(4 * ((size_t)K / 5 + 1000)); size_t np = 0;
    for (uint32_t i = 2; i <= K; i++) {
        if (!spf[i]) { spf[i] = i; pr[np++] = i; }
        for (size_t j = 0; j < np && pr[j] <= spf[i] && (uint64_t)i * pr[j] <= K; j++) spf[(size_t)i * pr[j]] = pr[j];
    }
    gk = malloc(4 * 500000); gpig = malloc(4 * 500000);
    for (uint32_t k = 2; k <= K; k++) {
        uint32_t m = k; uint64_t tau = 1; uint32_t gp = 1;
        while (m > 1) { uint32_t p = spf[m]; uint32_t e = 0; while (m % p == 0) { m /= p; e++; } tau *= e + 1; if (p > gp) gp = p; }
        if (k % (2 * tau) == 0) {
            size_t a = 0, b = np; while (a < b) { size_t mid = (a + b) / 2; if (pr[mid] <= gp) a = mid + 1; else b = mid; }
            gk[ng] = k; gpig[ng] = (uint32_t)a; ng++;
        }
    }
    free(spf);
    fprintf(stderr, "generators: %zu\n", ng);
    /* grid */
    FILE *f = fopen(argv[1], "r"); char line[512];
    uint64_t X[200], Lcsv[200], Xn[200], Pcsv[200]; int nr = 0;
    if (!fgets(line, sizeof line, f)) return 1;
    while (fgets(line, sizeof line, f)) {
        unsigned long long i, x, L, xn, pn;
        if (sscanf(line, "%llu,%llu,%llu,%llu,%llu", &i, &x, &L, &xn, &pn) == 5) { X[nr] = x; Lcsv[nr] = L; Xn[nr] = xn; Pcsv[nr] = pn; nr++; }
    }
    fclose(f);
    int nextra = argc - 4; uint64_t Xe[64];
    for (int i = 0; i < nextra; i++) Xe[i] = strtoull(argv[4 + i], 0, 10);
    int nL = nr + nextra;
    uint64_t *XL = malloc(8 * nL);
    for (int r = 0; r < nr; r++) XL[r] = X[r];
    for (int i = 0; i < nextra; i++) XL[nr + i] = Xe[i];
    /* base primes up to sqrt(YS), from pr[] if possible */
    uint64_t sq = isqrt64(YS) + 1;
    bprimes = malloc(4 * np); nbp = 0;
    for (size_t j = 0; j < np && pr[j] <= sq; j++) if (pr[j] >= 17) bprimes[nbp++] = pr[j];
    if (sq > K) { fprintf(stderr, "YS too large for base primes\n"); return 1; }
    /* queries */
    size_t cap = (size_t)nL * ng + 1000; Qs = malloc(sizeof(query) * cap); nqs = 0;
    FILE *fb = fopen("bigargs.txt", "w");
    long double *Lbig = calloc(nL, sizeof(long double)); (void)Lbig;
    int64_t *Lsum = calloc(nL, 8);
    uint64_t nbig = 0;
    double t_pc = 0;
    for (int r = 0; r < nL; r++) {
        for (size_t j = 0; j < ng; j++) {
            uint64_t y = XL[r] / gk[j];
            if (y <= YS) { Qs[nqs].y = y; Qs[nqs].row = r; Qs[nqs].pig = gpig[j]; nqs++; }
            else {
                int64_t p = primecount_pi((int64_t)y);
                fprintf(fb, "%llu %lld\n", (unsigned long long)y, (long long)p);
                int64_t t = p - (int64_t)gpig[j]; if (t > 0) Lsum[r] += t;
                nbig++;
            }
        }
    }
    uint64_t *pin = calloc(nr, 8); int *pin_isbig = calloc(nr, sizeof(int));
    for (int r = 0; r < nr; r++) {
        if (Xn[r] <= YS) { Qs[nqs].y = Xn[r]; Qs[nqs].row = r; Qs[nqs].pig = UINT32_MAX; nqs++; }
        else { pin[r] = (uint64_t)primecount_pi((int64_t)Xn[r]); pin_isbig[r] = 1; fprintf(fb, "%llu %llu\n", (unsigned long long)Xn[r], (unsigned long long)pin[r]); }
    }
    fclose(fb);
    (void)t_pc;
    fprintf(stderr, "queries: %zu sieve, %llu primecount\n", nqs, (unsigned long long)nbig);
    qsort(Qs, nqs, sizeof(query), cmpq);
    Qans = calloc(nqs, 8);
    /* chunks */
    nch = (YS + 1 + CH - 1) / CH;
    chcount = calloc(nch, 8); chq0 = calloc(nch, sizeof(size_t)); chq1 = calloc(nch, sizeof(size_t));
    { size_t q = 0; for (uint64_t c = 0; c < nch; c++) { chq0[c] = q; while (q < nqs && Qs[q].y < (c + 1) * CH) q++; chq1[c] = q; } }
    /* pattern: index i <-> n = 2i+1, zero if divisible by 3,5,7,11,13 */
    pattern = malloc(PER + SEGB + 16);
    for (size_t i = 0; i < PER + SEGB + 16; i++) { uint64_t n = 2 * (uint64_t)i + 1; pattern[i] = (n % 3 && n % 5 && n % 7 && n % 11 && n % 13) ? 1 : 0; }
    pthread_t tid[64];
    for (int i = 0; i < th; i++) pthread_create(&tid[i], 0, chunk_worker, 0);
    for (int i = 0; i < th; i++) pthread_join(tid[i], 0);
    /* global offsets */
    uint64_t off = 0; uint64_t *choff = malloc(8 * nch);
    for (uint64_t c = 0; c < nch; c++) { choff[c] = off; off += chcount[c]; }
    fprintf(stderr, "own sieve: pi(%llu) = %llu\n", (unsigned long long)YS, (unsigned long long)(off + (YS >= 2)));
    for (uint64_t c = 0; c < nch; c++) for (size_t q = chq0[c]; q < chq1[c]; q++) {
        uint64_t piy = choff[c] + Qans[q] + (Qs[q].y >= 2);
        if (Qs[q].pig == UINT32_MAX) pin[Qs[q].row] = piy;
        else { int64_t t = (int64_t)piy - (int64_t)Qs[q].pig; if (t > 0) Lsum[Qs[q].row] += t; }
        Qans[q] = piy;
    }
    /* sample of sieve answers for cross-check */
    FILE *fs = fopen("sieve_sample.txt", "w");
    for (size_t q = 0; q < nqs; q += 9973) fprintf(fs, "%llu %llu\n", (unsigned long long)Qs[q].y, (unsigned long long)Qans[q]);
    if (nqs) fprintf(fs, "%llu %llu\n", (unsigned long long)Qs[nqs - 1].y, (unsigned long long)Qans[nqs - 1]);
    fclose(fs);
    printf("row,x_i,L_ref,L_csv,x_next,pi_ref,pi_csv,pi_src,slack_ref,ok\n");
    int bad = 0;
    for (int r = 0; r < nr; r++) {
        int64_t slack = 2 * Lsum[r] - (int64_t)pin[r];
        int ok = (slack >= 0) && (r + 1 >= nr || Xn[r] == X[r + 1]);
        if (!ok || (uint64_t)Lsum[r] != Lcsv[r] || pin[r] != Pcsv[r]) bad++;
        printf("%d,%llu,%lld,%llu,%llu,%llu,%llu,%s,%lld,%s\n", r, (unsigned long long)X[r], (long long)Lsum[r], (unsigned long long)Lcsv[r],
               (unsigned long long)Xn[r], (unsigned long long)pin[r], (unsigned long long)Pcsv[r], pin_isbig[r] ? "primecount" : "own-sieve",
               (long long)slack, ok ? "ok" : "FAIL");
    }
    for (int i = 0; i < nextra; i++) printf("extra x=%llu L_ref=%lld\n", (unsigned long long)Xe[i], (long long)Lsum[nr + i]);
    printf("rows %d, rows with any discrepancy or failure: %d\n", nr, bad);
    return 0;
}
