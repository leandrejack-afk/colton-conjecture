// pointwise.c: exact check of Colton's inequality pi(n) <= 2 T(n) for every
// n <= N, where T(n) counts the tau numbers <= n (tau(m) | m).
//
// Every n is factored by a segmented sieve (trial division by the primes up to
// sqrt(N) inside each segment; the cofactor left is 1 or one prime), giving
// tau(n).  n is prime iff tau(n) == 2; n is a tau number iff n % tau(n) == 0.
//
// Usage: pointwise N threads [checkpoint ...]
// Output: T(c) and pi(c) for every power of 10 up to N and every extra
// checkpoint c <= N given on the command line; T(N), pi(N); the minimum of
// 2T(n) - pi(n) over all n in [1, N] and where it is attained.
// 2T(n) - pi(n) only decreases at primes, so it is evaluated at every prime
// (n = 1, a tau number and not prime, gives 2T - pi = 2).

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#define SEG (1u << 18)
#define MAXCP 4096

static uint64_t N;
static uint32_t *primes;
static size_t nprimes;
static uint64_t cps[MAXCP];   // sorted checkpoints
static int ncp = 0;

typedef struct {
  uint64_t tcount, pcount;   // tau numbers and primes in the segment
  int64_t minval;            // min over primes n in the segment of 2 t(n) - p(n), local counts up to n
  uint64_t minarg;           // first n attaining minval in the segment
  uint64_t nmin;             // number of n attaining minval in the segment
} SegInfo;

static SegInfo *info;
static uint64_t *cp_t, *cp_p;  // local counts in the checkpoint's segment, up to and including it
static uint64_t nseg;
static uint64_t next_seg = 0;
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

static void process(uint64_t s, uint64_t *rem, uint32_t *tau) {
  uint64_t lo = s * (uint64_t)SEG + 1;      // segment covers [lo, hi)
  uint64_t hi = lo + SEG;
  if (hi > N + 1) hi = N + 1;
  uint64_t len = hi - lo;
  for (uint64_t i = 0; i < len; i++) { rem[i] = lo + i; tau[i] = 1; }
  for (size_t t = 0; t < nprimes; t++) {
    uint64_t p = primes[t];
    if (p * p > hi - 1) break;
    uint64_t start = ((lo + p - 1) / p) * p;
    for (uint64_t m = start; m < hi; m += p) {
      uint64_t i = m - lo;
      uint64_t r = rem[i] / p;
      uint32_t e = 1;
      while (r % p == 0) { r /= p; e++; }
      rem[i] = r;
      tau[i] *= (e + 1);
    }
  }
  // checkpoints inside this segment
  int c0 = 0;
  while (c0 < ncp && cps[c0] < lo) c0++;
  SegInfo inf = {0, 0, INT64_MAX, 0, 0};
  for (uint64_t i = 0; i < len; i++) {
    uint64_t n = lo + i;
    if (rem[i] > 1) tau[i] *= 2;
    if (n % tau[i] == 0) inf.tcount++;
    if (tau[i] == 2) {
      inf.pcount++;
      int64_t v = 2 * (int64_t)inf.tcount - (int64_t)inf.pcount;
      if (v < inf.minval) { inf.minval = v; inf.minarg = n; inf.nmin = 1; }
      else if (v == inf.minval) inf.nmin++;
    }
    while (c0 < ncp && cps[c0] == n) { cp_t[c0] = inf.tcount; cp_p[c0] = inf.pcount; c0++; }
  }
  info[s] = inf;
}

static void *worker(void *arg) {
  (void)arg;
  uint64_t *rem = malloc(sizeof(uint64_t) * SEG);
  uint32_t *tau = malloc(sizeof(uint32_t) * SEG);
  for (;;) {
    pthread_mutex_lock(&lock);
    uint64_t s = next_seg++;
    pthread_mutex_unlock(&lock);
    if (s >= nseg) break;
    process(s, rem, tau);
  }
  free(rem); free(tau);
  return 0;
}

static int cmpu64(const void *a, const void *b) {
  uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
  return x < y ? -1 : x > y;
}

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: pointwise N threads [checkpoint ...]\n"); return 1; }
  N = strtoull(argv[1], 0, 10);
  int nth = atoi(argv[2]);
  for (uint64_t d = 1; d <= N && ncp < MAXCP; d *= 10) { cps[ncp++] = d; if (d > UINT64_MAX / 10) break; }
  for (int a = 3; a < argc && ncp < MAXCP; a++) {
    uint64_t c = strtoull(argv[a], 0, 10);
    if (c >= 1 && c <= N) cps[ncp++] = c;
  }
  qsort(cps, ncp, sizeof(uint64_t), cmpu64);
  int u = 0;
  for (int i = 0; i < ncp; i++) if (i == 0 || cps[i] != cps[i - 1]) cps[u++] = cps[i];
  ncp = u;
  cp_t = calloc(ncp, sizeof(uint64_t)); cp_p = calloc(ncp, sizeof(uint64_t));

  uint64_t r = (uint64_t)sqrtl((long double)N);
  while (r * r > N) r--;
  while ((r + 1) * (r + 1) <= N) r++;
  char *comp = calloc(r + 2, 1);
  primes = malloc(sizeof(uint32_t) * (r + 2));
  for (uint64_t i = 2; i <= r; i++) if (!comp[i]) { primes[nprimes++] = (uint32_t)i; for (uint64_t j = i * i; j <= r; j += i) comp[j] = 1; }
  nseg = (N + SEG - 1) / SEG;
  info = calloc(nseg, sizeof(SegInfo));
  pthread_t th[256];
  if (nth < 1) nth = 1;
  if (nth > 256) nth = 256;
  for (int i = 0; i < nth; i++) pthread_create(&th[i], 0, worker, 0);
  for (int i = 0; i < nth; i++) pthread_join(th[i], 0);

  // combine: prefix counts at segment starts
  uint64_t T = 0, P = 0;
  int64_t gmin = INT64_MAX; uint64_t gargs[64]; int ng = 0; uint64_t gcount = 0;
  int c = 0;
  for (uint64_t s = 0; s < nseg; s++) {
    SegInfo *in = &info[s];
    uint64_t lo = s * (uint64_t)SEG + 1, hi = lo + SEG; if (hi > N + 1) hi = N + 1;
    if (in->pcount > 0) {
      int64_t v = 2 * (int64_t)T - (int64_t)P + in->minval;
      if (v < gmin) { gmin = v; ng = 0; gcount = 0; }
      if (v == gmin) { if (ng < 64) gargs[ng++] = in->minarg; gcount += in->nmin; }
    }
    while (c < ncp && cps[c] < hi) {
      printf("checkpoint n=%llu T=%llu pi=%llu 2T-pi=%lld\n", (unsigned long long)cps[c],
             (unsigned long long)(T + cp_t[c]), (unsigned long long)(P + cp_p[c]),
             (long long)(2 * (int64_t)(T + cp_t[c]) - (int64_t)(P + cp_p[c])));
      c++;
    }
    T += in->tcount; P += in->pcount;
  }
  printf("N=%llu T(N)=%llu pi(N)=%llu\n", (unsigned long long)N, (unsigned long long)T, (unsigned long long)P);
  printf("min over 1<=n<=N of 2T(n)-pi(n) = %lld, attained at %llu prime(s) (first in each segment listed; between consecutive tau numbers the value cannot drop below its value at the preceding prime):",
         (long long)(gmin < 2 ? gmin : 2), (unsigned long long)gcount);
  for (int i = 0; i < ng; i++) printf(" %llu", (unsigned long long)gargs[i]);
  printf("\n");
  return 0;
}
