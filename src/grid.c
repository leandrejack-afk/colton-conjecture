// grid.c: cover [x0, X] by grid points x_0 < x_1 < ... with 2 L(x_i) >= pi(x_{i+1}),
// where L is the generator lower bound for T (see lbound.c).  Then for every n
// in [x_i, x_{i+1}]:  pi(n) <= pi(x_{i+1}) <= 2 L(x_i) <= 2 T(x_i) <= 2 T(n).
//
// x_{i+1} is chosen as the largest integer y in (x_i, 2 x_i] with
// pi(y) <= 2 L(x_i), found by binary search (pi is nondecreasing).
//
// Usage: grid gens.bin x0 X > grid.csv
// Columns: i, x_i, L(x_i), x_{i+1}, pi(x_{i+1}), 2L(x_i) - pi(x_{i+1}), ratio x_{i+1}/x_i

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <primecount.h>

// Parse a nonnegative integer given either as digits or in the form "<int>e<exp>" (e.g. 1e9, 2e13).
static uint64_t parse_u64(const char *s) {
  char *end;
  uint64_t m = strtoull(s, &end, 10);
  if (*end == 'e' || *end == 'E') {
    long e = strtol(end + 1, &end, 10);
    while (e-- > 0) m *= 10;
  }
  if (*end != 0) { fprintf(stderr, "bad number: %s\n", s); exit(1); }
  return m;
}

static uint64_t *rec;
static int64_t *pig;
static size_t n;

static int64_t Lx(uint64_t x) {
  int64_t L = 0;
  for (size_t i = 0; i < n; i++) {
    uint64_t k = rec[2 * i], g = rec[2 * i + 1];
    uint64_t y = x / k;
    if (y <= g) continue;
    int64_t t = primecount_pi((int64_t)y) - pig[i];
    if (t > 0) L += t;
  }
  return L;
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: grid gens.bin x0 X\n"); return 1; }
  FILE *f = fopen(argv[1], "rb");
  if (!f) { perror("fopen"); return 1; }
  fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
  n = sz / 16;
  rec = malloc(sz);
  if (fread(rec, 16, n, f) != n) { fprintf(stderr, "short read\n"); return 1; }
  fclose(f);
  pig = malloc(sizeof(int64_t) * n);
  for (size_t i = 0; i < n; i++) pig[i] = primecount_pi((int64_t)rec[2 * i + 1]);

  uint64_t x = parse_u64(argv[2]);
  uint64_t X = parse_u64(argv[3]);
  printf("i,x_i,L_x_i,x_next,pi_x_next,slack,ratio\n");
  for (int i = 0; ; i++) {
    int64_t L = Lx(x);
    int64_t target = 2 * L;
    // largest y in [x, 2x] with pi(y) <= target
    uint64_t lo = x, hi = 2 * x;
    if (primecount_pi((int64_t)lo) > target) {
      fprintf(stderr, "FAIL at x=%llu: pi(x) > 2L(x)\n", (unsigned long long)x);
      return 2;
    }
    if (primecount_pi((int64_t)hi) <= target) lo = hi;
    else {
      while (hi - lo > 1) {
        uint64_t mid = lo + (hi - lo) / 2;
        if (primecount_pi((int64_t)mid) <= target) lo = mid; else hi = mid;
      }
    }
    uint64_t y = lo;
    int64_t py = primecount_pi((int64_t)y);
    printf("%d,%llu,%lld,%llu,%lld,%lld,%.6f\n", i, (unsigned long long)x, (long long)L,
           (unsigned long long)y, (long long)py, (long long)(target - py), (double)y / (double)x);
    fflush(stdout);
    if (y <= x) { fprintf(stderr, "no progress at x=%llu\n", (unsigned long long)x); return 3; }
    if (y >= X) break;
    x = y;
  }
  return 0;
}
