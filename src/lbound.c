// lbound.c: lower bound for T(x), the number of tau numbers <= x.
//
// For every generator k (2*tau(k) | k) and every prime p > g(k) (largest prime
// factor of k), n = k*p is a tau number: tau(kp) = 2 tau(k) divides k*p.
// Distinct pairs (k, p) give distinct n, because p = g(n) appears to the first
// power in n and k = n / p.  Hence
//     T(x) >= L(x) = sum over generators k of max(0, pi(x/k) - pi(g(k))).
// The sum may be truncated to any subset of generators (all terms are >= 0).
//
// Usage: lbound gens.bin x1 [x2 ...]
// Prints, for each x: L(x), pi(x), 2L(x)/pi(x).

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

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: lbound gens.bin x...\n"); return 1; }
  FILE *f = fopen(argv[1], "rb");
  if (!f) { perror("fopen"); return 1; }
  fseek(f, 0, SEEK_END);
  long sz = ftell(f);
  fseek(f, 0, SEEK_SET);
  size_t n = sz / 16;
  uint64_t *rec = malloc(sz);
  if (fread(rec, 16, n, f) != n) { fprintf(stderr, "short read\n"); return 1; }
  fclose(f);
  // pi(g(k)) for each generator, computed once
  int64_t *pig = malloc(sizeof(int64_t) * n);
  for (size_t i = 0; i < n; i++) pig[i] = primecount_pi((int64_t)rec[2 * i + 1]);

  for (int a = 2; a < argc; a++) {
    uint64_t x = parse_u64(argv[a]);
    __int128 L = 0;
    size_t used = 0;
    for (size_t i = 0; i < n; i++) {
      uint64_t k = rec[2 * i], g = rec[2 * i + 1];
      uint64_t y = x / k;
      if (y <= g) continue;  // no prime p > g with k*p <= x
      int64_t t = primecount_pi((int64_t)y) - pig[i];
      if (t > 0) { L += t; used++; }
    }
    int64_t px = primecount_pi((int64_t)x);
    printf("x=%llu L=%lld pi=%lld 2L/pi=%.6f generators_used=%zu\n",
           (unsigned long long)x, (long long)L, (long long)px,
           2.0 * (double)L / (double)px, used);
    fflush(stdout);
  }
  return 0;
}
