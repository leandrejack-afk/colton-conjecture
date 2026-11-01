// gens.c: list every k <= K with 2*tau(k) | k (Zelinsky's "p-generators":
// for such k and every prime p not dividing k, k*p is a tau number), with
// g(k) = largest prime factor of k.
//
// Usage: gens K out.bin
// out.bin: records of (uint64 k, uint64 g) in increasing k.
// Prints: count, sum of 1/k, sum of pi(g(k)).
//
// Method: segmented trial division by the primes up to sqrt(K); after the
// small primes are removed the cofactor is 1 or a single prime.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <primecount.h>

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: gens K out.bin\n"); return 1; }
  uint64_t K = strtoull(argv[1], 0, 10);
  FILE *out = fopen(argv[2], "wb");
  if (!out) { perror("fopen"); return 1; }

  uint64_t r = (uint64_t)sqrtl((long double)K);
  while (r * r > K) r--;
  while ((r + 1) * (r + 1) <= K) r++;
  // primes up to r
  char *comp = calloc(r + 2, 1);
  uint64_t *pr = malloc(sizeof(uint64_t) * (r + 2));
  size_t np = 0;
  for (uint64_t i = 2; i <= r; i++) {
    if (!comp[i]) {
      pr[np++] = i;
      for (uint64_t j = i * i; j <= r; j += i) comp[j] = 1;
    }
  }

  const uint64_t S = 1 << 20;
  uint64_t *rem = malloc(sizeof(uint64_t) * S);
  uint32_t *tau = malloc(sizeof(uint32_t) * S);
  uint64_t *gpf = malloc(sizeof(uint64_t) * S);
  uint64_t count = 0;
  long double A = 0, Acomp = 0;  // Kahan sum of 1/k
  long double B = 0;              // sum of pi(g(k))

  for (uint64_t lo = 1; lo <= K; lo += S) {
    uint64_t hi = lo + S;  // exclusive
    if (hi > K + 1) hi = K + 1;
    uint64_t len = hi - lo;
    for (uint64_t i = 0; i < len; i++) { rem[i] = lo + i; tau[i] = 1; gpf[i] = 1; }
    for (size_t t = 0; t < np; t++) {
      uint64_t p = pr[t];
      if (p * p > hi - 1) break;
      uint64_t start = ((lo + p - 1) / p) * p;
      for (uint64_t m = start; m < hi; m += p) {
        uint64_t i = m - lo;
        uint32_t e = 0;
        while (rem[i] % p == 0) { rem[i] /= p; e++; }
        tau[i] *= (e + 1);
        gpf[i] = p;
      }
    }
    for (uint64_t i = 0; i < len; i++) {
      uint64_t k = lo + i;
      if (rem[i] > 1) { tau[i] *= 2; gpf[i] = rem[i]; }
      if (k % (2ULL * tau[i]) == 0) {
        uint64_t rec[2] = { k, gpf[i] };
        fwrite(rec, sizeof(uint64_t), 2, out);
        count++;
        long double y = 1.0L / (long double)k - Acomp;
        long double tsum = A + y;
        Acomp = (tsum - A) - y;
        A = tsum;
        B += (long double)primecount_pi((int64_t)gpf[i]);
      }
    }
  }
  fclose(out);
  printf("K=%llu generators=%llu sum_1/k=%.12Lf sum_pi(g)=%.0Lf\n",
         (unsigned long long)K, (unsigned long long)count, A, B);
  return 0;
}
