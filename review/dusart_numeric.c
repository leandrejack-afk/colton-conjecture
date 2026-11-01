/* Referee: numerical check of the two Dusart (6.6) inequalities on a finite range,
 * as statements about REAL y (worst cases are at primes):
 *   lower: pi(y) >= y/(ln y - 1) for real y >= 5393
 *          <=> for every prime q >= 5393: pi(q) - 1 >= q/(ln q - 1)   (left limit at q)
 *              and at y = 5393 itself (5393 is prime): pi(5393) >= 5393/(ln 5393 - 1)
 *   upper: pi(y) <= y/(ln y - 1.1) for real y >= 60184
 *          <=> pi(60184) <= 60184/(ln 60184 - 1.1) and pi(q) <= q/(ln q - 1.1) for primes q > 60184
 * Also Rosser-Schoenfeld (3.5): pi(y) > y/ln y for y >= 17 (left limits at primes).
 * Reports the minimum relative slack of each inequality and where; long double.
 * Also reports the largest prime q below the range where the lower bound fails
 * (to confirm the 5393 threshold).  Single-threaded segmented sieve, own code.
 * Usage: dusart_numeric Y
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv) {
    uint64_t Y = strtoull(argv[1], 0, 10);
    uint64_t sq = (uint64_t)sqrtl((long double)Y) + 2;
    uint8_t *small = calloc(sq + 1, 1);
    uint32_t *bp = malloc(4 * (sq / 2 + 10)); size_t nbp = 0;
    for (uint64_t i = 2; i <= sq; i++) if (!small[i]) { bp[nbp++] = (uint32_t)i; for (uint64_t j = i * i; j <= sq; j += i) small[j] = 1; }
    const uint64_t SEG = 1u << 22;
    uint8_t *seg = malloc(SEG);
    uint64_t pi = 0;
    long double minLo = 1e30L, minUp = 1e30L, minRS = 1e30L; uint64_t atLo = 0, atUp = 0, atRS = 0;
    uint64_t lastFailLo = 0, lastFailUp = 0, lastFailRS = 0;
    uint64_t *nextm = malloc(8 * nbp);
    for (size_t j = 0; j < nbp; j++) nextm[j] = (uint64_t)bp[j] * bp[j];
    int checked60184 = 0;
    for (uint64_t lo = 0; lo <= Y; lo += SEG) {
        uint64_t hi = lo + SEG; if (hi > Y + 1) hi = Y + 1;
        memset(seg, 1, hi - lo);
        for (size_t j = 0; j < nbp; j++) {
            uint64_t p = bp[j]; if (p * p >= hi) break;
            uint64_t m = nextm[j];
            for (; m < hi; m += p) seg[m - lo] = 0;
            nextm[j] = m;
        }
        for (uint64_t n = lo; n < hi; n++) {
            if (n == 60184 && !checked60184) {
                long double l = logl((long double)n), b = n / (l - 1.1L), s = (b - pi) / b;
                if (s < minUp) { minUp = s; atUp = n; }
                if (s < 0) lastFailUp = n;
                checked60184 = 1;
            }
            if (n < 2 || !seg[n - lo]) continue;
            pi++;
            long double q = (long double)n, l = logl(q);
            if (n >= 5393) {
                long double b = q / (l - 1.0L);
                long double lhs = (n == 5393) ? (long double)pi : (long double)(pi - 1);
                long double s = (lhs - b) / b;
                if (s < minLo) { minLo = s; atLo = n; }
                if (s < 0) lastFailLo = n;
            } else if (n > 3) {
                /* below threshold: record last failure of the left-limit inequality */
                long double b = q / (l - 1.0L);
                if ((long double)(pi - 1) < b) lastFailLo = n;
            }
            if (n > 60184) {
                long double b = q / (l - 1.1L), s = (b - pi) / b;
                if (s < minUp) { minUp = s; atUp = n; }
                if (s < 0) lastFailUp = n;
            } else if (n > 10) {
                long double b = q / (l - 1.1L);
                if ((long double)pi > b) lastFailUp = n;
            }
            if (n >= 17) {
                long double b = q / l;
                long double s = ((long double)(pi - 1) - b) / b;   /* left limit at q (for q > 17) */
                if (n == 17) s = ((long double)pi - b) / b;
                if (s < minRS) { minRS = s; atRS = n; }
                if (s <= 0) lastFailRS = n;
            } else if (n > 2) {
                if ((long double)(pi - 1) <= q / l) lastFailRS = n;
            }
        }
    }
    printf("Y=%llu pi(Y)=%llu\n", (unsigned long long)Y, (unsigned long long)pi);
    printf("lower  pi(y) >= y/(ln y - 1), y>=5393:   min rel slack %.3Le at prime %llu; last failing prime (any range) %llu\n", minLo, (unsigned long long)atLo, (unsigned long long)lastFailLo);
    printf("upper  pi(y) <= y/(ln y - 1.1), y>=60184: min rel slack %.3Le at %llu; last failing point (any range) %llu\n", minUp, (unsigned long long)atUp, (unsigned long long)lastFailUp);
    printf("RS3.5  pi(y) > y/ln y, y>=17:            min rel slack %.3Le at prime %llu; last failing prime %llu\n", minRS, (unsigned long long)atRS, (unsigned long long)lastFailRS);
    return 0;
}
