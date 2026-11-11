// greedy_qubo.c -- O(N^2) greedy for dense QUBO ground states (confidence-ordered greedy alignment + descent).
// See README.md in this folder for the algorithm step by step.
//
// Problem:  minimize E(x) = - sum_{i,j} q_ij x_i x_j,  x_i in {0,1},  q symmetric, zero diagonal.
//           The double sum runs over both orders of each pair. Reported energy density: e = E / N^{3/2}.
//
// Build:    cc -O3 -o greedy_qubo greedy_qubo.c -lm
// Run:      ./greedy_qubo instance.txt
//              instance.txt: first line N, then one line "i j q_ij" per nonzero pair (1-based, i < j).
//              Prints E, E/N^{3/2}, the time, and the 0/1 solution.
//           ./greedy_qubo -gen N ninst seed
//              Generates ninst random instances of the Fig. 6 ensemble of Boettcher, PRR 1, 033142 (2019)
//              (integers uniform in [-100,100], each pair nonzero with probability 0.1), prints the mean of
//              E/N^{3/2} +- standard error, and the reference value from that paper's data where available.
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
    int n;
    double *q;          // n*n, row-major, symmetric, zero diagonal
    double *h;          // emergent field h_i = sum_j q_ij
    signed char *s;     // spins: 0 = undecided, +1 (x=1) or -1 (x=0)
    double *f;          // local field f_i = h_i + sum_{decided j} q_ij s_j
} Qubo;

// Fix spin i along its field (once, final). Updates every field: one O(N) row.
static void assign(Qubo *P, int i) {
    int n = P->n; signed char v = P->f[i] >= 0 ? 1 : -1; P->s[i] = v;
    const double *r = P->q + (size_t)i * n;
    for (int j = 0; j < n; j++) P->f[j] += r[j] * v;
}
// Flip spin i (only called when it strictly lowers the energy). One O(N) row.
static void flip(Qubo *P, int i) {
    int n = P->n; double d = -2.0 * P->s[i]; const double *r = P->q + (size_t)i * n;
    for (int j = 0; j < n; j++) P->f[j] += d * r[j];
    P->s[i] = -P->s[i];
}

static void greedy(Qubo *P) {
    int n = P->n;
    for (int i = 0; i < n; i++) { P->s[i] = 0; P->f[i] = P->h[i]; }

    // R_i: remaining coupling strength, sum of q_ij^2 over unfixed j. c: small regulariser.
    double *R = malloc(n * sizeof *R), mean = 0;
    for (int i = 0; i < n; i++) { const double *r = P->q + (size_t)i * n; double a = 0;
        for (int j = 0; j < n; j++) a += r[j] * r[j]; R[i] = a; mean += a; }
    double c = 0.01 * mean / n + 1.0;

    for (int t = 0; t < n; t++) {
        // confidence: field squared over remaining coupling strength^(5/4); fix the most confident spin
        int b = -1; double bv = -1;
        for (int i = 0; i < n; i++) if (!P->s[i]) {
            double rc = R[i] + c, v = P->f[i] * P->f[i] / (rc * sqrt(sqrt(rc)));
            if (v > bv) { bv = v; b = i; }
        }
        assign(P, b);
        const double *rb = P->q + (size_t)b * n;
        for (int j = 0; j < n; j++) R[j] -= rb[j] * rb[j];
    }
    free(R);

    // steepest descent: flip the spin that lowers the energy most, until none does
    for (;;) {
        int b = -1; double bv = -1e-9;
        for (int i = 0; i < n; i++) { double v = P->s[i] * P->f[i]; if (v < bv) { bv = v; b = i; } }
        if (b < 0) break;
        flip(P, b);
    }
}

// ---------------------------------------------------------------------------------------------------------
static double energy(const Qubo *P) {  // E = -sum_{i,j} q_ij x_i x_j, from scratch
    int n = P->n; double E = 0;
    for (int i = 0; i < n; i++) if (P->s[i] > 0) { const double *r = P->q + (size_t)i * n;
        for (int j = 0; j < n; j++) if (P->s[j] > 0) E -= r[j]; }
    return E;
}
static void alloc(Qubo *P, int n) {
    P->n = n; P->q = calloc((size_t)n * n, sizeof(double)); P->h = malloc(n * sizeof(double));
    P->s = malloc(n); P->f = malloc(n * sizeof(double));
}
static void fields(Qubo *P) {
    for (int i = 0; i < P->n; i++) { const double *r = P->q + (size_t)i * P->n; double a = 0;
        for (int j = 0; j < P->n; j++) a += r[j]; P->h[i] = a; }
}
static uint64_t splitmix(uint64_t *s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31);
}
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

// <E/N^{3/2}> for the Fig. 6 ensemble from Boettcher's data (tau-EO, averaged over many instances)
static const double REF[][2] = { {31, -9.666594}, {44, -10.073529}, {63, -10.317722}, {80, -10.425666}, {100, -10.508364},
    {127, -10.563889}, {160, -10.611198}, {255, -10.678669}, {511, -10.749101}, {1023, -10.775857}, {4095, -10.790297} };

int main(int argc, char **argv) {
    if (argc >= 5 && !strcmp(argv[1], "-gen")) {
        int n = atoi(argv[2]), K = atoi(argv[3]); uint64_t seed0 = strtoull(argv[4], 0, 10);
        double a = 0, aa = 0, T = 0, N15 = pow(n, 1.5); Qubo P; alloc(&P, n);
        for (int k = 0; k < K; k++) {
            uint64_t s = seed0 * 1000003ULL + k; memset(P.q, 0, (size_t)n * n * sizeof(double));
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++)
                if ((splitmix(&s) >> 11) * 0x1.0p-53 < 0.1) {
                    double v = (double)((int)(splitmix(&s) % 201) - 100); P.q[(size_t)i * n + j] = P.q[(size_t)j * n + i] = v; }
            fields(&P);
            double t0 = now(); greedy(&P); T += now() - t0;
            double e = energy(&P) / N15; a += e; aa += e * e;
        }
        double mu = a / K, se = K > 1 ? sqrt((aa / K - mu * mu) / (K - 1)) : 0;
        printf("N=%d  instances=%d  <E/N^1.5> = %.4f +- %.4f   %.3f ms/instance\n", n, K, mu, se, 1e3 * T / K);
        for (size_t r = 0; r < sizeof REF / sizeof REF[0]; r++) if ((int)REF[r][0] == n)
            printf("reference (Boettcher, Fig. 6 data): %.4f   difference: %+.2f%%\n", REF[r][1], 100 * (mu - REF[r][1]) / -REF[r][1]);
        return 0;
    }
    if (argc < 2) { fprintf(stderr, "usage: %s instance.txt  |  %s -gen N ninst seed\n", argv[0], argv[0]); return 1; }
    FILE *fp = fopen(argv[1], "r"); if (!fp) { perror(argv[1]); return 1; }
    int n; if (fscanf(fp, "%d", &n) != 1 || n < 2) { fprintf(stderr, "bad header\n"); return 1; }
    Qubo P; alloc(&P, n); int i, j; double v;
    while (fscanf(fp, "%d %d %lf", &i, &j, &v) == 3) {
        if (i < 1 || j < 1 || i > n || j > n || i == j) { fprintf(stderr, "bad pair %d %d\n", i, j); return 1; }
        P.q[(size_t)(i - 1) * n + (j - 1)] = P.q[(size_t)(j - 1) * n + (i - 1)] = v;
    }
    fclose(fp); fields(&P);
    double t0 = now(); greedy(&P); double t = now() - t0, E = energy(&P);
    printf("N=%d  E=%.10g  E/N^1.5=%.6f  time=%.3f ms\nx=", n, E, E / pow(n, 1.5), 1e3 * t);
    for (int k = 0; k < n; k++) putchar(P.s[k] > 0 ? '1' : '0');
    putchar('\n'); return 0;
}
