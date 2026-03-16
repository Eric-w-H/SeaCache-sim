#include "util.h"
#include "headers.h"
#include <cmath>

// calculate log base 2 of x
int getlog(int x) {
    int ret = -1;
    while (x) {
        ret++;
        x >>= 1;
    }
    return ret;
}

// Hash function h1 and h2
double hash1(int x, double a, double b, int pmod) {
    // printf("%lf %lf\n", a*x+b, (double)(std::fmod(a * x + b, pmod)));
    return (double)(std::fmod(a * x + b, pmod)) / pmod;
}

double hash2(int x, double a, double b, int pmod) {
    return (double)(std::fmod(a * x + b, pmod)) / pmod;
}

// Function to generate random coefficients
double getRandomCoefficient() { return dis(gen); }

bool sampleP() { return dis(gen) < samplep; }