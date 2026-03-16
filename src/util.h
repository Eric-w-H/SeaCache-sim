#ifndef UTIL_H
#define UTIL_H

#include <random>
#include <vector>

int getlog(int x);

// Global variables
static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<> dis(0.0, 1.0);

double hash1(int x, double a, double b, int pmod);
double hash2(int x, double a, double b, int pmod);
double getRandomCoefficient();
bool sampleP();

#endif // UTIL_H
