#ifndef DATA_H
#define DATA_H

#include <queue>
#include <set>
#include <unordered_set>
#include <vector>
#include "headers.h"

static const int MAXN = 3000000;

// sparse matrices A, B and their transposes Ac, Bc
extern std::vector<u64> *A, *Ac;
extern std::vector<int> *sparchA;
extern std::vector<u64> *B, *Bc;

// store the offsets for A, Ac, B, Bc
extern u64 *offsetarrayA, *offsetarrayAc;
extern u64 *offsetarrayB, *offsetarrayBc;

// sample matrix
extern int SI, SK;
extern std::vector<u64> *SA, *SAc;
extern std::vector<u64> *SBc, *SB;

#endif // DATA_H
