#include "data.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "headers.h"


u64 **A   = NULL;  // A[r]  = [c_1, c_2, ..., c_u] are non-zero column indices in row r
u64 **Ac  = NULL; // Ac[c] = [r_1, r_2, ..., r_v] are non-zero row indices in column c
u64 **B   = NULL;
u64 **Bc  = NULL;

std::vector<int> *sparchA   = nullptr;

u64 *offsetarrayA     = nullptr;
u64 *offsetarrayAc    = nullptr;
u64 *offsetarrayB     = nullptr;
u64 *offsetarrayBc        = nullptr;

// Let row r be sampled from A, and it is the (r')th sampled row.
// Then for all indices [(r,c_1), (r,c_2), ..., (r,c_u)] in the row:
std::vector<u64> *SA = nullptr;  //  SA[r']  = [c_1, c_2, ..., c_u]
std::vector<u64> *SAc = nullptr; //  SAc[c]  = [hash1(r, ...)] x u (list of u copies of same hash)

// Let col c be sampled from B, and it is the (c')th sampled row.
// Then for all indices [(r_1,c), (r_2,c), ..., (r_v,c)] in the col:
std::vector<u64> *SBc = nullptr; // SBc[c']  = [r_1, r_2, ..., r_v]
std::vector<u64> *SB = nullptr;  // SB[r]    = [hash2(c, ...)] x v (list of v copies of same hash)

extern double ha1, hb1;
extern double ha2, hb2;
