#ifndef SIMULATOR_H
#define SIMULATOR_H

extern int cachesize;
extern int prefetchSize;
extern bool adaptive_prefetch;

extern int cacheScheme;

void configPartial(float partialA, float partialB, float partialC);
void reinitialize();
void runTile(int iii, int jjj, int kkk, long long tti, long long ttj, long long ttk);
void run();

extern int *currsizeB;
extern int *currsizeBc;
extern int *beginB;

extern int TI, TJ, TK;

// Tcnt are for counting non-zeros in each subtile
// initialize to 0 at each round
extern int Tcnt[2][2];
// store the now sum of sizejk in ttjsum*ttksum tiles
// store all 7 types
extern long long sizejksum[10];
// ttjsum, ttksum means the now ttj*ttk
// initialize to zero if jjj*kkk changes,
extern int tilesum;

void initialize_simulator();
void deinitialize_simulator();

#endif
