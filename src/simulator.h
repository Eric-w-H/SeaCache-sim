#ifndef SIMULATOR_H
#define SIMULATOR_H

struct config {
    enum DataFlow   dataflow;
    enum InterOrder interorder;
    enum Format     format;

    Coord   I, // num rows in A
            J, // num cols in A, equiv. num rows in B
            K; // num cols in B
    Coord   tti,
            ttj,
            ttk;
    // block (tile?) size
    Coord   iii,
            jjj,
            kkk;
};

struct simulator_state {
    struct config cfg;
};


extern struct Arena *global_persist;
extern struct Arena *global_temp;
extern struct simulator_state sim;

extern int prefetchSize;
extern bool adaptive_prefetch;

void configPartial(f32 partialA, f32 partialB, f32 partialC);
void reinitialize();
void runTile(int kkk);
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

struct simulator_state initialize_simulator(const struct config *cfg);

#endif
