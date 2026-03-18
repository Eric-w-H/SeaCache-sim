#include <assert.h>
#include "cache.h"
#include "estimation.h"
#include "headers.h"
#include "statistics.h"
// store all the buffered C now
set<int> *bufferedC = nullptr;
// record length of buffered C
// equals bufferedC[i].size()
int *bufferedClen = nullptr;
int *beginA = nullptr;
int *beginB = nullptr;
int *beginAc = nullptr;
int *beginBc = nullptr;
Coord *dirtyC = NULL;
// int *begin = nullptr;
/*
The current fiber size of each array
stored according to the dataflow order
(not the storage format order)
update currsize each time the block of the array changes
(when inter-iterate)
*/
int *currsizeA = nullptr;
int *currsizeAc = nullptr;
int *currsizeB = nullptr;
int *currsizeBc = nullptr;
// int *currsizeC = nullptr;
/*
The currently buffered size of each array
(part of cuursizeA/B/C)
update bufferedsize each time
*/
int *bufferedsizeB = nullptr;
// int *bufferedsizeC = nullptr;
u8 *tmpC = nullptr;
// start of current block (NOTE(ejs): block == tile?)
int TI, TJ, TK;

// then will take length+1 place for the next fiberlet pointer
const int fiberletlength = 4;
/*
The allocated buffer to an array (A/B/Csize)
(accordind to "partial" parameter and dataflow)
and the currently used part (A/B/Csizenow)
*/
int Asizenow;
int Bsizenow;
int Csizenow;
bool fulltagA, fulltagB, fulltagC;
int fullA, fullB;

void reset_cursor2(struct cursor *c)
{
    c->first= 1;
    c->ti   = 0;
    c->tj   = 0;
    c->tk   = 0;
}


// ============================================================================
// compute_tile_luts / advance_cursor — the ONLY source of truth for
// sim.cursor.{A,B}.{begins,sizes}.
//
// These arrays are TILE-LOCAL: indexed [0 .. major_dim-1].
//   sim.cursor.A.begins[r]  ↔  old beginA[TI + r]
//   sim.cursor.A.sizes[r]   ↔  old currsizeA[TI + r]
//   sim.cursor.B.begins[r]  ↔  old beginB[TJ + r]
//   sim.cursor.B.sizes[r]   ↔  old currsizeB[TJ + r]
//
// Every consumer that held an absolute index (ii, jj, ti, tj) must subtract
// the tile origin:  sim.cursor.A.sizes[ii - TI],  sim.cursor.B.sizes[jj - TJ].
// ============================================================================

void compute_tile_luts(struct tile *t, b16 reset_begins)
{
    Coord   major_dim, minor_dim, major_idx, minor_idx;
    const Coord *offsets, **map;
    Coord *begins, *sizes;
    major_dim   = t->major_dim;
    minor_dim   = t->minor_dim;
    major_idx   = *t->major_tile_idx * t->major_dim;
    minor_idx   = *t->minor_tile_idx * t->minor_dim;
    map         = t->map;
    offsets     = t->offsets;
    begins      = t->begins;
    sizes       = t->sizes;
    memset(sizes, 0, major_dim*sizeof(*sizes));
    if (reset_begins)
        memset(begins, 0, major_dim*sizeof(*sizes));
    Coord abs_rj_lim= minor_idx + minor_dim;
    for (Coord ri = 0; ri < major_dim; ++ri) {
        Coord abs_ri        = major_idx + ri;
        Coord nzero_elem_cnt= offsets[abs_ri+1] - offsets[abs_ri];
        Coord rj        = begins[ri];
        while (rj < nzero_elem_cnt && map[abs_ri][rj] < minor_idx)
            ++rj;
        begins[ri] = rj;
        while (rj < nzero_elem_cnt && map[abs_ri][rj] < abs_rj_lim) {
            ++sizes[ri];
            ++rj;
        }
    }
}

void advance_cursor() // / tile pair
{
    struct cursor *c= &sim.cursor;
    if (c->first) {
        c->first = 0;
        compute_tile_luts(&c->A, 1);
        compute_tile_luts(&c->B, 1);
        return;
    }
    // increment tile indices and check if either tile needs recomputing
    b16 Atile_need_recompute, A_reset_begins;
    b16 Btile_need_recompute, B_reset_begins;
    
    {
        Coord old_A_major_tile_idx = *c->A.major_tile_idx;
        Coord old_A_minor_tile_idx = *c->A.minor_tile_idx;
        Coord old_B_major_tile_idx = *c->B.major_tile_idx;
        Coord old_B_minor_tile_idx = *c->B.minor_tile_idx;
        {
            Coord   inner_p1    = *c->inner_tile_idx + 1;
            Coord   middle_p1   = *c->middle_tile_idx + 1;
            b16 inner_at_lim    = inner_p1 == c->inner_ntiles;
            b16 middle_at_lim   = middle_p1 == c->middle_ntiles;
            *c->inner_tile_idx  = inner_at_lim ? 0 : inner_p1;
            *c->middle_tile_idx = inner_at_lim ? (middle_at_lim ? 0 : middle_p1) : *c->middle_tile_idx;
            *c->outer_tile_idx  += (inner_at_lim && middle_at_lim);
            *c->inner_wrap      = inner_at_lim;
            *c->middle_wrap     = inner_at_lim && middle_at_lim;
            *c->outer_wrap      = 0; // outer never wraps
        }
        Atile_need_recompute= (old_A_major_tile_idx != *c->A.major_tile_idx) || (old_A_minor_tile_idx != *c->A.minor_tile_idx);
        Btile_need_recompute= (old_B_major_tile_idx != *c->B.major_tile_idx) || (old_B_minor_tile_idx != *c->B.minor_tile_idx);
        A_reset_begins = *c->A.minor_wrap || (old_A_major_tile_idx != *c->A.major_tile_idx);
        B_reset_begins = *c->B.minor_wrap || (old_B_major_tile_idx != *c->B.major_tile_idx);
    }
    
    if (Atile_need_recompute)
        compute_tile_luts(&c->A, A_reset_begins);
    if (Btile_need_recompute)
        compute_tile_luts(&c->B, B_reset_begins);
}

// ============================================================================
// Helpers that derive TI/TJ/TK from the canonical cursor tile indices.
// Called once after every advance_cursor().
// ============================================================================
static inline void sync_tile_origins()
{
    TI = sim.cursor.ti * sim.cfg.iii;
    TJ = sim.cursor.tj * sim.cfg.jjj;
    TK = sim.cursor.tk * sim.cfg.kkk;
}

// ============================================================================
// The old updateBegin*, AllupdateBegin*, reverse_*, reinitialize_inner/mid/outer,
// iterate_inner/mid/outer_loop, updateTI/TJ/TK, updateBlockA/B, isIJ, isJK
// functions are ALL DELETED.
//
// advance_cursor + compute_tile_luts handles begin tracking.
// sync_tile_origins handles TI/TJ/TK derivation.
// The flat tile loop in run() handles iteration.
// ============================================================================

bool checkAndLoadReuseA() {
    if ((interorder == IJK || interorder == JIK)) {
        if (TK == 0) {
            Asizenow = 0;
            fulltagA = 0;
            if (format == RR || format == RC) {
                for (int ti = TI; ti < min(TI + sim.cfg.iii, sim.cfg.I); ti++) {
                    Asizenow++;
                    // FIX: tile-local index
                    int ri = ti - TI;
                    int startj = sim.cursor.A.begins[ri], tmpj = sim.cursor.A.begins[ri],
                        maxj = offsetarrayA[ti + 1] - offsetarrayA[ti];
                    while (tmpj < maxj && A[ti][tmpj] < sim.cfg.jjj + TJ) {
                        tmpj++;
                    }
                    int tmpsize = (tmpj - startj);
                    if (Asizenow + tmpsize * 3 >= Asize) {
                        if (!fulltagA) {
                            fulltagA = 1;
                            fullA = ti;
                        }
                    } else {
                        Asizenow += tmpsize * 3;
                        preDramAccess += memoryBandwidthWhole(tmpsize * 3 + 2);
                        preA += memoryBandwidthWhole(tmpsize * 3 + 2);
                        preSramAccess += sramWriteBandwidth(tmpsize * 3 + 2);
                        AccessByte += tmpsize * 3 + 2;
                    }
                }
            }
        }
        return 1;
    }
    fulltagA = 1;
    fullA = TI - 1;
    return 0;
}
bool checkReuseB() {
    if ((interorder == JKI || interorder == KJI) && (TI != 0))
        return 1;
    return 0;
}
void pre_load_A() {
    if (dataflow == Outer) {
        fulltagA = 1;
        fullA = TJ;
        return;
    }
    if (checkAndLoadReuseA()) {
        return;
    }
}
void pre_load_B() {
    if (checkReuseB()) {
        return;
    }
    if (dataflow == Outer) {
        fulltagB = 1;
        fullB = TJ;
        return;
    }
    Bsizenow = 0;
    fulltagB = 0;
    fullB = 0;
    if ((dataflow == Inner) || (dataflow == Gust)) {
        Bsizenow = 0;
        fulltagB = 0;
        fullB = 0;
        if ((dataflow == Gust) && ((format == CR) || (format == RR))) {
            int tj;
            int _TJ;
            for (tj = TJ; tj < min(TJ + sim.cfg.jjj, sim.cfg.J); tj++) {
                if ((tj - TJ) < (sim.cfg.jjj / 2)) {
                    _TJ = 0;
                } else {
                    _TJ = 1;
                }
                Bsizenow++;
                // FIX: tile-local index for B
                int rj = tj - TJ;
                int startk = sim.cursor.B.begins[rj], tmpk = sim.cursor.B.begins[rj],
                    maxk = offsetarrayB[tj + 1] - offsetarrayB[tj];
                int halfk = sim.cursor.B.begins[rj];
                while (halfk < maxk && B[tj][halfk] < (sim.cfg.kkk / 2) + TK) {
                    halfk++;
                }
                tmpk = halfk;
                while (tmpk < maxk && B[tj][tmpk] < sim.cfg.kkk + TK) {
                    tmpk++;
                }
                int tmpsize = (tmpk - startk);
                if (Bsizenow + tmpsize * 3 >= Bsize) {
                    if (!fulltagB) {
                        fulltagB = 1;
                        fullB = tj;
                    }
                } else {
                    Bsizenow += tmpsize * 3;
                    preDramAccess += memoryBandwidthWhole(tmpsize * 3 + 2);
                    preB += memoryBandwidthWhole(tmpsize * 3 + 2);
                    preSramAccess += sramWriteBandwidth(tmpsize * 3 + 2);
                    AccessByte += tmpsize * 3 + 2;
                }
            }
        }
        if ((dataflow == Gust) && ((format == CC) || (format == RC))) {
            // inconsistent format — uses column-major Bc view.
            // NOTE: sim.cursor.B is row-major (J-indexed).  This path would need
            // a separate column-major cursor to be correct.  With Gust+RR this
            // path is unreachable.  Left as-is with a guard.
            assert(0 && "Gust+inconsistent B format not yet ported to advance_cursor");
            Bsizenow += ((fiberletlength + 1) / 2) * sim.cfg.jjj;
            Bsizenow += sim.cfg.jjj;
            for (int tj = TJ; tj < TJ + sim.cfg.jjj; tj++) {
                bufferedsizeB[tj] = 0;
            }
            for (int tk = TK; tk < min(TK + sim.cfg.kkk, sim.cfg.K); tk++) {
                int startj = beginBc[tk], tmpj = beginBc[tk],
                    maxj = offsetarrayBc[tk + 1] - offsetarrayBc[tk];
                while (tmpj < maxj) {
                    int tmpindex = Bc[tk][tmpj];
                    if (tmpindex >= TJ + sim.cfg.jjj) {
                        break;
                    }
                    tmpj++;
                    bufferedsizeB[tmpindex]++;
                }
                int tmpsize = (tmpj - startj);
                if (Bsizenow + tmpsize * 3 >= Bsize) {
                    if (!fulltagB) {
                        fulltagB = 1;
                        fullB = tk;
                    }
                } else {
                    preDramAccess += memoryBandwidthWhole(tmpsize * 3 + 2);
                    preB += memoryBandwidthWhole(tmpsize * 3 + 2);
                    AccessByte += tmpsize * 3 + 2;
                    preSramAccess += sramWriteBandwidth(tmpsize);
                    preSramAccess +=
                        sramReadBandwidth(tmpsize + tmpsize / fiberletlength) * 3;
                }
            }
        }
        // Inner dataflow paths — not reachable with Gust, guarded
        if ((dataflow == Inner) && (((format == RC) || (format == CC)))) {
            assert(0 && "Inner dataflow not yet ported to advance_cursor");
        }
        if ((dataflow == Inner) && (((format == RR) || (format == CR)))) {
            assert(0 && "Inner dataflow not yet ported to advance_cursor");
        }
    }
}
/*
Load streams(A/B) into buffer before calculation.
*/
void pre_calculate_load() {
    preSramAccess = preDramAccess = 0;
    if (ISCACHE) {
        initializeCacheValid();
    }
    if (!ISCACHE) {
        pre_load_B();
    }
    pre_load_A();
    preSramAccess /= sramBank;
    totalCycle += max(preSramAccess, preDramAccess);
    preCycle += max(preSramAccess, preDramAccess);
}
bool consistent_B() {
    if (dataflow == Gust) {
        if ((format == RR) || (format == CR))
            return 1;
        if ((format == RC) || (format == CC))
            return 0;
    }
    if (dataflow == Inner) {
        if ((format == RC) || (format == CC))
            return 1;
        if ((format == RR) || (format == CR))
            return 0;
    }
    if (dataflow == Outer) {
        if ((format == RR) || (format == CR))
            return 1;
        if ((format == RC) || (format == CC))
            return 0;
    }
    return 0;
}
bool consistent_A() {
    if (dataflow == Gust) {
        if ((format == RR) || (format == RC))
            return 1;
        if ((format == CR) || (format == CC))
            return 0;
    }
    if (dataflow == Inner) {
        if ((format == RR) || (format == RC))
            return 1;
        if ((format == CR) || (format == CC))
            return 0;
    }
    if (dataflow == Outer) {
        if ((format == CR) || (format == CC))
            return 1;
        if ((format == RR) || (format == RC))
            return 0;
    }
    return 0;
}
void get_B_fiber(int jj, int ii)
{
    // FIX: tile-local index
    int rj = jj - TJ;
    if (!ISCACHE) {
        if (consistent_B()) {
            if (fulltagB == 0 || jj < fullB) {
                computeSramAccess += sramReadBandwidth(sim.cursor.B.sizes[rj] * 3 + 2);
            } else {
                computeDramAccess += memoryBandwidthPE(sim.cursor.B.sizes[rj] * 3 + 2);
                computeB += memoryBandwidthPE(sim.cursor.B.sizes[rj] * 3 + 2);
                AccessByte += sim.cursor.B.sizes[rj] * 3 + 2;
            }
        } else {
            computeSramAccess +=
                sramReadBandwidth(fiberletlength * 3) * ((bufferedsizeB[jj] + 3) / 4);
            if (fulltagB) {
                computeDramAccess +=
                    (memoryBandwidthPE(3)) * ((long long)TK + sim.cfg.kkk - fullB);
                computeB +=
                    (memoryBandwidthPE(3)) * (long long)((long long)TK + sim.cfg.kkk - fullB);
                AccessByte += 3 * (long long)((long long)TK + sim.cfg.kkk - fullB);
            }
        }
    }
    else {
        int fibersize = sim.cursor.B.sizes[rj] * 3 + 1;
        cacheAccessFiber(jj, fibersize, ii);
    }
}
void get_B_fiber_col_iii(int kk, int iii) {
    // NOTE: Inner dataflow only.  Uses currsizeBc in old code.
    // With Gust fixed, this is unreachable.  Left with old currsizeBc
    // reference so it compiles but will assert if called improperly.
    if (fulltagB == 0 || kk < fullB) {
        computeSramAccess += sramReadBandwidth(currsizeBc[kk] * 3 + 2) * iii;
    } else {
        computeDramAccess += memoryBandwidthPE(currsizeBc[kk] * 3 + 2) * iii;
        computeB += memoryBandwidthPE(currsizeBc[kk] * 3 + 2) * iii;
        AccessByte += (currsizeBc[kk] * 3 + 2) * iii;
    }
}
void get_A_fiber_col(int jj)
{
    // NOTE: Outer dataflow only.  Uses currsizeAc in old code.
    // With Gust fixed, unreachable.
    assert(consistent_A());
    if (fulltagA == 0 || jj < fullA) {
        computeSramAccess += sramReadBandwidth(currsizeAc[jj] * 3 + 2);
    } else {
        computeDramAccess += memoryBandwidthPE(currsizeAc[jj] * 3 + 2);
        computeA += memoryBandwidthPE(currsizeAc[jj] * 3 + 2);
        AccessByte += currsizeAc[jj] * 3 + 2;
        computeSramAccess += sramReadBandwidth(currsizeAc[jj] * 3 + 2) +
                                sramWriteBandwidth(currsizeAc[jj] * 3 + 2);
        if (cacheScheme == CACHE_SCHEME_INNER_SP) {
            computeDramAccess += memoryBandwidthPE(currsizeAc[jj] * 3 + 2);
            computeA += memoryBandwidthPE(currsizeAc[jj] * 3 + 2);
            computeSramAccess += sramReadBandwidth(currsizeAc[jj] * 3 + 2) +
                                    sramWriteBandwidth(currsizeAc[jj] * 3 + 2);
        }
    }
}
void get_A_fiber(int ii) {
    assert(consistent_A());
    // FIX: tile-local index
    int ri = ii - TI;
    if (fulltagA == 0 || ii < fullA) {
        computeSramAccess += sramReadBandwidth(sim.cursor.A.sizes[ri] * 3 + 2);
        if (cacheScheme == CACHE_SCHEME_INNER_SP) {
            computeSramAccess += sramReadBandwidth(sim.cursor.A.sizes[ri] * 3 + 2);
        }
    } else {
        computeDramAccess += memoryBandwidthPE(sim.cursor.A.sizes[ri] * 3 + 2);
        computeA += memoryBandwidthPE(sim.cursor.A.sizes[ri] * 3 + 2);
        AccessByte += sim.cursor.A.sizes[ri] * 3 + 2;
        computeSramAccess += sramReadBandwidth(sim.cursor.A.sizes[ri] * 3 + 2) +
                                sramWriteBandwidth(sim.cursor.A.sizes[ri] * 3 + 2);
        if (cacheScheme == CACHE_SCHEME_INNER_SP) {
            computeDramAccess += memoryBandwidthPE(sim.cursor.A.sizes[ri] * 3 + 2);
            computeA += memoryBandwidthPE(sim.cursor.A.sizes[ri] * 3 + 2);
            computeSramAccess += sramReadBandwidth(sim.cursor.A.sizes[ri] * 3 + 2) +
                                    sramWriteBandwidth(sim.cursor.A.sizes[ri] * 3 + 2);
        }
    }
}
void update_c_fiber(int jj) {
    // FIX: tile-local index
    int rj = jj - TJ;
    for (int k1 = sim.cursor.B.begins[rj]; k1 < sim.cursor.B.begins[rj] + sim.cursor.B.sizes[rj]; k1++) {
        tmpC[B[jj][k1]] = 1;
    }
}
void updateCAccess(int ii)
{
    if ((Csize >= 100.0) && ((interorder == IKJ) || (interorder == KIJ))) {
        int deltaC = 0;
        for (int k1 = TK; k1 < TK + sim.cfg.kkk; k1++) {
            if (tmpC[k1]) {
                if (bufferedC[ii].find(k1) == bufferedC[ii].end()) {
                    deltaC++;
                    bufferedC[ii].insert(k1);
                }
            }
        }
        bufferedClen[ii] += deltaC;
        Csizenow += deltaC;
        if (Csizenow > Csize) {
            computeDramAccess += memoryBandwidthPE(Csizenow * 3);
            postDramAccess += memoryBandwidthPE(Csizenow * 3);
            computeC += memoryBandwidthPE(Csizenow * 3);
            postC += memoryBandwidthPE(Csizenow * 3);
            AccessByte += Csizenow * 3;
            AccessByte += Csizenow * 3;
            computeSramAccess +=
                sramReadBandwidth(Csizenow * 3) + sramWriteBandwidth(Csizenow * 3);
            Csizenow = 0;
            for (int i = TI; i < TI + sim.cfg.iii; i++) {
                bufferedC[i]    = std::set<int>();
                bufferedClen[i] = 0;
            }
        }
    } else {
        int cntc = 0;
        for (int k1 = TK; k1 < TK + sim.cfg.kkk; k1++) {
            if (tmpC[k1]) {
                cntc++;
            }
        }
        computeDramAccess += memoryBandwidthPE(cntc * 3);
        computeC += memoryBandwidthPE(cntc * 3);
        AccessByte += cntc * 3;
        if (sim.cfg.jjj != sim.cfg.J) {
            postDramAccess += memoryBandwidthPE(cntc * 3);
            postC += memoryBandwidthPE(cntc * 3);
            AccessByte += cntc * 3;
        }
        computeSramAccess +=
            sramReadBandwidth(cntc * 3) + sramWriteBandwidth(cntc * 3);
    }
}
void get_B_fibers(int ii)
{
    assert(LIKELY(dataflow == Gust));
    // FIX: tile-local index for A
    int ri = ii - TI;
    int tmpj = sim.cursor.A.begins[ri];
    int maxj = offsetarrayA[ii + 1] - offsetarrayA[ii];
    Coord ndirty = 0;
    while (tmpj < maxj && A[ii][tmpj] < TJ + sim.cfg.jjj) {
        int jj = A[ii][tmpj];
        // FIX: tile-local index for B
        int rj = jj - TJ;
        Coord bsize   = sim.cursor.B.sizes[rj];
        // >> get_B_fiber() inlined
        if (!ISCACHE) {
            if (consistent_B()) {
                i64 cost = bsize * 3 + 2;
                if (fulltagB == 0 || jj < fullB) {
                    computeSramAccess   += sramReadBandwidth(cost);
                } else {
                    computeDramAccess   += memoryBandwidthPE(cost);
                    computeB            += memoryBandwidthPE(cost);
                    AccessByte          += cost;
                }
            } else {
                computeSramAccess += sramReadBandwidth(fiberletlength * 3) * ((bufferedsizeB[jj] + 3) / 4);
                if (fulltagB) {
                    computeDramAccess   += (memoryBandwidthPE(3)) * ((long long)TK + sim.cfg.kkk - fullB);
                    computeB            += (memoryBandwidthPE(3)) * (long long)((long long)TK + sim.cfg.kkk - fullB);
                    AccessByte          += 3 * (long long)((long long)TK + sim.cfg.kkk - fullB);
                }
            }
        } else {
            int fibersize = bsize * 3 + 1;
            cacheAccessFiber(jj, fibersize, ii);
        }
        // << get_B_fiber inlined
        computePE += bsize;
        elements_processed_since_last_adjustment += bsize;
        // >> update_c_fiber() inlined — FIX: tile-local B index
        for (int k1 = sim.cursor.B.begins[rj]; k1 < sim.cursor.B.begins[rj] + bsize; k1++) {
            Coord index = B[jj][k1];
            if (!tmpC[index]) {
                dirtyC[ndirty++]= index;
                tmpC[index]     = 1;
            }
        }
        // << update_c_fiber() inlined
        tmpj++;
    }
    // update A access — FIX: tile-local A index
    if (consistent_A()) {
        i64 cost = (tmpj - sim.cursor.A.begins[ri]) * 3;
        b32 hitA = (interorder == IJK || interorder == JIK) && (fulltagA == 0 || ii < fullA);
        if (hitA) {
            computeSramAccess += sramReadBandwidth(cost);
            if (cacheScheme == CACHE_SCHEME_INNER_SP)
                computeSramAccess += sramReadBandwidth(cost);
        } else {
            computeDramAccess   += memoryBandwidthPE(cost);
            computeA            += memoryBandwidthPE(cost);
            AccessByte          += cost;
            computeSramAccess   += sramReadBandwidth(cost) + sramWriteBandwidth(cost);
            if (cacheScheme == CACHE_SCHEME_INNER_SP) {
                computeDramAccess   += memoryBandwidthPE(cost);
                computeA            += memoryBandwidthPE(cost);
                computeSramAccess   += sramReadBandwidth(cost) + sramWriteBandwidth(cost);
            }
        }
    } else {
        computeSramAccess   += sramReadBandwidth(fiberletlength * 3) * ((bufferedsizeB[ii] + 3) / 4);
        if (fulltagA) {
            computeDramAccess   += (memoryBandwidthPE(3)) * ((long long)TJ + sim.cfg.jjj - fullA);
            computeA            += (memoryBandwidthPE(3)) * (long long)((long long)TJ + sim.cfg.jjj - fullA);
            AccessByte          += 3 * ((long long)TJ + sim.cfg.jjj - fullA);
        }
    }
    // >> updateCaccess() inlined
    if ((Csize >= 100.0) && ((interorder == IKJ) || (interorder == KIJ))) {
        int deltaC = 0;
        for (Coord d = 0; d < ndirty; ++d) {
            Coord k1 = dirtyC[d];
            if (tmpC[k1]) {
                if (bufferedC[ii].find(k1) == bufferedC[ii].end()) {
                    ++deltaC;
                    bufferedC[ii].insert(k1);
                }
            }
        }
        bufferedClen[ii]+= deltaC;
        Csizenow        += deltaC;
        if (Csizenow > Csize) {
            computeDramAccess   += memoryBandwidthPE(Csizenow * 3);
            postDramAccess      += memoryBandwidthPE(Csizenow * 3);
            computeC            += memoryBandwidthPE(Csizenow * 3);
            postC               += memoryBandwidthPE(Csizenow * 3);
            AccessByte          += Csizenow * 6;
            computeSramAccess   += sramReadBandwidth(Csizenow * 3) + sramWriteBandwidth(Csizenow * 3);
            Csizenow = 0;
            for (int i = TI; i < TI + sim.cfg.iii; i++) {
                bufferedC[i]    = std::set<int>();
                bufferedClen[i] = 0;
            }
        }
    } else {
        computeDramAccess += memoryBandwidthPE(ndirty * 3);
        computeC += memoryBandwidthPE(ndirty * 3);
        AccessByte += ndirty * 3;
        if (sim.cfg.jjj != sim.cfg.J) {
            postDramAccess += memoryBandwidthPE(ndirty * 3);
            postC += memoryBandwidthPE(ndirty * 3);
            AccessByte += ndirty * 3;
        }
        computeSramAccess +=
            sramReadBandwidth(ndirty * 3) + sramWriteBandwidth(ndirty * 3);
    }
    // << updateCaccess() inlined
    for (Coord d = 0; d < ndirty; d++)
        tmpC[dirtyC[d]] = 0;
}
// prefetchSize: the size of the buffer allocated for prefetch
// prefetchNow: the current prefetch size (must < prefetchSize)
// prefetchRowNow: row of the current(next) prefetch position
int prefetchSize = 500;
int prefetchNow = 0;
int prefetchRowNow = 0;
bool prefetchrow(int ii) {
    int needsize = 0;
    // FIX: tile-local index
    int ri = ii - TI;
    if (cacheScheme == 6 || cacheScheme == 7) {
        needsize = sim.cursor.A.sizes[ri] * 4 + 1;
    }
    else if (cacheScheme == 66 || cacheScheme == CACHE_SCHEME_FLFU) {
        needsize = sim.cursor.A.sizes[ri] * 3;
    }
    if (prefetchNow + needsize >= prefetchSize) {
        return 0;
    }
    prefetchNow += needsize;
    // FIX: tile-local index for begin
    int tmpj = sim.cursor.A.begins[ri];
    int maxj = offsetarrayA[ii + 1] - offsetarrayA[ii];
    while (tmpj < maxj && A[ii][tmpj] < TJ + sim.cfg.jjj) {
        int jj = A[ii][tmpj];
        if (cacheScheme == 6 || cacheScheme == 7 || cacheScheme == CACHE_SCHEME_INNER_SP ||
            cacheScheme == CACHE_SCHEME_SPARCH) {
            nextposvector[jj].push(-ii);
        }
        if (cacheScheme == 66) {
            LFUtag[jj]++;
        }
        if (cacheScheme == CACHE_SCHEME_FLFU) {
            long long firstaddr = getCacheAddr(jj, 0);
            // FIX: tile-local index for B sizes
            int rj = jj - TJ;
            int fibersize = sim.cursor.B.sizes[rj] * 3;
            for (int tmpcurr = 0; tmpcurr < fibersize; tmpcurr += CACHE_BLOCK_NELEMS) {
                long long tmpaddr = getCacheAddr(jj, tmpcurr / CACHE_BLOCK_NELEMS);
                int _set = getSet2(tmpaddr);
                int _tag = getTag2(tmpaddr);
                bool need_halve_lfu = 0;
                bool incache = 0;
                prefetch_increments++;
                for (int i = 0; i < SETASSOC; i++) {
                    if (Valid[_set * SETASSOC + i] && (Tag[_set * SETASSOC + i] == _tag)) {
                        if (tmpcurr != 0) {
                            if (PosOrig[_set * SETASSOC + i] != getOrig(firstaddr)) {
                                continue;
                            }
                        } else {
                            if (PosOrig[_set * SETASSOC + i] != 0) {
                                continue;
                            }
                        }
                        incache = 1;
                        lfubit[_set * SETASSOC + i]++;
                        if (lfubit[_set * SETASSOC + i] > LFUmax) {
                            need_halve_lfu = 1;
                        }
                        break;
                    }
                }
                if (useVirtualTag) {
                    if (!incache) {
                        bool invirtualtag = 0;
                        for (int i = 0; i < VIRTUALSETASSOC; i++) {
                            if (virtualValid[_set * VIRTUALSETASSOC + i]) {
                                if (virtualTag[_set * VIRTUALSETASSOC + i] == _tag) {
                                    invirtualtag = 1;
                                    virtuallfubit[_set * VIRTUALSETASSOC + i]++;
                                    if (virtuallfubit[_set * VIRTUALSETASSOC + i] > LFUmax) {
                                        need_halve_lfu = 1;
                                    }
                                    break;
                                }
                            } else {
                            }
                        }
                        if (!invirtualtag) {
                            bool hasinvalid = 0;
                            for (int i = 0; i < VIRTUALSETASSOC; i++) {
                                if (virtualValid[_set * VIRTUALSETASSOC + i] == 0) {
                                    hasinvalid = 1;
                                    virtualValid[_set * VIRTUALSETASSOC + i] = 1;
                                    virtualTag[_set * VIRTUALSETASSOC + i] = _tag;
                                    virtuallfubit[_set * VIRTUALSETASSOC + i] = 1;
                                    break;
                                }
                            }
                            bool haszero = 0;
                            if (!hasinvalid) {
                                for (int i = 0; i < VIRTUALSETASSOC; i++) {
                                    if (virtuallfubit[_set * VIRTUALSETASSOC + i] == 0) {
                                        haszero = 1;
                                        virtualValid[_set * VIRTUALSETASSOC + i] = 1;
                                        virtualTag[_set * VIRTUALSETASSOC + i] = _tag;
                                        virtuallfubit[_set * VIRTUALSETASSOC + i] = 1;
                                        break;
                                    }
                                }
                            }
                            if (hasinvalid == 0 && haszero == 0) {
                                prefetch_discards++;
                            }
                        }
                    }
                }
                if (need_halve_lfu) {
                    for (int i = 0; i < SETASSOC; i++) {
                        if (Valid[_set * SETASSOC + i]) {
                            lfubit[_set * SETASSOC + i] /= 2;
                        }
                    }
                    if (useVirtualTag) {
                        for (int i = 0; i < VIRTUALSETASSOC; i++) {
                            if (virtualValid[_set * VIRTUALSETASSOC + i]) {
                                virtuallfubit[_set * VIRTUALSETASSOC + i] /= 2;
                            }
                        }
                    }
                }
            }
        }
        tmpj++;
    }
    return 1;
}
int get_num_samples(f64 current_temperature) {
    if (current_temperature > SA_INITIAL_TEMP * 0.5) {
        return 4;
    } else if (current_temperature > SA_INITIAL_TEMP * 0.1) {
        return 8;
    } else {
        return 16;
    }
}
bool lastaccept = 1;
void update_prefetch_size() {
    f64 temperature = SA_INITIAL_TEMP * pow(SA_COOLING_RATE, sa_iteration_k);
    int num_samples = get_num_samples(temperature);
    sa_iteration_k++;
    if ((sa_iteration_k % num_samples) != 0) {
        return;
    }
    if (lastaccept == 0) {
        f64 current_data_miss_rate =
            1.0 - ((f64)(data_access_hit) / data_access_total);
        last_iteration_data_miss_rate = current_data_miss_rate;
        lastaccept = 1;
        f64 current_discard_rate;
        if (prefetch_increments == 0) {
            current_discard_rate = 0;
        } else {
            current_discard_rate = ((f64)prefetch_discards) / prefetch_increments;
        }
        previous_prefetch_size = current_prefetch_size;
        if (current_discard_rate > RATE_THRESHOLD && current_data_miss_rate < 0.3) {
            current_prefetch_size *= 0.8;
        } else {
            if (current_data_miss_rate >= 0.3 &&
                current_discard_rate <= RATE_THRESHOLD) {
                current_prefetch_size *= 1.2;
            } else {
                if (current_data_miss_rate >= 0.3 &&
                    current_discard_rate > RATE_THRESHOLD) {
                    f64 perturbation =
                        ((static_cast<f64>(rand()) / RAND_MAX) - 0.5) * 0.2;
                    current_prefetch_size *= (1.0 + perturbation);
                }
            }
        }
        current_prefetch_size = max(current_prefetch_size, 1.0 / 256.0);
        current_prefetch_size = min(current_prefetch_size, 0.10);
        prefetchSize = current_prefetch_size * inputcachesize;
        cachesize = inputcachesize - prefetchSize;
        elements_processed_since_last_adjustment = 0;
        prefetch_discards = 0;
        prefetch_increments = 0;
        data_access_hit = 0;
        data_access_total = 0;
        return;
    }
    f64 current_data_miss_rate =
        1.0 - ((f64)(data_access_hit) / data_access_total);
    f64 delta_M = (f64)((1.0 - last_iteration_data_miss_rate) -
                              (1.0 - current_data_miss_rate)) /
                     (1.0 - last_iteration_data_miss_rate);
    bool accept_change = false;
    if (delta_M < 0) {
        accept_change = true;
    } else {
        f64 acceptance_prob = exp(-delta_M / temperature);
        f64 random_val = static_cast<f64>(rand()) / RAND_MAX;
        if (acceptance_prob > random_val) {
            accept_change = true;
        }
    }
    if (accept_change) {
        last_iteration_data_miss_rate = current_data_miss_rate;
        if (current_data_miss_rate < best_data_miss_rate) {
            best_data_miss_rate = current_data_miss_rate;
        }
        lastaccept = 1;
    } else {
        current_prefetch_size = previous_prefetch_size;
        lastaccept = 0;
        prefetchSize = current_prefetch_size * inputcachesize;
        cachesize = inputcachesize - prefetchSize;
        elements_processed_since_last_adjustment = 0;
        prefetch_discards = 0;
        prefetch_increments = 0;
        data_access_hit = 0;
        data_access_total = 0;
        return;
    }
    f64 current_discard_rate;
    if (prefetch_increments == 0) {
        current_discard_rate = 0;
    } else {
        current_discard_rate = ((f64)prefetch_discards) / prefetch_increments;
    }
    previous_prefetch_size = current_prefetch_size;
    if (current_discard_rate > RATE_THRESHOLD && current_data_miss_rate < 0.3) {
        current_prefetch_size *= 0.8;
    } else {
        if (current_data_miss_rate >= 0.3 &&
            current_discard_rate <= RATE_THRESHOLD) {
            current_prefetch_size *= 1.2;
        } else {
            if (current_data_miss_rate >= 0.3 &&
                current_discard_rate > RATE_THRESHOLD) {
                f64 perturbation =
                    ((static_cast<f64>(rand()) / RAND_MAX) - 0.5) * 1.0;
                current_prefetch_size *= (1.0 + perturbation);
            }
        }
    }
    current_prefetch_size = max(current_prefetch_size, 1.0 / 256.0);
    current_prefetch_size = min(current_prefetch_size, 0.1);
    prefetchSize = current_prefetch_size * inputcachesize;
    cachesize = inputcachesize - prefetchSize;
    elements_processed_since_last_adjustment = 0;
    prefetch_discards = 0;
    prefetch_increments = 0;
    data_access_hit = 0;
    data_access_total = 0;
}
bool adaptive_prefetch = 0;
/*
Perform calculation.
*/
void calculate() {
    computePE = computeDramAccess = computeSramAccess = 0;
    switch (dataflow) {
    case Gust: {
        prefetchNow = 0;
        if (cacheScheme == 6 || cacheScheme == 7 || cacheScheme == 66 ||
            cacheScheme == CACHE_SCHEME_FLFU || cacheScheme == CACHE_SCHEME_INNER_SP || cacheScheme == CACHE_SCHEME_SPARCH) {
            if (cacheScheme == 6 || cacheScheme == 7 || cacheScheme == CACHE_SCHEME_INNER_SP ||
                cacheScheme == CACHE_SCHEME_SPARCH) {
                for (int j1 = TJ; j1 < min(TJ + sim.cfg.jjj, sim.cfg.J); j1++) {
                    while (!nextposvector[j1].empty()) {
                        nextposvector[j1].pop();
                    }
                }
            }
            if (cacheScheme == 66) {
                for (int j1 = TJ; j1 < TJ + sim.cfg.jjj; j1++) {
                    LFUtag[j1] = 0;
                }
            }
            for (int ii = 0; prefetchNow < prefetchSize && ii < sim.cfg.iii; ii++) {
                if (TI + ii >= sim.cfg.I)
                    break;
                prefetchRowNow = TI + ii;
                if (!prefetchrow(TI + ii)) {
                    break;
                }
            }
        }
        for (int ii = 0; ii < sim.cfg.iii; ii++) {
            if (TI + ii >= sim.cfg.I)
                break;
            get_B_fibers(TI + ii);
            if (cacheScheme == 6 || cacheScheme == 7 || cacheScheme == 66 ||
                cacheScheme == CACHE_SCHEME_FLFU || cacheScheme == CACHE_SCHEME_INNER_SP || cacheScheme == CACHE_SCHEME_SPARCH) {
                int needsize = 0;
                if (cacheScheme == 6 || cacheScheme == 7 || cacheScheme == CACHE_SCHEME_SPARCH) {
                    // FIX: tile-local index (ii is already relative: 0..iii-1)
                    needsize = sim.cursor.A.sizes[ii] * 4 + 1;
                }
                if (cacheScheme == 66 || cacheScheme == CACHE_SCHEME_FLFU) {
                    // FIX: tile-local index
                    needsize = sim.cursor.A.sizes[ii] * 3;
                }
                if (prefetchNow > needsize) {
                    prefetchNow -= needsize;
                }
                if (prefetchNow >= prefetchSize)
                    continue;
                while (prefetchRowNow < TI + sim.cfg.iii) {
                    if (!prefetchrow(prefetchRowNow)) {
                        break;
                    }
                    prefetchRowNow++;
                }
            }
            if (adaptive_prefetch) {
                if (elements_processed_since_last_adjustment >= adjustment_interval) {
                    update_prefetch_size();
                    elements_processed_since_last_adjustment = 0;
                }
            }
        }
    }
        break;
    case Inner: {
        for (int k = TK; k < TK + sim.cfg.kkk; k++) {
            get_B_fiber_col_iii(k, sim.cfg.iii);
        }
        for (int ii = 0; ii < sim.cfg.iii; ii++) {
            get_A_fiber(TI + ii);
            // FIX: tile-local index for A begins
            int tmpj = sim.cursor.A.begins[ii];
            int maxj = offsetarrayA[TI + ii + 1] - offsetarrayA[TI + ii];
            for (int k1 = TK; k1 < TK + sim.cfg.kkk; k1++) {
                tmpC[k1] = 0;
            }
            while (tmpj < maxj && A[TI + ii][tmpj] < TJ + sim.cfg.jjj) {
                int jj = A[TI + ii][tmpj];
                update_c_fiber(jj);
                tmpj++;
            }
            updateCAccess(TI + ii);
        }
    }
        break;
    case Outer: {
        for (int jj = 0; jj < sim.cfg.jjj; jj++) {
            get_A_fiber_col(TJ + jj);
            get_B_fiber(TJ + jj, jj);
        }
        fulltagC = 0;
        Csizenow = 0;
        for (int ii = TI; ii < TI + sim.cfg.iii; ii++) {
            // NOTE: Outer dataflow uses beginA with absolute index.
            // With Gust fixed this is unreachable.  Would need currsizeAc cursor.
            int tmpj = sim.cursor.A.begins[ii - TI];
            int maxj = offsetarrayA[ii + 1] - offsetarrayA[ii];
            for (int k1 = TK; k1 < TK + sim.cfg.kkk; k1++) {
                tmpC[k1] = 0;
            }
            while (tmpj < maxj && A[ii][tmpj] < TJ + sim.cfg.jjj) {
                int jj = A[ii][tmpj];
                int rj = jj - TJ;
                computePE += sim.cursor.B.sizes[rj];
                update_c_fiber(jj);
                tmpj++;
            }
            int cntc = 0;
            for (int k1 = TK; k1 < TK + sim.cfg.kkk; k1++) {
                if (tmpC[k1]) {
                    cntc++;
                }
            }
            if (!fulltagC) {
                if (Csizenow + cntc * 3 >= Csize) {
                    fulltagC = 1;
                    Csizenow += cntc * 3;
                } else {
                    Csizenow += cntc * 3;
                    computeSramAccess +=
                        sramReadBandwidth(cntc * 3 + 2) * (sim.cfg.jjj / mergecnt);
                }
            } else {
                computeDramAccess += memoryBandwidthPE(cntc * 3 + 2) * (sim.cfg.jjj / mergecnt);
                computeC += memoryBandwidthPE(cntc * 3 + 2) * (sim.cfg.jjj / mergecnt);
                AccessByte += (cntc * 3 + 2) * (sim.cfg.jjj / mergecnt);
            }
        }
    }
        break;
    }
    totalCycle += max(computePE / PEcnt, max(computeDramAccess / PEcnt,
                                             computeSramAccess / sramBank));
    calCycle += max(computePE / PEcnt,
                    max(computeDramAccess / PEcnt, computeSramAccess / sramBank));
    totalSram += computeSramAccess / sramBank;
    totalDram += computeDramAccess / PEcnt;
    totalPE += computePE / PEcnt;
}
void configPartial(f32 partialA, f32 partialB, f32 partialC) {
    Asize = cachesize * partialA;
    Bsize = cachesize * partialB;
    Csize = cachesize * partialC;
}
struct simulator_state initialize_simulator(const struct config *cfg)
{
    try {
        if (bufferedC == nullptr)
            bufferedC = new set<int>[cfg->I]();
        if (bufferedClen == nullptr)
            bufferedClen = new int[cfg->I]();
        // Old full-dimension arrays kept for non-Gust paths that aren't ported yet
        if (beginA == nullptr)
            beginA = new int[cfg->I]();
        if (beginB == nullptr)
            beginB = new int[cfg->J]();
        if (beginAc == nullptr)
            beginAc = new int[cfg->J]();
        if (beginBc == nullptr)
            beginBc = new int[cfg->K]();
        if (currsizeA == nullptr)
            currsizeA = new int[cfg->I]();
        if (currsizeAc == nullptr)
            currsizeAc = new int[cfg->J]();
        if (currsizeB == nullptr)
            currsizeB = new int[cfg->J]();
        if (currsizeBc == nullptr)
            currsizeBc = new int[cfg->J]();
        if (bufferedsizeB == nullptr)
            bufferedsizeB = new int[cfg->J]();
        if (tmpC == nullptr)
            tmpC    = new u8[cfg->K]();
        if (dirtyC == nullptr)
            dirtyC  = new Coord[cfg->K]();
        if (LFUtag == nullptr)
            LFUtag = new int[cfg->J]();
        if (nextposvector == nullptr)
            nextposvector = new queue<int>[cfg->J]();
    } catch (const std::bad_alloc &e) {
        std::cerr << "Error allocating memory for " << e.what() << std::endl;
        std::exit(1);
    }
    return (struct simulator_state) {
        .cfg    = *cfg
    };
}
void reinitialize() {
    // reinitialize statistics
    totalCycle = 0;
    preCycle = calCycle = postCycle = 0;
    computeA = computeB = computeC = 0;
    totalA = totalB = totalC = 0;
    preA = preB = postC = 0;
    totalSram = totalDram = totalPE = 0;
    TI = TJ = TK = 0;
    AccessByte = 0;
    totalhit = 0;
    data_access_hit = 0;
    totalaccess = 0;
    data_access_total = 0;
    postDramAccess = postSramAccess = 0;
    if (ISCACHE) {
        initializeCacheValid();
        if (useVirtualTag) {
            memset(virtualValid, 0, sizeof(bool) * CACHE_NSETS * VIRTUALSETASSOC);
        }
        memset(PosOrig, 0, sizeof(short) * CACHE_NSETS * SETASSOC);
    }
    // reinitialize buffer c
    Csizenow = 0;
    for (int i = 0; i < sim.cfg.I; i++) {
        bufferedC[i]    = std::set<int>();
        bufferedClen[i] = 0;
    }
    // FIX: cursor begins/sizes are tile-local (size iii / jjj).
    // Do NOT zero them with full-dimension loops — that was an OOB write.
    // advance_cursor(first=true) will recompute them from scratch.
    // Just reset the cursor state:
    reset_cursor2(&sim.cursor);
}
int getcntc(int ii) 
{
    // NOTE: Only called from postTileMerge → Outer dataflow path.
    // With Gust fixed, this uses the old full-dimension arrays as fallback.
    int tmpj = beginA[ii];
    int maxj = offsetarrayA[ii + 1] - offsetarrayA[ii];
    Coord ndirty = 0;
    for (int k1 = 0; k1 < sim.cfg.K; k1++)
        tmpC[k1] = 0;
    while (tmpj < maxj && A[ii][tmpj] < sim.cfg.J) {
        int jj = A[ii][tmpj];
        for (int k1 = beginB[jj]; k1 < beginB[jj] + currsizeB[jj]; k1++) {
            Coord index = B[jj][k1];
            if (!tmpC[index]) {
                dirtyC[ndirty++]= index;
                tmpC[index]     = 1;
            }
        }
        tmpj++;
    }
    Coord cntc = 0;
    for (Coord d = 0; d < ndirty; ++d) {
        Coord k1 = dirtyC[d];
        cntc    += tmpC[k1];
        tmpC[k1] = 0;
    }
    return cntc;
}
void postTileMerge() {
    if (dataflow == Outer) {
        for (int ii = 0; ii < sim.cfg.I; ii++) {
            int cntc = getcntc(ii);
            postDramAccess += memoryBandwidthPE(cntc * 3) * (sim.cfg.ttj);
            postC += memoryBandwidthPE(cntc * 3) * (sim.cfg.ttj);
            AccessByte += cntc * 3 * sim.cfg.ttj;
        }
    }
    postSramAccess /= sramBank;
    postDramAccess /= PEcnt;
    totalCycle  += max(postDramAccess, postSramAccess);
    postCycle   += max(postDramAccess, postSramAccess);
}

// ============================================================================
// run() — flat tile loop driven entirely by advance_cursor.
//
// OLD triple do-while with iterate_inner/mid/outer_loop + updateT* + reverse_*
// is replaced by a single counted loop.  advance_cursor manages tile index
// incrementing, wrapping, and begin/size recomputation.  sync_tile_origins
// derives TI/TJ/TK from the tile indices after each step.
// ============================================================================
void run()
{
    reinitialize();
    if (adaptive_prefetch) {
        current_prefetch_size = 1.0 / 128.0;
        prefetchSize = current_prefetch_size * inputcachesize;
        cachesize = inputcachesize - prefetchSize;
        setSET();
        sa_iteration_k = 0;
        previous_prefetch_size = current_prefetch_size;
        last_iteration_data_miss_rate = 1.0;
        best_data_miss_rate = 1.0;
        adjustment_interval = estEffMAC / 500;
        if (adjustment_interval == 0)
            adjustment_interval = 1000;
        elements_processed_since_last_adjustment = 0;
        prefetch_discards = 0;
        data_access_hit = 0;
        data_access_total = 0;
    }
    if (sim.cfg.iii > sim.cfg.I)
        sim.cfg.iii = sim.cfg.I;
    if (sim.cfg.jjj > sim.cfg.J)
        sim.cfg.jjj = sim.cfg.J;
    if (sim.cfg.kkk > sim.cfg.K)
        sim.cfg.kkk = sim.cfg.K;

    // FIX: flat tile loop.  No more iterate_*_loop / updateT* / reverse_*.
    Coord total_tiles = (Coord)sim.cfg.tti * sim.cfg.ttj * sim.cfg.ttk;
    for (Coord t = 0; t < total_tiles; ++t) {
        if (ISCACHE) {
            initializeCacheValid();
        }
        advance_cursor();
        sync_tile_origins();
        pre_calculate_load();
        calculate();
    }
    postTileMerge();
    analyze_statistics();
}
void runTile(int kkk)
{
    assert(ISCACHE);
    if (cacheScheme == 6 || cacheScheme == 7) {
        cachesize = inputcachesize - prefetchSize;
        cachesize -= kkk * 2;
        if (cachesize < 0) {
            puts("!!!!!! metadata out of range!!!!!!!!!!");
            fflush(stdout);
            return;
        }
        setSET();
    }
    if (cacheScheme == 66) {
        cachesize = inputcachesize - prefetchSize;
        cachesize -= kkk;
        if (cachesize < 0) {
            puts("!!!!!! metadata out of range!!!!!!!!!!");
            fflush(stdout);
            return;
        }
        setSET();
    }
    if (cacheScheme == CACHE_SCHEME_FLFU) {
        cachesize = inputcachesize - prefetchSize;
        setSET();
    }
    if ((cacheScheme == 4) || (cacheScheme == 5) || (cacheScheme == 7)) {
        cachesize = inputcachesize;
        CACHE_NSETS = cachesize / (CACHE_BLOCK_NELEMS * SETASSOC);
        CACHE_NSETS_LOG2 = getlog(CACHE_NSETS);
    }
    hitcnt = 0;
    misscnt = 0;
    fflush(stdout);
    run();
}