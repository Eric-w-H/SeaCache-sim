#include "headers.h"
#include "arena.h"
#include "util.h"

long long estEffMAC;

struct pair_hash {
    template <typename T1, typename T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};

void getParameter()
{
    struct Arena_Mark mark = arena_snap(global_temp);
    b8  *seen;
    u64 *dirty;
    seen = (b8 *)arena_push(global_temp, sim.cfg.K*sizeof(*seen), __alignof__(*seen), 1);
    dirty= (u64*)arena_push(global_temp, sim.cfg.K*sizeof(*dirty), __alignof__(*dirty), 0);

    estEffMAC   = 0;
    u64 estnnzC = 0;

    for (u64 i = 0; i < sim.cfg.I; ++i) {
        u64 lenA = offsetarrayA[i+1] - offsetarrayA[i];
        u64 ndirty = 0;

        for (u64 j = 0; j < lenA; ++j) {
            u64 tmpx = A[i][j];
            u64 lenB = offsetarrayB[tmpx+1] - offsetarrayB[tmpx];
            estEffMAC += lenB;

            for (u64 k = 0; k < lenB; ++k) {
                u64 tmpk = B[tmpx][k];
                if (!seen[tmpk]) {
                    seen[tmpk]      = 1;
                    dirty[ndirty++] = tmpk;
                }
            }
        }

        estnnzC += ndirty;
        for (u64 d = 0; d < ndirty; d++)
            seen[dirty[d]] = 0;
    }

    printf("MAC: %lld   nnzC: %lld\n", estEffMAC, estnnzC);
    arena_rewind(mark);
}