/*
NOTE: For regression checking, check against branch "main_with_indexing_fix"
(specifically, commit: e8f2d9bac4bbf0ad0619930dbdd54d089d9cb351).

This branch fixes a "read one past end of array allocation" bug for offsetarrayA/Ac/B/Bc.
*/
#include "arena.h"
#include <stdlib.h>
#include "cache.h"
#include "estimation.h"
#include "headers.h"
#include "json.hpp"
#include "simulator.h"
#include <cstdlib>
#include <fstream>

struct Arena *global_persist;
struct Arena *global_temp;
struct simulator_state sim;
struct cache cache;

using json = nlohmann::json;

struct matrix {
    b16 transpose;

    u64 nrows;
    u64 ncols;
    u64 nzM;    // number of non-zero elements

    u64 *M_backing;
    u64 *Mc_backing;

    u64 **M,
        **Mc;
    u64
        *offsetarrayM,
        *offsetarrayMc;
};

i32 cmp_u64(const void *a, const void *b)
{
    u64 x = *(const u64 *)a;
    u64 y = *(const u64 *)b;
    return (x > y) - (x < y);
}

void parse_matrix(FILE *f, struct matrix *x)
{
    struct Arena_Mark mark = arena_snap(global_temp);

    b16 transpose = x->transpose;

    // Lines limited to 1024 characters by spec
    #define BUFFER_NBYTES (u64)1024
    char readbuffer[BUFFER_NBYTES];

    u64 *raw_rows   = (u64 *)arena_push(global_temp, x->nzM*sizeof(*raw_rows), __alignof__(*raw_rows), 0);
    u64 *raw_cols   = (u64 *)arena_push(global_temp, x->nzM*sizeof(*raw_cols), __alignof__(*raw_cols), 0);
    u64 *row_lens   = (u64 *)arena_push(global_temp, x->nrows*sizeof(*row_lens), __alignof__(*row_lens), 1); // must be zeroed
    u64 *col_lens   = (u64 *)arena_push(global_temp, x->ncols*sizeof(*col_lens), __alignof__(*row_lens), 1); // must be zeroed
    u64 *offsets;

    x->M_backing    = (u64 *)arena_push(global_persist, x->nzM*sizeof(*x->M_backing), __alignof__(*x->M_backing), 0);
    x->Mc_backing   = (u64 *)arena_push(global_persist, x->nzM*sizeof(*x->Mc_backing), __alignof__(*x->Mc_backing), 0);
    x->M            = (u64 **)arena_push(global_persist, (x->nrows+1)*sizeof(*x->M), __alignof__(*x->M), 0);
    x->Mc           = (u64 **)arena_push(global_persist, (x->ncols+1)*sizeof(*x->Mc), __alignof__(*x->Mc), 0);
    x->offsetarrayM = (u64 *)arena_push(global_persist, (x->nrows+1)*sizeof(x->offsetarrayM[0]), __alignof__(x->offsetarrayM[0]), 0);
    x->offsetarrayMc= (u64 *)arena_push(global_persist, (x->ncols+1)*sizeof(x->offsetarrayMc[0]), __alignof__(x->offsetarrayMc[0]), 0);

    for (u64 i = 0; i < x->nzM; ++i) {
        /* NOTE(ejs): Each mtx row is a matrix entry.
        xx, yy [zz, lala]
            - xx  = row index
            - yy  = col index
            - zz  = real component of value? (ignored)
            - lala= imag. component of value? (ignored)
        */
        assert(fgets(readbuffer, BUFFER_NBYTES, f));
        u64 xx, yy;
        // f64 zz, lala;
        // assert(sscanf(readbuffer, "%llu %llu %lf %lf", &xx, &yy, &zz, &lala) >= 2);
        assert(sscanf(readbuffer, "%llu %llu", &xx, &yy) == 2);

        // WARNING(ejs): mtx indices are stored 1-based; it converts 0-based representation in code
        u64 row_index = xx-1;
        u64 col_index = yy-1;
        if (transpose)
            swap(row_index, col_index);
        raw_rows[i] = row_index;
        raw_cols[i] = col_index;
        ++row_lens[row_index];
        ++col_lens[col_index];
    }

    // csr
    // prefix sum
    offsets = x->offsetarrayM;
    offsets[0] = 0;
    for (u64 i = 0; i < x->nrows; ++i)
        offsets[i+1] = offsets[i] + row_lens[i];

    // scatter cols into row-major order (reuse row_lens as cursor)
    memset(row_lens, 0, x->nrows * sizeof(*row_lens));
    for (u64 i = 0; i < x->nzM; ++i) {
        u64 r = raw_rows[i];
        u64 pos = offsets[r] + row_lens[r]++;
        x->M_backing[pos] = raw_cols[i];
    }

    // sort within each row by column
    for (u64 r = 0; r < x->nrows; ++r) {
        u64 *base = x->M_backing + offsets[r];
        qsort(base, row_lens[r], sizeof(*base), cmp_u64);
    }
    
    // build csr pointers
    for (u64 i = 0; i < x->nrows; ++i)
        x->M[i] = x->M_backing + offsets[i];
    x->M[x->nrows] = x->M_backing + x->offsetarrayM[x->nrows];

    // csc
    offsets = x->offsetarrayMc;
    offsets[0] = 0;
    for (u64 i = 0; i < x->ncols; ++i)
        offsets[i+1] = offsets[i] + col_lens[i];

    memset(col_lens, 0, x->ncols * sizeof(*col_lens));
    for (u64 i = 0; i < x->nzM; ++i) {
        u64 c = raw_cols[i];
        u64 pos = offsets[c] + col_lens[c]++;
        x->Mc_backing[pos] = raw_rows[i];
    }
    for (u64 c = 0; c < x->ncols; ++c) {
        u64 *base = x->Mc_backing + offsets[c];
        qsort(base, col_lens[c], sizeof(*base), cmp_u64);
    }
    for (u64 i = 0; i < x->ncols; ++i)
        x->Mc[i] = x->Mc_backing + offsets[i];
    x->Mc[x->ncols] = x->Mc_backing + x->offsetarrayMc[x->ncols];

    arena_rewind(mark);
}

int main(int argc, char *argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <matrix1 name> <matrix2 name> <config filepath>\n"
                  << "config_filepath is a fully qualified path to the .json config for the run (likely config/config.json).\n"
                  << std::endl;
        return 1;
    }

    global_persist  = arena_alloc(16*GB, MB);
    global_temp     = arena_alloc(16*GB, MB);
    cache.backing   = arena_alloc(16*GB, MB);

    std::string matrix1_name    = argv[1];
    std::string matrix2_name    = argv[2];
    std::string config_filepath = argv[3];

    std::ifstream file(config_filepath);
    if (!file.is_open()) {
        std::cerr << "Error opening config file." << std::endl;
        return 1;
    }

    json config;
    file >> config;

    int transpose           = config["transpose"].get<int>();
    f32 tmpsram             = config["cachesize"].get<f32>();
    f32 tmpbandw            = config["memorybandwidth"].get<f32>();
    int baselinetest        = config["baselinetest"].get<int>();
    bool condensedOP        = config["condensedOP"].get<bool>();
    std::string tile_dir    = config["tileDir"].get<std::string>();
    std::string output_dir  = config["outputDir"].get<std::string>();

    cachesize = tmpsram * 262144 * 0.9;
    inputcachesize = cachesize;
    HBMbandwidth = (tmpbandw / 4.0) * 0.6;
    int tmpPE = config["PEcnt"].get<int>();
    PEcnt = tmpPE;
    mergecnt = tmpPE;
    HBMbandwidthperPE = HBMbandwidth / PEcnt;
    int tmpbank = config["srambank"].get<int>();
    sramBank = tmpbank;

    std::string matrix1_filepath = "data/" + matrix1_name + ".mtx";
    std::string matrix2_filepath = "data/" + matrix2_name + ".mtx";
    std::string output_filepath = output_dir
        +   (1 ? "C" : "_") // +   (ISCACHE ? "C" : "_")
        +   printDataFlow[dataflow]
        +   (baselinetest ? "Base_" : "570Cache_")
        +   std::to_string(tmpsram)
        +   "MB_" + std::to_string(tmpbandw)
        +   "GBs_" + std::to_string(tmpPE)
        +   "PEs_" + std::to_string(tmpbank) + "sbanks_"
        +   "_" + matrix1_name + "_" + matrix2_name + "_"
        +   printFormat[format] + "_" + (transpose ? "1" : "0") + ".txt";

    FILE *matrix1_file  = fopen(matrix1_filepath.c_str(), "r");
    FILE *matrix2_file  = fopen(matrix2_filepath.c_str(), "r");
    assert(matrix1_file);
    assert(matrix2_file);
    assert(freopen(output_filepath.c_str(), "w", stdout));

    struct matrix matA = {0};
    struct matrix matB = {0};
    {
        #define TEMP_BUFFER_NBYTES 1024
        char buf[TEMP_BUFFER_NBYTES];

        // read and ignore annotation '%' lines
        while (fgets(buf, TEMP_BUFFER_NBYTES, matrix1_file)) {
            if (buf[0] != '%')
                break;
        }
        sscanf(buf, "%llu%llu%llu", &matA.nrows, &matA.ncols, &matA.nzM);

        // read and ignore annotation '%' lines
        while (fgets(buf, TEMP_BUFFER_NBYTES, matrix2_file)) {
            if (buf[0] != '%')
                break;
        }
        sscanf(buf, "%llu%llu%llu", &matB.nrows, &matB.ncols, &matB.nzM);
    }

    matA.transpose = (b16)transpose;
    matB.transpose = (b16)((matB.nrows == matB.ncols) ? transpose : !transpose);
    if (matA.transpose)
        swap(matA.ncols, matA.nrows);
    if (matB.transpose)
        swap(matB.ncols, matB.nrows);
    assert(matA.ncols == matB.nrows);

    FILE *tile_file = fopen((tile_dir + matrix1_name).c_str(), "r");
    assert(tile_file);
    u64 t_i, t_j, t_k;
    {
        #define TEMP_BUFFER_NBYTES 1024
        char buf[TEMP_BUFFER_NBYTES];

        assert(fgets(buf, TEMP_BUFFER_NBYTES, tile_file));
        assert(sscanf(buf, "%llu%llu%llu", &t_i, &t_j, &t_k) == 3);
    }

    const struct config cfg = {
        .dataflow   = Gust,
        .interorder = IJK,
        .format     = RR,

        .I          = matA.nrows,
        .J          = matA.ncols,
        .K          = matB.ncols,
        .iii        = t_i,
        .jjj        = t_j,
        .kkk        = t_k,
        .tti        = div_rup(matA.nrows, t_i),
        .ttj        = div_rup(matA.ncols, t_j),
        .ttk        = div_rup(matB.ncols, t_k),
    };
    sim = initialize_simulator(&cfg);

    parse_matrix(matrix1_file, &matA);
    A   = matA.M;
    Ac  = matA.Mc;
    offsetarrayA = matA.offsetarrayM;
    offsetarrayAc= matA.offsetarrayMc;

    if (condensedOP) {
        assert(0); // this path is not maintained

        // // memory management for sparchA
        // sparchA = new std::vector<int>[sim.cfg.J]();
        // if (sparchA == nullptr) {
        //     if (sparchA != nullptr)
        //         delete[] sparchA;
        //     std::cerr << "Error allocating memory for sparchA" << std::endl;
        //     std::exit(1);
        // }

        // // if use the condensed OP dataflow, need to preprocess the A matrix into
        // // the condensed format first. first put the data into sparchA[], then put
        // // it back to A[], and call gust dataflow
        // for (int j = 0; j < sim.cfg.J; j++) {
        //     for (int i = 0; i < sim.cfg.I; i++) {
        //         if (static_cast<int>(A[i].size()) > j) {
        //             sparchA[j].push_back(A[i][j]);
        //         }
        //     }
        // }

        // for (int j = 0; j < sim.cfg.J; j++) {
        //     A[j].clear();
        //     for (int i = 0; i < static_cast<int>(sparchA[j].size()); i++) {
        //         A[j].push_back(sparchA[j][i]);
        //     }
        // }

        // delete[] sparchA;
    }

    long long totalempty = 0;

    long long totalincache = 0;

    long long totaltagmatch48 = 0;
    long long totaltagmatch16 = 0;

    for (int i = 1; i < sim.cfg.I; i++) {
        u64 size_im1 = offsetarrayA[i] - offsetarrayA[i-1];
        if (size_im1 < 48) {
            totalempty += (48 - size_im1);
        }
        totalincache    += min(48, size_im1);
        totaltagmatch48 += div_rup(size_im1, 48);
        totaltagmatch16 += div_rup(size_im1, 16);
    }

    /////////////////////////// input B /////////////////////////////////
    parse_matrix(matrix2_file, &matB);
    B   = matB.M;
    Bc  = matB.Mc;
    offsetarrayB = matB.offsetarrayM;
    offsetarrayBc= matB.offsetarrayMc;

    /******************Config************************************/

    printf("Matrix A: %llu x %llu, number of non-zeros = %llu\n", matA.nrows, matA.ncols, matA.nzM);
    printf("*** ratio of empty %lf, ratio of not empty %lf\n", totalempty / (sim.cfg.I * 48.0), 1 - (totalempty / (sim.cfg.I * 48.0)));
    printf("*** ratio of in cache %lf\n", totalincache / ((f64)matA.nzM));
    printf("** ratio tag access 48 %lf\n", sim.cfg.I / ((f64)sim.cfg.I + totaltagmatch48));
    printf("** ratio tag access 16 %lf\n", sim.cfg.I / ((f64)sim.cfg.I + totaltagmatch16));
    printf("Matrix B: %llu x %llu, number of non-zeros = %llu\n", matB.nrows, matB.ncols, matB.nzM);
    printf("transpose: %d\n", matB.transpose);
    printf("I = %llu, K = %llu, J = %llu\n", sim.cfg.I, sim.cfg.K, sim.cfg.J);
    /************************************************************/

    getParameter(); // sets estEffMAC

    configPartial(0.05, 0.5, 0.45);

    /////////////// Baseline configurations

    if (baselinetest) {
        // EWH
        // Incorporate SeaCache into baseline
        puts("***************** SeaCache *******************");
        printf("nnzB:%llu  K:%llu  J/TJ:%llu  nzlB:%llu\n",
            matB.nzM, sim.cfg.K, (sim.cfg.J + sim.cfg.jjj - 1) / sim.cfg.jjj,
            matB.nzM / (sim.cfg.K * ((sim.cfg.J + sim.cfg.jjj - 1) / sim.cfg.jjj))
        );

        ISCACHE = 1;
        cachesize               = 262144;
        CACHE_BLOCK_NELEMS      = 16;
        CACHE_BLOCK_NELEMS_LOG2 = getlog(CACHE_BLOCK_NELEMS);
        setSET();

        cache.cfg = {
            .block_nelems       = 1,
            .block_nelems_log2  = 1,
            .scheme             = CACHE_SCHEME_FLFU,
        };

        adaptive_prefetch = 1;
        useVirtualTag = 1;
        cacheScheme = CACHE_SCHEME_FLFU;
        cachesize = inputcachesize;

        runTile(sim.cfg.kkk);

        adaptive_prefetch = 0;
        useVirtualTag = 0;

        adaptive_prefetch = 0;

        ////////////  InnserSP
        // static FLRU + 16 words scheme0
        puts("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   test InnerSP   "
             "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        prefetchSize = inputcachesize / 6;
        cacheScheme = CACHE_SCHEME_INNER_SP;
        cachesize = inputcachesize;
        CACHE_BLOCK_NELEMS = 16;
        CACHE_BLOCK_NELEMS_LOG2 = 4;
        setSET();
        runTile(sim.cfg.kkk);

        ////////////  Sparch
        // dynamic FLRU + 128KB prefetch size + 144 words scheme0
        puts("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   test Sparch   "
             "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        ISCACHE = 1;
        cacheScheme = CACHE_SCHEME_SPARCH;
        prefetchSize = inputcachesize / 6;
        cachesize = inputcachesize - prefetchSize;
        CACHE_BLOCK_NELEMS = 144;
        CACHE_BLOCK_NELEMS_LOG2 = 8;
        setSET();
        // calculate metadata overhead.
        // if metadata overflow, choose smaller tile
        int newkkk = sim.cfg.kkk;
        int newttk = sim.cfg.ttk;
        // if can keep, just use current kkk
        if (cachesize > sim.cfg.kkk * 2) {
            cachesize -= sim.cfg.kkk * 2;
        } else {
            // if can't keep, use smaller kkk
            // (make kkk*2 to be half cachesize)
            newkkk = cachesize / 4;
            newttk = (sim.cfg.K + sim.cfg.kkk - 1) / sim.cfg.kkk;
            cachesize -= sim.cfg.kkk * 2;
        }
        runTile(newkkk);
        // return to the default setting
        CACHE_BLOCK_NELEMS = 16;
        CACHE_BLOCK_NELEMS_LOG2 = 4;
        cachesize = inputcachesize;
        setSET();

        ////////////  X-cache
        puts("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   test X-cache   "
             "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        // LRU + 4 words scheme0
        // just same as using scheme0 with cacheline = 4
        ISCACHE = 1;
        cacheScheme = CACHE_SCHEME_BASE;
        cachesize = inputcachesize;
        CACHE_BLOCK_NELEMS = 4;
        CACHE_BLOCK_NELEMS_LOG2 = 2;
        setSET();
        runTile(sim.cfg.kkk);

        // return to the default setting
        CACHE_BLOCK_NELEMS = 16;
        CACHE_BLOCK_NELEMS_LOG2 = 4;
        setSET();

        puts("!!!!!!!!!!!!!!!!!!!!  Scratchpad   !!!!!!!!!!!!!!!!!!!!!!!");
        ISCACHE = 0;

        configPartial(0.05, 0.9, 0.05);

        reinitialize();

        run();
    }

    if (!baselinetest) {
        puts("\n!!!!!!!!!!!!!!!!!!!! EECS570 !!!!!!!!!!!!!!!!!!!!");

        /*****************************************
        adaptive_prefetch = 1;
        useVirtualTag = 2;
        cacheScheme;
        cachesize = inputcachesize;

        runTile(kkk);
        adaptive_prefetch = 0;
        useVirtualTag = 0;
        *****************************************/
    }

    bool ablationtest = 0;
    if (ablationtest) {

        adaptive_prefetch = 0;

        /////////////// ablation test

        puts("\n!!!!!!!!!!!!!!!!!!!!!!!!!! scheme0 (base)   "
             "!!!!!!!!!!!!!!!!!!!!!!!!");
        puts("CacheScheme 0");
        ISCACHE = 1;
        cacheScheme = CACHE_SCHEME_BASE;
        cachesize = inputcachesize;
        setSET();
        runTile(sim.cfg.kkk);

        puts("\n!!!!!!!!!!!!!!!!!!!!!!!!!! scheme1 (mapping)   "
             "!!!!!!!!!!!!!!!!!!!!!!!!");

        puts("CacheScheme 1");
        ISCACHE = 1;
        cacheScheme = CACHE_SCHEME_MAPPING;
        cachesize = inputcachesize;
        setSET();
        runTile(sim.cfg.kkk);

        puts("\n!!!!!!!!!!!!!!!!!!!!!!!!!! scheme88 without virtue   "
             "!!!!!!!!!!!!!!!!!!!!!!!!");

        useVirtualTag = 0;
        cacheScheme = CACHE_SCHEME_FLFU;
        cachesize = inputcachesize;
        prefetchSize = cachesize / 6;
        runTile(sim.cfg.kkk);

        puts("\n!!!!!!!!!!!!!!!!!!!!!!!!!! scheme88 with virtue   "
             "!!!!!!!!!!!!!!!!!!!!!!!!");

        puts("CacheScheme 88 practical FLFU  with virtual tag 1/6");
        useVirtualTag = 1;
        cacheScheme = CACHE_SCHEME_FLFU;
        cachesize = inputcachesize;
        prefetchSize = cachesize / 6;
        runTile(sim.cfg.kkk);
        useVirtualTag = 0;

        puts("CacheScheme 88 practical FLFU  with virtual tag 1/16");
        useVirtualTag = 1;
        cacheScheme = CACHE_SCHEME_FLFU;
        cachesize = inputcachesize;
        prefetchSize = cachesize / 16;
        runTile(sim.cfg.kkk);
        useVirtualTag = 0;
    }

    return 0;
}
