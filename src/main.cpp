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
struct simulator_state sim;

using json = nlohmann::json;

struct matrix {
    const char *mat_name;
    i16 transpose;

    u64 nrows;
    u64 ncols;
    u64 nzM;    // number of non-zero elements

    u64 *M_backing;
    u64 *Mc_backing;

    std::vector<u64>
        *M,
        *Mc; 
    u64
        *offsetarrayM,
        *offsetarrayMc;
};

struct kv {
    u64 key;
    u64 val;
};

i32 compare_kv(const void *p, const void *q)
{
    u64 x = ((struct kv *)p)->key;
    u64 y = ((struct kv *)q)->key;
    return (x > y) - (x < y);
}

void parse_matrix(FILE *f, struct matrix *x)
{
    i16 transpose = x->transpose;

    // Lines limited to 1024 characters by spec
    #define BUFFER_NBYTES (u64)1024
    char readbuffer[BUFFER_NBYTES];

    // void *arena_push(struct Arena *a, uint64_t size, uint64_t align, uint8_t zero)
    // u64 total_nbytes    = sizeof(u64)*x->nzM;
    u64 *raw_rows, *raw_cols, *flat_indices;
    struct kv *kvs;
    raw_rows    = (u64 *)arena_push(global_persist, x->nzM*sizeof(*raw_rows), __alignof__(*raw_rows), 0);
    raw_cols    = (u64 *)arena_push(global_persist, x->nzM*sizeof(*raw_cols), __alignof__(*raw_cols), 0);
    flat_indices= (u64 *)arena_push(global_persist, x->nzM*sizeof(*flat_indices), __alignof__(*flat_indices), 0);
    kvs         = (struct kv *)arena_push(global_persist, x->nzM*sizeof(kvs[0]), __alignof__(kvs[0]), 0);

    x->M_backing    = (u64 *)arena_push(global_persist, x->nzM*sizeof(*x->M_backing), __alignof__(*x->M_backing), 0);
    x->Mc_backing   = (u64 *)arena_push(global_persist, x->nzM*sizeof(*x->Mc_backing), __alignof__(*x->Mc_backing), 0);
    // x->offsetarrayM = (u64 *)arena_push(global_persist, x->nrows*sizeof(x->offsetarrayM[0]), __alignof__(x->offsetarrayM[0]), 0);
    // x->offsetarrayMc= (u64 *)arena_push(global_persist, x->ncols*sizeof(x->offsetarrayMc[0]), __alignof__(x->offsetarrayMc[0]), 0);
    u64 *row_lens   = (u64 *)arena_push(global_persist, x->nrows*sizeof(*row_lens), __alignof__(*row_lens), 1); // must be zeroed
    u64 *col_lens   = (u64 *)arena_push(global_persist, x->ncols*sizeof(*col_lens), __alignof__(*row_lens), 1); // must be zeroed

    try {
        if (x->M == nullptr)
            x->M = new std::vector<u64>[x->nrows]();
        if (x->Mc == nullptr)
            x->Mc = new std::vector<u64>[x->ncols]();
        if (x->offsetarrayM == nullptr)
            x->offsetarrayM = new u64[x->nrows]();
        if (x->offsetarrayMc == nullptr)
            x->offsetarrayMc = new u64[x->ncols]();
    } catch (const std::bad_alloc &e) {
        std::cerr << "Error allocating memory for " << e.what() << std::endl;
        std::exit(1);
    }

    std::string input;

    for (u64 i = 0; i < x->nzM; ++i) {
        assert(fgets(readbuffer, BUFFER_NBYTES, f));
        input = readbuffer;

        std::istringstream iss(input);
        std::vector<std::string> tokens;
        std::string token;

        // Splits the file input on whitespace
        while (iss >> token) {
            tokens.push_back(token);
        }

        int xx, yy;
        // double zz, lala;

        /* NOTE(ejs): Each mtx row is a matrix entry.
        xx, yy [zz, lala]
            - xx  = row index
            - yy  = col index
            - zz  = real component of value? (ignored)
            - lala= imag. component of value? (ignored)
        */
        if (tokens.size() == 2) { // pattern (nonzero values ommitted)

            std::istringstream(tokens[0]) >> xx;
            std::istringstream(tokens[1]) >> yy;
            // std::cout << "values: " << xx << ", " << yy << std::endl;
        } else if (tokens.size() == 3) { // real or integer matrix

            std::istringstream(tokens[0]) >> xx;
            std::istringstream(tokens[1]) >> yy;
            // std::istringstream(tokens[2]) >> zz;
            // std::cout << "values: " << xx << ", " << yy << ", " << zz << std::endl;
        } else if (tokens.size() == 4) { // complex matrix (we only take the real part, unfortunately)
            std::istringstream(tokens[0]) >> xx;
            std::istringstream(tokens[1]) >> yy;
            // std::istringstream(tokens[2]) >> zz;
            // std::istringstream(tokens[3]) >> lala;
            // std::cout << "values: " << xx << ", " << yy << ", " << zz << std::endl;
        } else {

            std::cout << "Format Incorrect! " << std::endl;
            cout << tokens.size() << endl;
            cout << i << endl
                 << input << endl;
            return;
        }

        u64 row_index = xx-1;
        u64 col_index = yy-1;
        if (transpose)
            swap(row_index, col_index);
        raw_rows[i] = row_index;
        raw_cols[i] = col_index;

        // WMRNING(ejs): mtx indices are stored 1-based; it converts 0-based representation in code
        if (transpose) {
            x->Mc[xx - 1].push_back(yy - 1);
            x->M[yy - 1].push_back(xx - 1);

            // x->M_backing[col_major_flat_index]  = row_index; 
            // x->Mc_backing[row_major_flat_index] = col_index;
            // ++row_lens[col_index];
            // ++col_lens[row_index];
        } else {
            x->M[xx - 1].push_back(yy - 1);
            x->Mc[yy - 1].push_back(xx - 1);

            // x->M_backing[row_major_flat_index]  = col_index;
            // x->Mc_backing[col_major_flat_index] = row_index;
            // ++row_lens[row_index];
            // ++col_lens[col_index];
        }
    }

    {
        // row major sort
        for (u64 i = 0; i < x->nzM; ++i) {
            kvs[i] = {
                .key = raw_rows[i]*x->ncols + raw_cols[i],
                .val = i
            };
        }

        qsort(kvs, x->nzM, sizeof(*kvs), compare_kv);

        for (u64 i = 0; i < x->nzM; ++i) {
            u64 rdi = kvs[i].val;
            u64 row = raw_rows[rdi];
            u64 col = raw_cols[rdi];

            ++row_lens[row];
            x->M_backing[i] = col;
        }

        // col major sort
        for (u64 i = 0; i < x->nzM; ++i) {
            kvs[i] = {
                .key = raw_cols[i]*x->nrows + raw_rows[i],
                .val = i
            };
        }

        qsort(kvs, x->nzM, sizeof(*kvs), compare_kv);

        for (u64 i = 0; i < x->nzM; ++i) {
            u64 rdi = kvs[i].val;
            u64 row = raw_rows[rdi];
            u64 col = raw_cols[rdi];

            ++col_lens[col];
            x->Mc_backing[i] = row;
        }
    }

    for (u64 i = 0; i < x->nrows; ++i)
        x->offsetarrayM[i]  = x->offsetarrayM[i-1] + row_lens[i-1];
    for (u64 i = 0; i < x->ncols; ++i)
        x->offsetarrayMc[i] = x->offsetarrayMc[i-1]+ col_lens[i-1];

    for (int i = 0; i < x->nrows; i++) {
        sort(x->M[i].begin(), x->M[i].end());
    }
    for (int j = 0; j < x->ncols; j++) {
        sort(x->Mc[j].begin(), x->Mc[j].end());
    }

}

int main(int argc, char *argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <matrix1 name> <matrix2 name> <config filepath>\n"
                  << "config_filepath is a fully qualified path to the .json config for the run (likely config/config.json).\n"
                  << std::endl;
        return 1;
    }

    global_persist = arena_alloc(16*GB, MB);

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
    float tmpsram           = config["cachesize"].get<float>();
    float tmpbandw          = config["memorybandwidth"].get<float>();
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
    ISCACHE = 1;

    std::string matrix1_filepath = "data/" + matrix1_name + ".mtx";
    std::string matrix2_filepath = "data/" + matrix2_name + ".mtx";
    FILE *matrix1_file = fopen(matrix1_filepath.c_str(), "r");
    FILE *matrix2_file = fopen(matrix2_filepath.c_str(), "r");
    assert(matrix1_file);
    assert(matrix2_file);

    struct matrix matA = {.mat_name="A"};
    struct matrix matB = {.mat_name="B"};
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

    matA.transpose = (i16)transpose;
    matB.transpose = (i16)((matB.nrows == matB.ncols) ? transpose : !transpose);
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


    if (!freopen((output_dir + (ISCACHE ? "C" : "_") + printDataFlow[dataflow] +
                  (baselinetest ? "Base_" : "570Cache_") +
                  std::to_string(tmpsram) + "MB_" + std::to_string(tmpbandw) +
                  "GBs_" + std::to_string(tmpPE) + "PEs_" +
                  std::to_string(tmpbank) + "sbanks_" + "_" + matrix1_name + "_" +
                  matrix2_name + "_" + printFormat[format] + "_" +
                  (transpose ? "1" : "0") + ".txt")
                     .c_str(),
                 "w", stdout)) {
        std::cerr << "Error opening output folder." << std::endl;
        return 1;
    }

    parse_matrix(matrix1_file, &matA);
    A   = matA.M;
    Ac  = matA.Mc;
    offsetarrayA = matA.offsetarrayM;
    offsetarrayAc= matA.offsetarrayMc;

    if (condensedOP) {
        assert(0); // this path is not maintained

        // memory management for sparchA
        sparchA = new std::vector<int>[sim.cfg.J]();
        if (sparchA == nullptr) {
            if (sparchA != nullptr)
                delete[] sparchA;
            std::cerr << "Error allocating memory for sparchA" << std::endl;
            std::exit(1);
        }

        // if use the condensed OP dataflow, need to preprocess the A matrix into
        // the condensed format first. first put the data into sparchA[], then put
        // it back to A[], and call gust dataflow
        for (int j = 0; j < sim.cfg.J; j++) {
            for (int i = 0; i < sim.cfg.I; i++) {
                if (static_cast<int>(A[i].size()) > j) {
                    sparchA[j].push_back(A[i][j]);
                }
            }
        }

        for (int j = 0; j < sim.cfg.J; j++) {
            A[j].clear();
            for (int i = 0; i < static_cast<int>(sparchA[j].size()); i++) {
                A[j].push_back(sparchA[j][i]);
            }
        }

        delete[] sparchA;
    }

    long long totalempty = 0;

    long long totalincache = 0;

    long long totaltagmatch48 = 0;
    long long totaltagmatch16 = 0;

    for (int i = 1; i < sim.cfg.I; i++) {
        if (A[i - 1].size() < 48) {
            totalempty += (48 - A[i - 1].size());
        }
        totalincache += min(48, (int)A[i - 1].size());
        totaltagmatch48 += ((int)A[i - 1].size() + 47) / 48;
        totaltagmatch16 += ((int)A[i - 1].size() + 15) / 16;
    }

    /////////////////////////// input B /////////////////////////////////
    parse_matrix(matrix2_file, &matB);
    B   = matB.M;
    Bc  = matB.Mc;
    offsetarrayB = matB.offsetarrayM;
    offsetarrayBc= matB.offsetarrayMc;

    if (ISCACHE == 1) {
        setSET();
    }

    /******************Config************************************/

    printf("Matrix A: %llu x %llu, number of non-zeros = %llu\n", matA.nrows, matA.ncols, matA.nzM);
    printf("*** ratio of empty %lf, ratio of not empty %lf\n", totalempty / (sim.cfg.I * 48.0), 1 - (totalempty / (sim.cfg.I * 48.0)));
    printf("*** ratio of in cache %lf\n", totalincache / ((double)matA.nzM));
    printf("** ratio tag access 48 %lf\n", sim.cfg.I / ((double)sim.cfg.I + totaltagmatch48));
    printf("** ratio tag access 16 %lf\n", sim.cfg.I / ((double)sim.cfg.I + totaltagmatch16));
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
        CACHEBLOCK = 16;
        CACHEBLOCKLOG = 4;
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
        CACHEBLOCK = 144;
        CACHEBLOCKLOG = 8;
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
        CACHEBLOCK = 16;
        CACHEBLOCKLOG = 4;
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
        CACHEBLOCK = 4;
        CACHEBLOCKLOG = 2;
        setSET();
        runTile(sim.cfg.kkk);

        // return to the default setting
        CACHEBLOCK = 16;
        CACHEBLOCKLOG = 4;
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
