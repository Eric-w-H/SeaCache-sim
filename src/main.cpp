#include "cache.h"
#include "estimation.h"
#include "headers.h"
#include "json.hpp"
#include "simulator.h"
#include <cstdlib>
#include <fstream>

struct simulator_state sim;

using json = nlohmann::json;

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <matrix1 name> <matrix2 name> <config filepath>\n"
                  << "config_filepath is a fully qualified path to the .json config for the run (likely config/config.json).\n"
                  << std::endl;
        return 1;
    }

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

    // dataflow = Gust;
    // format = RR;
    int transpose = config["transpose"].get<int>();
    float tmpsram = config["cachesize"].get<float>();
    cachesize = tmpsram * 262144 * 0.9;
    inputcachesize = cachesize;
    float tmpbandw = config["memorybandwidth"].get<float>();
    HBMbandwidth = (tmpbandw / 4.0) * 0.6;
    int tmpPE = config["PEcnt"].get<int>();
    PEcnt = tmpPE;
    mergecnt = tmpPE;
    HBMbandwidthperPE = HBMbandwidth / PEcnt;
    int tmpbank = config["srambank"].get<int>();
    sramBank = tmpbank;
    ISCACHE = 1;
    int baselinetest = config["baselinetest"].get<int>();
    bool condensedOP = config["condensedOP"].get<bool>();
    std::string tile_dir = config["tileDir"].get<std::string>();
    std::string output_dir = config["outputDir"].get<std::string>();

    if (!freopen(("data/" + matrix1_name + ".mtx").c_str(), "r", stdin)) {
        std::cerr << "Error opening input file." << std::endl;
        return 1;
    }

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

    // Lines limited to 1024 characters by spec
    const std::size_t BUFFERSIZE = 1024;
    char readbuffer[BUFFERSIZE];

    // read and ignore annotation '%' lines
    while (std::cin.getline(readbuffer, BUFFERSIZE)) {
        if (readbuffer[0] != '%') {
            break;
        }
    }

    std::sscanf(readbuffer, "%d%d%d", &N, &M, &nzA);

    if (transpose) {
        swap(N, M);
    }

    printf("Matrix A: %d x %d, number of non-zeros = %d\n", N, M, nzA);
    fflush(stdout);

    samplek = 100;
    samplep = 0.1;

    sim.cfg.I = N;
    sim.cfg.J = M;

    initialize_data_A();

    string input;

    fflush(stdout);

    for (int i = 1; i <= nzA; i++) {

        std::getline(std::cin, input);

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
            return 0;
        }

        // WARNING(ejs): mtx indices are stored 1-based; it converts 0-based representation in code
        if (transpose) {
            Ac[xx - 1].push_back(yy - 1);
            A[yy - 1].push_back(xx - 1);
        } else {
            A[xx - 1].push_back(yy - 1);
            Ac[yy - 1].push_back(xx - 1);
        }
    }

    for (int i = 0; i < sim.cfg.I; i++) {
        sort(A[i].begin(), A[i].end());
    }
    for (int j = 0; j < sim.cfg.J; j++) {
        sort(Ac[j].begin(), Ac[j].end());
    }

    if (condensedOP) {
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
        offsetarrayA[i] = offsetarrayA[i - 1] + A[i - 1].size();
        if (A[i - 1].size() < 48) {
            totalempty += (48 - A[i - 1].size());
        }
        totalincache += min(48, (int)A[i - 1].size());
        totaltagmatch48 += ((int)A[i - 1].size() + 47) / 48;
        totaltagmatch16 += ((int)A[i - 1].size() + 15) / 16;
    }

    printf("*** ratio of empty %lf, ratio of not empty %lf\n",
           totalempty / (sim.cfg.I * 48.0), 1 - (totalempty / (sim.cfg.I * 48.0)));
    printf("*** ratio of in cache %lf\n", totalincache / ((double)nzA));

    printf("** ratio tag access 48 %lf\n", sim.cfg.I / ((double)sim.cfg.I + totaltagmatch48));
    printf("** ratio tag access 16 %lf\n", sim.cfg.I / ((double)sim.cfg.I + totaltagmatch16));

    for (int i = 1; i < sim.cfg.J + 2; i++) {
        offsetarrayAc[i] = offsetarrayAc[i - 1] + Ac[i - 1].size();
    }

    initsample();

    sampleA();

    fclose(stdin);

    /////////////////////////// input B /////////////////////////////////

    cin.clear();
    if (!freopen(("data/" + matrix2_name + ".mtx").c_str(), "r", stdin)) {
        std::cerr << "Error opening input file." << std::endl;
        return 1;
    }

    // read and ignore annotation '%' lines
    while (std::cin.getline(readbuffer, BUFFERSIZE)) {
        if (readbuffer[0] != '%') {
            break;
        }
    }

    std::sscanf(readbuffer, "%d%d%d", &N, &M, &nzB);

    printf("Matrix B: %d x %d, number of non-zeros = %d\n", N, M, nzB);
    fflush(stdout);

    // FIXME(ejs): This is extremely confusing. It should not automatically try to invert the matrix
    // if (and only if !!?) it is non-square. The user should be responsible for storing a separate transpose
    // version of the matrix (or add some intermediate helper that does the transposition).
    if (N != M)
        transpose ^= 1; // when transposeA = 0 -> transposeB = 1; when tranposeA=
                        // 1-> transposeB = 0

    if (transpose) {
        swap(N, M);
    }

    printf("transpose: %d\n", transpose);

    if (sim.cfg.J != N) {
        printf("Mismatch J!\n");
        return 0;
    }

    sim.cfg.K = M;

    initialize_data_B();

    // std::getline(std::cin, input);

    for (int i = 1; i <= nzB; i++) {

        std::getline(std::cin, input);

        std::istringstream iss(input);
        std::vector<std::string> tokens;
        std::string token;

        // 将输入行分割为单词
        while (iss >> token) {
            tokens.push_back(token);
        }

        int xx, yy;
        double zz, lala;

        if (tokens.size() == 2) {

            std::istringstream(tokens[0]) >> xx;
            std::istringstream(tokens[1]) >> yy;
            // std::cout << "values: " << xx << ", " << yy << std::endl;
        } else if (tokens.size() == 3) {

            std::istringstream(tokens[0]) >> xx;
            std::istringstream(tokens[1]) >> yy;
            std::istringstream(tokens[2]) >> zz;
            // std::cout << "values: " << xx << ", " << yy << ", " << zz << std::endl;
        } else if (tokens.size() == 4) {

            std::istringstream(tokens[0]) >> xx;
            std::istringstream(tokens[1]) >> yy;
            std::istringstream(tokens[2]) >> zz;
            std::istringstream(tokens[2]) >> lala;
            // std::cout << "values: " << xx << ", " << yy << ", " << zz << std::endl;
        } else {
            std::cout << "Format Incorrect! " << std::endl;
            cout << i << endl
                 << input << endl;
            return 0;
        }

        if (transpose) {
            Bc[xx - 1].push_back(yy - 1);
            B[yy - 1].push_back(xx - 1);
        } else {
            B[xx - 1].push_back(yy - 1);
            Bc[yy - 1].push_back(xx - 1);
        }
    }

    // cout << N << endl<<M <<endl<< nz << endl << nz/N <<endl;

    for (int j = 0; j < sim.cfg.J; j++) {
        sort(B[j].begin(), B[j].end());
    }
    for (int k = 0; k < sim.cfg.K; k++) {
        sort(Bc[k].begin(), Bc[k].end());
    }

    for (int j = 1; j < sim.cfg.J; j++) {
        offsetarrayB[j] = offsetarrayB[j - 1] + B[j - 1].size();
    }
    // two problem:
    // 1) this calculate way just calculate the minimum
    // 2) the + J will change is tiling J  -> but actually long will alos change
    // -> counteract? but the above calculate seems don't consider the emptys
    // (larger than real) so maybe counteract

    // move this to above for the weights (1 -> )
    // shortpart += J/(CACHEBLOCKSHORT);

    for (int k = 1; k < sim.cfg.K; k++) {
        offsetarrayBc[k] = offsetarrayBc[k - 1] + Bc[k - 1].size();
    }

    sampleB();

    if (ISCACHE == 1) {
        setSET();
    }

    /******************Config************************************/

    // notation of J and K in the code is swapped as in the paper
    // use the paper's notation as print output
    printf("I = %d, K = %d, J = %d\n", sim.cfg.I, sim.cfg.K, sim.cfg.J);
    /************************************************************/

    // getParameterSample();
    getParameter();

    configPartial(0.05, 0.5, 0.45);

    int t_i, t_j, t_k;

    if (!freopen((tile_dir + matrix1_name).c_str(), "r", stdin)) {
        std::cerr << "Error opening " << (tile_dir + matrix1_name) << std::endl;
        return 1;
    }

    if (std::scanf("%d%d%d", &t_i, &t_j, &t_k) != 3) {
        std::cerr << "Error reading " << (tile_dir + matrix1_name) << ", expected three integers." << std::endl;
        return 1;
    }
    fclose(stdin);

    const struct config cfg = {
        .dataflow   = Gust,
        .interorder = IJK,
        .format     = RR,

        .I          = (u64)N,
        .J          = (u64)M,
        .K          = (u64)M, 
        .iii        = (u64)t_i,
        .jjj        = (u64)t_j,
        .kkk        = (u64)t_k,
        .tti        = (u64)div_rup(sim.cfg.I, t_i),
        .ttj        = (u64)div_rup(sim.cfg.J, t_j),
        .ttk        = (u64)div_rup(sim.cfg.K, t_k),
    };

    sim = initialize_simulator(&cfg);

    /////////////// Baseline configurations

    if (baselinetest) {
        // EWH
        // Incorporate SeaCache into baseline
        puts("***************** SeaCache *******************");
        printf("nnzB:%d  K:%d  J/TJ:%d  nzlB:%d\n", nzB, sim.cfg.K, (sim.cfg.J + sim.cfg.jjj - 1) / sim.cfg.jjj,
               nzB / (sim.cfg.K * ((sim.cfg.J + sim.cfg.jjj - 1) / sim.cfg.jjj)));

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

        fflush(stdout);

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

        fflush(stdout);

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

        fflush(stdout);

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
