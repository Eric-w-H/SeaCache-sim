#include "headers.h"
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
    // These only live in this function
    std::unique_ptr<map<int, bool>[]> estC;
    std::unique_ptr<int[]> endA, endB;

    unordered_set<double> hashset;
    priority_queue<double> hashqu;

    unordered_set<pair<double, int>, pair_hash> hashsetr;
    priority_queue<pair<double, int>> hashqur;

    unordered_set<double> Fset;
    unordered_set<pair<double, int>, pair_hash> Fsetr;

    unordered_set<double> hashsetTJ[256];
    priority_queue<double> hashquTJ[256];

    long long estnnzC, nnzCTk[33];

    vector<double> vectorTK[8195];

    try {
        estC = std::make_unique<map<int, bool>[]>(max(sim.cfg.I, sim.cfg.J));
        endA = std::make_unique<int[]>(sim.cfg.I);
        endB = std::make_unique<int[]>(sim.cfg.J);
    } catch (const std::bad_alloc &e) {
        std::cerr << "Allocation failed: " << e.what() << std::endl;
        exit(1);
    }

    // get parameters in force

    // calculate estEffMAC and nnzC
    estEffMAC = 0;
    estnnzC = 0;

    for (int i = 0; i < sim.cfg.I; i++)
        endA[i] = A[i].size();

    for (int j = 0; j < sim.cfg.J; j++)
        endB[j] = B[j].size();

    for (int i = 0; i < sim.cfg.I; i++) {
        int lenA = endA[i];
        for (int j = 0; j < lenA; j++) {
            int tmpx = A[i][j];

            int lenB = endB[tmpx];
            estEffMAC += lenB;

            for (int k = 0; k < lenB; k++) {

                int tmpk = B[tmpx][k];
                if (!estC[i][tmpk]) {
                    estC[i][tmpk] = 1;
                    estnnzC++;
                }
            }
        }
    }

    printf("MAC: %lld   nnzC: %lld\n", estEffMAC, estnnzC);
}
