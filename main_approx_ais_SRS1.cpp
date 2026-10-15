#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <numeric>
#include <tuple>
#include <cmath>
#include <climits>
#include <iomanip>

#include "read_data.cpp"
#include "motif_id.cpp"

using namespace std;


// =============================================================
// Data structures
// =============================================================

struct hyperwedge {
    int a;
    int b;
    int C_ab;
};

struct EdgeInfo {
    int a;
    int b;
    int C_ab;
    int q;
    int bin;
};


// =============================================================
// Original helper function
// =============================================================

inline long long convert_id(int hyperedge_a, int hyperedge_b)
{
    return hyperedge_a * (1LL << 31) + hyperedge_b;
}




inline double elapsed_ms(clock_t start, clock_t end)
{
    return (double)(end - start) / CLOCKS_PER_SEC * 1000.0;
}




int main(int argc, char* argv[])
{


    if (argc < 2) {
        cerr << "Usage: " << argv[0]
             << " <total_samples> [graph_file]" << endl;
        cerr << endl;
        cerr << "Example:" << endl;
        cerr << "  " << argv[0]
             << " 10000 dblp_graph.txt" << endl;
        return 1;
    }


    // =========================================================
    // Sampling budget
    // =========================================================

    long long total_samples;

    try {
        total_samples = stoll(argv[1]);
    }
    catch (...) {
        cerr << "Error: invalid total_samples: "
             << argv[1] << endl;
        return 1;
    }



    random_device rd;

    auto now =
        chrono::high_resolution_clock::now()
            .time_since_epoch()
            .count();

    seed_seq seed_source {
        rd(),
        rd(),
        static_cast<unsigned int>(now),
        static_cast<unsigned int>(
            static_cast<unsigned long long>(now) >> 32
        )
    };

    unsigned int seed;

    seed_source.generate(
        &seed,
        &seed + 1
    );
    // =========================================================
    //g++ -O3 -std=c++11 main_approx_ais_SRS1.cpp -o approx_ais_SRS1;
    //./approx_ais_SRS1 10000
    string graphFile = "email-Eu-full.txt";

    if (argc >= 3) {
        graphFile = argv[2];
    }


    // =========================================================
    // Validate sampling budget
    // =========================================================

    if (total_samples <= 0) {
        cerr << "Error: total_samples must be positive." << endl;
        return 1;
    }


    cout << fixed << setprecision(3);

    cout << "==============================================" << endl;
    cout << "SRS Approximate h-Motif Counting" << endl;
    cout << "Dataset       : " << graphFile << endl;
    cout << "Sampling size : " << total_samples << endl;
    cout << "Random seed   : " << seed << endl;
    cout << "==============================================" << endl << endl;


    // =========================================================
    // Total runtime starts BEFORE dataset loading
    // =========================================================

    clock_t global_start = clock();


    // =========================================================
    // 1. Read data
    // =========================================================

    clock_t start = clock();

    vector<vector<int>> node2hyperedge;
    vector<vector<int>> hyperedge2node;
    vector<unordered_set<int>> hyperedge2node_set;

    read_data(
        graphFile,
        node2hyperedge,
        hyperedge2node,
        hyperedge2node_set
    );

    clock_t read_end = clock();

    int V = node2hyperedge.size();
    int E = hyperedge2node.size();

    cout << "# of nodes: " << V << '\n';
    cout << "# of hyperedges: " << E << '\n';

    cout << "Reading data done: "
         << elapsed_ms(start, read_end)
         << " ms" << endl;

    cout << "------------------------------------------"
         << endl << endl;

    if (E == 0) {
        cerr << "Error: the input hypergraph contains no hyperedges."
             << endl;
        return 1;
    }


    // =========================================================
    // 2. Adjacency-list construction
    // =========================================================

    start = clock();

    hyperedge2node.resize(E);
    hyperedge2node_set.resize(E);

    vector<vector<pair<int, int>>> hyperedge_adj;
    vector<unordered_map<int, int>> hyperedge_inter;

    hyperedge_adj.resize(E);
    hyperedge_inter.resize(E);

    vector<hyperwedge> W;
    vector<long long> upd_time(E, -1LL);


    // ---------------------------------------------------------
    // Original adjacency-list construction
    // ---------------------------------------------------------

    for (int hyperedge_a = 0; hyperedge_a < E; hyperedge_a++) {
        long long l_hyperedge_a = (long long)hyperedge_a;

        for (const int& node : hyperedge2node[hyperedge_a]) {
            for (const int& hyperedge_b : node2hyperedge[node]) {
                if (hyperedge_b == hyperedge_a) {
                    continue;
                }

                if ((upd_time[hyperedge_b] >> 31) ^ hyperedge_a) {
                    upd_time[hyperedge_b] =
                        (l_hyperedge_a << 31)
                        + (long long)hyperedge_adj[hyperedge_b].size();

                    hyperedge_adj[hyperedge_b].push_back({
                        hyperedge_a,
                        0
                    });
                }

                hyperedge_adj[hyperedge_b]
                    [(int)(upd_time[hyperedge_b] & 0x7FFFFFFFLL)]
                    .second++;
            }
        }
    }


    // ---------------------------------------------------------
    // Complete adjacency-list construction
    // ---------------------------------------------------------

    for (int hyperedge_a = 0; hyperedge_a < E; hyperedge_a++) {
        int deg_a = hyperedge_adj[hyperedge_a].size();

        hyperedge_inter[hyperedge_a].rehash(deg_a);

        for (int i = 0; i < deg_a; i++) {
            int hyperedge_b = hyperedge_adj[hyperedge_a][i].first;
            int C_ab = hyperedge_adj[hyperedge_a][i].second;

            if (hyperedge_a < hyperedge_b) {
                W.push_back(hyperwedge{
                    hyperedge_a,
                    hyperedge_b,
                    C_ab
                });
            }

            hyperedge_inter[hyperedge_a].insert({
                hyperedge_b,
                C_ab
            });
        }
    }


    clock_t projection_end = clock();

    long long total_pairs = (long long)W.size();

    cout << "# of hyperwedges: "
         << total_pairs << "\n";

    cout << "Adjacency list construction done: "
         << elapsed_ms(start, projection_end)
         << " ms" << endl;

    cout << "------------------------------------------"
         << endl << endl;

    if (total_pairs == 0) {
        cerr << "Error: no hyperwedges were found." << endl;
        return 1;
    }

    if (total_samples > total_pairs) {
        cerr << "Error: total_samples = "
             << total_samples
             << " exceeds |Lambda| = "
             << total_pairs
             << "."
             << endl;
        return 1;
    }


    // =========================================================
    // 3. Degree-bin assignment
    // =========================================================

    start = clock();

    vector<int> degree_bins = {
        0, 10, 40, 160, 320, 640, INT_MAX
    };

    int num_degree_bins = degree_bins.size() - 1;

    vector<int> hyperedge_bin(E, 0);
    vector<int> hyperedge_size(E, 0);

    int max_hyperedge_size = 0;

    for (int a = 0; a < E; ++a) {
        hyperedge_size[a] = hyperedge2node[a].size();

        max_hyperedge_size = max(
            max_hyperedge_size,
            hyperedge_size[a]
        );

        int deg = hyperedge_adj[a].size();
        int assigned_bin = num_degree_bins - 1;

        if (deg == 0) {
            assigned_bin = 0;
        }
        else {
            for (int i = 0; i < num_degree_bins; ++i) {
                if (deg > degree_bins[i]
                    && deg <= degree_bins[i + 1]) {

                    assigned_bin = i;
                    break;
                }
            }
        }

        hyperedge_bin[a] = assigned_bin;
    }

    cout << "Degree bins K: "
         << num_degree_bins
         << endl;


    // =========================================================
    // 4. First scan:
    // Count hyperwedges in each stratum
    // =========================================================

    unordered_map<long long, long long> layer_count;

    for (const hyperwedge& hw : W) {
        int a = hw.a;
        int b = hw.b;
        int q = hw.C_ab;

        int bin = max(
            hyperedge_bin[a],
            hyperedge_bin[b]
        );

        /*
         * s(Lambda_ab)
         *
         * = K * |e_a intersection e_b|
         *
         * + max(bin(e_a), bin(e_b))
         */

        long long key =
            (long long)q * num_degree_bins + bin;

        layer_count[key]++;
    }


    clock_t stratification_end = clock();

    long long H = (long long)layer_count.size();

    cout << "# nonempty strata H: "
         << H
         << endl;

    cout << "Stratification done: "
         << elapsed_ms(start, stratification_end)
         << " ms" << endl;

    cout << "------------------------------------------"
         << endl << endl;


    // =========================================================
    // Positive allocation requires r >= H
    // =========================================================

    if (total_samples < H) {
        cerr << "Error: sampling budget is too small." << endl;
        cerr << "total_samples = "
             << total_samples << endl;
        cerr << "number of nonempty strata H = "
             << H << endl;
        cerr << "Positive allocation requires r >= H."
             << endl;
        return 1;
    }


    // =========================================================
    // 5. Positive sample allocation
    // =========================================================

    start = clock();

    vector<long long> layer_keys;
    vector<long long> N_list;
    vector<long long> alloc_list;
    vector<long double> remainders;

    layer_keys.reserve(layer_count.size());
    N_list.reserve(layer_count.size());
    alloc_list.reserve(layer_count.size());
    remainders.reserve(layer_count.size());


    // ---------------------------------------------------------
    // One sample for every nonempty stratum
    // ---------------------------------------------------------

    for (const auto& kv : layer_count) {
        layer_keys.push_back(kv.first);
        N_list.push_back(kv.second);
        alloc_list.push_back(1);
        remainders.push_back(0.0L);
    }


    /*
     * R = r - H
     */

    long long remaining_samples =
        total_samples - H;


    /*
     * Residual population = |Lambda| - H
     */

    long long remaining_capacity =
        total_pairs - H;


    // ---------------------------------------------------------
    // Residual proportional allocation
    // ---------------------------------------------------------

    if (remaining_samples > 0) {
        if (remaining_capacity <= 0) {
            cerr << "Error: invalid residual sampling capacity."
                 << endl;
            return 1;
        }

        for (size_t i = 0; i < layer_keys.size(); ++i) {
            long long capacity_s = N_list[i] - 1;

            /*
             *                  (r-H)(N_s-1)
             * a_tilde_s = ----------------------
             *                    |Lambda|-H
             */

            long double exact =
                (long double)remaining_samples
                * (long double)capacity_s
                / (long double)remaining_capacity;

            long long base =
                (long long)floor(exact);

            if (base > capacity_s) {
                base = capacity_s;
            }

            alloc_list[i] += base;
            remainders[i] =
                exact - (long double)base;
        }
    }


    // ---------------------------------------------------------
    // Samples allocated after flooring
    // ---------------------------------------------------------

    long long sum_alloc =
        accumulate(
            alloc_list.begin(),
            alloc_list.end(),
            0LL
        );


    // ---------------------------------------------------------
    // Largest remainder method
    // ---------------------------------------------------------

    if (sum_alloc < total_samples) {
        vector<size_t> idx(alloc_list.size());

        iota(idx.begin(), idx.end(), 0);

        stable_sort(
            idx.begin(),
            idx.end(),
            [&](size_t i, size_t j) {
                if (remainders[i] != remainders[j]) {
                    return remainders[i] > remainders[j];
                }

                return layer_keys[i] < layer_keys[j];
            }
        );

        for (size_t i : idx) {
            if (sum_alloc >= total_samples) {
                break;
            }

            if (alloc_list[i] < N_list[i]) {
                alloc_list[i]++;
                sum_alloc++;
            }
        }
    }


    // =========================================================
    // Allocation consistency check
    // =========================================================

    if (sum_alloc != total_samples) {
        cerr << "Error: sample allocation failed." << endl;
        cerr << "Allocated = "
             << sum_alloc
             << ", requested = "
             << total_samples
             << endl;
        return 1;
    }


    long long zero_alloc_layers = 0;
    long long min_alloc = LLONG_MAX;
    long long max_alloc = 0;

    for (size_t i = 0; i < layer_keys.size(); ++i) {
        if (alloc_list[i] <= 0) {
            zero_alloc_layers++;
        }

        if (alloc_list[i] > N_list[i]) {
            cerr << "Error: allocation exceeds "
                 << "stratum population."
                 << endl;

            cerr << "key = "
                 << layer_keys[i]
                 << ", N_s = "
                 << N_list[i]
                 << ", r_s = "
                 << alloc_list[i]
                 << endl;

            return 1;
        }

        min_alloc = min(min_alloc, alloc_list[i]);
        max_alloc = max(max_alloc, alloc_list[i]);
    }


    if (zero_alloc_layers != 0) {
        cerr << "Error: positive allocation failed." << endl;
        return 1;
    }


    // =========================================================
    // Build allocation map
    // =========================================================

    unordered_map<long long, long long> alloc_map;

    alloc_map.rehash(layer_count.size());

    for (size_t i = 0; i < layer_keys.size(); ++i) {
        alloc_map.insert({
            layer_keys[i],
            alloc_list[i]
        });
    }


    clock_t allocation_end = clock();

    cout << "Sample allocation done: "
         << elapsed_ms(start, allocation_end)
         << " ms" << endl;

    cout << "# nonempty strata        : "
         << H << endl;

    cout << "# zero-allocation strata : "
         << zero_alloc_layers << endl;

    cout << "Minimum r_s              : "
         << min_alloc << endl;

    cout << "Maximum r_s              : "
         << max_alloc << endl;

    cout << "Allocated samples        : "
         << sum_alloc << endl;

    cout << "[OK] Every nonempty stratum satisfies "
         << "1 <= r_s <= N_s."
         << endl;

    cout << "------------------------------------------"
         << endl << endl;


    // =========================================================
    // 6. Stratified reservoir sampling
    // =========================================================

    start = clock();


    // =========================================================
    // IMPORTANT:
    //
    // A new random seed is automatically generated
    // for every independent run.
    //
    // The generated seed is printed in the output
    // for reproducibility.
    // =========================================================

    mt19937 gen(seed);

    unordered_map<long long, vector<EdgeInfo>> reservoirs;
    unordered_map<long long, long long> seen;

    reservoirs.rehash(layer_count.size());
    seen.rehash(layer_count.size());


    // ---------------------------------------------------------
    // Second scan over all hyperwedges
    // ---------------------------------------------------------

    for (const hyperwedge& hw : W) {
        int a = hw.a;
        int b = hw.b;
        int C_ab = hw.C_ab;

        int q = C_ab;

        int bin = max(
            hyperedge_bin[a],
            hyperedge_bin[b]
        );

        long long key =
            (long long)q * num_degree_bins + bin;

        auto allocation_it =
            alloc_map.find(key);

        if (allocation_it == alloc_map.end()) {
            cerr << "Error: missing allocation "
                 << "for stratum "
                 << key
                 << endl;
            return 1;
        }

        long long n_sample =
            allocation_it->second;

        vector<EdgeInfo>& res =
            reservoirs[key];

        long long& cnt =
            seen[key];


        // -----------------------------------------------------
        // Insert the first r_s hyperwedges
        // -----------------------------------------------------

        if (cnt < n_sample) {
            res.push_back(
                EdgeInfo{
                    a,
                    b,
                    C_ab,
                    q,
                    bin
                }
            );
        }
        else {
            /*
             * Standard reservoir sampling.
             *
             * Current item number = cnt + 1.
             *
             * Draw j uniformly from:
             *
             * {0,...,cnt}
             */

            uniform_int_distribution<long long> dist(0, cnt);

            long long j = dist(gen);

            if (j < n_sample) {
                res[(size_t)j] =
                    EdgeInfo{
                        a,
                        b,
                        C_ab,
                        q,
                        bin
                    };
            }
        }

        cnt++;
    }


    // =========================================================
    // Collect sampled hyperwedges and their weights
    // =========================================================

    vector<tuple<int, int, int, long double>> samples;

    samples.reserve((size_t)total_samples);

    for (auto& kv : reservoirs) {
        long long key = kv.first;

        long long n_sample =
            alloc_map[key];

        long long N =
            layer_count[key];

        /*
         * Inverse inclusion-probability weight:
         *
         * w_s = N_s / r_s
         */

        long double weight =
            (long double)N
            / (long double)n_sample;

        for (const EdgeInfo& e : kv.second) {
            samples.emplace_back(
                e.a,
                e.b,
                e.C_ab,
                weight
            );
        }
    }


    clock_t sampling_end = clock();

    cout << "Actual sampled hyperwedges: "
         << samples.size()
         << endl;

    cout << "Reservoir sampling done: "
         << elapsed_ms(start, sampling_end)
         << " ms" << endl;

    cout << "------------------------------------------"
         << endl << endl;


    if ((long long)samples.size() != total_samples) {
        cerr << "Error: actual sample size = "
             << samples.size()
             << ", requested sample size = "
             << total_samples
             << endl;
        return 1;
    }


    // =========================================================
    // 7. Weighted h-motif counting
    // =========================================================

    start = clock();

    vector<long double> h_motif(30, 0.0L);

    vector<int> intersection(
        max(E, max_hyperedge_size),
        0
    );

    fill(
        upd_time.begin(),
        upd_time.end(),
        -1LL
    );


    for (size_t sample = 0;
         sample < samples.size();
         ++sample) {

        int hyperedge_a;
        int hyperedge_b;
        int C_ab;

        long double weight;

        tie(
            hyperedge_a,
            hyperedge_b,
            C_ab,
            weight
        ) = samples[sample];


        int size_a =
            hyperedge2node[hyperedge_a].size();

        int size_b =
            hyperedge2node[hyperedge_b].size();

        int deg_a =
            hyperedge_adj[hyperedge_a].size();

        int deg_b =
            hyperedge_adj[hyperedge_b].size();

        long long current_sample =
            (long long)sample;

        upd_time[hyperedge_a] = current_sample;
        upd_time[hyperedge_b] = current_sample;


        // -----------------------------------------------------
        // Compute e_a intersection e_b
        // -----------------------------------------------------

        int min_ab = hyperedge_a;
        int max_ab = hyperedge_b;

        if (size_a > size_b) {
            min_ab = hyperedge_b;
            max_ab = hyperedge_a;
        }

        const auto& nodes_ab =
            hyperedge2node_set[max_ab];

        auto it_end_ab =
            nodes_ab.end();

        int cnt = 0;

        for (const int& node :
             hyperedge2node[min_ab]) {

            if (nodes_ab.find(node) != it_end_ab) {
                intersection[cnt++] = node;
            }
        }


        // -----------------------------------------------------
        // Candidates adjacent to hyperedge_b
        // -----------------------------------------------------

        for (int i = 0; i < deg_b; i++) {
            int hyperedge_c =
                hyperedge_adj[hyperedge_b][i].first;

            int C_bc =
                hyperedge_adj[hyperedge_b][i].second;

            if (upd_time[hyperedge_c] ^ current_sample) {
                upd_time[hyperedge_c] = current_sample;

                int size_c =
                    hyperedge2node[hyperedge_c].size();

                int C_ca = 0;
                int g_abc = 0;

                auto it_ca =
                    hyperedge_inter[hyperedge_a].find(
                        hyperedge_c
                    );

                if (it_ca !=
                    hyperedge_inter[hyperedge_a].end()) {

                    C_ca = it_ca->second;
                }


                // ---------------------------------------------
                // Triple intersection
                // ---------------------------------------------

                const auto& nodes_c =
                    hyperedge2node_set[hyperedge_c];

                auto it_end_c =
                    nodes_c.end();

                for (int k = 0; k < C_ab; k++) {
                    if (nodes_c.find(intersection[k])
                        != it_end_c) {

                        g_abc++;
                    }
                }


                int motif_index =
                    get_motif_index_new(
                        size_a,
                        size_b,
                        size_c,
                        C_ab,
                        C_bc,
                        C_ca,
                        g_abc
                    );


                /*
                 * Number of constituent hyperwedges:
                 *
                 * open motif   -> 2
                 * closed motif -> 3
                 */

                int nonzero =
                    (C_ab > 0)
                    + (C_bc > 0)
                    + (C_ca > 0);


                if (nonzero >= 2) {
                    h_motif[motif_index] +=
                        weight
                        / (long double)nonzero;
                }
            }
        }


        // -----------------------------------------------------
        // Candidates adjacent only to hyperedge_a
        // -----------------------------------------------------

        for (int i = 0; i < deg_a; i++) {
            int hyperedge_c =
                hyperedge_adj[hyperedge_a][i].first;

            int C_ca =
                hyperedge_adj[hyperedge_a][i].second;

            if (upd_time[hyperedge_c] ^ current_sample) {
                upd_time[hyperedge_c] = current_sample;

                int size_c =
                    hyperedge2node[hyperedge_c].size();


                /*
                 * c was not found in b's neighborhood,
                 * therefore C_bc = 0.
                 */

                int C_bc = 0;


                /*
                 * If C_bc = 0,
                 * the triple intersection must also be 0.
                 */

                int g_abc = 0;


                int motif_index =
                    get_motif_index_new(
                        size_a,
                        size_b,
                        size_c,
                        C_ab,
                        C_bc,
                        C_ca,
                        g_abc
                    );


                int nonzero =
                    (C_ab > 0)
                    + (C_bc > 0)
                    + (C_ca > 0);


                if (nonzero >= 2) {
                    h_motif[motif_index] +=
                        weight
                        / (long double)nonzero;
                }
            }
        }
    }


    clock_t counting_end = clock();


    // =========================================================
    // 8. End total runtime
    // =========================================================

    clock_t global_end = clock();


    // =========================================================
    // 9. Output motif estimates
    // =========================================================

    cout << "Estimated h-motif counts:"
         << endl;

    int index = 0;

    for (int i = 0; i < 30; i++) {
        if (i == 0
            || i == 1
            || i == 4
            || i == 6) {

            continue;
        }

        cout << "h-motif "
             << ++index
             << ": "
             << h_motif[i]
             << endl;
    }


    // =========================================================
    // 10. Runtime summary
    // =========================================================

    double loading_ms =
        elapsed_ms(
            global_start,
            read_end
        );

    double projection_ms =
        elapsed_ms(
            read_end,
            projection_end
        );

    double stratification_ms =
        elapsed_ms(
            projection_end,
            stratification_end
        );

    double allocation_ms =
        elapsed_ms(
            stratification_end,
            allocation_end
        );

    double sampling_ms =
        elapsed_ms(
            allocation_end,
            sampling_end
        );

    double counting_ms =
        elapsed_ms(
            sampling_end,
            counting_end
        );

    double total_runtime_ms =
        elapsed_ms(
            global_start,
            global_end
        );


    cout << "\n"
         << "=============================================="
         << endl;

    cout << "Runtime summary" << endl;

    cout << "=============================================="
         << endl;

    cout << "Dataset loading        : "
         << loading_ms << " ms" << endl;

    cout << "Adjacency preprocessing: "
         << projection_ms << " ms" << endl;

    cout << "Stratification         : "
         << stratification_ms << " ms" << endl;

    cout << "Sample allocation      : "
         << allocation_ms << " ms" << endl;

    cout << "Reservoir sampling     : "
         << sampling_ms << " ms" << endl;

    cout << "Weighted motif counting: "
         << counting_ms << " ms" << endl;

    cout << "TOTAL RUNNING TIME     : "
         << total_runtime_ms << " ms" << endl;

    cout << "=============================================="
         << endl;


    // =========================================================
    // Machine-readable summary
    // =========================================================

    cout << "\nSUMMARY"
         << ",dataset=" << graphFile
         << ",samples=" << total_samples
         << ",seed=" << seed
         << ",K=" << num_degree_bins
         << ",H=" << H
         << ",zero_alloc=" << zero_alloc_layers
         << ",min_rs=" << min_alloc
         << ",max_rs=" << max_alloc
         << ",actual_samples=" << samples.size()
         << ",total_ms=" << total_runtime_ms
         << endl;


    // =========================================================
    // Clear large structures
    // =========================================================

    W.clear();
    node2hyperedge.clear();
    hyperedge2node.clear();
    hyperedge2node_set.clear();
    hyperedge_adj.clear();
    hyperedge_inter.clear();

    return 0;
}
