#include "parlay/parallel.h"
#include <fstream>
#include <cinttypes>
#include "../../dyncomm/TemporalWeightedMinhashGadget.cpp"


using namespace std;

class DynCommDriver {
public:
    static uint64_t GROUPS_TO_REPORT;
    static bool ONE_BASED_TIMELINE;
    static string INPUT_FILENAME;
    static bool PR_USE_PRUNING;
    static bool PR_USE_GROUPS;
    static bool PR_ONLY_PRUNING;
    static bool PR_COMPUTE_ALL_EIGS;
    static double PR_ESTIMATE;
    static bool PR_CALC_ESTIMATE;
    static uint64_t PR_ESTIMATION_HASH_TRIES;
    static uint64_t PR_ESTIMATION_INTERVALS;
    static uint64_t HASH_BANDS;
    static uint64_t HASH_ROWS_PER_BAND;
    static bool HASH_USE_LOG_SCALES;
    static uint64_t HASH_LOG_SCALE_MULTIPLIER;
    static bool RW_USE_TEMPORAL_RW;
    static double RW_ALPHA;
    static uint64_t RW_BINS_TO_EXAMINE;

    static void multiscale_driver(const string &fileName);

    static double lowerKey(const map<double, set<TemporalNode>> &map, double hashScore);
};

uint64_t DynCommDriver::GROUPS_TO_REPORT = 10;
bool DynCommDriver::ONE_BASED_TIMELINE = false;
bool DynCommDriver::PR_USE_PRUNING = true;
bool DynCommDriver::PR_USE_GROUPS = true;
bool DynCommDriver::PR_ONLY_PRUNING = false;
bool DynCommDriver::PR_COMPUTE_ALL_EIGS = false;
double DynCommDriver::PR_ESTIMATE = 0.016;
bool DynCommDriver::PR_CALC_ESTIMATE = true;
uint64_t DynCommDriver::PR_ESTIMATION_HASH_TRIES = 20;
uint64_t DynCommDriver::PR_ESTIMATION_INTERVALS = 10;
uint64_t DynCommDriver::HASH_BANDS = 1;
uint64_t DynCommDriver::HASH_ROWS_PER_BAND = 2;
bool DynCommDriver::HASH_USE_LOG_SCALES = true;
uint64_t DynCommDriver::HASH_LOG_SCALE_MULTIPLIER = 5;
bool DynCommDriver::RW_USE_TEMPORAL_RW = false;
double DynCommDriver::RW_ALPHA = 0.15;
uint64_t DynCommDriver::RW_BINS_TO_EXAMINE = 20;

double DynCommDriver::lowerKey(const map<double, set<TemporalNode>> &map, double hashScore) {

    auto it = map.lower_bound(hashScore);
    if (it != map.begin()) {
        it--;
    }
    return it->first;
}


void DynCommDriver::multiscale_driver(const string &fileName) {
    typedef chrono::duration<uint64_t> Duration;
    typedef chrono::high_resolution_clock Instant;

    AdjListWeight list;
    string basePath = "../";
    set<int64_t> possible_pivots;

    cout << "\t% Begin PHASR:" << endl;
    cout << "\t%%% Begin Import:" << endl;

    list.addFromFile("../" + fileName);

    cout << "\t%%%%%% Finished Import:" << "ms\n";


    uint64_t big = list.getLastVertex();
    uint64_t timeline = list.getMaxTime() + 1;
    uint64_t maxScale = list.getMaxTime() + 1;
    RealPruneContainer pc;

    if (DynCommDriver::PR_USE_PRUNING) {

        cout << "\t%%% Begin Pruning:" << "ms\n";

        string eigstr = DynCommDriver::PR_COMPUTE_ALL_EIGS ? "ALL eigenvalues" : "selected eigenvalues";
        string grstr = DynCommDriver::PR_USE_GROUPS ? "will" : "will NOT";
        cout << "Pruning will compute " + eigstr + " and " + grstr + " try to prune groups of intervals.";

        RealPruneContainer::FAST_ESTIMATE = DynCommDriver::PR_CALC_ESTIMATE;
        RealPruneContainer::HASH_ATTEMPTS = DynCommDriver::PR_ESTIMATION_HASH_TRIES;
        RealPruneContainer::INTERVALS_TO_CHECK = DynCommDriver::PR_ESTIMATION_INTERVALS;

        pc = RealPruneContainer(list, timeline, basePath, fileName, DynCommDriver::ONE_BASED_TIMELINE,
                                DynCommDriver::PR_ESTIMATE,
                                DynCommDriver::PR_USE_GROUPS, DynCommDriver::PR_COMPUTE_ALL_EIGS);


        cout << "\t%%%%%% Finished precomputation:" << "s\n";


        pc.improveEstimate(DynCommDriver::PR_ESTIMATE);
        pc.updatePruning();


        cout << "\t%%%%%% Finished w/ PruneContainer: " << endl;


        cout << "\t%%% Finished pruning: " << endl;


        pc.printPruningStats();

        int i;
        if (HASH_USE_LOG_SCALES) {
            for (i = HASH_LOG_SCALE_MULTIPLIER; i < maxScale; i *= HASH_LOG_SCALE_MULTIPLIER) {
                if (!(pc).isPrunedDuration(i - 1)) {
                    possible_pivots.emplace(2 * maxScale / i);
                }
            }
        } else {
            for (i = maxScale; i > 1; --i) {
                if (!(pc).isPrunedDuration(i - 1)) {
                    possible_pivots.emplace(2 * maxScale / i);
                }
            }
        }

        if (!(pc).isPrunedDuration(maxScale - 1)) {
            possible_pivots.emplace(1);
        }
    } else {}

    if (PR_ONLY_PRUNING) {
        cout << "" << endl;
    } else {
        cout << "Max scale: " << maxScale << endl;
        for (auto const &elem: possible_pivots) {
            cout << "Possible pivots: " << elem << endl;
        }

        vector<double> bestCond;
        vector<double> bestCondHashScore;
        vector<int64_t> bestSetRank;
        vector<set<int64_t>> bestSet;
        MatrixXd bestTimes;

        int64_t timebits;
        for (timebits = 0; timebits < GROUPS_TO_REPORT; ++timebits) {
            bestCond.emplace_back(1.0);
            bestCondHashScore.emplace_back(0.0);
            bestSetRank.emplace_back(RW_BINS_TO_EXAMINE + 1);
            bestSet.emplace_back();
            bestTimes.resize(GROUPS_TO_REPORT, 2);
        }

        auto var15 = possible_pivots.rbegin();

        while (var15 != possible_pivots.rend()) {

            timebits = *var15;
            cout << "\tWorking with " << timebits << " pivots." << endl;
            auto partitionStart = Instant::now();


            cout << "\t%%% Begin single-scale PHASR (" << timebits << " pivots): " << endl;


            int target = 2 * (list.getMaxTime() - list.getMinTime() + 1) / timebits;
            if (target >= timeline) {
                target = timeline - 1;
            }
            cout << "\t%%%%%% Begin hashing: " << endl;
            TemporalWeightedMinhashGadget hasher(big + 1, HASH_BANDS, HASH_ROWS_PER_BAND, timebits, target);
            cout << "\t%%%%%%%%% Finished RNG: " << endl;
            map<double, set<TemporalNode>> hashResult = hasher.performHashing(list, target, pc, RW_USE_TEMPORAL_RW);
            cout << "\t%%%%%% Finished hashing: " << endl;

            if (hashResult.empty()) {
                cout << "\t\tNO HASHING BINS AT THIS SCALE." << endl;
            } else {
                    cout << "\t%%%%%% Begin random walks: " << endl;

                auto hashScoreTemp = hashResult.rbegin();
                double hashScore = hashScoreTemp->first;
                int counter = 0;

                cout << "\t%%%%%% Finished RW setup: " << endl;

                for (bool pruningCheck = false;
                     !isnan(hashScore) && counter < RW_BINS_TO_EXAMINE; hashScore = lowerKey(hashResult, hashScore)) {
                    ++counter;
                    set<TemporalNode> seeds = hashResult.find(hashScore)->second;
                    bool skip = false;
                    if (pruningCheck) {

                        int range = seeds.rbegin()->time - seeds.begin()->time;
                        auto var27 = TemporalNode::extractTimes(seeds).begin();

                        while (var27 != TemporalNode::extractTimes(seeds).end()) {
                            int time = *var27;
                            if (pc.isPruned(time, range)) {
                                skip = true;
                                break;
                            }
                        }
                    }

                    if (!skip) {
                        set<int> verts = {};
                        vector<int> times = TemporalNode::span(seeds);//times the Node neighborhood lasts
                        set<TemporalNode> swipeResult = {};
                        vector<double> conductance = {};
                        map<double, int> tempMap = {};


                        tempMap = list.standardRWR(TemporalNode::extractVertices(seeds), RW_ALPHA, true, false,
                                                   seeds.begin()->time, seeds.rbegin()->time + 1);
                        verts = list.standardSweeps(tempMap, conductance, true, seeds.begin()->time,
                                                    seeds.rbegin()->time + 1);
                        cout << "" << endl;

                    }
                }
            }
            ++var15;
        }
    }
}

