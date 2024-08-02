#include "parlay/parallel.h"
#include "../adj_listweight/adj_listweight.h"
#include "../adj_listweight/adj_listweight.cpp"
#include <fstream>
#include <cinttypes>
#include <cmath>
#include <../../eigen/Eigen/Core>
#include <../../eigen/Eigen/SparseCore>
#include <../../spectra/include/Spectra/MatOp/SparseSymMatProd.h>
#include <../../spectra/include/Spectra/SymEigsShiftSolver.h>
#include <../../spectra/include/Spectra/SymEigsSolver.h>
#include <bulk/bulk.hpp>
#include <bulk/backends/thread/thread.hpp>

using namespace std;
using namespace Eigen;

class Prune {

public:
    static bool NORMALIZE;
    static double NORM_CONSTANT;

    map<pair<int64_t, int64_t>, vector<double>> preComputedVolumes;
    map<pair<int64_t, int64_t>, double> halfLambdas;
    map<pair<int64_t, int64_t>, SparseMatrix<double>> aggregations;

    AdjListWeight timegraph;
    NodeMap nodegraph;

    int64_t maxNode = 0;
    int64_t totalTime = 0;
    int64_t startTime = 0;
    int64_t endTime  = 0;


     explicit Prune(const string &path) {
        NORMALIZE = true;
        NORM_CONSTANT = 0.0;
        readGraph(path);
    }

    void readGraph(const string &path) {
        maxNode = numeric_limits<uint64_t>::min();
        timegraph.addFromFile(path);
        startTime = timegraph.getMinTime();
        endTime = timegraph.getMaxTime();
        totalTime = endTime - startTime;
        maxNode = timegraph.getMaxVertex();
    };

    void preaggregate();

    double getBounds(int64_t start, int64_t end);

    static double normalizeConductance(double phi, int64_t num_timestamps);

    static double unnormalizeConductance(double phi, int64_t num_timestamps);

    void getNodeWeights();

    vector<double> computeVolumes(int64_t start, int64_t end);

    SparseMatrix<double> getAggregation(int64_t start, int64_t end);

    vector<pair<int64_t,int64_t>> getSpans(int64_t start,int64_t end);

    double getIntervalWeight(int64_t subStart, int64_t subEnd, int64_t wholeStart, int64_t wholeEnd);

    double getCachedDegree(int64_t start, int64_t end, int64_t node, const vector<pair<int64_t,int64_t>>& spans);
};

bool Prune::NORMALIZE = true;
double Prune::NORM_CONSTANT = 0.0; //0.0


void Prune::preaggregate() {

    //in the first step form the Laplacian Graph for every timestamp denoted with {i,i}
    int64_t size = 0;
    //timespan
    pair<int64_t, int64_t> wholeSpan = {};
    bulk::thread::environment env;
    auto start = this->startTime;
    auto end = this->endTime;
    auto tGraph = this->timegraph;
    auto maxNode2 = this->maxNode;
    pair<int64_t, int64_t> wholeSpan2 = {};
    auto aggregations2 = this->aggregations;
    auto totalTime1 = totalTime;

    env.spawn(1, [&start, &end, &tGraph, &maxNode2, &wholeSpan2, &aggregations2](bulk::world &world) {
        int s = world.rank();
        int p = world.active_processors();
        int size = 0;

        auto numbers = bulk::queue<int>(world);
        if (s == 0) {
            for (size = start; size <= end; ++size) {
                numbers(std::hash<double>{}(size) % p).send(size);
            }
        }
        world.sync();

        for (auto time: numbers) {
            SparseMatrix<double> lsm(maxNode2 + 1, maxNode2 + 1);
            Edge tempEdges = tGraph.getEdges().find(time);
            auto it4 = tempEdges.lock_table();
            typedef Triplet<double> T;
            vector<T> tripletList;
            tripletList.reserve(2);

            for (const auto &innerTbl: it4) {
                for (const auto &edge: innerTbl.second) {
                    auto n1 = innerTbl.first;
                    auto n2 = edge.first;
                    auto weight = edge.second;
                    tripletList.emplace_back(n1, n2, (lsm.coeff(n1, n2) - weight));
                    tripletList.emplace_back(n1, n1, (lsm.coeff(n1, n1) + weight));
                }
            }

            lsm.setFromTriplets(tripletList.begin(), tripletList.end());
            wholeSpan2 = pair(time, time);
            aggregations2.insert_or_assign(wholeSpan2, lsm);
            cout << "Whole Span: " << wholeSpan2.first <<","<< wholeSpan2.second << endl;
        }
        world.sync();
    });

    env.spawn(1, [&start, &end, &totalTime1, &tGraph, &maxNode2, &wholeSpan2, &aggregations2](bulk::world &world) {//2
        int s = world.rank();
        int p = world.active_processors();
        int size = 0;

        auto numbers = bulk::queue<int, int>(world);
        if (s == 0) {

            for (size = 2; size < totalTime1; size *= 2) {
                for (int64_t i = start; i + (size - 1) <= end; i += size) {
                    numbers(std::hash<int>{}(size) % p).send(i, size);
                }
            }
        }
        world.sync();


        for (auto [i, size2]: numbers) {

            pair<int64_t, int64_t> wholeSpan = pair(i, i + (size2 - 1));
            pair<int64_t, int64_t> firstSpan = pair(i, i + (size2 / 2 - 1));
            pair<int64_t, int64_t> secondSpan = pair(i + size2 / 2, i + (size2 - 1));

//            cout << "Whole Span: " << wholeSpan.first << "," << wholeSpan.second << endl;
//            cout << "First Span: " << firstSpan.first << "," << firstSpan.second << endl;
//            cout << "Second Span: " << secondSpan.first << "," << secondSpan.second << endl;
//            cout << "---------------------------" << endl;

            SparseMatrix<double> firstM(0, 0);
            SparseMatrix<double> secondM(0, 0);
            firstM = aggregations2.find(firstSpan)->second;
            secondM = aggregations2.find(secondSpan)->second;
            SparseMatrix<double> wholeM = firstM + secondM;
            aggregations2.insert_or_assign(wholeSpan, wholeM);
        }

        world.sync();

    });
    this->aggregations = aggregations2;
}


/**
 * get lower bounds for lambda given a timespan
 * takes the aggregated graph->calculates the normalized laplacian -> computes lower bound lambda(EV) for it
 * @param start
 * @param end
 * @return
 */
double Prune::getBounds(int64_t start, int64_t end) {
//steps: aggregated graphs (Laplacian Graph) ->  Normalized Laplacian of aggregated Graphs -> compute EV

    pair<int64_t, int64_t> span = {start, end};

    //check if lambda (eigenvalues) is already calculated and in halfLambdas
    if (halfLambdas.find(span)!= halfLambdas.end()) {
        double lowerBounds = halfLambdas.find(span)->second;
        return lowerBounds;

    } else {

        /**
         * find the aggregated graph in aggregations
         * getAggregation creates intervals based on Bits and iterates over them, looks for laplacians in aggregations map
         * and aggregates them.
         */
         //check aggregations-> check precompute and precompute construction pattern of intervals
        SparseMatrix<double> lsm = getAggregation(start,end);


        //check if the vols are already in preComputedVolumes (the diagonal degree of the laplacian matrix);
        //inserts it in preCmputedVolumes
        if (preComputedVolumes.find(span) == preComputedVolumes.end()) {

            vector<double> vols(maxNode+1);
            //fill vols array with all the vols
            for (int64_t i = 0; i <= maxNode; ++i) {
                vols[i] = lsm.coeff(i, i);
            }
            //insert the vols into the precomputed array
            preComputedVolumes.insert_or_assign(span, vols);
        }

        //Calculate the normalized Laplacian
        SparseMatrix<double> normalizedLaplacian(maxNode + 1, maxNode + 1);

        while (true) {

            double val = 0;
            double normalizedVal = 0;
            //iterate over the nonzero elements
            for (int k = 0; k < lsm.outerSize(); ++k) {
                for (SparseMatrix<double>::InnerIterator it(lsm, k); it; ++it) {
                    int64_t c = it.col();
                    int64_t r = it.row();
                    val = it.value();
                    //this calculates the normalizedLaplacian Volumes
                    if (lsm.coeff(c, c) == 0.0 && lsm.coeff(r, r) == 0.0) {
                        normalizedLaplacian.coeffRef(r, c) = 0.0;
                    } else {
                        normalizedVal = val / (sqrt(lsm.coeff(c, c)) * sqrt(lsm.coeff(r, r)));
                        normalizedLaplacian.coeffRef(r, c) = normalizedVal;//maybe coeffRef here
                    }
                }
            }

            // Construct matrix operation object using the wrapper class
            Spectra::SparseSymMatProd<double> op(normalizedLaplacian);
            Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op,
                                              (int)log((double)maxNode),maxNode/2);

            eigs.init();
            eigs.compute(Spectra::SortRule::SmallestMagn);
            //retrieve results
            VectorXd evals;
            if (eigs.info() == Spectra::CompInfo::Successful) evals = eigs.eigenvalues();
            //cout<<evals<<endl;
            val = evals[evals.size()-2]; //retrieve smallest eigenvalue
            normalizedVal = val/2.0;
            halfLambdas.insert_or_assign(span, normalizedVal);
            return normalizedVal;
        }
    }
}
//normalize the conductance; currently 0 means /1 so no normalization
double Prune::normalizeConductance(double phi, int64_t num_timestamps) {
    return !NORMALIZE ? phi : phi / pow((double) num_timestamps, NORM_CONSTANT);
}
double Prune::unnormalizeConductance(double phi, int64_t num_timestamps) {
    return !NORMALIZE ? phi : phi * pow((double) num_timestamps, NORM_CONSTANT);
}

void Prune::getNodeWeights() {
    nodegraph = timegraph.nodes;
}
/**
 * creates Spans based on highestBit and lowestBit of the Integer
 * example output for (0,44) -> (0,31) (32,39) (40,43) (44,44)
 * @param start
 * @param end
 * @return
 */
vector<pair<int64_t, int64_t>> Prune::getSpans(int64_t start, int64_t end) {

    vector<pair<int64_t, int64_t>> result;
    int64_t adv;

    for (result = {}; start <= end; start += adv) {
        //finds the lowest one bit
        adv = start & -start;

        if (adv == 0) {
            //finds highest one bit
            adv = 1 << (31 - __builtin_clz(end + 1));  // Find the highest power of 2 less than or equal to end + 1
        }

        while (start + adv > end + 1) adv /= 2;
        pair<int64_t, int64_t> p={start, start + adv - 1};
        result.emplace_back(p);
    }
    return result;
}

/**
 * gets aggregation, uses getSpans to return lsm aggregation combinations that preaggregate doesnt have.
 * @param start
 * @param end
 * @return
 */
SparseMatrix<double> Prune::getAggregation(int64_t start, int64_t end) {

    SparseMatrix<double> result(maxNode +1, maxNode +1);
    result.setZero();
    //creates intervals based on Bits
    vector<pair<int64_t,int64_t>> intervals = getSpans(start,end);
    //looks if intervals exist in aggregations that were previously precomputed

//    for(auto span:intervals){
//        cout<<"SpansfromAggregation:"<<span.first<<","<<span.second<<endl;
//    }

    for (auto &interval: intervals) {
        if (aggregations.find(interval) == aggregations.end()) {
            cout << "Missing aggregation:" << interval.first << interval.second << endl;
        }

        SparseMatrix<double> curM = aggregations.find(interval)->second;
        //aggregates the found intervals
        result += curM;
    }
    return result;
}
/**
 * returns the volumes of the aggregated lsm
 * @param start
 * @param end
 * @return
 */

vector<double> Prune::computeVolumes(int64_t start, int64_t end) {
    //gets the aggregations of lsm between the start, end based on getSpan intervals
    //aggregations are laplacian matrices
    SparseMatrix<double> lsm = getAggregation(start,end);
    vector<double> vols(maxNode+1);
    //fill vols array with all the vols
    for (int64_t i = 0; i <= maxNode; ++i) {
        vols[i] = lsm.coeff(i, i);
    }
    return vols;
}

double Prune::getIntervalWeight(int64_t subStart, int64_t subEnd, int64_t wholeStart, int64_t wholeEnd) {

    double weight = numeric_limits<double>::max();
    /**
     * spans: current interval that is checked (1,1)
     * wholespans:  (1,1) (2,3) (4,4)
     */
    vector<pair<int64_t, int64_t>> spans = getSpans(subStart, subEnd);
    vector<pair<int64_t, int64_t>> wholeSpans = getSpans(wholeStart, wholeEnd);

    for (int64_t node = 0; node <= maxNode; ++node) {
        //gets vol (1,1) of node in position [node]
        double top = getCachedDegree(subStart, subEnd, node, spans);

        if (top != 0.0) {
            //
            double bottom = getCachedDegree(wholeStart, wholeEnd, node, wholeSpans);
            //
            double ratio = top / bottom;
            weight = ratio < weight ? ratio : weight;
        }
    }

    return weight;
}

double Prune::getCachedDegree(int64_t start, int64_t end, int64_t node, const vector<pair<int64_t, int64_t>>& spans) {
    //current span: (1,1)
    pair<int64_t, int64_t> span(start, end);

    vector<double> vols;

    if (preComputedVolumes.find(span)!=preComputedVolumes.end()) {
        //get vols of the current span/interval (1,1); reminder vols are computed for every interval
        vols = preComputedVolumes.find(span)->second;
        if (vols[node] > 0.0) {
            return vols[node];
        }
    } else {
        vols.resize(maxNode + 1);
    }

    double result = 0.0;

    for (auto p: spans) {
        if (preComputedVolumes.find(p)!=preComputedVolumes.end()) {
            vector<double> v = preComputedVolumes.find(p)->second;
            result += v[node];
        }
    }
    vols[node] = result;
    return result;
}



