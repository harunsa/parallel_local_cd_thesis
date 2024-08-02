#include <cinttypes>
#include <iostream>
#include <random>
#include <../../eigen/Eigen/Core>
#include <map>
#include <set>
#include "../dyncomm/RealPruneContainer.cpp"
#include <bulk/bulk.hpp>
#include <bulk/backends/thread/thread.hpp>

using namespace Eigen;
using namespace std;



class TemporalWeightedMinhashGadget{


public:
    vector<MatrixXd> r;

    vector<MatrixXd> c;

    vector<MatrixXd> beta;

    MatrixXd psi;

    random_device rg;

    vector<int64_t> mapping;

    int64_t size;

    int64_t bands;

    int64_t rows;

    int64_t timebits = 0;

    int64_t target;

    map<double,set<TemporalNode>> performHashing(AdjListWeight &g, int64_t scale,RealPruneContainer& pc,  bool split);

    TemporalWeightedMinhashGadget(int64_t size, int64_t bands, int64_t rows, int64_t timebits);
    TemporalWeightedMinhashGadget(int64_t size, int64_t bands, int64_t rows, int64_t timebits,int64_t target);

    //map<int,int> tl;
    static string pad(long k, int digits);
    static string turnHex(int t);

    double scoreBin(const set<TemporalNode> &bin, int scale);
    double scoreBin(const set<TemporalNode> &bin);
};
double TemporalWeightedMinhashGadget::scoreBin(const set<TemporalNode> &bin, int scale) {

    double full = (double)(TemporalNode::extractVertices(bin).size() * scale);
    return (double)bin.size() / full;
}

double TemporalWeightedMinhashGadget::scoreBin(const set<TemporalNode> &bin) {

    set<int> nodes = {};
    set<int> times = {};
    auto var5 = bin.begin();

    while(var5!=bin.end()) {
        TemporalNode tn = *var5;
        nodes.emplace(tn.getNode());
        times.emplace(tn.getTime());
    }

    return (double)(nodes.size() * times.size() * bin.size()) * 1.0;
}

string TemporalWeightedMinhashGadget::pad(long k, int digits){
    std::ostringstream ss;
    ss << std::setw( digits ) << std::setfill( '0' ) << k;
    string result = ss.str();
    return result;
}

string TemporalWeightedMinhashGadget::turnHex(int t){
    std::stringstream sstream;
    sstream << std::hex << t;
    std::string result = sstream.str();
    return result;
}

TemporalWeightedMinhashGadget::TemporalWeightedMinhashGadget(int64_t size, int64_t bands, int64_t rows,
                                                             int64_t timebits,int64_t target) : TemporalWeightedMinhashGadget (size,bands,rows,timebits) {
    this->target=target;
}

TemporalWeightedMinhashGadget::TemporalWeightedMinhashGadget(int64_t size, int64_t bands, int64_t rows,int64_t timebits) {
    r.resize(bands);
    for (auto &elem: r) {
    elem.resize(rows,size);
    }

    c.resize(bands);
    for (auto &elem: c) {
        elem.resize(rows,size);
    }

    beta.resize(bands);
    for (auto &elem: beta) {
        elem.resize(rows,size);
    }

    psi.resize(bands,timebits);
    this->mapping = {};
    this->size = size;
    this->bands = bands;
    this->rows = rows;
    mt19937 gen(rg());
    uniform_real_distribution<> distrib(0.0,1.0);

    for(int i = 0; i < bands; ++i) {
        int p;
        for(p = 0; p < rows; ++p) {
            for(int k = 0; k < size; ++k) {//every node gets a rng number
                r[i](p,k) = -1.0 * (log(distrib(gen)) + log(distrib(gen)));
                c[i](p,k) = -1.0 * (log(distrib(gen)) + log(distrib(gen)));
                beta[i](p,k) = distrib(gen);
            }
        }
        for(p = 0; p < timebits; ++p) {
           this->psi(i,p) = distrib(gen);
        }
    }
    this->target = 200 / timebits;
}

map<double, set<TemporalNode>> TemporalWeightedMinhashGadget::performHashing(AdjListWeight &g, int64_t scale,RealPruneContainer& pc, bool split) {

    map<double, set<TemporalNode>> result = {};
    set<set<TemporalNode>> bins = {};
    cout << "Generating timebits... "<<endl;
    MatrixXd pivots;
    pivots.resize(this->size,this->timebits);
    pivots.setZero();
    int64_t duration = g.getMaxTime() - g.getMinTime();

    int64_t tl;
    int digits;

    for(tl = 0; tl < this->bands; ++tl) {
        for(digits = 0; digits < this->timebits; ++digits) {
            pivots(tl,digits) = (double)duration * this->psi(tl,digits) + (double)g.getMinTime();
        }
    }
    cout << "done. "<<endl;

    tl = (int)(log((double)this->timebits) / log(16.0) + 1.0);

    if (this->timebits < 16) {
        tl = 1;
    }
    digits = (int)(log((double)this->size) / log(16.0) + 1.0);

    int64_t numNeighborhoods = 0;

    //care here
    if (this->mapping.empty() ) {

        for (int i = 0; i < this->bands; ++i) {

            libcuckoo::cuckoohash_map<string, set<TemporalNode>> table = {};
            for (int curT = g.getMinTime(); curT < g.getMaxTime(); ++curT) {
                label110:

                if (!pc.isPruned(curT, this->target - 1)) {

                    //cout<<curT<<endl;
                    numNeighborhoods += g.tNodes(curT).size();
                    set<TemporalNode> neighborhood = g.tNodes(curT);
                    auto var17 = neighborhood.begin();

                    while (true) {//get neighborhood from node

                        TemporalNode tn;
                        std::map<long,double> wMap;
                        do {
                            if (var17 == neighborhood.end()) {
                                curT++;
                                goto label110;
                            }

                            tn = *var17;
                            wMap = g.getWeightedNeighborhood(tn);
                            ++var17;
                        } while (wMap.empty());

                        string sig;
                        string mhglob;

                        auto r2 = this->r;
                        auto c2 = this->c;
                        auto beta2 = this->beta;

                        for (int j = 0; j < this->rows; ++j) {//hash node neighborhood row times


                            bulk::thread::environment env;

                            env.spawn(1, [&mhglob,wMap,r2,beta2,c2,i,j,digits](
                                    bulk::world &world) {

                                int s = world.rank();
                                int p = world.active_processors();
                                string mh;

                                double maxW = std::numeric_limits<double>::min();
                                double minA = std::numeric_limits<double>::infinity();
                                auto numbers = bulk::queue<uint64_t, double>(world);

                                if (s == 0) {
                                    for (auto [k, w]: wMap) {
                                        numbers(std::hash<int>{}(k)%p).send(k, w);
                                    }
                                }

                                world.sync();

                                for (auto [k, w]: numbers) {

                                    if (!isnan(w)) {
                                        if (w > maxW) {
                                            maxW = w;
                                        }
                                        int t = (int) floor(log(w) / r2[i](j, k) + beta2[i](j, k));
                                        double y = exp(r2[i](j, k) * ((double) t - beta2[i](j, k)));
                                        double a = c2[i](j, k) / (y * exp(r2[i](j, k)));
                                        if (a < minA) {
                                            minA = a;
                                            string pad1 = pad((long) k, digits);
                                            string pad2 = turnHex(t);
                                            mh = pad1 + pad2;
                                        }
                                    }
                                }
                                //cout<<mh<<endl;
                                auto report = bulk::queue<double,string>(world);
                                report(0).send(minA,mh);
                                world.sync();

                                map<double,string> temp = {};
                                for (const auto&[a,string]:report) {
                                    temp.insert({a,string});
                                }
                                mh = temp.begin()->second;
                                mhglob = mh;
                            });
                            cout<<sig<<endl;
                            sig = sig.append(mhglob);
                        }

                        long timeHash = 0L;

                        for (int p = 0; p < this->timebits; ++p) {
                            timeHash *= 2L;
                            if ((double) tn.getTime() >= pivots(i, p)) {
                                ++timeHash;
                            }
                        }

                        sig = pad(timeHash, tl).append(sig);
                        set<TemporalNode> nodes;
                        nodes.emplace(tn);

                        if(table.contains(sig)){
                            auto it = table.find(sig);
                            for(auto const &node:it){
                                nodes.emplace(node);
                            }
                        }

                        table.insert_or_assign(sig,nodes);//crashes on std::map too much memory
                    }
                }
            }


            cout<<"Total neighborhood calculations: " << numNeighborhoods<<endl;
            cout<<"Splitting bins..."<<endl;

            auto tableTemp = table;
            auto temp = tableTemp.lock_table();
            auto var36 = temp.begin();

            while(true){
                while(true){
                    string key;
                do{

                    if(var36==temp.end()){
                        cout<<"done!"<<endl;
                        cout<<"Scoring bins... "<<endl;
                        auto var37 = bins.begin();

                        //split bin into a coarray

                        while(var37!=bins.end()){
                            auto bin = *var37;
                            if(scale==0){
                                result.insert_or_assign(scoreBin(bin),bin);
                                ++var37;
                            }else{
                                result.insert_or_assign(scoreBin(bin, scale),bin);
                                ++var37;
                            }
                        }



                        cout<<"done!"<<endl;
                        //i++;
                        goto label139;
                    }
                    auto temp2 = *var36;
                    key= temp2.first;
                    ++var36;

                }while(table.find(key).size()==1);
                      set<TemporalNode> temp3 =  table.find(key);
                      bins.emplace(temp3);
                      if(var36!=temp.end())++var36;
                }
            }
            label139:
            continue;
        }
    }
    return result;
}
