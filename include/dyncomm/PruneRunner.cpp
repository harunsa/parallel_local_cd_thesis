#include "parlay/parallel.h"
//#include "../adj_listweight/adj_listweight.h"
#include <fstream>
#include <cinttypes>
#include "Prune.cpp"
#include <../../eigen/Eigen/Core>

using namespace std;
using namespace Eigen;

class PruneRunner {

public:

    std::array<int64_t,35> PREC_SIZES = {1, 2, 4, 5, 8, 9, 14, 15, 20, 21, 27, 28, 35, 36, 44, 45, 54, 55, 65, 66, 77,
                                        78, 90, 91, 103, 104, 120, 137, 155, 174, 194, 215, 237, 260, 284};
    Prune loadFile(const string& path);
    void precompute(Prune& p, MatrixXd& bounds, bool printEigs);
    MatrixXd getNewMixedBounds(Prune& p ,double estimate, MatrixXd& bounds, bool group);
    void getNewMixedBoundsHelper1(Prune &p,MatrixXd& bounds);
    void getNewMixedBoundsHelper2(Prune &p,double estimate, MatrixXd& bounds,bool group);
    void calcNewVolumes(Prune& p);
    vector<int64_t> findNewIntervals(int64_t start, int64_t end, int64_t maxSize);
    set<set<int64_t>> findNewIntervalsHelp(int64_t start, int64_t end, int64_t maxSize);
    double getPrecomputedBounds(Prune& prune, int64_t start, int64_t end, MatrixXd& bounds);

};

Prune PruneRunner::loadFile(const string& path) {
    return Prune(path);
}

/**
 * precomputes bounds for a given Prune, Lengths are taken from PREC_SIZES
 * 1.Step fills the initial column with the bounds of (i,i)
 * 2.fills the bounds matrix with bounds from (t,end) into index (t,PREC_SIZES-1)
 *
 * Precomputes the composite bound at scales lexp(i) from PREC_SIZES
 *
 * @param p
 * @param bounds
 * @param printEigs
 */
void PruneRunner::precompute(Prune& p, MatrixXd& bounds, bool printEigs) {

    //creates and aggregates Laplacian Graphs taken from the Adjacency Matrix
    p.preaggregate();
    //gets the eigenvalue of the aggregated normalized laplacian -> eigenvalue serves as lowerBound for conductance (spectralBound)
    double lowerBound = p.getBounds(p.startTime,p.endTime);
    //initialize matrix were precomputed bounds are saved
    bounds.resize(p.endTime+1,p.endTime+1);
    //initialize all elements to 0
    bounds.setZero();
    cout<<bounds.rows()<<","<<bounds.cols()<<endl;

    //1.(0,10)
    bounds(p.startTime,p.endTime-p.startTime) = Prune::normalizeConductance(lowerBound, p.endTime - p.startTime + 1);
    //cout<<bounds.rows()<<","<<bounds.cols()<<endl;

    int64_t size = 0;

    bulk::thread::environment env;
    env.spawn(2, [&p, &bounds](bulk::world &world) {
        int s = world.rank();
        int ap = world.active_processors();
        int size = 0;

        auto numbers = bulk::queue<int>(world);
        if (s == 0) {
            for (size = p.startTime; size < p.endTime; ++size) {
                numbers(std::hash<double>{}(size) % ap).send(size);
            }
        }
        world.sync();

        for (auto time: numbers) {
            double lowerBound = p.getBounds(time, time);//lowerBounds for (1,1) (2,2) (3,3)

            //fill initial column 0 with bounds of (i,i)
            bounds(time, 0) = Prune::normalizeConductance(lowerBound, 1); //(1,0) (2,0) (3,0)
             cout<<bounds.rows()<<","<<bounds.cols()<<endl;
             cout << time << "," << 0 << " normCond:" << bounds(time, 0) << " of boundInterval " << time << "," << time<< endl;

        }
        world.sync();
    });

    //compute lambda, scheme follows with increasing lengths according to PREC_SIZES
    auto var11 = PREC_SIZES;
    int64_t var10 = var11.size();

    env.spawn(2, [&var10, &var11, &p, &bounds](bulk::world &world) {
        int s = world.rank();
        int ap = world.active_processors();
        int size = 0;


        auto numbers = bulk::queue<int>(world);


        if (s == 0) {
            for (int var9 = 0; var9 < var10; ++var9) {
                size = var11[var9];
                numbers(std::hash<double>{}(size) % ap).send(size);
            }
        }
        world.sync();

        for (auto time: numbers) {
            if (time != 1) {
                for (int64_t t = p.startTime; t + time - 1 <= p.endTime; t += time) {
                    int64_t end = t + time - 1;
                    double lowerBound = p.getBounds(t, end);
                    //fills cols col number dictates scale/PREC_SIZE
                    bounds(t, time - 1) = Prune::normalizeConductance(lowerBound, time);
                    //cout << t <<","<<size-1 <<endl;
                   cout << t << "," << time - 1 << " normCond:" << bounds(t, time - 1) << " of boundInterval " << t
                         << "," << end << endl;
                }
            }

        }
        world.sync();
    });

    cout<<bounds<<endl;
}

void PruneRunner::getNewMixedBoundsHelper1(Prune &p,MatrixXd& bounds) {
    int64_t totalTime = p.totalTime;
    int64_t maxSize = totalTime;
    int64_t t, full_size,i, half_size, temp, myEnd, range;
    double compositeBound, answer;
    //cout<<bounds<<endl;
      cout<<"getNewMixedBoundsHelper1:"<<endl;

    for (t = p.startTime; t <= p.endTime; ++t) {

        for (i = 0; i < PREC_SIZES.size(); ++i) {
            full_size = PREC_SIZES[i];
            //cout<<"PREC_SIZE/FULL_SIZE:"<<full_size<<endl;
            //cout <<"IF Statement bound = 0, check:"<<t<<","<< full_size - 1<<endl;
            if (full_size != 1 && t + (full_size - 1) <= p.endTime && bounds(t, full_size - 1) <= 0.0) {
                half_size = t + full_size - 1;
                //cout<<"HALF_SIZE:"<<half_size<<endl;
                vector<int64_t> myIntervals = findNewIntervals(t, half_size, maxSize);
                //cout<<"myIntervals"<<endl;
                //for (auto elem:myIntervals) {
                //    cout << elem<<endl;
                //}
                double compositeBound2 = 0.0;
                compositeBound = 0;
                temp = 0;
                int64_t interval;
                //sum compositeBound based on intervals
                for (auto it = myIntervals.begin(); it != myIntervals.end(); temp += interval) {
                    interval = *it;
                    myEnd = t + temp;
                    range = myEnd + (interval - 1);
                    compositeBound2 = getPrecomputedBounds(p, myEnd, range, bounds);
                    answer = p.getIntervalWeight(myEnd, range, t, half_size);
                    //cout << "getPrecompBounds:" << myEnd << "," << range << endl;

                    //cout << "getIntervalWeight:" << myEnd << "," << range << "/" << t << "," << half_size<< endl;

                    compositeBound += answer * compositeBound2;
                    ++it;
                }
                compositeBound = Prune::normalizeConductance(compositeBound, full_size);
                bounds(t, full_size - 1) = compositeBound;
 //               cout << "normalizeCond:" << compositeBound << "," << full_size  << endl;
 //               cout << "insert into:" << t << "," << full_size - 1 << endl;
                //cout << "-------------------END OF INSERT------------------------" << endl;
            }
        }
    }
    //cout<<bounds<<endl;
}

void PruneRunner::getNewMixedBoundsHelper2(Prune &p,double estimate, MatrixXd& bounds,bool group){
    int64_t totalTime = p.totalTime;
    int64_t maxSize = totalTime;
    int64_t t2 = 0;
    int64_t full_size,half_size,myEnd,temp;
    double compositeBound, answer2;
    cout<<"getNewMixedBoundsHelper2:"<<endl;

    for (int64_t t = p.startTime; t < p.endTime; ++t) {
        int64_t s, interval, subEnd, i,range, tmp, dur;
        double ratio;
        //cout<<"Time:"<<t<<endl;

        for (i = 1; i < PREC_SIZES.size(); ++i) {
            full_size = PREC_SIZES[i];
            //cout<<"PREC_SIZE/FULL_SIZE:"<<full_size<<endl;

            if (t + (full_size - 1) > p.endTime) {
                //cout<<"IF-statement: t + (full_size - 1)>p.endTime:"<<t + (full_size - 1)<<endl;
                //cout<<"BREAK OUT OF LOOP"<<endl;
                break;
            }

            half_size = PREC_SIZES[i - 1];
            //cout<<"HALF_SIZE:"<<half_size<<endl;

            if (half_size != 1) {
                //cout<<"IF HALF_SIZE!=1:"<<half_size<<endl;
                double answer = 0.0;
                pair<int64_t, int64_t> sp2 = {};
                range = t + full_size - 1;
                myEnd = t + half_size - 1;
                vector<int64_t> myIntervals = {};
                double r;
                //calculates answer2 & begins at 9,14
                //goes in here if intervals are atleast 4 , check next statement

                /**-------------- CALC GroupBound PREFIX--------------  */
                if (group && full_size - half_size >= 4) {
                    //cout<<"range,myEnd:"<<range<<","<<myEnd<<endl;
                    //cout<<"time:"<<t<<endl;

                    //cout << "IF - Statement: full_size - half_size>=4:"<<full_size<<","<<half_size<<",,"<< full_size - half_size<<endl;
                    sp2 = {t, range};
                    vector<double> fv = {};

                    if (p.preComputedVolumes.find(sp2) != p.preComputedVolumes.end()) {
                        fv = p.preComputedVolumes.find(sp2)->second;
                    } else {
                        cout << "ERROR! Missing precomputed volume in getMixedBounds!" << endl;
                        cout << "\tsp2 = " << sp2.first << "," << sp2.second << endl;
                        exit(-3);
                    }

                    vector<pair<int64_t, int64_t>> intervals = {};
                    myIntervals = findNewIntervals(t, myEnd, maxSize);
                    tmp = 0;

                    for (auto it = myIntervals.begin(); it != myIntervals.end(); tmp += s) {
                        s = *it;
                        interval = t + tmp;
                        subEnd = interval + (s - 1);
                        pair<int64_t, int64_t> p1 = {interval, subEnd};
                        intervals.emplace_back(p1);
                        ++it;
                    }
                    double ratio_thirdLoop;
                    double lamb;

                    for (auto it = intervals.begin(); it != intervals.end(); ++it) {
                        pair<int64_t, int64_t> sp = *it;
                        vector<double> iv = {};
                        if (p.preComputedVolumes.find(sp) != p.preComputedVolumes.end()) {
                            iv = p.preComputedVolumes.find(sp)->second;
                        } else {
                            cout << "ERROR! Missing inner precomputed volume in getMixedBounds!" << endl;
                            cout << "\tsp2 = " << sp.first << "," << sp.second << endl;
                            exit(-3);
                        }
                        lamb = 0.0;
                        if (p.halfLambdas.find(sp) != p.halfLambdas.end()) {
                            lamb = p.halfLambdas.find(sp)->second; //spectralBound
                            //cout << "getPrecompBounds:" << sp.first << "," << sp.second << endl;
                        } else {
                            cout << "ERROR! Missing precomputed volume in getMixedBounds!" << endl;
                            cout << "\tsp2 = " << sp.first << "," << sp.second << endl;
                            exit(-3);
                        }

                        ratio_thirdLoop = 1.0;

                        for (int64_t v = 0; v <= p.maxNode; ++v) {
                            r = iv[v] / fv[v];
                            if (r < ratio_thirdLoop) {
                                ratio_thirdLoop = r;
                            }
                        }
                        //cout << "getIntervalWeight:" << sp.first << "," << sp.second << "/" << sp2.first << ","<< sp2.second<< endl;
                        answer += ratio_thirdLoop * lamb;
                    }
                    answer = Prune::normalizeConductance(answer, full_size);
                    //cout << "normalizeCond:" << answer << "," << full_size << endl;
                    //cout << "answer :" << answer << endl;
                }
                /** ---------------INSERT INTO PRUNING GROUP-------------------  */
                if (estimate < answer) {//for 0 entries this is false
                    //cout << "IF STATEMENT estimate < answer: " << estimate << "," << answer << endl;

                    for (range = half_size; range < full_size - 1; ++range) {
                        if (bounds(t, range) < answer) {
                            bounds(t, range) = answer;
                            //cout << t << "," << range << ":" << bounds(t, range) << endl;
                        }
                    }
                } else {
                /** --------------------------COMPBOUND FOR SINGLE INTERVAL -------------------------  */

                    //cout << "ELSE STATEMENT estimate < answer: " << estimate << "," << answer << endl;

                    for (range = half_size; range < full_size - 1; ++range) {
                        //cout << "Range:" <<range<< endl;
                        if (bounds(t, range) <= 0.0) {
                           // cout << "IF-Statement bounds(t, range) <= 0.0:" <<bounds(t, range)<< endl;
                            myEnd = t + range;
                            myIntervals = findNewIntervals(t, myEnd, maxSize);
                            answer = 0.0;
                            ratio = 0;
                            dur = 0;

                            for (auto it = myIntervals.begin(); it != myIntervals.end(); dur += interval) {
                                interval = *it;
                                int64_t subStart = t + dur;
                                subEnd = subStart + (interval - 1);
                                ratio = getPrecomputedBounds(p, subStart, subEnd, bounds);
                                r = p.getIntervalWeight(subStart, subEnd, t, myEnd);
                                //cout << "getPrecomputedBounds:" << subStart << "," << subEnd << endl;
                                //cout << "getIntervalWeight:" << subStart << "," << subEnd << "/" << t << ","<< myEnd << endl;
                                answer += r * ratio;
                                ++it;
                            }
                            answer = Prune::normalizeConductance(answer, range + 1);
                            bounds(t, range) = answer;
                            //cout << "normalizeCond:" << answer << "," << range + 1 << endl;
                            //cout << "insert into:" << t << "," << range << endl;
                            //cout << "END OF INSERT----------" << endl;
                        }
                    }
                }
            }
            //cout<<"endofloop"<<endl;
        }
//      cout<<bounds<<endl;
        //cout<<"LAST PART OF LOOP"<<endl;
        half_size = PREC_SIZES[i-1];
        int64_t mid = t+half_size -1;
        temp = p.endTime;
        full_size = temp -t +1;
        bool final_group = false;
        double lamb;
        //cout<<"halfsize:"<<half_size<<endl;
        //cout << "IF-Statement temp - mid >3:" <<temp<<"-"<<mid<<"="<< temp-mid  <<";mid=("<<t<<"+"<<half_size<<"-1"<<")"<<endl;


        /** ------------------CALC PREFIX + PRUNING GROUP------------------------   */

        if(group &&temp - mid >3) {
            pair<int64_t, int64_t> sp2 = {t, temp};
            vector<double> fv = {};
            if (p.preComputedVolumes.find(sp2) != p.preComputedVolumes.end()) {
                fv = p.preComputedVolumes.find(sp2)->second;
            } else {
                cout << "ERROR! Missing precomputed volume in getMixedBounds!" << endl;
                cout << "\tsp2 = " << sp2.first << "," << sp2.second << endl;
                exit(-3);
            }

            vector<pair<int64_t, int64_t>> intervals = {};
            vector<int64_t> intervalSizes = findNewIntervals(t, mid, maxSize);
            int64_t tmp_secondLoop = 0;

            for (auto it = intervalSizes.begin(); it != intervalSizes.end(); ) {
                tmp = *it;
                dur = t + tmp_secondLoop;//tmp vs tmp_secondLoop
                interval = dur + (tmp - 1);
                pair<int64_t, int64_t> p1 = {dur, interval};
                intervals.emplace_back(p1);
                tmp_secondLoop+= tmp;
                ++it;
            }

            answer2 = 0.0;
            auto var57 = intervals.begin();

            while (true) {
                //cout<<"While-true loop:"<<endl;
                if (var57 == intervals.end()) {
                   // cout<<"var57=intervals.end()"<<endl;
                    answer2 = Prune::normalizeConductance(answer2, full_size);
                   // cout << "normalizeCond:" << answer2 << "," << full_size << endl;

                    if (estimate >= answer2) {
                      //  cout<<"IF-statement: estimate >= answer2:"<<estimate<<endl;
                       // cout<<"BREAK OUT OF LOOP"<<endl;
                        break;
                    }

                    final_group = true;
                    dur = half_size;

                    while(true){
                       // cout <<"In 2nd While-true loop:" <<endl;
                        if(t+dur>p.endTime){
                      //      cout << "break out of first While-true " << endl;
                            goto label202;
                        }
                        ++t2;
                        if(bounds(t,dur)<answer2){
                            //++full_size;
                            bounds(t,dur) = answer2;
                       //     cout << "normalizeCond:" << answer2 << "," << full_size << endl;
                       //     cout << "insert into:" << t << "," << dur << endl;
                        }
                        ++dur;
                    }
                }

                pair<int64_t,int64_t> sp = *var57;

                vector<double> iv = {};

                if(p.preComputedVolumes.find(sp)!=p.preComputedVolumes.end()){
                    iv = p.preComputedVolumes.find(sp)->second;
                }else{
                    cout << "ERROR! Missing inner precomputed volume in getMixedBounds!" << endl;
                    cout << "\tsp2 = " << sp.first << "," << sp.second << endl;
                    exit(-3);
                }

                lamb = 0.0;
                if(p.halfLambdas.find(sp)!=p.halfLambdas.end()){
                    lamb = p.halfLambdas.find(sp)->second;
                  //  cout << "getPrecomputedBounds:" << sp.first << "," << sp.second << endl;

                }else{
                    cout << "ERROR! Missing precomputed volume in getMixedBounds!" << endl;
                    cout << "\tsp2 = " << sp.first << "," << sp.second << endl;
                    exit(-3);
                }

                ratio = 1.0;

                for(int v = 0; v <= p.maxNode; ++v) {
                    double r = iv[v] / fv[v];
                    if (r < ratio) {
                        ratio = r;
                    }
                }
                //cout << "getIntervalWeight:" << sp.first << "," << sp.second << "/" << sp2.first << "," << sp2.second;
                answer2 += ratio * lamb;
                ++var57;
            }
            label202:;
        }
        /** --------------------------COMPBOUND FOR SINGLE INTERVAL -------------------------  */
        if(!final_group){
         //   cout<<"FINAL_GROUP"<<endl;
          //  cout << "IF Statement: Final_Group false"<< endl;

            for(range = half_size; t + range<=p.endTime;++range){
                if(bounds(t, range) <= 0.0){
                    myEnd = t +range;
                    vector<int64_t> myIntervals = findNewIntervals(t,myEnd,maxSize);
                    compositeBound = 0.0;
                    tmp = 0;

                    for(auto it = myIntervals.begin();it!=myIntervals.end();tmp+=s){
                        s=*it;
                        interval = t+tmp;
                        subEnd = interval + (s-1);
                        lamb = getPrecomputedBounds(p,interval,subEnd,bounds);
                        ratio = p.getIntervalWeight(interval,subEnd,t,myEnd);
             //           cout << "getPrecomputedBounds:" << interval << "," << subEnd << endl;
              //          cout << "getIntervalWeight:" << interval << "," << subEnd << "/" << t << ","<< myEnd << endl;
                        compositeBound += ratio *lamb;
                        ++it;
                    }
                    compositeBound = Prune::normalizeConductance(compositeBound,range+1);
                    bounds(t,range) = compositeBound;
            //        cout << "normalizeCond:" << compositeBound << "," << range+1  << endl;
           //         cout << "insert into:" << t << "," << range << endl;
            //        cout << "---------------------------END OF INSERT" << endl;
                }
            }
        }
    }
}

MatrixXd PruneRunner::getNewMixedBounds(Prune &p, double estimate, MatrixXd &bounds, bool group) {

    int64_t totalTime = p.totalTime;
    int64_t maxSize = totalTime;
    p.getNodeWeights();
    cout << "Calculating Volumes... " << endl;
    calcNewVolumes(p);
    cout << "DONE!" << endl;

    getNewMixedBoundsHelper1(p,bounds);
    getNewMixedBoundsHelper2(p,estimate,bounds,group);

    cout<<"Done with getMixedBounds, grouping pruned " <<"t_firstLoop"<< " intervals, changed "<<"full_size"<<" bounds."<<endl;
    return bounds;
}


//Bounds Matrix is not complete during precomputation, fills gaps, no matrix computation
//calculates newVols based on precomputedVols, computes everyVol.

void PruneRunner::calcNewVolumes(Prune &p) {


//    cout<<"calcVolumes for Bounds"<<endl;
    vector<double> endVols = {};
    int64_t i = 0;
    pair<int64_t, int64_t> span1 = {};

    for (i = p.startTime; i <= p.endTime; ++i) {
        span1 = {i, i};
        //compute vols of single timestamps (i,i)
        vector<double> vols = p.computeVolumes(i, i);
        if (i == p.endTime) endVols = vols;
        p.preComputedVolumes.insert_or_assign(span1, vols);
    }

    //iterate over prec_sizes, it determines the distance in the interval
    for (i = 0; i < PREC_SIZES.size(); ++i) {
        //curr size
        int64_t full_size = PREC_SIZES[i];
        if (full_size != 1) {
            int64_t half_size = PREC_SIZES[i - 1];

            for (int64_t t = p.startTime; t + full_size - 1 <= p.endTime; ++t) {
                int64_t end = t + full_size - 1;
                pair<int64_t, int64_t> span = {t, end};//0,1
                pair<int64_t, int64_t> s1 = {t, t + half_size - 1};//0,0 will always exist, given via pattern
                pair<int64_t, int64_t> s2 = {t + half_size, end};//1,1 instances where it exists but not always
//                cout << "Whole Span newVols: " << span.first <<","<< span.second << endl;//(0,7)
//                cout << "First Span newVols: " << s1.first <<","<< s1.second << endl;//(0,4)
//               cout << "Second Span newVols: " << s2.first <<","<< s2.second << endl;//(5,7)
                vector<double> v1 = p.preComputedVolumes.find(s1)->second;//(0,4)
                vector<double> v2 = {};

                if (p.preComputedVolumes.find(s2) !=p.preComputedVolumes.end()) {
                    v2 = p.preComputedVolumes.find(s2)->second;//(5,7) doesnt exist
                } else {
                    v2.resize(v1.size());//(5,7) init. vector

                    int64_t st = t + half_size;

                    while (true) {
                        //1.st=5 end=7 ;  2.st=10 end=7
                        if (st > end) {
                            p.preComputedVolumes.insert_or_assign(s2, v2);//insert vol(5,9) into (5,7)
//                            cout << "insert s2: " << s2.first <<","<< s2.second  <<endl;
                            //cout<<"-------------------------"<<endl;

                            break;
                        }
                        vector<double> v3 (v1.size());//init v3
                        int64_t thisSize = 1;
                        int64_t m;
                        for (m = i - 1; m >= 0; --m) {
                            thisSize = PREC_SIZES[m];

                            pair<int64_t, int64_t> s3{st, st + thisSize - 1};//(5,9) fill v3 with it
                            if (p.preComputedVolumes.find(s3) != p.preComputedVolumes.end()) {
                                v3 = p.preComputedVolumes.find(s3)->second;//(5,9) exists
 //                               cout<<"INSERT IF doesnt exist:"<<s3.first<<","<<s3.second<<endl;
                                break;
                            }
                        }

                        for (m = 0; m < v1.size(); ++m) {
                            v2[m] += v3[m];//add v3(5,9) to v2(5,7)
                        }
                        st += thisSize;
                    }
                }
                vector<double> sum(v1.size());
                for (int j = 0; j < v1.size(); ++j) {
                    sum[j] = v1[j] + v2[j];//(0,4)+(5,9) vol
                }
                //
                p.preComputedVolumes.insert_or_assign(span, sum);
//                cout << "insert summ: " << span.first <<","<< span.second << endl;
//                cout<<"-------------------------"<<endl;

            }
        }
    }
    //fill last col
    for (i = p.endTime - 1; i >= p.startTime; --i) {
        span1 = {i, p.endTime};
        pair<int64_t, int64_t> now = {i, i};
        vector<double> v1 = p.preComputedVolumes.find(now)->second;
        vector<double> sum(v1.size());
        for (int j = 0; j < v1.size(); ++j) {
            sum[j] = v1[j] + endVols[j];
        }
        p.preComputedVolumes.insert_or_assign(span1, sum);
//        cout << "insert sum: " << span1.first <<","<< span1.second << endl;
//        cout<<"-------------------------"<<endl;

        endVols = sum;
    }
}

vector<int64_t> PruneRunner::findNewIntervals(int64_t start, int64_t end, int64_t maxSize) {

    vector<int64_t> intervals = {};
    set<set<int64_t>> result = findNewIntervalsHelp(start,end, maxSize);

    for(const auto& elem:result){
        int64_t a = *elem.begin();
        int64_t b = *elem.rbegin();
        intervals.emplace_back(b - a + 1);
    }
    return intervals;
}

set<set<int64_t>> PruneRunner::findNewIntervalsHelp(int64_t start, int64_t end, int64_t maxSize) {

    set<set<int64_t>> ranges = {};
    int64_t cover_start = -1;
    int64_t cover_end = -1;

    int64_t i;
    for (i = PREC_SIZES.size() - 1; i >= 0; --i) {
        int64_t size = PREC_SIZES[i];
        if (size <= maxSize) {
            for (int64_t t = start; t <= end; ++t) {
                if (t % size == 0) {
                    int64_t st = t;

                    for (int64_t et = t + (size - 1); et <= end; et += size) {
                        if (cover_start == -1) {
                            cover_start = t;
                        }
                        cover_end = et;
                        set<int64_t> pair = {};
                        pair.insert(st);
                        pair.insert(et);
                        ranges.insert(pair);
                        st += size;
                    }

                    if (cover_start != -1) {
                        break;
                    }
                }
            }

            if (cover_start != -1) {
                break;
            }
        }
    }
    if (cover_start != -1 && i > 0 && cover_start > start) {
        auto temp = findNewIntervalsHelp(start, cover_start - 1, PREC_SIZES[i - 1]);
        for (const auto &elem: temp) {
            ranges.insert(elem);
        }
    }

    if (cover_end != -1 && i > 0 && cover_end < end) {
        auto temp = findNewIntervalsHelp(cover_end + 1, end, PREC_SIZES[i - 1]);
        for (const auto &elem: temp) {
            ranges.insert(elem);
        }
    }
    return ranges;
}

double PruneRunner::getPrecomputedBounds(Prune& prune, int64_t start, int64_t end, MatrixXd& bounds) {
    if(bounds(start,end-start)>0.0){
        return bounds(start,end-start);
    }else{
        double b = prune.getBounds(start,end);
        bounds(start,end-start) = Prune::normalizeConductance(b,end-start+1);
        return b;
    }
}


