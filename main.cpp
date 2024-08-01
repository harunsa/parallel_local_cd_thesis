#include "../../dyncomm/DynCommDriver.cpp"
using namespace Spectra;
using namespace std;
#include <bulk/bulk.hpp>
#include <bulk/backends/thread/thread.hpp>


int main() {
    //synthV1000T30
    //""
    DynCommDriver::multiscale_driver("synthV1000T30");
    return 0;
}

//    bulk::thread::environment env;
//
//    /**
//     * take a vector
//     * create a que
//     */
//     vector<int> test = {1,2,3,4,5,6,7,8,9};
//     vector<int> newtest = {};
//
//
//    env.spawn(3, [&test,&newtest](bulk::world &world) {
//        int s = world.rank();
//        int p = world.active_processors();
//
//        // We create a queue that holds the individual words of the text
//        auto numbers = bulk::queue<int>(world);
//
//        if (s == 0) {
//
//            for (int i = 0; i < test.size(); ++i) {
//                //std::hash<int>{}(test[i])%p
//                numbers(std::hash<int>{}(test[i])%p).send(test[i]);
//            }
//        }
//
//        world.sync();
//
//        //everyone creates their own little vector, where the values are increased by 1
//        //increment numbers in the que by 1
//        auto ownVec = vector<int>{};
//        for(auto number:numbers){
//            ownVec.emplace_back(number+1);
//            world.log("%d/,owned by,%d",number,s);
//
//        }
//
//        // We send the results back to the master
//        auto report = bulk::queue<int>(world);
//        for (auto number : ownVec) {
//            report(0).send(number);
//        }
//        world.sync();
//
//
//                  for (auto elem : report) {
//                      newtest.emplace_back(elem);
//                  }
//
//        }
//    );
//
//
//
//    for (int i : test) {
//        std::cout << i << " ";
//    }
//    std::cout << "-------------" <<endl;
//
//    for (int i : newtest) {
//        std::cout << i << " ";
//    }


//                auto count=0;
//                for (auto [n2,weight]:pairs) {
//                    world.log("(%d/,%d-%f),owned by,%d",innerTbl.first,n2,weight,s);
//                    count++;
//                }
//                world.log("%d,owned by,%d",count,s);

//                world.sync();



//////////////////////////////////////////////////


//        auto ownVec = vector<int>{};
//        for(auto number:numbers){
//            ownVec.emplace_back(number);
//            world.log("%d,owns, %d",s,number);
//        }



/**
*

        bulk::thread::environment env;
        int end = this->bands;
        auto target1 = this->target;

        env.spawn(6, [&g,&end,&target1,&pc,&numNeighborhoods](bulk::world &world) {
            int s = world.rank();
            int p = world.active_processors();


            auto numbers = bulk::queue<int,int>(world);
            if (s == 0) {
                for (int i=0; i< end; ++i) {
                    for (int curT = g.getMinTime(); curT < g.getMaxTime(); ++curT) {
                        if (!pc.isPruned(curT, target1 - 1)) {
                            numbers(std::hash<int>{}(i) % p).send(i,curT);
                        }
                    }

                }
            }
            world.sync();

            for (auto [i,curT]: numbers) {

                numNeighborhoods += g.tNodes(curT).size();
                set<TemporalNode> neighborhood = g.tNodes(curT);
                auto var17 = neighborhood.begin();
                while (true) {//get neighborhood from node
                    TemporalNode tn;
                    std::map<long,double> wMap;
                    do {
                        if (var17 == neighborhood.end()) {
                            curT++;
                        }

                        tn = *var17;
                        wMap = g.getWeightedNeighborhood(tn);
                        ++var17;
                    } while (wMap.empty());

            }

        });





 *
 *
 *
 *
*/