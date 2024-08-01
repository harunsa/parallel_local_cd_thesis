#include "adj_listweight.h"
#include "parlay/parallel.h"
#include <cmath>
#include <fstream>
#include <cinttypes>

typedef libcuckoo::cuckoohash_map<uint64_t, std::vector<std::pair<uint64_t,double>>> Edge;
typedef libcuckoo::cuckoohash_map<uint64_t, Edge> NestedMap;

/**
 * Checks if the given edge exists in the graph.
 * @param source node of the edge
 * @param destination node of the edge
 * @param time timestamp of the edge
 * @return true if the edge was found
 */
bool AdjListWeight::findEdge(uint64_t source, uint64_t destination, uint64_t time) {

    bool flag = false;
    edges.find_fn(time,
                  [&destination, &flag, &source](Edge &e) {
                      e.find_fn(source,
                                [&flag, &destination](std::vector<std::pair<uint64_t, double>> &d) {
                                    flag = (std::find_if(d.begin(), d.end(),
                                                         [&destination](std::pair<uint64_t, double> &pair) {
                                                             return pair.first == destination;
                                                         }) != d.end());
                                });
                  });
    return flag;

}

//has a race condition issue when insertions/deletions are occurring
/**
 * Checks if the an edge between @p source and @p destination exists in the given range of timestamps.
 * @param source node of the edge
 * @param destination node of the edge
 * @param start of the range inclusive
 * @param end end of the range exclusive
 * @return true if at least one edge was found
 * @overload
 */
bool AdjListWeight::findEdge(uint64_t source, uint64_t destination, uint64_t start, uint64_t end){
    bool flag = false;
    auto uniqueTimesMap = genUniqueTimeMap(start, end);
    for (auto time : uniqueTimesMap) {
        flag = findEdge(source, destination, time.second);
        if (flag) break;
    }
    return flag;
}

/**
 * Inserts an edge into the given nested cuckoo map @p map.
 * @param source node of the edge
 * @param destination node of the edge
 * @param time timestamp of the edge
 * @param map
 */
void AdjListWeight::insertEdgeDirected(uint64_t source, std::pair<uint64_t,double> destination, uint64_t time, NestedMap &map) {

    if (map.contains(time)) {
        map.update_fn(time,
                      [&source, &destination](Edge &e) {
                          e.upsert(source,
                                   [&destination](std::vector<std::pair<uint64_t,double>> &d) { d.push_back(destination);},
                                   std::vector<std::pair<uint64_t,double>>{destination});
                      });
    } else {
        Edge e;
        e.insert(source, std::vector<std::pair<uint64_t,double>>{destination});
        map.insert(time, e);
    }
}

/**
 * Checks if the given edge to be inserted is already in the graph. If not, calls insertEdgeDirected twice to insert the
 * the edge in both directions (@p source -> @p destination and @p destination -> @p source).
 * @see AdjList::insertEdgeDirected
 * @param source node of the edge
 * @param destination node of the edge
 * @param time timestamp of the edge
 */
void AdjListWeight::insertEdgeUndirected(uint64_t source, std::pair<uint64_t,double> destination, uint64_t time) {
    //filters out duplicates
    if (findEdge(source, destination.first, time)) return;

    //insert edges from source
    insertEdgeDirected(source, destination, time, edges);
    //insert edges from destination
    insertEdgeDirected(destination.first, std::pair<uint64_t ,double>{source,destination.second}, time, edges);
}

/**
 * Deletes the given edge from the graph. Functions similar to insertEdgeDirected.
 * Also removes keys if their values are empty.
 * @param source node of the edge
 * @param destination node of the edge
 * @param time timestamp of the edge
 */
void AdjListWeight::deleteEdgeDirected(uint64_t source, uint64_t destination, uint64_t time) {
    bool isDestinationEmpty = false;
    bool isEdgeEmpty = false;

    if (edges.contains(time)) {

        edges.update_fn(time, [&isDestinationEmpty, &source, &destination](Edge &e) {
            e.update_fn(source, [&isDestinationEmpty, &destination](std::vector<std::pair<uint64_t,double>> &d) {
                d.erase(std::find_if(d.begin(), d.end(), [&destination](std::pair<uint64_t, double> &pair){
                    return pair.first==destination;
                }));
                if (d.empty()) isDestinationEmpty = true;
            });
        });
        //delete source node if it has no edges (destinations)
        if (isDestinationEmpty) {
            edges.find_fn(time, [&isEdgeEmpty, &source](Edge &e) {
                e.erase(source);
                //causes performance issues
                if (e.empty()) isEdgeEmpty = true;
            });
        }
        //delete timestamp if edges is empty
        if (isEdgeEmpty) {
            edges.erase(time);
            uniqueTimestamps.erase(time);
        }
    }
}

/**
 * Works similar to insertEdgeUndirected.
 * @param source node of the edge
 * @param destination node of the edge
 * @param time timestamp of the edge
 */
void AdjListWeight::deleteEdgeUndirected(uint64_t source, uint64_t destination, uint64_t time) {
    //check if edge to be deleted exists
    if (!findEdge(source, destination, time)) return;

    //delete edges from source
    deleteEdgeDirected(source, destination, time);
    //delete edges from destination
    deleteEdgeDirected(destination, source, time);
}

/**
 * Prints all timestamps and the edge contained in them.
 */
void AdjListWeight::printGraph() {
    std::cout << "Printing all edges:" << std::endl << std::endl;
    uint64_t count = 0;

    auto lt = edges.lock_table();

    for (const auto &innerTbl: lt) {
        Edge edgeData = innerTbl.second;
        auto lt2 = edgeData.lock_table();
        printf("Time %" PRIu64 " contains edges\n", innerTbl.first);

        for (const auto &vector: lt2) {
            for (auto &edge: vector.second) {
                printf("    - between: %" PRIu64 " and %" PRIu64 " with weight %d"  " at time %" PRIu64 "\n", vector.first, edge.first, (int)edge.second, innerTbl.first);
                count++;
            }
        }
        std::cout << std::endl;
    }
    std::cout << "Total number of edges: " << count << std::endl;
    printf("%" PRIu64 "\n", edges.size());
}

void AdjListWeight::printNodeGraph() {
    std::cout << "Printing all Nodes with volumes:" << std::endl << std::endl;
    auto lt = nodes.lock_table();

    for (const auto &innerTbl: lt) {
        auto edgeData = innerTbl.second;
        auto lt2 = edgeData.lock_table();
        printf("Time %" PRIu64 " contains Nodes\n", innerTbl.first);

        for (const auto &vector: lt2) {
            printf("    - Node: %" PRIu64  " with Volume %f"  " at time %" PRIu64 "\n", vector.first,vector.second,innerTbl.first);

        }
        std::cout << std::endl;
    }
}

/**
 * Helper function to organize the read data.
 * @param sourceVector container for all read sources with the same command.
 * @param destinationVector container for all read destinations with the same command.
 * @param timeVector container for all read timestamps with the same command.
 * @param uniqueTimes container for all read unique timestamps with the same command.
 * @param source value to be inserted
 * @param destination value to be inserted
 * @param time value to be inserted
 */
void AdjListWeight::fileReaderHelper(std::vector<uint64_t> &sourceVector, std::vector<uint64_t> &destinationVector, std::vector<double> &weightVector,
                               std::vector<uint64_t> &timeVector, std::set<uint64_t> &uniqueTimes,
                               uint64_t source, uint64_t destination, double weight, uint64_t time){
    sourceVector.push_back(source);
    destinationVector.push_back(destination);
    weightVector.push_back(weight);
    timeVector.push_back(time);
    uniqueTimes.insert(time);
}

/**
 * Helper function to organise unique timestamps. Inserts or removes them from uniqueTimestamps.
 * @param uniqueTimesMap container for unique times mapped to a key that iterates up from 0
 * @param uniqueTimes values to be inserted
 * @param insert determines whether to insert or delete from uniqueTimestamps
 */
void AdjListWeight::uniqueTimesHelper(std::unordered_map<uint64_t, uint64_t> &uniqueTimesMap, std::set<uint64_t> &uniqueTimes, bool insert){
    int i = 0;
    for (const uint64_t &time: uniqueTimes) {
        uniqueTimesMap[i] = time;
        if (insert) uniqueTimestamps.insert(time);
        i++;
    }
}


/**
 * Reads and extracts data from the file and calls functions to use the data on the graph.
 * @param path input file
 */
void AdjListWeight::addFromFile(const std::string &path) {
    std::ifstream file(path);
    if (file.is_open()) {
        std::string command;
        uint64_t source, destination,time;
        double weight;
        std::vector<uint64_t> sourceAdds{}, destinationAdds{}, timeAdds{};
        std::vector<double> weightAdds{};
        std::vector<uint64_t> sourceDels{}, destinationDels{}, timeDels{};
        std::vector<double> weightDels{};
        std::set<uint64_t> uniqueTimesAdd, uniqueTimesDel;
        std::unordered_map<uint64_t, uint64_t> uniqueTimesAddMap, uniqueTimesDelMap;

        while (file >> command >> source >> destination >> time >> weight) {
            if (command == "add") {
                fileReaderHelper(sourceAdds, destinationAdds, weightAdds, timeAdds, uniqueTimesAdd, source, destination, weight, time);
            }
            if (command == "delete") {
                fileReaderHelper(sourceDels, destinationDels, weightDels, timeDels, uniqueTimesDel, source, destination, weight, time);
            }
        }
        file.close();

        uniqueTimesHelper(uniqueTimesAddMap, uniqueTimesAdd, true);
        uniqueTimesHelper(uniqueTimesDelMap, uniqueTimesDel, false);


        //Create new hash map, keys are timestamps,values are Edges (source, <destination>).
        //This is then filled by sortBatch function.
        libcuckoo::cuckoohash_map<uint64_t, Edge> groupedDataAdds, groupedDataDels;

        sortBatch(sourceAdds, destinationAdds, weightAdds, timeAdds,  groupedDataAdds);
        batchOperationParlay(true, groupedDataAdds, uniqueTimesAddMap);

        sortBatch(sourceDels, destinationDels, weightDels, timeDels, groupedDataDels);
        batchOperationParlay(false, groupedDataDels, uniqueTimesDelMap);

        auto f = [](uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
            //printf("    - RangeQueryTest between: %" PRIu64 " and %" PRIu64 " at time %" PRIu64 "\n", b, c, a);
        };

        //rangeQueryToDestParlay(2013, 2019, f);
    }
    nodes = getNodeGraph();
}

/**
 * Iterated through @p groupedData and calls insertEdgeUndirected or deleteEdgeUndirected accordingly.
 * @param insert dictates whether to insert or delete the given data
 * @param groupedData Nested map of edges that are to be inserted/deleted
 */
void AdjListWeight::batchOperation(bool insert, NestedMap &groupedData) {
    auto t1 = std::chrono::high_resolution_clock::now();
    auto lt = groupedData.lock_table();

    for (const auto &innerTbl: lt) {
        Edge edgeData = innerTbl.second;
        auto lt2 = edgeData.lock_table();

        for (const auto &vector: lt2) {
            for (auto &edge: vector.second) {
                if (insert) insertEdgeUndirected(vector.first, edge, innerTbl.first);
                else deleteEdgeUndirected(vector.first, edge.first, innerTbl.first);
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "batchOperationCuckoo has taken " << ms_int.count() << "ms\n";
}

/**
 * Works similar to batchOperation. Iterates in parallel.
 * @param insert dictates whether to insert or delete the given data
 * @param groupedData Nested cuckoo map of edges that are to be inserted/deleted
 * @param uniqueTimesMap helper parameter to iterate through @p groupedData
 */
void AdjListWeight::batchOperationParlay(bool insert, NestedMap &groupedData, std::unordered_map<uint64_t, uint64_t> uniqueTimesMap) {
    auto t1 = std::chrono::high_resolution_clock::now();
    auto lt = groupedData.lock_table();

    parlay::parallel_for(0, uniqueTimesMap.size(), [&](uint64_t i) {
        uint64_t time = uniqueTimesMap[i];

        Edge innerTbl = lt.find(time)->second;
        auto lt2 = innerTbl.lock_table();

        for (const auto &vector: lt2) {
            for (auto &edge: vector.second) {
                if (insert) insertEdgeUndirected(vector.first, edge, time);
                else deleteEdgeUndirected(vector.first, edge.first, time);
            }
        }
    });
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "addBatchCuckooParlay has taken " << ms_int.count() << "ms\n";
}

/**
 * Organises the read data and calls insertEdgeDirected to insert them into @p groupedData.
 * @param sourceAdds list of source nodes
 * @param destinationAdds list of destination nodes
 * @param timeAdds list of timestamps
 * @param groupedData Nested cuckoo map and container for the organised data
 */
void AdjListWeight::sortBatch(const std::vector<uint64_t> &sourceAdds, const std::vector<uint64_t> &destinationAdds,
                       const std::vector<double> &weightAdds, const std::vector<uint64_t> &timeAdds, NestedMap &groupedData) {
    auto t1 = std::chrono::high_resolution_clock::now();

    // Determine the number of iterations
    size_t numIterations = timeAdds.size();

    // Group edges by time using hash map
    for (size_t i = 0; i < numIterations; ++i) {
        uint64_t source = sourceAdds[i];
        uint64_t destination = destinationAdds[i];
        double weight = weightAdds[i];
        uint64_t time = timeAdds[i];

        insertEdgeDirected(source, std::pair<uint64_t ,double> {destination,weight}, time, groupedData);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "sortBatch has taken " << ms_int.count() << "ms\n";
}

/**
 * Prints the edges in @p groupedData similar to printGraph.
 * @param groupedData
 */
void AdjListWeight::printGroupedData(libcuckoo::cuckoohash_map<uint64_t, Edge> &groupedData) {
    std::cout << "Printing grouped data:" << std::endl << std::endl;
    for (auto &it: groupedData.lock_table()) {
        uint64_t source = it.first;
        Edge &DestTime = it.second;
        auto lt2 = DestTime.lock_table();


        // Print the source and its associated edges
        printf("Time %" PRIu64 " contains %" PRIu64 " edges:\n", source, DestTime.size());
        for (const auto &edge: lt2) {
            for (auto &v: edge.second) {
                printf("    - between %" PRIu64 " and %" PRIu64 " with weight %f" PRIu64 "\n", edge.first, v.first,v.second);
            }
        }
        std::cout << std::endl;
    }
    std::cout << "-------------------------------------" << std::endl;
}

/**
 * Applies the given function @p func to all edges within the given range.
 * @param start of the range inclusive
 * @param end of the range exclusive
 * @param func
 */
void AdjListWeight::rangeQuery(uint64_t start, uint64_t end, const std::function<void(uint64_t, uint64_t, uint64_t, double)> &func) {
    auto t1 = std::chrono::high_resolution_clock::now();
    auto uniqueTimesMap = genUniqueTimeMap(start, end);

    for (auto &time: uniqueTimesMap) {
        Edge e = edges.find(time.second);

        for (const auto &vector: e.lock_table()) {
            for (auto &edge: vector.second) {
                func(time.second, vector.first, edge.first, edge.second);
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "rangeQuery has taken " << ms_int.count() << "ms\n";
}

/**
 * finds edges via uniqueTimesMap and applies given function to those edges
 * @tparam F
 * @param start
 * @param end
 * @param f
 */
template <typename F>
void AdjListWeight::rangeQueryToTimeParlay(uint64_t start, uint64_t end, F&& f) {
    //auto t1 = std::chrono::high_resolution_clock::now();
    auto uniqueTimesMap = genUniqueTimeMap(start, end);
    //auto lt = edges.lock_table();

    parlay::parallel_for(0, uniqueTimesMap.size(), [&](uint64_t i) {
        uint64_t time = uniqueTimesMap[i];
        Edge innerTbl = edges.find(time);
        //f applied to time layer
        f(time, innerTbl);
    });
    /*
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "rangeQueryToTimeParlay has taken " << ms_int.count() << "ms\n";
     */
}

/**
 * Applies the given function @p func to all inner map entries within the given range.
 * @param start of the range inclusive
 * @param end of the range exclusive
 * @param func
 */
template <typename F>
void AdjListWeight::rangeQueryToSourceParlay(uint64_t start, uint64_t end, F&& f) {
    auto innerF = [&f](uint64_t time, Edge edgeMap){
        auto lt = edgeMap.lock_table();
        //f applied to source Layer
        for (auto &edge: lt) {
            f(time,edge);
        }
    };
    rangeQueryToTimeParlay(start, end, innerF);
}

/**
 * Calls rangeQueryToSourceParlay in order to apply the given function @p func to all edges within the given range.
 * @param start of the range inclusive
 * @param end of the range exclusive
 * @param func
 */
template <typename F>
void AdjListWeight::rangeQueryToDestParlay(uint64_t start, uint64_t end, F&& f) {
    auto innerF = [&f](uint64_t time, const std::pair<uint64_t, std::vector<std::pair<uint64_t,double>>>& v){
        uint64_t source = v.first;
        for (std::pair<uint64_t ,uint64_t> destination: v.second) {
            //f applied to dest layer
            f(time, source, destination.first,destination.second);
        }
    };
    rangeQueryToSourceParlay(start, end, innerF);
}

/**
 *
 * @param start of the range inclusive
 * @param end of the range exclusive
 * @return map containing all unique vertices that have an edge within the time range
 */
libcuckoo::cuckoohash_map<uint64_t, std::pair<uint64_t, bool>> AdjListWeight::getVertices(uint64_t start, uint64_t end){
    //cuckoomap because it's threadsafe, tried parlay::sequence which was not threadsafe in my tests
    libcuckoo::cuckoohash_map<uint64_t, std::pair<uint64_t, bool>> map;
    auto f = [&map](uint64_t time, const std::pair<uint64_t, std::vector<std::pair<uint64_t,double>>>& sourceMap){
        map.insert(sourceMap.first, std::pair<uint64_t ,bool>{0,false});
    };
    rangeQueryToSourceParlay(start, end, f);
    return map;
}

/**
 *
 * @param start of the range inclusive
 * @param end of the range exclusive
 * @param source node whose neighbours are to be found
 * @return Cuckoomap with key = timestamp (within the range), value = vector of neighbours
 */
Edge AdjListWeight::getNeighboursOld(uint64_t start, uint64_t end, uint64_t source){
    Edge map;
    auto f = [&map, &source](uint64_t time, const std::pair<uint64_t, std::vector<std::pair<uint64_t,double>>>& v){
        if (v.first == source) map.insert(time, v.second);
    };
    rangeQueryToSourceParlay(start, end, f);
    return map;
}


void AdjListWeight::getNeighboursHelper(uint64_t start, uint64_t end, uint64_t source, libcuckoo::cuckoohash_map<uint64_t, std::pair<uint64_t ,bool>> &map){
    std::set<uint64_t> set;

    auto f = [&](uint64_t time, const Edge& innerTbl){
        if (innerTbl.contains(source)){
           auto destinationVector = innerTbl.find(source);
            for (auto destination : destinationVector) {

                if (!map.contains(destination.first)){
                    set.insert(destination.first);
                    map.insert(destination.first,std::pair<uint64_t,bool>{destination.second,false});
                }
            }
        }
    };
    rangeQueryToTimeParlay(start, end, f);
    for (uint64_t nextSource:set) {
        getNeighboursHelper(start, end, nextSource, map);
    }
}

uint64_t AdjListWeight::memoryConsumption() {

    auto lt = edges.lock_table();
    uint64_t memory=0;

    for(auto &it:lt){
        uint64_t key = it.first;
        Edge edgeData = it.second;
        auto lt2 = edgeData.lock_table();

        for(const auto &vector: lt2){
            for(auto &edge: vector.second){
                memory+=sizeof (key) + sizeof (vector.first) + sizeof (edge);
            }
        }
    }
    std::cout << "Memory consumption in Bytes:" << memory << std::endl;
    return memory;
}

/**
 * Generates a map of timestamps that are in the graph and within the given range.
 * @param start of the range inclusive
 * @param end of the range exclusive
 * @return the generated map
 */
std::map<uint64_t, uint64_t> AdjListWeight::genUniqueTimeMap(uint64_t start, uint64_t end) {
    std::map<uint64_t, uint64_t> map;
    int i = 0;
    for (uint64_t time : uniqueTimestamps) {
        if (start <= time && time < end){
            map[i] = time;
            i++;
        }
    }
    return map;
}

Edge AdjListWeight::computeComponents(uint64_t start, uint64_t end){
    auto t1 = std::chrono::high_resolution_clock::now();
    //no weights
    libcuckoo::cuckoohash_map<uint64_t, std::pair<uint64_t ,bool>> allNodes = getVertices(start, end);
    libcuckoo::cuckoohash_map<uint64_t, std::pair<uint64_t ,bool>> map;

    Edge components;
    uint64_t cKey = 0;

    auto lt = allNodes.lock_table();

    for (auto sources: lt) {
        uint64_t source = sources.first;

        if (!map.contains(source)){
            map.insert(source, std::pair<uint64_t ,bool>{0,false});
            //here gets weights
            getNeighboursHelper(start, end, source, map);
            libcuckoo::cuckoohash_map<uint64_t, std::pair<uint64_t ,bool>>  tempMap = map;

            components.insert(cKey, std::vector<std::pair<uint64_t ,double >>{});

            auto lt2 = tempMap.lock_table();
            for (auto &connectedNodes:lt2) {
                uint64_t connectedNode = connectedNodes.first;
                if (!map.find_fn(connectedNode,[&](std::pair<uint64_t,bool>& pair){
                    return pair.first;
                })) {
                    map.update(connectedNode, std::pair<uint64_t ,bool>{connectedNodes.second.first,true});
                    components.update_fn(cKey, [&connectedNodes](std::vector<std::pair<uint64_t ,double >> &vector) {
                        vector.emplace_back(std::pair<uint64_t,uint64_t>{connectedNodes.first,connectedNodes.second.first});
                    });
                }
            }
            cKey++;

            /*
            components.insert(cKey, std::vector<uint64_t>{v.first});
            computeComponentsHelper(start, end, allNodes, v.first, components, cKey);
             */
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "computeComponents has taken " << ms_int.count() << "ms\n";
    return components;
}

void AdjListWeight::computeComponentsHelper(uint64_t start, uint64_t end,libcuckoo::cuckoohash_map<uint64_t, std::pair<uint64_t ,bool>> &allNodes,
                                            uint64_t key, Edge &components, uint64_t cKey) {
    allNodes.update_fn(key, [&](std::pair<uint64_t,bool>& innerPair){
        innerPair.second=true;
    });
    //printf("key %" PRIu64 "\n", key);
    //auto lt = allNodesIt.lock_table();
    Edge neighbours = getNeighboursOld(start, end, key);
    //std::vector<uint64_t> v = neighbours.find(0);
    //printf("size %" PRIu64 "\n", v.size());
    auto lt = neighbours.lock_table();

    for (auto &k: lt) {
        for (auto &node: k.second) {
            //printf("node %" PRIu64 "\n", node);

            if (!allNodes.find(node.first).first){
                components.update_fn(cKey, [&node](std::vector<std::pair<uint64_t ,double>> &vector){vector.push_back(node);});
                computeComponentsHelper(start, end, allNodes, node.first, components, cKey);
            }
        }
    }

    /*


for (auto &v: neighbours.lock_table()) {

for (auto node: v.second) {
    if (!allNodes.find(node)){
        components.update_fn(cKey, [&node](std::vector<uint64_t> &vector){vector.push_back(node);});
        computeComponentsHelper(start, end, allNodes, node, components, cKey);
    }

}
     */
}

size_t AdjListWeight::getSize() {
    return edges.size();
}

uint64_t AdjListWeight::getEdgeCount(uint64_t timestamp){
    uint64_t count = 0;
    edges.find_fn(timestamp,[&count](Edge &e){
        auto lt = e.lock_table();
        for (const auto &vector: lt) {
            for (auto &edge: vector.second) {
                count++;
            }
        }
    });
    return count;
}

uint64_t AdjListWeight::getInnerTblCount(uint64_t timestamp){
    uint64_t count = 0;
    edges.find_fn(timestamp,[&count](Edge &e){
        auto lt = e.lock_table();
        for (const auto &vector: lt) {
            count++;
        }
    });
    return count;
}

uint64_t AdjListWeight::getDestSize(uint64_t timestamp, uint64_t source){
    uint64_t destSize;
    edges.find_fn(timestamp,[&source,&destSize](Edge &e){
        e.find_fn(source,[&destSize](std::vector<std::pair<uint64_t, double>>&destinationsWithWeight){
            destSize = destinationsWithWeight.size();
        });
    });
    return destSize;
}

NestedMap AdjListWeight::getEdges(){
    return edges;
}

int64_t AdjListWeight::getMaxTime(){
    auto it = uniqueTimestamps.rbegin();
    uint64_t last = 0;
    if (it != uniqueTimestamps.rend()) last = *it;
    return last;
}

int64_t AdjListWeight::getMinTime(){
    auto it = uniqueTimestamps.begin();
    int64_t first = 0;
    if (it != uniqueTimestamps.end()) first = *it;
    return first;
}

int64_t AdjListWeight::getLastVertex() {
    return this->lastVertex;
}

int64_t AdjListWeight::getMaxVertex(){
    return this->maxVertex;
}

NodeMap AdjListWeight::getNodeGraph() {
//use nodegraph for where i use get Vertex?
    NodeMap nodeGraph = {};
    double sumNodeWeight = 0;
    this->maxVertex = 0;
    libcuckoo::cuckoohash_map<uint64_t, double> temp = {};

    auto it = edges.lock_table();

    //iterate through times
    for (const auto &innerTbl: it) {
        Edge edgeData = innerTbl.second;
        temp = {};
        auto lt2 = edgeData.lock_table();

        //iterate through sources
        for (const auto &vector: lt2) {
            //iterate through destination/weight pairs
            for (auto &edge: vector.second) {
                sumNodeWeight += edge.second;
            }
            temp.insert(vector.first, sumNodeWeight);
            sumNodeWeight = 0;
            this->maxVertex = vector.first > this->maxVertex ? vector.first : this->maxVertex;
            this->lastVertex = vector.first;
        }
        nodeGraph.insert(innerTbl.first, temp);
    }
    return nodeGraph;
}

std::set<TemporalNode> AdjListWeight::tNodes(int64_t t) const {

    std::set<TemporalNode> result = {};

   Edge tempEdges = edges.find(t);
   auto var4 = tempEdges.lock_table();

   for(const auto&innerTbl :var4){
       TemporalNode tn = {t,(int64_t) innerTbl.first};
       result.emplace(tn);
   }
    return result;
}

std::map<long,double> AdjListWeight::getWeightedNeighborhood(TemporalNode tn) const {

    int64_t v = tn.getNode();
    int64_t t = tn.getTime();
    std::map<long,double> neighborhoodMap;
    std::vector<std::pair<uint64_t,double>> neighborhood = edges.find(t).find(v);

    for(auto node:neighborhood){
        neighborhoodMap.emplace(node.first,node.second);
    }
    return neighborhoodMap;
}

map<double,int> AdjListWeight::standardRWR(const set<int>& sources, double alpha, bool byVertex, bool normalizeByDegree, int t1, int t2) {
    if (!sources.empty()) {
        map<int, double> tempMap = {};
        auto var10 = sources.begin();
        /** ASSIGN WEIGHT 1.0 TO ALL SOURCE NODES*/
        while(var10!=sources.end()) {
            int v = *var10;
            tempMap.insert_or_assign(v, 1.0);
            ++var10;
        }

        return this->standardRWR(tempMap, alpha, byVertex, normalizeByDegree, t1, t2);
    } else {
        return {};
    }
}

map<double,int>AdjListWeight::standardRWR(const map<int,double>& sources, double alpha, bool byVertex, bool normalizeByDegree, int t1, int t2) {

    bool sanity = true;
    int STEPS = true;

    if(!sources.empty()){
        if(!(alpha<0.0) && !(alpha>1.0) !=0){

            double sumSources = 0.0;
            int v;

            /** NUMBER OF SOURCE NODES -> Stored in sumSources*/
            for(auto var13 = sources.begin(); var13!=sources.end(); sumSources += (double)sources.find(v)->second) {
                auto temp = *var13;
                v = temp.first;
                ++var13;
            }

            map<int, double> mass = {};

            auto var14 = sources.begin();

            /**  MASS  (NODE, 1/sumSOURCES) */
            while(var14!=sources.end()) {
                auto temp = *var14;
                int v2 = temp.first;
                mass.insert_or_assign(v2, sources.find(v2)->second / sumSources);
                ++var14;
            }

            map<int,double> wSums ={}; //wsum = hold source node and its mass (t1-t2)
            double score;
            map<int,double> temp;

            /** START */
            for(int step = 0; step < 4; ++step) {
                temp = {};
                auto var19 = mass.begin();
                double value;
                int from;

                /** iterate through mass */
                while(var19!=mass.end()) {
                    auto temp2 = *var19;
                    from = temp2.first;
                    score = 0.0;
                    int t;

                    /** SET SCORE FROM WSUMS */
                    if (wSums.find(from)!=wSums.end()) {
                        score = (double)wSums.find(from)->second;

                    } else {
                        t = t1;

                        while(true) {
                            if (t >= t2) {
                              /** INSERT NNeighborhood score to from */
                                wSums.insert_or_assign(from, score);
                                break;
                            }
                            double w;

                            /** iterate through timespan of NNeighborhood and SUM Nodeweights -> Mass (1,10) + (1,11) till t2  */
                            auto list = this->edges.find(t).find(from);
                            for(auto var24 = list.begin();var24!=list.end();score +=w) {
                                auto temp3 = *var24;
                                w = temp3.second;
                                ++var24;
                            }
                            ++t;
                        }
                    }

                    /**  grade all NN(from) INTO TEMP OVER ALL T-T2 */
                    for(t = t1; t < t2; ++t) {
                        auto wMap = this->edges.find(t).find(from);
                        auto var25 = wMap.begin();

                        while(var25!=wMap.end()) {//goes through wMap at time t & source ->(get destinations,weights)-> insert destination + calc value into temp
                            auto temp4 = *var25;
                            int to = temp4.first;
                            value = temp.find(to)!=temp.end() ? (double)temp.find(to)->second : 0.0;

                            value += (double)std::find_if(wMap.begin(), wMap.end(),[to](const std::pair<uint64_t , double>& element){
                                return element.first == to;})->second * (1.0 - alpha) * (1.0 / score) * (double)mass.find(from)->second;
                            temp.insert_or_assign(to, value);
                            ++var25;
                        }
                    }
                    ++var19;
                }

                auto var20 = sources.begin();
                /**ADD SOURCE NODE TO TEMP*/
                while(var20!=sources.end()) {
                    auto temp4 = *var20;
                    from = temp4.first;
                    value = temp.find(from)!=temp.end() ? (double)temp.find(from)->second : 0.0;
                    value += alpha * (double)sources.find(from)->second / sumSources;
                    temp.insert_or_assign(from, value);
                    ++var20;
                }

                mass = temp;

                /**SANITY CHECK */
                if (sanity) {//every nodes degree adds up to 1 checks
                    double sum = 0.0;

                    int v2;
                    for(auto var21 = temp.begin(); var21!=temp.end(); sum += (double)mass.find(v2)->second) {
                        auto temp5 = *var21;
                        v2 = temp5.first;
                        ++var21;
                    }

                    cout<<"Step " << step << ": " << sum<<endl;
                }
            }
            cout<<""<<endl;
            if(byVertex){

                map<double, int> result;
                auto var19 = mass.begin();

                while(var19!=mass.end()){
                    auto e = *var19;
                    result.insert_or_assign(e.second,e.first);
                    ++var19;
                }
                return result;
            }
        }
    }
    return {};
}

set<int>AdjListWeight::standardSweeps(const map<double, int>& map,  vector<double>& conductance, bool useTC, int t1, int t2){


    if(!map.empty()){

        set<int> result;
        set<int64_t> current;
        auto temp = map.rbegin();
        double curKey = temp->first;
        result.emplace(map.find(curKey)->second);
        current.emplace(map.find(curKey)->second);
        double cut  = 0.0;

        double w = 0;

        for(int t = t1; t<t2;++t){
            auto list = this->edges.find(t).find(map.find(curKey)->second);
            for(auto var13=list.begin();var13!=list.end();cut+=w){
                w = var13->second;
                ++var13;
            }
        }

        double internal = cut;
        double external = 0.0;

        NestedMapWeight subMap;

        for (int i = t1; i < t2; ++i) {
            subMap.insert_or_assign(i,this->edges.find(i));
        }

        auto lt = subMap.lock_table();

        for (const auto &innerTbl: lt) {
            Edge edgeData = innerTbl.second;
            auto lt2 = edgeData.lock_table();
            for (const auto &vector: lt2) {
                for (auto &edge: vector.second) {
                  external += edge.second;
                }
            }
        }
        external -= cut;
        conductance.emplace_back(1.0);


        while(true){

            curKey= lowerKey(map,curKey);

            int curNode = map.find(curKey)->second;

            NestedMapWeight subMap2;

            for (int i = t1; i < t2; ++i) {
                subMap2.insert_or_assign(i, this->edges.find(i));
            }

            auto lt1 = subMap2.lock_table();

            for (const auto &innerTbl: lt) {
                Edge edgeData = innerTbl.second;
                auto lt2 = edgeData.lock_table();
                for (const auto &vector: lt2) {
                    if(vector.first==curNode) {
                        for (auto &edge: vector.second) {
                            double w3 = edge.second;
                            internal += w3;
                            external -= w3;
                            if (current.find(edge.first) != current.end()) {
                                cut -= w3;
                            } else {
                                cut += w3;
                            }
                        }
                    }

                }
            }

            current.emplace(curNode);
            double curCond = cut / min(external, internal);

            if (curCond < conductance[0]) {
                conductance[0] = curCond;
                for (auto elem: current) {
                    result.emplace(elem);
                }
            }
            if (current.size() * 2 > map.size()) {
                break;
            }
        }
            if(useTC){
                conductance[0] = this->normalizeConductance(conductance[0], t2 - t1);
            }
            return result;
        }else{
        return {};
    }
}

double AdjListWeight::lowerKey(const map<double, int>& map, double hashScore) {

    auto it  = map.lower_bound(hashScore);
    if (it != map.begin()) {
        it--;
    }
    return it->first;
}

double AdjListWeight::normalizeConductance(double phi, int duration) {
    return phi / exp(0.0 * (double)duration);
}
