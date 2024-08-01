#ifndef TEMPUS_ADJ_LISTWEIGHT_H
#define TEMPUS_ADJ_LISTWEIGHT_H


#include <set>
#include <map>
#include <cstdint>
#include <unordered_map>
#include "../libcuckoo/libcuckoo/cuckoohash_map.hh"
#include "../dyncomm/TemporalNode.cpp"

typedef libcuckoo::cuckoohash_map<uint64_t, libcuckoo::cuckoohash_map<uint64_t, std::vector<std::pair<uint64_t,double>>>> NestedMapWeight;
typedef libcuckoo::cuckoohash_map<uint64_t, libcuckoo::cuckoohash_map<uint64_t, double>> NodeMap;

class AdjListWeight{
public:
    //time < source < list of destinations>>
    NestedMapWeight edges;
    //nodes mapped with their volumes time -> node -> volume
    NodeMap nodes;
    void addFromFile(const std::string& path);
    void printGraph();
    void printNodeGraph();
    size_t getSize();
    bool findEdge(uint64_t source, uint64_t destination, uint64_t time);
    bool findEdge(uint64_t source, uint64_t destination, uint64_t start, uint64_t end);
    void batchOperation(bool insert, NestedMapWeight &groupedData);
    void batchOperationParlay(bool insert, NestedMapWeight &groupedData, std::unordered_map<uint64_t, uint64_t> uniqueTimesMap);
    void rangeQuery(uint64_t start, uint64_t end, const std::function<void(uint64_t,uint64_t,uint64_t, double)> &func);
    uint64_t memoryConsumption();
    size_t getEdgeCount(uint64_t timestamp);
    uint64_t getInnerTblCount(uint64_t timestamp);
    uint64_t getDestSize(uint64_t timestamp, uint64_t source);
    libcuckoo::cuckoohash_map<uint64_t,std::pair<uint64_t ,bool>> getVertices(uint64_t start, uint64_t end);
    libcuckoo::cuckoohash_map<uint64_t, std::vector<std::pair<uint64_t,double>>>
    getNeighboursOld(uint64_t start, uint64_t end, uint64_t source);
    libcuckoo::cuckoohash_map<uint64_t, std::vector<std::pair<uint64_t ,double>>> computeComponents(uint64_t start, uint64_t end);
    void addFromGraph(const NestedMapWeight& map, uint64_t timestamp);
    void insertEdgeUndirected(uint64_t source, std::pair<uint64_t,double> destination, uint64_t time);
    void deleteEdgeUndirected(uint64_t source, uint64_t destination, uint64_t time);
    NestedMapWeight getEdges();
    void getNodeCount(uint64_t timestamp);
    std::set<uint64_t> getUinqueTimestamps();
    int64_t getLastVertex();
    int64_t getMaxTime();
    int64_t getMinTime();
    int64_t getMaxVertex();
    NodeMap getNodeGraph();
    std::set<TemporalNode> tNodes(int64_t t) const;
    std::map<long,double> getWeightedNeighborhood(TemporalNode tn) const;
    map<double,int>
    standardRWR(const set<int> &sources, double alpha, bool byVertex, bool normalizeByDegree, int t1, int t2);
    map<double,int>
    standardRWR(const map<int,double> &sources, double alpha, bool byVertex, bool normalizeByDegree, int t1, int t2);
    set<int> standardSweeps(const map<double, int>& map,  vector<double>& conductance, bool useTC, int t1, int t2);
    double lowerKey(const map<double, int> &map, double hashScore);
    double normalizeConductance(double phi, int duration);

private:


    std::set<int64_t> uniqueTimestamps;
    //the last vertex in the map
    int64_t lastVertex;
    //the maximum Node
    int64_t maxVertex;

    //TODO: std::unorderedmap<uint64_t, libcuckoo::cuckoohash_map<uint64_t, std::vector<uint64_t>>>
    //TODO: std::map<uint64_t, libcuckoo::cuckoohash_map<uint64_t, std::vector<uint64_t>>>
    static void insertEdgeDirected(uint64_t source, std::pair<uint64_t,double> destination, uint64_t time, NestedMapWeight &map);

    void deleteEdgeDirected(uint64_t source, uint64_t destination, uint64_t time);

    static void sortBatch(const std::vector<uint64_t>& sourceAdds, const std::vector<uint64_t>& destinationAdds,
                          const std::vector<double>& weightAdds, const std::vector<uint64_t>& timeAdds, NestedMapWeight &groupedData);
    static void printGroupedData(NestedMapWeight &groupedData);
    static void fileReaderHelper(std::vector<uint64_t> &sourceVector, std::vector<uint64_t> &destinationVector, std::vector<double> &weightVector,
                                 std::vector<uint64_t> &timeVector, std::set<uint64_t> &uniqueTimes, uint64_t source,
                                 uint64_t destination, double weight, uint64_t time);
    void uniqueTimesHelper(std::unordered_map<uint64_t, uint64_t> &uniqueTimesMap, std::set<uint64_t> &uniqueTimes, bool insert);
    std::map<uint64_t, uint64_t> genUniqueTimeMap(uint64_t start, uint64_t end);
    template<typename F>
    void rangeQueryToSourceParlay(uint64_t start, uint64_t end, F &&f);
    template<typename F>
    void rangeQueryToDestParlay(uint64_t start, uint64_t end, F &&f);


    void
    getNeighboursHelper(uint64_t start, uint64_t end, uint64_t source, libcuckoo::cuckoohash_map<uint64_t, std::pair<uint64_t ,bool>> &map);

    libcuckoo::cuckoohash_map<uint64_t, std::pair<uint64_t ,bool>> getNeighbours(uint64_t start, uint64_t end, uint64_t source);

//    void
//    computeComponentsHelper(uint64_t start, uint64_t end, libcuckoo::cuckoohash_map<uint64_t, bool> &allNodes,
//                            uint64_t key,
//                            libcuckoo::cuckoohash_map<uint64_t, std::vector<uint64_t>> &components, uint64_t cKey,
//                            libcuckoo::cuckoohash_map<uint64_t, bool> &map);
//
//    void
//    computeComponentsHelper(uint64_t start, uint64_t end,
//                            libcuckoo::cuckoohash_map<uint64_t, std::pair<uint64_t ,std::pair<uint64_t,bool>>> &allNodes, uint64_t key,
//                            libcuckoo::cuckoohash_map<uint64_t, std::vector<std::pair<uint64_t,uint64_t>>> &components, uint64_t cKey);
    template<typename F>
    void rangeQueryToTimeParlay(uint64_t start, uint64_t end, F &&f);

    void computeComponentsHelper(uint64_t start, uint64_t end,
                                 libcuckoo::cuckoohash_map<uint64_t, std::pair<uint64_t, bool>> &allNodes, uint64_t key,
                                 libcuckoo::cuckoohash_map<uint64_t, std::vector<std::pair<uint64_t, double>>> &components,
                                 uint64_t cKey);

    //return allVertices over all timestamps, not sure overall or over 1 timestamp only
    libcuckoo::cuckoohash_map<uint64_t, std::vector<std::pair<uint64_t, double>>> getVertices();

};




#endif //TEMPUS_ADJ_LISTWEIGHT_H
