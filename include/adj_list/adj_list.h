#ifndef TEMPUS_ADJ_LIST_H
#define TEMPUS_ADJ_LIST_H

#include <set>
#include <map>
#include <cstdint>
#include "libcuckoo/cuckoohash_map.hh"

typedef libcuckoo::cuckoohash_map<uint64_t, libcuckoo::cuckoohash_map<uint64_t, std::vector<uint64_t>>> NestedMap;

class AdjList{
public:
    void addFromFile(const std::string& path);
    void printGraph();
    size_t getSize();
    bool findEdge(uint64_t source, uint64_t destination, uint64_t time);
    bool findEdge(uint64_t source, uint64_t destination, uint64_t start, uint64_t end);
    void batchOperation(bool insert, NestedMap &groupedData);
    void batchOperationParlay(bool insert, NestedMap &groupedData, std::unordered_map<uint64_t, uint64_t> uniqueTimesMap);
    void rangeQuery(uint64_t start, uint64_t end, const std::function<void(uint64_t,uint64_t,uint64_t)> &func);
    uint64_t memoryConsumption();
    size_t getEdgeCount(uint64_t timestamp);
    uint64_t getInnerTblCount(uint64_t timestamp);
    uint64_t getDestSize(uint64_t timestamp, uint64_t source);
    libcuckoo::cuckoohash_map<uint64_t, bool> getVertices(uint64_t start, uint64_t end);
    libcuckoo::cuckoohash_map<uint64_t, std::vector<uint64_t>>
    getNeighboursOld(uint64_t start, uint64_t end, uint64_t source);
    libcuckoo::cuckoohash_map<uint64_t, std::vector<uint64_t>> computeComponents(uint64_t start, uint64_t end);
    void addFromGraph(const NestedMap& map, uint64_t timestamp);
    void insertEdgeUndirected(uint64_t source, uint64_t destination, uint64_t time);
    void deleteEdgeUndirected(uint64_t source, uint64_t destination, uint64_t time);
    NestedMap getEdges();
    void getNodeCount(uint64_t timestamp);

private:
    //time < source < list of destinations>>
    NestedMap edges;
    std::set<uint64_t> uniqueTimestamps;

    //TODO: std::unorderedmap<uint64_t, libcuckoo::cuckoohash_map<uint64_t, std::vector<uint64_t>>>
    //TODO: std::map<uint64_t, libcuckoo::cuckoohash_map<uint64_t, std::vector<uint64_t>>>
    static void insertEdgeDirected(uint64_t source, uint64_t destination, uint64_t time, NestedMap &map);

    void deleteEdgeDirected(uint64_t source, uint64_t destination, uint64_t time);

    static void sortBatch(const std::vector<uint64_t>& sourceAdds, const std::vector<uint64_t>& destinationAdds,
                          const std::vector<uint64_t>& timeAdds, NestedMap &groupedData);
    static void printGroupedData(NestedMap &groupedData);
    static void fileReaderHelper(std::vector<uint64_t> &sourceVector, std::vector<uint64_t> &destinationVector,
                          std::vector<uint64_t> &timeVector, std::set<uint64_t> &uniqueTimes, uint64_t source,
                          uint64_t destination, uint64_t time);
    void uniqueTimesHelper(std::unordered_map<uint64_t, uint64_t> &uniqueTimesMap, std::set<uint64_t> &uniqueTimes, bool insert);
    std::map<uint64_t, uint64_t> genUniqueTimeMap(uint64_t start, uint64_t end);
    template<typename F>
    void rangeQueryToSourceParlay(uint64_t start, uint64_t end, F &&f);
    template<typename F>
    void rangeQueryToDestParlay(uint64_t start, uint64_t end, F &&f);


    void
    getNeighboursHelper(uint64_t start, uint64_t end, uint64_t source, libcuckoo::cuckoohash_map<uint64_t, bool> &map);

    libcuckoo::cuckoohash_map<uint64_t, bool> getNeighbours(uint64_t start, uint64_t end, uint64_t source);

    void
    computeComponentsHelper(uint64_t start, uint64_t end, libcuckoo::cuckoohash_map<uint64_t, bool> &allNodes,
                            uint64_t key,
                            libcuckoo::cuckoohash_map<uint64_t, std::vector<uint64_t>> &components, uint64_t cKey,
                            libcuckoo::cuckoohash_map<uint64_t, bool> &map);

    void
    computeComponentsHelper(uint64_t start, uint64_t end, libcuckoo::cuckoohash_map<uint64_t, bool> &allNodes,
                            uint64_t key,
                            libcuckoo::cuckoohash_map<uint64_t, std::vector<uint64_t>> &components, uint64_t cKey);
    template<typename F>
    void rangeQueryToTimeParlay(uint64_t start, uint64_t end, F &&f);
};

#endif //TEMPUS_ADJ_LIST_H
