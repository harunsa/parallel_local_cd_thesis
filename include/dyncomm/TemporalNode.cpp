#include <fstream>
#include <cinttypes>
#include <set>
#include <vector>
#include <limits>


using namespace std;

class TemporalNode{

public:
    int64_t time;
    int64_t node;
    TemporalNode(){
    };

    TemporalNode(int64_t time, int64_t node){
        this->time = time;
        this->node = node;
    }
    bool operator<(const TemporalNode& other) const {
        if (time != other.time) {
            return time < other.time;
        }
        return node < other.node;
    }

    int64_t getNode(){
        return this->node;
    }

    int64_t getTime(){
        return this->time;
    }

    static set<int> extractVertices(const set<TemporalNode>& nodes) {
        set<int> result = {};
        auto var3 = nodes.begin();

        while(var3!=nodes.end()) {
            TemporalNode tn = *var3;
            result.emplace(tn.getNode());
            ++var3;
        }
        return result;
    }

    static set<int> extractTimes(const set<TemporalNode>& nodes) {
        set<int> result = {};
        auto var3 = nodes.begin();

        while(var3!=nodes.end()) {
            TemporalNode tn = *var3;
            result.emplace(tn.getTime());
            ++var3;
        }
        return result;
    }

    static vector<int> span(const set<TemporalNode>& nodes) {
        vector<int> result = {std::numeric_limits<int>::max(), std::numeric_limits<int>::min()};
        auto var3 = nodes.begin();

        while(var3!=nodes.end()) {
            TemporalNode tn = *var3;
            int i = tn.getTime();
            if (i < result[0]) {
                result[0] = i;
            }

            if (i > result[1]) {
                result[1] = i;
            }
            ++var3;
        }
        return result;
    }
};