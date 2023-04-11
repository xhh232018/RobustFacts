#ifndef GRAPHMANAGER_H
#define GRAPHMANAGER_H


#include <cstdio>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <sys/time.h>


inline
unsigned long long getTime() {
    struct timeval tv;

    gettimeofday(&tv, NULL);

    unsigned long long ret = tv.tv_usec;

    /* Adds the seconds after converting them to microseconds (10^-6) */
    ret += (tv.tv_sec * 1000 * 1000);

    return ret;
};


inline
double getInterval(unsigned long long start, unsigned long long stop) {
    return (double) (stop - start) / 1000.0;
};


__inline__
uint32_t extractHigh32bits(uint64_t num) {
    return (uint32_t) ((num >> 32) & 0xFFFFFFFF);
}

__inline__
uint32_t extractLow32bits(uint64_t num) {
    return (uint32_t) (num & 0xFFFFFFFF);
}

__inline__
uint64_t combine2u32(uint32_t high, uint32_t low) {
    uint64_t res = 0;
    res = res | high;
    res = res << 32;
    res = res | low;
    return res;
}

//template <typename T1, typename T2, typename T3>
struct priority_entry {
    uint32_t nid;
    double received_prob; //probability received from source
    double toSource_prob;
    double score;// probability assigned. Bidirectional
    std::string patternSig;
    int globalPathPatternSignificance;
    std::unordered_set<uint32_t> parents; //to avoid going immediately back


    priority_entry(uint32_t n, double rp, double s, std::string p, int rpc = 1, double toSourceP = 0) {
        nid = n;
        received_prob = rp;
        toSource_prob = toSourceP;
        score = s; // -1 means not calculated yet
        patternSig = p;
        globalPathPatternSignificance = rpc;
    }

    bool operator==(const priority_entry &pe2) const {
        return nid == pe2.nid && patternSig == pe2.patternSig;
    }
};


namespace std {
    template<>
    struct hash<priority_entry> {
        typedef priority_entry argument_type;
        typedef std::size_t result_type;

        result_type operator()(argument_type const &pe) const noexcept {
            result_type const h(std::hash<uint32_t>{}(pe.nid));
            return h;
        }
    };
}

class pe_comparator {
public:
    bool operator()(const priority_entry &p1, const priority_entry &p2) {
        /// based on total score
//        return p1.score < p2.score;
        ///based on received prob
//        return p1.received_prob < p2.received_prob;
        ///linear scoring on score and significance.
        if (p1.score < p2.score) return true;
        if (p1.score > p2.score) return false;
        return p1.globalPathPatternSignificance < p2.globalPathPatternSignificance;

    }
};

class pe_comparator_pairNode {
public:
    bool operator()(const priority_entry &p1, const priority_entry &p2) {
        /// based on total score
//        return p1.score < p2.score;
        ///based on received prob
//        return p1.received_prob < p2.received_prob;
        ///linear scoring on score and significance.
        if (p1.score < p2.score) return true;
        if (p1.score > p2.score) return false;
        return p1.globalPathPatternSignificance < p2.globalPathPatternSignificance;

    }
};

class pe_comparator_pairNodeTopkHeap {
public:
    bool operator()(const priority_entry &p1, const priority_entry &p2) {
        /// based on total score
//        return p1.score > p2.score;
        ///based on received prob
//        return p1.received_prob > p2.received_prob;
        ///linear scoring on score and significance.
        if (p1.score > p2.score) return true;
        if (p1.score < p2.score) return false;
        return p1.globalPathPatternSignificance > p2.globalPathPatternSignificance;

    }
};

class GraphManager {
public:
    /**Graph Adjacency Lists.
     * Type specific store.*/
    size_t node_sz = 0; // 1 larger than the actual node number
    size_t edge_sz = 0;

    unsigned int *nodes = nullptr;
    unsigned long *edges = nullptr; //high 32-bit for edge type, low 32-bit for edge id

    /***Node type maps*/

    std::unordered_map<uint32_t,std::unordered_set<std::string>> nid2types;
    std::unordered_map<std::string, std::unordered_set<uint32_t>> type2nid;
    std::vector<std::string> nid2IRI;
    std::vector<std::unordered_map<uint32_t,double>> degStats;

    /**Edge Type to start node.*/

    std::vector<std::string> typeId2Name;
    std::unordered_map<std::string, uint32_t> IRI2nid;
    std::vector<std::unordered_set<uint32_t>> typeId2Count; // edge type to source node

    GraphManager() = default;

    ~GraphManager() {
        delete[] nodes;
        delete[] edges;
    }

    //load the graph from edge file
    void readEdges(const char *filename);

    //load node types
    void readNodeTypes(const char * filename);

    // serialization implementation
    void serializeToDisk(const char *filename);

    void deserializeFromDisk(const char *filename);

    /**Given an edge type, find the first edge position*/
    int firstEdgePos(uint32_t start_node, uint32_t eid);

    /** Set edge type to edge pos */
    std::unordered_map<uint64_t, int> snidEid2pos;
    void setEdgeType2Pos();

    std::unordered_map<std::string, uint32_t> iri2Nid;
    void setIRI2NidMap();
};


#endif //PATHPATTERNMINING_GRAPHMANAGER_H
