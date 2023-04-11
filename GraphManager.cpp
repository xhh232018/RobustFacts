#include "GraphManager.h"
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <queue>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/vector.hpp>

using namespace std;

int GraphManager::firstEdgePos(uint32_t start_node, uint32_t eid) {
    uint64_t edgeValue = combine2u32(eid, 0);
    int _min = nodes[start_node];
    int _max = nodes[start_node + 1];
    while (_min < _max) {
        int _d = (_min + _max) / 2;
        if (edges[_d] > edgeValue) {
            _max = _d;
        } else {
            _min = _d + 1;
        }
    }
    if (extractHigh32bits(edges[_min]) == eid)
        return _min;
    return -1;

}

void GraphManager::readEdges(const char *filename) {

    unordered_map<string, int> etypes2id;
    
    vector<pair<uint64_t, string>> tmp_edges;
    ifstream ifs(filename);
    if (!ifs.is_open()) {
        cout << "Invalid filepath\n";
        return;
    }

    string line;

    while (getline(ifs, line)) {
        istringstream iss(line);
        uint32_t sid, eid;
        string token;
        getline(iss, token, ',');
        if (IRI2nid.count(token))
            sid = IRI2nid[token];
        else {
            sid = (uint32_t) IRI2nid.size();
            IRI2nid[token] = sid;
        }

        getline(iss, token, ',');
        if (IRI2nid.count(token))
            eid = IRI2nid[token];
        else {
            eid = (uint32_t) IRI2nid.size();
            IRI2nid[token] = eid;
        }

        if (sid == eid) continue;

        getline(iss, token, ',');
        token.erase(std::remove(token.begin(),token.end(),'\r'),token.end());
        string reversedType = "r-" + token;
        if (!etypes2id.count(token)) {
            int sz = etypes2id.size();
            etypes2id[token] = sz;
            etypes2id[reversedType] = sz + 1;
        }
        tmp_edges.emplace_back(make_pair(combine2u32(sid, eid), token));
        tmp_edges.emplace_back(make_pair(combine2u32(eid, sid), reversedType));
    }

    ifs.close();

    sort(tmp_edges.begin(), tmp_edges.end());
    //remove repeatitions
    auto iter = std::unique(tmp_edges.begin(), tmp_edges.end());

    tmp_edges.resize(std::distance(tmp_edges.begin(), iter));

    node_sz = IRI2nid.size();
    edge_sz = tmp_edges.size();
    nodes = new uint32_t[node_sz + 1];
    edges = new uint64_t[edge_sz];

    uint32_t iter_j = 0;
    for (int i = 0; i < node_sz; ++i) {
        uint32_t startNode = extractHigh32bits(tmp_edges[iter_j].first);
        if (i <= startNode) {
            nodes[i] = iter_j;
            continue;
        }

        if (iter_j == edge_sz) {
            nodes[i] = (uint32_t) edge_sz;
            continue;
        }
        while (iter_j < edge_sz && extractHigh32bits(tmp_edges[iter_j].first) ==
                                   startNode) {
            iter_j++;
        }
        nodes[i] = iter_j;

    }
    nodes[node_sz] = (uint32_t) edge_sz;

    /**Set edges type*/
    for (int i = 0; i < edge_sz; ++i) {
        int tid = etypes2id[tmp_edges[i].second];
        uint32_t eid = extractLow32bits(tmp_edges[i].first);
        edges[i] = combine2u32((uint32_t) tid, eid);
    }

    typeId2Count.resize(etypes2id.size());
    degStats.resize(node_sz);
    for (uint32_t i = 0; i < node_sz; ++i) {
        //re-sort neighbors by edge types;
        int adj_sz = nodes[i + 1] - nodes[i];
        if (adj_sz == 0) continue;

        sort(edges + nodes[i], edges + nodes[i] + adj_sz);
        unordered_map<uint32_t, double> tid2count;
        for (int j = 0; j < adj_sz; ++j) {
            auto e = edges[nodes[i] + j];
            uint32_t tid = extractHigh32bits(e);

            if(degStats[i].count(tid)) degStats[i][tid]++;
            else degStats[i][tid] = 1;

            typeId2Count[tid].insert(i);
            if (tid2count.count(tid))
                tid2count[tid] += 1;
            else
                tid2count[tid] = 1;
        }
    }

    typeId2Name.resize(etypes2id.size());
    for (auto &kv_pair : etypes2id)
        typeId2Name[kv_pair.second] = kv_pair.first;

    nid2IRI.resize(IRI2nid.size());
    for (auto &kv_pair : IRI2nid)
        nid2IRI[kv_pair.second] = kv_pair.first;

    printf("reading finished!\n");
}

void GraphManager::readNodeTypes(const char *filename) {
    ifstream ifs(filename);
    if (!ifs.is_open()) {
        cout << "Invalid File Path" << endl;
    }
    string line;
    uint32_t currentNode;

    while (getline(ifs, line)) {
        istringstream iss(line);
        string iri,token,tToken;
        bool first = true;
        while(getline(iss, token, ',')){
            token.erase(std::remove(token.begin(),token.end(),'\r'),token.end());
            if(first && !IRI2nid.count(token)) break;
            if(first && IRI2nid.count(token)){
                currentNode = IRI2nid[token];
                first = false;
                continue;
            }

            nid2types[currentNode].insert(token);
        }
    }

    cout<<"Ntype finished!"<<endl;
}

void GraphManager::serializeToDisk(const char *filename) {
    ofstream ofs(filename);
    if (!ofs.is_open()) {
        cout << "Invalid File path\n";
        return;
    }

    ofs.write((char *) &node_sz, sizeof(size_t));
    ofs.write((char *) &edge_sz, sizeof(size_t));

    ofs.write((char *) nodes, sizeof(uint32_t) * (node_sz + 1));
    ofs.write((char *) edges, sizeof(uint64_t) * edge_sz);


    boost::archive::binary_oarchive oa(ofs);
    oa << typeId2Name;
    oa << typeId2Count;
    oa << nid2IRI;
    oa << nid2types;
    oa << type2nid;
    oa << snidEid2pos;
    oa << IRI2nid;
    oa << degStats;

    ofs.close();
}

void GraphManager::deserializeFromDisk(const char *filename) {
    ifstream ifs(filename);
    if (!ifs.is_open()) {
        cout << "Invalid File path\n";
        return;
    }

    ifs.read((char *) &node_sz, sizeof(size_t));
    ifs.read((char *) &edge_sz, sizeof(size_t));

    nodes = new uint32_t[node_sz + 1];
    edges = new uint64_t[edge_sz];


    ifs.read((char *) nodes, sizeof(uint32_t) * (node_sz + 1));
    ifs.read((char *) edges, sizeof(uint64_t) * edge_sz);

    boost::archive::binary_iarchive ia(ifs);
    ia >> typeId2Name;
    ia >> typeId2Count;
    ia >> nid2IRI;
    ia >> nid2types;
    ia >> type2nid;
    ia >> snidEid2pos;
    ia >> IRI2nid;
    ia >> degStats;

    ifs.close();
}

void GraphManager::setEdgeType2Pos() {

    ifstream ifs("snidEid2pos.txt");
    if (ifs.is_open()) {
        boost::archive::binary_iarchive ia(ifs);
        ia >> snidEid2pos;
        ifs.close();
        return;
    }
    ifs.close();
    for (uint32_t i = 0; i < node_sz; ++i) {
        int edge_start = nodes[i];
        int adj_sz = nodes[i + 1] - edge_start;
        for (int j = 0; j < adj_sz; ++j) {
            auto edge = edges[edge_start + j];
            auto eid = extractHigh32bits(edge);
            auto snidEid = combine2u32(i, eid);
            if (snidEid2pos.count(snidEid)) {
                continue;
            }
            snidEid2pos[snidEid] = edge_start + j;
        }
    }
    ofstream ofs("snidEid2pos.txt");
    boost::archive::binary_oarchive oa(ofs);
    oa << snidEid2pos;
    ofs.close();
}

void GraphManager::setIRI2NidMap() {
    ifstream ifs("iri2Nid.txt");
    if (ifs.is_open()) {
        boost::archive::binary_iarchive ia(ifs);
        ia >> iri2Nid;
        ifs.close();
        return;
    }
    ifs.close();
    string prefix = "https://www.wikidata.org/wiki/";
    for (uint32_t i = 0; i < node_sz; ++i) {
        auto iri = nid2IRI[i].substr(prefix.size() + 1);
        iri2Nid[iri] = i;
    }
    ofstream ofs("iri2Nid.txt");
    boost::archive::binary_oarchive oa(ofs);
    oa << iri2Nid;
    ofs.close();


}
