#ifndef QUERYPERTURBATION_H
#define QUERYPERTURBATION_H

#include <vector>
#include "GraphManager.h"
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <cmath>
#include <queue>
#include <iostream>

class QueryPerturbation{
public:
    std::unordered_set<uint32_t> pSpace;
    double givenScore;

public:

    QueryPerturbation() = default;
    ~QueryPerturbation() = default;

    void neiChecking(GraphManager &gm,uint32_t source,std::unordered_map<uint32_t,std::vector<uint32_t>> &srcNeighMap,double &A_Coeff); 

    double proxCal(GraphManager &gm,uint32_t ptbNode,std::unordered_map<uint32_t,std::vector<uint32_t>> cxtNeighMap,double &A_Coeff); 
    
    double queryPtbSSCal(GraphManager &gm,uint32_t ptb,std::vector<uint32_t> backEdges,uint32_t attr, uint32_t attrVal); 

    void pSpaceFromAttr(GraphManager &gm,std::vector<uint32_t> oriEdges,uint32_t attr, uint32_t attrVal,uint32_t cxtNode);

    void bySampleSizeBase(GraphManager &gm, std::string &pattern, double oriScore, double timeLimit, uint32_t attr, 
                            uint32_t attrVal, uint32_t cxtNode, int index, std::string opFileName, std::string prefix);

    void printBase(GraphManager &gm,int index,double trueScore,double oriScore,uint32_t caGlb,double proxGlb, double sScoreGlb,
                    double fullTime,std::string opFileName, std::string prefix);

    void bySampleSizeExp(GraphManager &gm, std::string &pattern, double oriScore, double timeLimit, double sampleRatio, uint32_t attr, 
                            uint32_t attrVal, uint32_t cxtNode, int index, std::string opFileName, std::string prefix);

    void printSampleExp(GraphManager &gm,double estScore, double oriScore, uint32_t ca,
                             double prox, double sScoreSP,double sampleTime,
                             int index,std::string opFileName, std::string prefix,int sampleSize, double sampleRatio);
};

#endif