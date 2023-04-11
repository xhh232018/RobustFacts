#ifndef DATAPERTURBATION_H
#define DATAPERTURBATION_H

#include <vector>
#include "GraphManager.h"
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <cmath>
#include <queue>
#include <iostream>
#include <tuple>

class DataPerturbation{
public:
    std::vector<std::unordered_set<uint32_t>> oriLevel,N_bar,N_star;
    std::vector<std::vector<uint32_t>> oriLevelVec;
    std::unordered_map<uint32_t,double> oriFreqTable;
    std::vector<std::unordered_map<uint32_t,std::unordered_set<uint32_t>>> NViVec;
    std::vector<uint32_t> V0Bar;
    double oriCount = 0,fullTime = 0, ESS = 0, normConst = 0;

public:

    DataPerturbation() = default;
    ~DataPerturbation() = default;

    void neiChecking(GraphManager &gm,uint32_t source,std::unordered_map<uint32_t,std::vector<uint32_t>> &srcNeighMap,double &A_Coeff); 

    void oriFreqTab(GraphManager &gm,uint32_t attr);

    void oriVec();

    bool fullMatch(GraphManager &gm,uint32_t srcNode,int currLevel,std::vector<uint32_t> backEdges,int maxLen,uint32_t attr);

    void prepSpace();

    void findChildren(GraphManager &gm,uint32_t srcNode,int currLevel,std::vector<uint32_t> backEdges,std::unordered_set<uint32_t> &children);

    bool findPeers(GraphManager &gm,uint32_t srcNode,uint32_t attr,uint32_t edgePrev);

    void NeiVi(GraphManager &gm);

    void NBarInitImp(GraphManager &gm, int level,std::unordered_set<uint32_t> &candidates,uint32_t attr,uint32_t edgeNext);

    void NBarMid(GraphManager &gm, int level,std::unordered_set<uint32_t> &candidates,uint32_t edgePrevFwd,uint32_t edgePrev,uint32_t edgeNext);

    void partialMatchGenerator(GraphManager &gm,std::vector<uint32_t> backEdges,std::vector<uint32_t> oriEdges,uint32_t cxtNode,uint32_t attr);

    void HRCGen(GraphManager &gm,int level,std::vector<uint32_t> backEdges);

    double proxLeftJCD(GraphManager &gm,std::vector<uint32_t> sourceNeiNi,uint32_t u_i,
                        std::unordered_map<uint32_t,std::unordered_map<uint32_t,double>> &levelRelMapCurr);

    double proxRightBase(GraphManager &gm,uint32_t u_i,uint32_t source,uint32_t edgePrev,
                            std::unordered_map<uint32_t,std::vector<uint32_t>> srcNeighMap,
                            std::unordered_map<uint32_t,std::unordered_map<uint32_t,double>> &levelRelMap,
                            int level,double A_coeff); 

    double dataPtbSSCal(GraphManager &gm,std::vector<uint32_t> backEdges,uint32_t validEnd,int currLevel,uint32_t attr,uint32_t attrVal,double oriScore);

    void exact(GraphManager &gm, std::string &pattern, double oriScore, double timeLimit, uint32_t attr, 
                uint32_t attrVal, uint32_t cxtNode, int index, std::string opFileName, std::string prefix);

    void sample(GraphManager &gm, std::string &pattern, double oriScore, double timeLimit, double sampleRatio, uint32_t attr, 
                uint32_t attrVal, uint32_t cxtNode, int index, std::string opFileName, std::string prefix);

    void printExact(GraphManager &gm,double trueScore, double oriScore,
                    uint32_t caSourceGlb,uint32_t caEndGlb,double proxGlb, double sScoreGlb,int caGlblevel,
                    double fullTime,int index, std::string prefix,std::string opFilename,int preSize,int aftSize);

    void printSample(GraphManager &gm,double estScore,double oriScore,
                        uint32_t caSourceSp,uint32_t caEndSp,double prox, double sScoreSP,int caSPlevel,double sampleTime,
                        int index,std::string prefix,std::string opFilename,int sampleSize,double sampleRatio,int aftSize);
};
#endif