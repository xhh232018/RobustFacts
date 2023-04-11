#include "QueryPerturbation.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <random>
#include <numeric>
#include <math.h>

using namespace std;

void QueryPerturbation::neiChecking(GraphManager &gm,uint32_t source,std::unordered_map<uint32_t,std::vector<uint32_t>> &srcNeighMap,double &A_Coeff){
    double acoeff = 0;
    auto edgeStart = gm.nodes[source];
    int adj_sz = gm.nodes[source + 1] - gm.nodes[source];
    for(int i = 0;i < adj_sz;++i){
        auto edge = gm.edges[edgeStart + i];
        auto etype = extractHigh32bits(edge);
        auto end_node = extractLow32bits(edge);
        //srcNeiTmp[etype].insert(end_node);
        srcNeighMap[etype].push_back(end_node);
        uint32_t eidRvs;
        if(etype % 2 == 1 ) eidRvs = etype - 1;
        else eidRvs = etype + 1;
        if(gm.degStats[end_node][eidRvs] > 1) acoeff += 2.0/gm.degStats[end_node][eidRvs];
    }

    for(auto &kv_pair:srcNeighMap){
        sort(kv_pair.second.begin(),kv_pair.second.end());
    }
    A_Coeff = acoeff;
}

double QueryPerturbation::proxCal(GraphManager &gm,uint32_t ptbNode,std::unordered_map<uint32_t,std::vector<uint32_t>> cxtNeighMap,double &A_Coeff){
    //prune nodes(CN <= 2)
    double proximity,B_coeff,M_coeff = 0;
    int cnt = 0,th = 2;
    std::unordered_map<uint32_t,std::vector<uint32_t>> sampleNodeNeighMap;
    neiChecking(gm,ptbNode,sampleNodeNeighMap,B_coeff);
    for(auto kv_pair:cxtNeighMap){
        if(!sampleNodeNeighMap.count(kv_pair.first)) continue;
        std::vector<uint32_t> neiShare;
        std::set_intersection(kv_pair.second.begin(), kv_pair.second.end(),
                                      sampleNodeNeighMap[kv_pair.first].begin(), sampleNodeNeighMap[kv_pair.first].end(),
                                      std::back_inserter(neiShare));
        if(neiShare.empty()) continue;
        uint32_t eidRvs;
        if(kv_pair.first % 2 == 1 ) eidRvs = kv_pair.first - 1;
        else eidRvs = kv_pair.first + 1;
        for(auto nodeS:neiShare){
            M_coeff += 2.0/gm.degStats[nodeS][eidRvs];
            cnt++;
        }
    }
    if(cnt <= th) proximity = 0;
    else proximity =  M_coeff / (A_Coeff + B_coeff - M_coeff);
    
    return proximity;
}

double QueryPerturbation::queryPtbSSCal(GraphManager &gm,uint32_t ptb,std::vector<uint32_t> backEdges,uint32_t attr, uint32_t attrVal){
    std::unordered_set<uint32_t> levelProp{ptb},peers;
    double newScore;
    for(int i = 0;i < backEdges.size();++i){
        if(i == backEdges.size() - 1){
            uint32_t lastRvs;
            if(backEdges[i] % 2 == 1) lastRvs = backEdges[i] -1;
            else lastRvs = backEdges[i] + 1;
            bool same = (attr == lastRvs);
            for(auto node:levelProp){
                auto key = combine2u32(node,backEdges[i]);
                auto ePos = gm.snidEid2pos[key];
                for(int k = ePos;k < gm.nodes[node+ 1];++k){
                    uint32_t eid = extractHigh32bits(gm.edges[k]);
                    if (eid != backEdges[i]) break;
                    uint32_t end_node = extractLow32bits(gm.edges[k]);
                    if((!same && gm.degStats[end_node].count(attr))||(same && gm.degStats[end_node][attr] > 1)) peers.insert(end_node);
                }
            }
        }
        else{
            uint32_t lastRvs;
            if(backEdges[i + 1] % 2 == 1) lastRvs = backEdges[i + 1] -1;
            else lastRvs = backEdges[i + 1] + 1;
            bool same = (backEdges[i] == lastRvs);
            std::unordered_set<uint32_t> lvpptmp;
            for(auto node:levelProp){
                auto key = combine2u32(node,backEdges[i]);
                auto ePos = gm.snidEid2pos[key];
                for(int k = ePos;k < gm.nodes[node+ 1];++k){
                    uint32_t eid = extractHigh32bits(gm.edges[k]);
                    if (eid != backEdges[i]) break;
                    uint32_t end_node = extractLow32bits(gm.edges[k]);
                    if((!same && gm.degStats[end_node].count(backEdges[i+1]))||(same && gm.degStats[end_node][backEdges[i+1]] > 1))lvpptmp.insert(end_node);
                }
            }
            levelProp = lvpptmp;
        }
    }

    std::unordered_map<uint32_t,double> freqTable;
    int peerCnt = 0;
    double total = 0;
    for(auto peer:peers){
        auto key = combine2u32(peer,attr);
        auto ePos = gm.snidEid2pos[key];
        for(int j = ePos;j < gm.nodes[peer + 1];++j){
            auto edge = gm.edges[j];
            uint32_t etype = extractHigh32bits(edge);
            if(etype != attr) break;
            uint32_t end_node = extractLow32bits(edge);
            if(freqTable.count(end_node)) freqTable[end_node]++;
            else freqTable[end_node] = 1;
            total++;
        }
        peerCnt++;
    }

    double valid = 0;
    for(auto kv_pair:freqTable){
        if(kv_pair.second > freqTable[attrVal]) valid += kv_pair.second;
    }

    if(peerCnt <= 10) newScore = -1;
    else newScore = valid/total;

    return newScore;
}

void QueryPerturbation::pSpaceFromAttr(GraphManager &gm,std::vector<uint32_t> oriEdges,uint32_t attr, uint32_t attrVal,uint32_t cxtNode){
    std::unordered_set<uint32_t> levelProp;
    // attrV -> PeerCandidate
    uint32_t attrRvs;
    if(attr % 2 == 1) attrRvs = attr - 1;
    else attrRvs = attr + 1;
    bool same = (attr == oriEdges.front());
    auto key = combine2u32(attrVal,attrRvs);
    auto ePos = gm.snidEid2pos[key];
    for(int k = ePos;k < gm.nodes[attrVal+ 1];++k){
        uint32_t eid = extractHigh32bits(gm.edges[k]);
        if (eid != attrRvs) break;
        uint32_t end_node = extractLow32bits(gm.edges[k]);
        if((!same && gm.degStats[end_node].count(oriEdges.front()))||(same && gm.degStats[end_node][oriEdges.front()] > 1)){
            levelProp.insert(end_node);
        } 
    }

    // PeerCandidate -> ptbSpace
    for(int i = 0;i < oriEdges.size();++i){
        if(i == oriEdges.size() - 1){
            
            for(auto node:levelProp){
                auto key = combine2u32(node,oriEdges[i]);
                auto ePos = gm.snidEid2pos[key];
                for(int k = ePos;k < gm.nodes[node+ 1];++k){
                    uint32_t eid = extractHigh32bits(gm.edges[k]);
                    if (eid != oriEdges[i]) break;
                    uint32_t end_node = extractLow32bits(gm.edges[k]);
                    if(end_node != cxtNode) this->pSpace.insert(end_node);
                }
            }
        }
        else{
            uint32_t currRvs;
            if(oriEdges[i] % 2 == 1) currRvs = oriEdges[i] -1;
            else currRvs = oriEdges[i] + 1;
            bool same = (oriEdges[i + 1] == currRvs);
            std::unordered_set<uint32_t> lvTmp;
            for(auto node:levelProp){
                auto key = combine2u32(node,oriEdges[i]);
                auto ePos = gm.snidEid2pos[key];
                for(int k = ePos;k < gm.nodes[node+ 1];++k){
                    uint32_t eid = extractHigh32bits(gm.edges[k]);
                    if (eid != oriEdges[i]) break;
                    uint32_t end_node = extractLow32bits(gm.edges[k]);
                    if((!same && gm.degStats[end_node].count(oriEdges[i+1]))||(same && gm.degStats[end_node][oriEdges[i+1]] > 1)){
                        lvTmp.insert(end_node);
                    }
                }
            }
            levelProp = lvTmp;
        }
    }
}

void QueryPerturbation::bySampleSizeBase(GraphManager &gm, std::string &pattern, double oriScore, double timeLimit, uint32_t attr, 
                                            uint32_t attrVal, uint32_t cxtNode, int index, std::string opFileName, std::string prefix){
    std::istringstream iss(pattern);
    std::string token;
    std::vector<std::string> tmp;
    std::vector<uint32_t> oriEdges,backEdges;
    std::unordered_map<uint32_t,std::vector<uint32_t>> cxtNeighMap;
    double fullTime = 0;
    double A_coeff;

    while(iss >> token){
        tmp.push_back(token);
    }
    for(int i = 0;i<tmp.size();++i){
        if(i % 3 == 1) oriEdges.push_back(stoul(tmp[i]));
    }
    for(int i = oriEdges.size() -1 ;i >= 0; --i){
        if(oriEdges[i] % 2 == 1) backEdges.push_back(oriEdges[i]-1);
        else backEdges.push_back(oriEdges[i] + 1);
    }

    neiChecking(gm,cxtNode,cxtNeighMap,A_coeff);

    auto s_time = getTime();
    pSpaceFromAttr(gm,oriEdges,attr,attrVal,cxtNode);
    auto e_time = getTime();
    cout << "pspace size: " << this -> pSpace.size() << endl;
    printf("pspace time = %.3lf ms.\n", getInterval(s_time, e_time));

    fullTime += getInterval(s_time, e_time);

    if(this -> pSpace.empty()) return;

    std::unordered_map<uint32_t, double> ssMapSp;
    std::vector<double> weightVec;
    std::vector<uint32_t> ptbSpaceVec;

    s_time = getTime();
    for(auto ptbNode:this->pSpace){
        if(ptbNode == cxtNode) continue;
        double proximity  = proxCal(gm,ptbNode,cxtNeighMap,A_coeff);
        if(proximity == 0) continue;
        ptbSpaceVec.push_back(ptbNode);
        weightVec.push_back(proximity);
    }
    e_time = getTime();
    printf("relevance time = %.3lf ms.\n", getInterval(s_time, e_time));
    if(weightVec.empty()) return;
    fullTime += getInterval(s_time, e_time);

    this->givenScore = oriScore;

    double meanBase = 0;
    double caGlb = -1,caLhGlb,caScoreGlb;
    std::vector<double> estScoreVec,proxVecOrd;
    double totalWeightFull = 0;
    uint32_t ptbMaxGlb;


    //exact
    double totalValid = 0;
    for(int i = 0;i<ptbSpaceVec.size();++i){
        auto ss = getTime();
        uint32_t ptbNode = ptbSpaceVec[i];

        double ssTmp = queryPtbSSCal(gm,ptbNode,backEdges,attr,attrVal);
        if(ssTmp < 0) continue;
        double diff = this->givenScore - ssTmp;
        double ptbScoreTmp = diff * weightVec[i];
        meanBase += ssTmp * weightVec[i];
        totalWeightFull += weightVec[i];

        if(ptbScoreTmp > caGlb){
            caGlb = ptbScoreTmp;
            ptbMaxGlb = ptbNode;
            caLhGlb = weightVec[i];
            caScoreGlb = ssTmp;
        }
        auto ee = getTime();
        totalValid++;
        fullTime += getInterval(ss, ee);
        if(fullTime > timeLimit) break;
    }

    meanBase /= totalWeightFull;
    printf("full time = %.3lf ms.\n", fullTime);

    //print
    printBase(gm,index, meanBase, oriScore, ptbMaxGlb, caLhGlb, caScoreGlb, fullTime,opFileName,prefix);
}

void QueryPerturbation::printBase(GraphManager &gm,int index,double trueScore,double oriScore,uint32_t caGlb,double proxGlb, double sScoreGlb,
                                    double fullTime,std::string opFileName, std::string prefix){
    ofstream myout;
    myout.open(opFileName,std::ofstream::app);

    double cherryPickMeas = (oriScore - trueScore) / trueScore;

    string caGlbIriOri = gm.nid2IRI[caGlb];
    string caGlbIriQ = caGlbIriOri.substr(31);

    myout <<prefix <<"," << index << ","<< trueScore << "," << oriScore << ","
          << cherryPickMeas << ","<< caGlbIriQ << "," << proxGlb << "," << sScoreGlb << "," << fullTime << "," << this->pSpace.size() <<endl;

    myout.close();
}

void QueryPerturbation::bySampleSizeExp(GraphManager &gm, std::string &pattern, double oriScore, double timeLimit, double sampleRatio, uint32_t attr, 
                                        uint32_t attrVal, uint32_t cxtNode, int index, std::string opFileName, std::string prefix){
    
    std::istringstream iss(pattern);
    std::string token;
    std::vector<std::string> tmp;
    std::vector<uint32_t> oriEdges,backEdges;
    std::unordered_map<uint32_t,std::vector<uint32_t>> cxtNeighMap;
    double A_coeff, sampleTime = 0;
    this->givenScore = oriScore;

    bool sample = true;

    while(iss >> token){
        tmp.push_back(token);
    }
    for(int i = 0;i<tmp.size();++i){
        if(i % 3 == 1) oriEdges.push_back(stoul(tmp[i]));
    }
    for(int i = oriEdges.size() -1 ;i >= 0; --i){
        if(oriEdges[i] % 2 == 1) backEdges.push_back(oriEdges[i]-1);
        else backEdges.push_back(oriEdges[i] + 1);
    }

    neiChecking(gm,cxtNode,cxtNeighMap,A_coeff);

    auto s_time = getTime();
    pSpaceFromAttr(gm,oriEdges,attr,attrVal,cxtNode);
    auto e_time = getTime();
    cout << "pspace size: " << this -> pSpace.size() << endl;
    printf("pspace time = %.3lf ms.\n", getInterval(s_time, e_time));
    sampleTime += getInterval(s_time, e_time);

    if(this -> pSpace.empty()) return;

    std::unordered_map<uint32_t, double> ssMapSp;
    std::vector<double> weightVec;
    std::vector<uint32_t> ptbSpaceVec;

    s_time = getTime();
    for(auto ptbNode:this->pSpace){
        if(ptbNode == cxtNode) continue;
        double proximity  = proxCal(gm,ptbNode,cxtNeighMap,A_coeff);
        if(proximity == 0) continue;
        ptbSpaceVec.push_back(ptbNode);
        weightVec.push_back(proximity);
    }
    e_time = getTime();
    printf("relevance time = %.3lf ms.\n", getInterval(s_time, e_time));
    sampleTime += getInterval(s_time, e_time);

    std::discrete_distribution<> d(weightVec.begin(),weightVec.end());
    int batch = sampleRatio * weightVec.size();
    if(batch < 30) batch = 30; 

    double caSp = -1,caLhSp,caScoreSp;
    double meanOrd = 0;
    double totalWeightOrd = 0;
    uint32_t ptbMaxOrd;
    std::random_device rd;
    std::mt19937 gen(rd());

    auto s_start = getTime();
    
    for(int i = 0; i < batch; ++i){
        auto s_sample = getTime();
        int sample = d(gen);
        uint32_t ptbNode = ptbSpaceVec[sample];

        double ssTmp;
        if(ssMapSp.count(ptbNode)) ssTmp = ssMapSp[ptbNode];
        else{
            ssTmp = queryPtbSSCal(gm,ptbNode,backEdges,attr,attrVal);
            ssMapSp[ptbNode] = ssTmp;
        }
        if(ssTmp < 0) continue;    

        double diff = this->givenScore - ssTmp;

        double caTmpOrd = diff * weightVec[sample];

        totalWeightOrd += weightVec[sample];
        meanOrd += ssTmp * weightVec[sample];

        if(caTmpOrd > caSp){
            caSp = caTmpOrd;
            ptbMaxOrd = ptbNode;
            caLhSp = weightVec[sample];
            caScoreSp = ssTmp;
        }
        auto e_sample = getTime();
        sampleTime += getInterval(s_sample,e_sample);
        if(sampleTime > timeLimit) break; 
    }
    printf("sample time = %.3lf ms.\n", sampleTime);
    meanOrd /= totalWeightOrd;

    printSampleExp(gm,meanOrd,oriScore,ptbMaxOrd,caLhSp,caScoreSp,sampleTime,index,opFileName,prefix,batch,sampleRatio);
}

void QueryPerturbation::printSampleExp(GraphManager &gm,double estScore, double oriScore, uint32_t ca,
                                        double prox, double sScoreSP,double sampleTime,
                                        int index,std::string opFileName, std::string prefix,int sampleSize, double sampleRatio){
    ofstream myout;
    myout.open(opFileName,std::ofstream::app);

    double cherryPickMeas = (oriScore - estScore) / estScore;

    string caSPIriOri = gm.nid2IRI[ca];
    string caSPIriQ = caSPIriOri.substr(31);

    myout <<prefix <<"," << index << ","<< estScore << "," << oriScore << ","
          << cherryPickMeas << ","<< caSPIriQ << "," << prox << "," << sScoreSP << "," << sampleTime <<"," 
          << sampleSize << "," << sampleRatio <<endl;

    myout.close();
}