#include "DataPerturbation.h"
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <stack>
#include <random>
#include <numeric>
#include <math.h>
#include <tuple>

using namespace std;

void DataPerturbation::neiChecking(GraphManager &gm,uint32_t source,std::unordered_map<uint32_t,std::vector<uint32_t>> &srcNeighMap,double &A_Coeff){
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

void DataPerturbation::NeiVi(GraphManager &gm){
    for(int i = 0; i < this->oriLevel.size() - 1;++i){
        for(auto &nodeLastLevel:this->oriLevel[i]){
            int adj_sz = gm.nodes[nodeLastLevel + 1] - gm.nodes[nodeLastLevel];
            int edge_start = gm.nodes[nodeLastLevel];
            for(int j = 0;j<adj_sz;++j){
                //1-hop nei
                uint64_t edge = gm.edges[edge_start + j];
                uint32_t eid = extractHigh32bits(edge);
                uint32_t end_node = extractLow32bits(edge);
                this->NViVec[i][eid].insert(end_node);
            }
        }
    }
}

bool DataPerturbation::findPeers(GraphManager &gm,uint32_t srcNode,uint32_t attr,uint32_t edgePrev){
    bool valid = false;
    bool sameFlag = (attr == edgePrev);
    if((!sameFlag&&gm.degStats[srcNode].count(attr)) || (sameFlag && gm.degStats[srcNode][attr] > 1 )) valid = true;
    return valid;
}

void DataPerturbation::findChildren(GraphManager &gm,uint32_t srcNode,int currLevel,std::vector<uint32_t> backEdges,std::unordered_set<uint32_t> &children){
    //Only Next Level
    if(currLevel < backEdges.size() - 1){
        uint32_t edgeCurr = backEdges[currLevel],edgeNext;
        if(backEdges[currLevel + 1] % 2 ==1) edgeNext = backEdges[currLevel + 1] - 1 ;
        else edgeNext = backEdges[currLevel + 1] + 1 ;
        bool sameFlag = (edgeCurr == edgeNext);

        auto key = combine2u32(srcNode,edgeCurr);
        int ePos = gm.snidEid2pos[key];
        for(int j = ePos;j < gm.nodes[srcNode + 1];++j){
            auto eid = extractHigh32bits(gm.edges[j]);
            if (eid != edgeCurr) break;
            auto end_node = extractLow32bits(gm.edges[j]);
            if((!sameFlag&&gm.degStats[end_node].count(backEdges[currLevel + 1])) || (sameFlag && gm.degStats[end_node][backEdges[currLevel + 1]] > 1 )) children.insert(end_node);
        }
    }
    else{
        //last level
        uint32_t edgeCurr = backEdges[currLevel];
        auto key = combine2u32(srcNode,edgeCurr);
        int ePos = gm.snidEid2pos[key];
        for(int j = ePos;j < gm.nodes[srcNode + 1];++j){
            auto eid = extractHigh32bits(gm.edges[j]);
            if (eid != edgeCurr) break;
            auto end_node = extractLow32bits(gm.edges[j]);
            children.insert(end_node);
        }
    }
}

bool DataPerturbation::fullMatch(GraphManager &gm,uint32_t srcNode,int currLevel,std::vector<uint32_t> backEdges,int maxLen,uint32_t attr){
    bool valid = false;
    if(currLevel == maxLen){
        uint32_t edgePrev;
        if(backEdges.back()%2 == 1) edgePrev = backEdges.back() - 1;
        else edgePrev = backEdges.back() + 1;
        valid = findPeers(gm,srcNode,attr,edgePrev);
    }
    else{
        std::unordered_set<uint32_t> children;
        findChildren(gm,srcNode,currLevel,backEdges,children);
        if(children.empty()) valid = false;
        else{
            for(uint32_t child:children){
                if(this->oriLevel[currLevel + 1].count(child)) valid = true;
                else valid = valid | fullMatch(gm,child,currLevel + 1,backEdges,maxLen,attr);
            }
        }
    }
    if(valid) this->oriLevel[currLevel].insert(srcNode);
    return valid;
}

void DataPerturbation::oriFreqTab(GraphManager &gm,uint32_t attr){
    for(auto peer:this->oriLevel.back()){
        auto key = combine2u32(peer,attr);
        int ePos = gm.snidEid2pos[key];
        for(int j = ePos;j < gm.nodes[peer + 1];++j){
            auto eid = extractHigh32bits(gm.edges[j]);
            if (eid != attr) break;
            auto end_node = extractLow32bits(gm.edges[j]);
            if(this->oriFreqTable.count(end_node)) this->oriFreqTable[end_node]++;
            else this->oriFreqTable[end_node] = 1;
            this->oriCount++;
        }
    }
}

void DataPerturbation::oriVec(){
    for(int i = 0;i < this->oriLevel.size();++i){
        std::vector<uint32_t> tmpVec;
        for(auto node:this->oriLevel[i]){
            tmpVec.push_back(node);
        }
        sort(tmpVec.begin(),tmpVec.end());
        this->oriLevelVec.push_back(tmpVec);
    }
}

void DataPerturbation::prepSpace(){
    this->NViVec.resize(this->oriLevel.size());
    this->N_bar.resize(this->oriLevel.size());
    this->N_star.resize(this->oriLevel.size());
    /*this->N_Comp.resize(this->oriLevel.size());
    this->N_tilde.resize(this->oriLevel.size());
    this->N_hat.resize(this->oriLevel.size());
    this->N_CompVec.resize(this->oriLevel.size());
    this->N_tildeVec.resize(this->oriLevel.size());
    this->N_hatVec.resize(this->oriLevel.size());
    this->N_barVec.resize(this->oriLevel.size());*/
}

void DataPerturbation::NBarInitImp(GraphManager &gm, int level,std::unordered_set<uint32_t> &candidates,uint32_t attr,uint32_t edgeNext){
    bool sameFlag = (attr == edgeNext);
    for(auto &candidate: gm.typeId2Count[attr]){
        if((!sameFlag&&gm.degStats[candidate].count(edgeNext)) || (sameFlag && gm.degStats[candidate][edgeNext] > 1 )){
            if(!this->oriLevel.back().count(candidate)) this->V0Bar.push_back(candidate);
        }
    }
}

void DataPerturbation::NBarMid(GraphManager &gm, int level,std::unordered_set<uint32_t> &candidates,uint32_t edgePrevFwd,
                                uint32_t edgePrev,uint32_t edgeNext){
    bool same = (edgePrev == edgeNext);
    if(level == this->oriLevel.size() - 2) {
        for(auto &canSrc:this->V0Bar){
            auto key = combine2u32(canSrc,edgePrevFwd);
            int ePos =  gm.snidEid2pos[key];
            for(int k = ePos;k < gm.nodes[canSrc+ 1];++k){
                auto edge = gm.edges[k];
                uint32_t eid = extractHigh32bits(edge);
                if (eid != edgePrevFwd) break;
                auto end_node = extractLow32bits(edge);
                if((!same &&gm.degStats[end_node].count(edgeNext)) || (same &&gm.degStats[end_node][edgeNext] > 1)) candidates.insert(end_node);
            }
        }
    }
    else{
        for(auto &canSrc:this->N_bar[level + 1]){
            auto key = combine2u32(canSrc,edgePrevFwd);
            int ePos =  gm.snidEid2pos[key];
            for(int k = ePos;k < gm.nodes[canSrc+ 1];++k){
                auto edge = gm.edges[k];
                uint32_t eid = extractHigh32bits(edge);
                if (eid != edgePrevFwd) break;
                auto end_node = extractLow32bits(edge);
                if((!same &&gm.degStats[end_node].count(edgeNext)) || (same &&gm.degStats[end_node][edgeNext] > 1)) candidates.insert(end_node);
            }
        }
    }
}

void DataPerturbation::partialMatchGenerator(GraphManager &gm,std::vector<uint32_t> backEdges,std::vector<uint32_t> oriEdges,uint32_t cxtNode,uint32_t attr){
    for(int j = this->oriLevel.size() - 1; j > 0; j--){ //v0.1: j = orilevel-2, v1.0: j = orilevel - 1
        std::unordered_set<uint32_t> tmp;
        
        if(j == this->oriLevel.size() - 1){
            int backIdx = j;
            //int fwdIdx = this->oriLevel.size() - 1 - j;
            //NBarInit(gm,j,tmp,backEdges[backIdx],oriEdges[fwdIdx],attr);
            NBarInitImp(gm,j,tmp,attr,oriEdges.front());
            this->N_bar[j] = tmp;
        }

        else{
            int backIdx = j;
            int fwdIdx = this->oriLevel.size() - 1 - j;
            NBarMid(gm,j,tmp,oriEdges[fwdIdx - 1],backEdges[backIdx],oriEdges[fwdIdx]);
            this->N_bar[j] = tmp;
        }
    }
}

void DataPerturbation::HRCGen(GraphManager &gm,int level,std::vector<uint32_t> backEdges){
    std::unordered_set<uint32_t> &ViBar = this -> N_bar[level];
    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> &N_i = this->NViVec[level];
    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> &N_i1 = this->NViVec[level - 1];
    std::unordered_set<uint32_t> Vi1_hat,Vi_tilde,nodeSetEi;

    //Vi_tilde
    for(auto &v:ViBar){
        //head rel
        bool keepFlag = false,tail = false;
        for(auto &kv_pair:gm.degStats[v]){
            //commonNei
            if(N_i.count(kv_pair.first)){
                auto key = combine2u32(v,kv_pair.first);
                auto ePos = gm.snidEid2pos[key];
                for(int i = ePos;i < gm.nodes[v + 1];++i){
                    auto edge = gm.edges[i];
                    uint32_t etypeCurr = extractHigh32bits(edge);
                    if(etypeCurr != kv_pair.first) break;
                    uint32_t oneNeigh = extractLow32bits(edge);
                    if(N_i[kv_pair.first].count(oneNeigh)){
                        keepFlag = true;
                        break;
                    } 
                }
            }
            if(keepFlag){
                Vi_tilde.insert(v);
            } 
        }
    }

    uint32_t E_i;
    if(backEdges[level - 1]% 2 == 1)E_i = backEdges[level - 1] - 1;
    else E_i = backEdges[level - 1] + 1;

    for(auto &node:Vi_tilde){
        auto key = combine2u32(node,E_i);
        auto ePos = gm.snidEid2pos[key];
        for(int i = ePos;i < gm.nodes[node + 1];++i){
            auto edge = gm.edges[i];
            uint32_t etypeCurr = extractHigh32bits(edge);
            if(etypeCurr != E_i) break;
            uint32_t oneNeigh = extractLow32bits(edge);
            nodeSetEi.insert(oneNeigh);
        }
    }

    for(auto &node:nodeSetEi){
        bool keepFlag = false;
        for(auto kv_pair:gm.degStats[node]){
            //commonNei
            if(N_i1.count(kv_pair.first)){
                auto key = combine2u32(node,kv_pair.first);
                auto ePos = gm.snidEid2pos[key];
                for(int i = ePos;i < gm.nodes[node + 1];++i){
                    auto edge = gm.edges[i];
                    uint32_t etypeCurr = extractHigh32bits(edge);
                    if(etypeCurr != kv_pair.first) break;
                    uint32_t oneNeigh = extractLow32bits(edge);
                    if(N_i1[kv_pair.first].count(oneNeigh)){
                        keepFlag = true;
                        break;
                    } 
                }
            }
            if(keepFlag) break;
        }
        if(keepFlag) Vi1_hat.insert(node);
    }

    //cascading
    for(auto &v:Vi_tilde){
        
        auto key = combine2u32(v,E_i);
        auto ePos = gm.snidEid2pos[key];

        for(int i = ePos;i < gm.nodes[v + 1];++i){
            auto edge = gm.edges[i];
            uint32_t etypeCurr = extractHigh32bits(edge);
            if(etypeCurr != E_i) break;
            uint32_t vEi = extractLow32bits(edge);
            if(Vi1_hat.count(vEi)){
                this->N_star[level].insert(v);
                break;
            }
        }
    }
}

double DataPerturbation::proxLeftJCD(GraphManager &gm,std::vector<uint32_t> sourceNeiNi,uint32_t u_i,
                                    std::unordered_map<uint32_t,std::unordered_map<uint32_t,double>> &levelRelMapCurr){
    double leftLh  = 0,A_coeff;
    std::unordered_map<uint32_t,std::vector<uint32_t>> uiNei;
    neiChecking(gm,u_i,uiNei,A_coeff);

    for(auto n_i:sourceNeiNi){
        double tmpScore = 0;
        if(levelRelMapCurr[n_i].count(u_i)){
            leftLh += levelRelMapCurr[n_i][u_i];
            continue;
        }
        double B_coeff = 0,M_coeff = 0;
        std::unordered_map<uint32_t,std::vector<uint32_t>> sampleNodeNeighMap;
        neiChecking(gm,n_i,sampleNodeNeighMap,B_coeff);
        for(auto kv_pair:uiNei){
            if(!sampleNodeNeighMap.count(kv_pair.first)) continue;
            std::vector<uint32_t> neiShare;
            std::set_intersection(kv_pair.second.begin(), kv_pair.second.end(),
                                    sampleNodeNeighMap[kv_pair.first].begin(), sampleNodeNeighMap[kv_pair.first].end(),
                                    std::back_inserter(neiShare));
            uint32_t eidRvs;
            if(kv_pair.first % 2 == 1 ) eidRvs = kv_pair.first - 1;
            else eidRvs = kv_pair.first + 1;
            for(auto nodeS:neiShare){
                M_coeff += 2.0/gm.degStats[nodeS][eidRvs];
            }
        }

        tmpScore = M_coeff / (A_coeff + B_coeff - M_coeff);
        leftLh += tmpScore;
        levelRelMapCurr[n_i][u_i] = tmpScore;
    }
    return leftLh;

}

double DataPerturbation::proxRightBase(GraphManager &gm,uint32_t u_i,uint32_t source,uint32_t edgePrev,
                                        std::unordered_map<uint32_t,std::vector<uint32_t>> srcNeighMap,
                                        std::unordered_map<uint32_t,std::unordered_map<uint32_t,double>> &levelRelMap,
                                        int level,double A_coeff){
    double rightLh = 0;
    
    uint32_t edgePrevBack;
    if(edgePrev % 2 == 1) edgePrevBack = edgePrev - 1;
    else edgePrevBack = edgePrev + 1;
    std::vector<uint32_t> sampleSpace;


    auto key = combine2u32(u_i,edgePrevBack);
    int ePos =  gm.snidEid2pos[key];
    
    for(int k = ePos;k < gm.nodes[u_i+ 1];++k){
        auto edge = gm.edges[k];
        uint32_t eid = extractHigh32bits(edge);
        if (eid != edgePrevBack) break;
        auto end_node = extractLow32bits(edge);
        
        if(!levelRelMap.empty() && levelRelMap[source].count(end_node)) rightLh += levelRelMap[source][end_node];
        else sampleSpace.push_back(end_node);
    }

    //cout<<"Sample Space Size: "<< sampleSpace.size() <<endl;

    for(auto end_node:sampleSpace){
        double B_coeff,M_coeff = 0;
        std::unordered_map<uint32_t,std::vector<uint32_t>> sampleNodeNeighMap;
        neiChecking(gm,end_node,sampleNodeNeighMap,B_coeff);
        for(auto kv_pair:srcNeighMap){
            if(!sampleNodeNeighMap.count(kv_pair.first)) continue;
            std::vector<uint32_t> neiShare;
            std::set_intersection(kv_pair.second.begin(), kv_pair.second.end(),
                                      sampleNodeNeighMap[kv_pair.first].begin(), sampleNodeNeighMap[kv_pair.first].end(),
                                      std::back_inserter(neiShare));
            uint32_t eidRvs;
            if(kv_pair.first % 2 == 1 ) eidRvs = kv_pair.first - 1;
            else eidRvs = kv_pair.first + 1;
            for(auto nodeS:neiShare){
                M_coeff += 2.0/gm.degStats[nodeS][eidRvs];
            }
        }
        double tmpRel = M_coeff / (A_coeff + B_coeff - M_coeff);
        rightLh  +=  tmpRel;
    }
    return rightLh;
}

double DataPerturbation::dataPtbSSCal(GraphManager &gm,std::vector<uint32_t> backEdges,uint32_t validEnd,int endLevel,uint32_t attr,uint32_t attrVal,double oriScore){
    std::unordered_set<uint32_t> levelProp{validEnd},peersNew;
    double newScore;
    for(int i = endLevel;i < backEdges.size();++i){
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
                    if(this->oriLevel.back().count(end_node)) continue;
                    if((!same && gm.degStats[end_node].count(attr))||(same && gm.degStats[end_node][attr] > 1)) peersNew.insert(end_node);
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

    if(peersNew.empty()) newScore = oriScore;
    else{
        std::unordered_map<uint32_t,double> freqTable,oriFreqCp = this->oriFreqTable;
        double total = 0;
        for(auto peer:peersNew){
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
        }
        for(auto kv_pair:freqTable){
            if(oriFreqCp.count(kv_pair.first)) oriFreqCp[kv_pair.first]+= kv_pair.second;
            else oriFreqCp[kv_pair.first] = kv_pair.second;
        }
        double valid = 0;
        for(auto kv_pair:oriFreqCp){
            if(kv_pair.second > oriFreqCp[attrVal]) valid += kv_pair.second;
        }
        newScore = valid/(this->oriCount + total);
    }
    return newScore;
}

void DataPerturbation::exact(GraphManager &gm, std::string &pattern, double oriScore, double timeLimit, uint32_t attr, 
                                uint32_t attrVal, uint32_t cxtNode, int index, std::string opFileName, std::string prefix){

    auto s_time = getTime();
    std::istringstream iss(pattern);
    std::string token;
    std::vector<std::string> tmp;
    std::vector<uint32_t> oriEdges,backEdges;
    int spaceSize = 0;

    while(iss >> token){
        tmp.push_back(token);
    }
    for(int i = 0;i<tmp.size();++i){
        if(i % 3 == 1) oriEdges.push_back(stoul(tmp[i]));
    }
    bool sampling = true;
    for(int i = oriEdges.size() -1 ;i >= 0; --i){
        if(oriEdges[i] % 2 == 1) backEdges.push_back(oriEdges[i]-1);
        else backEdges.push_back(oriEdges[i] + 1);
    }

    this->oriLevel.resize(backEdges.size() + 1);
    
    bool valid = fullMatch(gm,cxtNode,0,backEdges,backEdges.size(),attr);
    auto e_time = getTime();
    this -> fullTime += getInterval(s_time,e_time);
    if(!valid) return;

    oriFreqTab(gm,attr);

    oriVec();

    prepSpace();

    NeiVi(gm);

    s_time = getTime();
    partialMatchGenerator(gm,backEdges,oriEdges,cxtNode,attr);

    for(int i = 1; i <= oriLevel.size()-2;++i){
        HRCGen(gm,i,backEdges);
    }
    e_time = getTime();
    this -> fullTime += getInterval(s_time,e_time);
    if(this -> N_star[this->oriLevel.size() - 2].empty()) return;

    std::unordered_map<uint32_t,std::unordered_map<uint32_t,double>> levelRelMap;
    double CAglb = -1,CAGlbSS,CAGlbLh;
    uint32_t CAGlbSource,CAGlbEnd;
    int CAGlbLevel;//counterargument
    
    for(int i = 1;i <= backEdges.size()-1;++i){
        if(this -> fullTime > timeLimit) break;
        std::unordered_map<uint32_t,std::unordered_map<uint32_t,double>> levelRelMapCurr;
        std::unordered_map<uint32_t,double> levelEndSSMap;
        
        for(auto source:oriLevel[i - 1]){
            
            std::vector<uint32_t> sourceNeiNi;
            std::unordered_set<uint32_t> N_i = this -> oriLevel[i];
            std::vector<uint32_t> N_ibarSrc;
            N_ibarSrc.insert(N_ibarSrc.end(),this->N_star[i].begin(), this->N_star[i].end());
            cout<<"Vi star size: " << N_ibarSrc.size() << endl;
            std::unordered_map<uint32_t,std::vector<uint32_t>> srcNeighMap;
            double A_coeff;
            neiChecking(gm,source,srcNeighMap,A_coeff);

            if(i == 1){ //first level
                
                sourceNeiNi = this->oriLevelVec[i];
                for(auto u_i:N_ibarSrc){
                    auto s_time = getTime();
                    double leftLh = proxLeftJCD(gm,sourceNeiNi,u_i,levelRelMapCurr); // head relevance
                    double rightLh = proxRightBase(gm,u_i,source,backEdges[i-1],srcNeighMap,levelRelMap,i,A_coeff); //tailRelevance

                    double endSS = dataPtbSSCal(gm,backEdges,u_i,i,attr,attrVal,oriScore);

                    double p_x = leftLh * rightLh;
                    this->ESS += endSS * p_x;
                    double CAGlbValue = (oriScore - endSS) * p_x;
                    this->normConst += p_x;

                    if(CAGlbValue > CAglb){
                        CAglb = CAGlbValue;
                        CAGlbSS = endSS;
                        CAGlbSource = source;
                        CAGlbEnd = u_i;
                        CAGlbLh = p_x;
                        CAGlbLevel = i;
                    }

                    auto e_time = getTime();
                    this->fullTime += getInterval(s_time,e_time);
                    if(this -> fullTime > timeLimit) break;
                }
                if(this -> fullTime > timeLimit) break;
                levelRelMap = levelRelMapCurr;
                levelRelMapCurr.clear();
            }
            else{
                //sourceNei in V_i
                
                auto key = combine2u32(source,backEdges[i - 1]);
                auto ePos = gm.snidEid2pos[key];
                for(int j = ePos;j<gm.nodes[source + 1];++j){
                    //1-hop nei
                    auto edge = gm.edges[j];
                    uint32_t etypeCurr = extractHigh32bits(edge);
                    if(etypeCurr != backEdges[i - 1]) break;
                    uint32_t oneNeigh = extractLow32bits(edge);
                    if(N_i.count(oneNeigh)) sourceNeiNi.push_back(oneNeigh);
                }

                for(auto u_i:N_ibarSrc){
                    if(this -> fullTime > timeLimit) break;
                    if(u_i == source) continue;

                    auto start = getTime();
                    double leftLh = proxLeftJCD(gm,sourceNeiNi,u_i,levelRelMapCurr);
                    auto end = getTime();
                    this -> fullTime += getInterval(start,end);
                    if(leftLh == 0) continue;

                    start = getTime();
                    double rightLh = proxRightBase(gm,u_i,source,backEdges[i-1],srcNeighMap,levelRelMap,i,A_coeff);
                    end = getTime();
                    this -> fullTime += getInterval(start,end);
                    if(rightLh == 0) continue;

                    start = getTime();
                    double endSS;
                    if(levelEndSSMap.count(u_i)) endSS = levelEndSSMap[u_i];
                    else{
                        endSS = dataPtbSSCal(gm,backEdges,u_i,i,attr,attrVal,oriScore);
                        levelEndSSMap[u_i] = endSS;
                    }

                    double p_x = leftLh * rightLh;
                    ESS += endSS * p_x;
                    double CAGlbValue = (oriScore - endSS) * p_x;
                    normConst += p_x;

                    if(CAGlbValue > CAglb){
                        CAglb = CAGlbValue;
                        CAGlbSS = endSS;
                        CAGlbSource = source;
                        CAGlbEnd = u_i;
                        CAGlbLh = p_x;
                        CAGlbLevel = i;
                    }
                    
                    end = getTime();

                    this->fullTime += getInterval(start,end);
                    if(this -> fullTime > timeLimit) break;
                }
                if(this -> fullTime > timeLimit) break;
            }
        }
        levelRelMap = levelRelMapCurr;
        levelRelMapCurr.clear();
        if(this -> fullTime > timeLimit) break;
    }

    this -> ESS /= this -> normConst;

    int nstar = 0;
    int nbar = 0;

    for(int i = 0;i < this->N_star.size() - 2;++i){
        nbar += this->oriLevel[i].size() * this->N_bar[i+1].size();
        nstar += this->oriLevel[i].size() * this->N_star[i+1].size();
    }

    //printExact
    printExact(gm,this->ESS,oriScore,CAGlbSource,CAGlbEnd,CAGlbLh,CAGlbSS,CAGlbLevel,this -> fullTime,index,prefix,opFileName,nbar,nstar);
}

void DataPerturbation::printExact(GraphManager &gm,double trueScore, double oriScore,
                                    uint32_t caSourceGlb,uint32_t caEndGlb,double proxGlb, double sScoreGlb,int caGlblevel,
                                    double fullTime,int index, std::string prefix,std::string opFilename,int preSize,int aftSize){
    ofstream myout;
    myout.open(opFilename,std::ofstream::app);

    double cherryPickMeas = (oriScore - trueScore) / trueScore;


    string caSourceGlbIriOri = gm.nid2IRI[caSourceGlb];
    string caSourceGlbIriQ = caSourceGlbIriOri.substr(31);

    string caEndGlbIriOri = gm.nid2IRI[caEndGlb];
    string caEndGlbIriQ = caEndGlbIriOri.substr(31);

    myout << prefix <<"," << index << "," << trueScore << "," << oriScore << "," << cherryPickMeas << "," << caSourceGlbIriQ << "," <<caEndGlbIriQ << 
          "," << proxGlb << "," << sScoreGlb  << "," << caGlblevel << "," << fullTime <<"," <<preSize <<","<<aftSize << endl;

    myout.close();
}

void DataPerturbation::printSample(GraphManager &gm,double estScore,double oriScore,
                                    uint32_t caSourceSp,uint32_t caEndSp,double prox, double sScoreSP,int caSPlevel,double sampleTime,
                                    int index,std::string prefix,std::string opFilename,int sampleSize,double sampleRatio,int aftSize){
    ofstream myout;
    myout.open(opFilename,std::ofstream::app);

    double cherryPickMeas = (oriScore - estScore) / estScore;
    
    string caSourceIriOri = gm.nid2IRI[caSourceSp];
    string caSourceIriQ = caSourceIriOri.substr(31);

    string caEndIriOri = gm.nid2IRI[caEndSp];
    string caEndIriQ = caEndIriOri.substr(31);

    myout <<prefix <<"," << index << "," <<estScore << ","<< oriScore << "," << cherryPickMeas <<","<< caSourceIriQ <<"," 
          << caEndIriQ << "," << prox << ","<< sScoreSP  << caSPlevel << ","<< sampleTime  <<","<< aftSize <<"," << sampleSize << "," << sampleRatio <<endl;

    myout.close();
}

void DataPerturbation::sample(GraphManager &gm, std::string &pattern, double oriScore, double timeLimit, double sampleRatio, uint32_t attr, 
                                uint32_t attrVal, uint32_t cxtNode, int index, std::string opFileName, std::string prefix){

    auto s_time = getTime();
    std::istringstream iss(pattern);
    std::string token;
    std::vector<std::string> tmp;
    std::vector<uint32_t> oriEdges,backEdges;
    int spaceSize = 0;

    while(iss >> token){
        tmp.push_back(token);
    }
    for(int i = 0;i<tmp.size();++i){
        if(i % 3 == 1) oriEdges.push_back(stoul(tmp[i]));
    }
    bool sampling = true;
    for(int i = oriEdges.size() -1 ;i >= 0; --i){
        if(oriEdges[i] % 2 == 1) backEdges.push_back(oriEdges[i]-1);
        else backEdges.push_back(oriEdges[i] + 1);
    }

    this->oriLevel.resize(backEdges.size() + 1);
    
    bool valid = fullMatch(gm,cxtNode,0,backEdges,backEdges.size(),attr);
    auto e_time = getTime();
    this -> fullTime += getInterval(s_time,e_time);
    if(!valid) return;

    oriFreqTab(gm,attr);

    oriVec();

    prepSpace();

    NeiVi(gm);

    s_time = getTime();

    partialMatchGenerator(gm,backEdges,oriEdges,cxtNode,attr);

    for(int i = 1; i <= oriLevel.size()-2;++i){
        HRCGen(gm,i,backEdges);
    }
    e_time = getTime();
    this -> fullTime += getInterval(s_time,e_time);
    if(this -> N_star[this->oriLevel.size() - 2].empty()) return;


    double nstar = 0;

    for(int i = 0;i < this->N_star.size() - 2;++i){
        nstar += this->oriLevel[i].size() * this->N_star[i+1].size();
    }

    int batch = sampleRatio * nstar;
    if(batch < 30) batch = 30; 

    std::default_random_engine rng;
    std::vector<uint32_t> sampleLevelVec;
    for(int i = 0;i < this->N_star.size()-1;++i){
        if(!N_star[i].empty()) sampleLevelVec.push_back(i);
    }
    std::random_device rd;
    std::mt19937 g(rd());

    std::vector<std::unordered_map<uint32_t,std::unordered_map<uint32_t,double>>> levelRelMapSPVec;
    std::vector<std::unordered_map<uint32_t,double>> spLevelSS;

    levelRelMapSPVec.resize(this->oriLevel.size());
    spLevelSS.resize(this->oriLevel.size());
    double totalWeightOrd = 0, meanOrd = 0;

    int CASpLevel;
    double caSP = -1,CASpSS,CASpLh;//counterargument
    uint32_t CASpSource,CASpEnd;

    for(int k = 0;k < batch;++k){
       
        std::shuffle(sampleLevelVec.begin(), sampleLevelVec.end(),g);
        std::vector<uint32_t> sourceNeiNi,sampleSp;
        int endLevel = sampleLevelVec.front();
        std::vector<uint32_t> N_i1Curr = this->oriLevelVec[endLevel - 1];
        std::uniform_int_distribution<int> unidSource(0,N_i1Curr.size() - 1);
        int sourceIdx = unidSource(rng);
        uint32_t source = this->oriLevelVec[endLevel - 1][sourceIdx];
        std::unordered_map<uint32_t,std::vector<uint32_t>> srcNeighMap;
        double A_coeff;
        neiChecking(gm,source,srcNeighMap,A_coeff);

        if(endLevel == 1){
            sourceNeiNi = this->oriLevelVec[endLevel];
            sampleSp.insert(sampleSp.end(),this->N_star[endLevel].begin(), this->N_star[endLevel].end());
        }
        else{
            auto key = combine2u32(source,backEdges[endLevel - 1]);
            auto ePos = gm.snidEid2pos[key];
            for(int j = ePos;j<gm.nodes[source + 1];++j){
                //1-hop nei
                auto edge = gm.edges[j];
                uint32_t etypeCurr = extractHigh32bits(edge);
                if(etypeCurr != backEdges[endLevel - 1]) break;
                uint32_t oneNeigh = extractLow32bits(edge);
                if(this->oriLevel[endLevel].count(oneNeigh)) sourceNeiNi.push_back(oneNeigh);
            }
            sampleSp.insert(sampleSp.end(),this->N_star[endLevel].begin(), this->N_star[endLevel].end());
        }

        std::uniform_int_distribution<int> unidEnd(0,sampleSp.size() - 1);
        int endIdx = unidEnd(g);
        uint32_t end = sampleSp[endIdx];   
        
        auto s_rlh =  getTime();
        double rightLh = proxRightBase(gm,end,source,backEdges[endLevel-1],srcNeighMap,levelRelMapSPVec[endLevel - 1],endLevel,A_coeff);
        auto e_rlh =  getTime();
        this ->fullTime += getInterval(s_rlh, e_rlh);
        if(rightLh == 0) continue;

        auto s_llh =  getTime();
        double leftLh = proxLeftJCD(gm,sourceNeiNi,end,levelRelMapSPVec[endLevel]);
        auto e_llh =  getTime();     
        this ->fullTime += getInterval(s_llh, e_llh);
        if(leftLh == 0) continue;  

        auto m_start =  getTime();
        double p_x = leftLh * rightLh;
        double q_x = 1.0/(sampleLevelVec.size() * 1.0 * (this->oriLevelVec[endLevel - 1].size() * 1.0) * (this->N_star[endLevel].size() * 1.0));
        double weight = p_x / q_x;
        totalWeightOrd += weight;       

        double endSS;
        if(spLevelSS[endLevel].count(end)) endSS = spLevelSS[endLevel][end];
        else{
                
            endSS = dataPtbSSCal(gm,backEdges,end,endLevel,attr,attrVal,oriScore);
            spLevelSS[endLevel][end] = endSS;
        }
        meanOrd += endSS * weight;
        double CAValue = (oriScore - endSS) * p_x;
        if(CAValue > caSP){
            caSP = CAValue;
            CASpSource = source;
            CASpEnd = end;
            CASpSS = endSS;
            CASpLh = p_x;
            CASpLevel = endLevel;
        }
        auto m_end = getTime();
        this -> fullTime += getInterval(m_start,m_end);

        if(this -> fullTime > timeLimit) break; 
    }
    meanOrd /= totalWeightOrd;

    //printSample
    printSample(gm,meanOrd,oriScore,CASpSource,CASpEnd,CASpLh,CASpSS,CASpLevel,this->fullTime,index,prefix,opFileName,batch,sampleRatio,nstar);
}