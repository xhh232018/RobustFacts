#include <iostream>
#include <dirent.h>
#include "QueryInterface.h"
#include "QueryPerturbation.h"
#include "DataPerturbation.h"

using namespace std;

void QueryInterface::perturbPrep(GraphManager &gm, const char *filename){
    ifstream ifs(filename);
    if (!ifs.is_open()) {
        cout << "Invalid filepath\n";
        return;
    }

    string line;
    while (getline(ifs, line)) {

        istringstream iss(line);
        string token,tokenInner,psig;
        getline(iss, token, ','); // target
        queryNidVec.push_back(stoul(token));

        getline(iss, token, ','); // patternSig
        patternVec.push_back(token);

        std::vector<std::string> tmp;
        std::istringstream iss2(token);
        std::string tk;
        while(iss2 >> tk){
            tmp.push_back(tk);
        }
        cxtNodeVec.push_back(stoul(tmp.back()));

        getline(iss, token, ','); // attrEid
        attrVec.push_back(stoul(token));

        getline(iss, token, ','); // attrNid
        attrValVec.push_back(gm.IRI2nid[token]);

        getline(iss, token, ','); //oriscore
        oriScoreVec.push_back(stod(token));
    }
}

void QueryInterface::processBaseQuery(GraphManager &gm, std::string nidQueryFile, std::string opFilename, std::string prefix){
    perturbPrep(gm, nidQueryFile.c_str());
    if (queryNidVec.empty()) {
        cout << "No query loaded!" << endl;
        return;
        
    }
    unordered_map<uint32_t,double> degmap;
    double avgDeg = 0;
    for (int i = 0; i < queryNidVec.size(); ++i) {

        uint32_t entityInterest = queryNidVec[i];
        if(!gm.nid2types.count(entityInterest)) continue;

        string pattern = patternVec[i];
        uint32_t attr = attrVec[i];
        uint32_t attrVal = attrValVec[i];

        cout << pattern <<","<< attr << "," << attrVal << endl;

        QueryPerturbation qp;

        qp.bySampleSizeBase(gm, pattern, this->oriScoreVec[i], this->timeLimit, attr, attrVal, this->cxtNodeVec[i], i, opFilename, prefix);

    }
}

void QueryInterface::processSPQuery(GraphManager &gm, std::string nidQueryFile,std::string opFilename, std::string prefix,double sampleRatio){

    perturbPrep(gm, nidQueryFile.c_str());
    if (queryNidVec.empty()) {
        cout << "No query loaded!" << endl;
        return;
    }
    for (int i = 0; i < queryNidVec.size(); ++i) {
        
        uint32_t entityInterest = queryNidVec[i];

        string pattern = patternVec[i];
        uint32_t attr = attrVec[i];
        uint32_t attrVal = attrValVec[i];

        cout << pattern <<","<< attr << "," << attrVal << endl;

        QueryPerturbation qp;

        qp.bySampleSizeExp(gm, pattern, this->oriScoreVec[i], this->timeLimit,sampleRatio,attr, attrVal, this->cxtNodeVec[i], i, opFilename, prefix);
    }
}

void QueryInterface::processBaseData(GraphManager &gm, std::string nidQueryFile, std::string opFilename, std::string prefix){
    perturbPrep(gm, nidQueryFile.c_str());
    if (queryNidVec.empty()) {
        cout << "No query loaded!" << endl;
        return;
        
    }
    unordered_map<uint32_t,double> degmap;
    double avgDeg = 0;
    for (int i = 0; i < queryNidVec.size(); ++i) {

        uint32_t entityInterest = queryNidVec[i];
        if(!gm.nid2types.count(entityInterest)) continue;

        string pattern = patternVec[i];
        uint32_t attr = attrVec[i];
        uint32_t attrVal = attrValVec[i];

        cout << pattern <<","<< attr << "," << attrVal << endl;

        DataPerturbation dp;

        dp.exact(gm, pattern, this->oriScoreVec[i], this->timeLimit, attr, attrVal, this->cxtNodeVec[i], i, opFilename, prefix);

    }
}

void QueryInterface::processSPData(GraphManager &gm, std::string nidQueryFile,std::string opFilename, std::string prefix,double sampleRatio){
    perturbPrep(gm, nidQueryFile.c_str());
    if (queryNidVec.empty()) {
        cout << "No query loaded!" << endl;
        return;
    }
    for (int i = 0; i < queryNidVec.size(); ++i) {
        
        uint32_t entityInterest = queryNidVec[i];

        string pattern = patternVec[i];
        uint32_t attr = attrVec[i];
        uint32_t attrVal = attrValVec[i];

        cout << pattern <<","<< attr << "," << attrVal << endl;

        DataPerturbation dp;

        dp.sample(gm, pattern, this->oriScoreVec[i], this->timeLimit,sampleRatio,attr, attrVal, this->cxtNodeVec[i], i, opFilename, prefix);
    }
}