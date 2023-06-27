#include <iostream>
#include <dirent.h>
#include "QueryInterface.h"
#include "EntityPerturbation.h"
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

void QueryInterface::processBaseEntity(GraphManager &gm, std::string nidQueryFile, std::string opFilename, std::string prefix){
    perturbPrep(gm, nidQueryFile.c_str());

    unordered_map<uint32_t,double> degmap;
    double avgDeg = 0;
    for (int i = 0; i < patternVec.size(); ++i) {

        string pattern = patternVec[i];
        uint32_t attr = attrVec[i];
        uint32_t attrVal = attrValVec[i];

        cout << pattern <<","<< attr << "," << attrVal << endl;

        EntityPerturbation ep;

        ep.bySampleSizeBase(gm, pattern, this->oriScoreVec[i], this->timeLimit, attr, attrVal, this->cxtNodeVec[i], i, opFilename, prefix);

    }
}

void QueryInterface::processSPEntity(GraphManager &gm, std::string nidQueryFile,std::string opFilename, std::string prefix,double sampleRatio){

    perturbPrep(gm, nidQueryFile.c_str());

    for (int i = 0; i < patternVec.size(); ++i) {

        string pattern = patternVec[i];
        uint32_t attr = attrVec[i];
        uint32_t attrVal = attrValVec[i];

        cout << pattern <<","<< attr << "," << attrVal << endl;

        EntityPerturbation ep;

        ep.bySampleSizeExp(gm, pattern, this->oriScoreVec[i], this->timeLimit,sampleRatio,attr, attrVal, this->cxtNodeVec[i], i, opFilename, prefix);
    }
}

void QueryInterface::processBaseData(GraphManager &gm, std::string nidQueryFile, std::string opFilename, std::string prefix){
    perturbPrep(gm, nidQueryFile.c_str());

    unordered_map<uint32_t,double> degmap;
    double avgDeg = 0;
    for (int i = 0; i < patternVec.size(); ++i) {

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

    for (int i = 0; i < patternVec.size(); ++i) {
        
        string pattern = patternVec[i];
        uint32_t attr = attrVec[i];
        uint32_t attrVal = attrValVec[i];

        cout << pattern <<","<< attr << "," << attrVal << endl;

        DataPerturbation dp;

        dp.sample(gm, pattern, this->oriScoreVec[i], this->timeLimit,sampleRatio,attr, attrVal, this->cxtNodeVec[i], i, opFilename, prefix);
    }
}