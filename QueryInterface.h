#ifndef PERTURBATION_QUERYINTERFACE_H
#define PERTURBATION_QUERYINTERFACE_H

#include "GraphManager.h"
#include "QueryPerturbation.h"
#include <fstream>

class QueryInterface{
public:

    std::vector<std::string> queryIRIVec,patternVec;
    std::vector<uint32_t> queryNidVec,attrVec,attrValVec,cxtNodeVec;
    std::vector<double> oriScoreVec;
    double timeLimit = 300000;
                           
    void perturbPrep(GraphManager &gm, const char *filename);

    void processBaseData(GraphManager &gm, std::string nidQueryFile, std::string opFilename, std::string prefix);

    void processSPData(GraphManager &gm, std::string nidQueryFile, std::string opFilename, std::string prefix, double sampleRatio);

    void processBaseQuery(GraphManager &gm, std::string nidQueryFile, std::string opFilename, std::string prefix);

    void processSPQuery(GraphManager &gm, std::string nidQueryFile, std::string opFilename, std::string prefix, double sampleRatio);

};

#endif