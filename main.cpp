#include <iostream>
#include <fstream>
#include <sstream>
#include "GraphManager.h"
#include "QueryInterface.h"
#include "DataPerturbation.h"
#include <unistd.h>

using namespace std;

int main(int argc,char *argv[]) {

    int opt;
    const char *optstring = "e:n:f:t:s:r:p:o:";
    //parameter setting
    string sampleEFile,sampleNFile;
    string factFile;
    string perturbationType;
    string sampleFlag;
    string factPrefix;
    string opFilename,srF;
    double sampleRatio;
    while ((opt=getopt(argc, argv, optstring)) != -1){
        switch(opt){
            case 'e':
                sampleEFile = optarg;
                break;
            case 'n':
                sampleNFile = optarg;
                break;
            case 'f':
                factFile = optarg;
                break;
            case 't':
                perturbationType = optarg;
                break;
            case 's':
                sampleFlag = optarg;
                break;
            case 'r':
                sampleRatio = stod(optarg);
                break;
            case 'p':
                factPrefix = optarg;
                break;
            case 'o':
                opFilename = optarg;
                break;
        }
    }
    
    //Graph Construction
    GraphManager gm;
    gm.readEdges(sampleEFile.c_str());
    gm.readNodeTypes(sampleNFile.c_str());
    gm.setEdgeType2Pos();

    //gm.deserializeFromDisk(serializeFile.c_str());
    printf("Node Num = %ld, edge num = %ld.\n", gm.node_sz, gm.edge_sz);
    cout<<"Graph reading finished!"<<endl;
    
    if(perturbationType == "data"){
        if(sampleFlag == "T" && sampleRatio > 0 && sampleRatio <= 1){
            cout << "data perturbation " << endl;
            QueryInterface qi;
            qi.processSPData(gm,factFile,opFilename,factPrefix,sampleRatio);
        }
        else if (sampleFlag == "F")
        {
            cout << "data perturbation exact" << endl;
            QueryInterface qi;
            qi.processBaseData(gm,factFile,opFilename,factPrefix);
        }
        else{
            cout << "Syntax Error - Incorrect Parameter Usage" << endl;
        }
    }
    else if(perturbationType == "query"){
        if(sampleFlag == "T" && sampleRatio > 0 && sampleRatio <= 1){
            cout << "query perturbation sampling" << endl;
            QueryInterface qi;
            qi.processSPQuery(gm,factFile,opFilename,factPrefix,sampleRatio);
        }
        else if (sampleFlag == "F")
        {
            cout << "query perturbation exact" << endl;
            QueryInterface qi;
            qi.processBaseQuery(gm,factFile,opFilename,factPrefix);
        }
        else{
            cout << "Syntax Error - Incorrect Parameter Usage" << endl;
        }
    }
    else{
        cout << "Syntax Error - Incorrect Parameter Usage" << endl;
    }

    cout<<"test finished"<<endl;
}
