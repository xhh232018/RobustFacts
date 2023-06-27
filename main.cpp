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
    const char *optstring = "d:f:t:s:r:p:o:l:";
    //parameter setting
    string datfile;
    string factFile;
    string perturbationType;
    string sampleFlag;
    string factPrefix;
    string opFilename,srF;
    double sampleRatio,tl;

    while ((opt=getopt(argc, argv, optstring)) != -1){
        switch(opt){
            case 'd':
                datfile = optarg;
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
            case 'l':
                tl = stod(optarg);
                break;
        }
    }
    tl *= 100;
    //Graph Construction
    GraphManager gm;
    gm.deserializeFromDisk(datfile.c_str());

    printf("Node Num = %ld, edge num = %ld.\n", gm.node_sz, gm.edge_sz);
    cout<<"Graph reading finished!"<<endl;
    
    if(perturbationType == "data"){
        if(sampleFlag == "T" && sampleRatio > 0 && sampleRatio <= 1){
            cout << "data perturbation " << endl;
            QueryInterface qi;
            qi.timeLimit = tl;
            qi.processSPData(gm,factFile,opFilename,factPrefix,sampleRatio);
        }
        else if (sampleFlag == "F")
        {
            cout << "data perturbation exact" << endl;
            QueryInterface qi;
            qi.timeLimit = tl;
            qi.processBaseData(gm,factFile,opFilename,factPrefix);
        }
        else{
            cout << "Syntax Error - Incorrect Parameter Usage" << endl;
        }
    }
    else if(perturbationType == "entity"){
        if(sampleFlag == "T" && sampleRatio > 0 && sampleRatio <= 1){
            cout << "entity perturbation sampling" << endl;
            QueryInterface qi;
            qi.timeLimit = tl;
            qi.processSPEntity(gm,factFile,opFilename,factPrefix,sampleRatio);
        }
        else if (sampleFlag == "F")
        {
            cout << "entity perturbation exact" << endl;
            QueryInterface qi;
            qi.timeLimit = tl;
            qi.processBaseEntity(gm,factFile,opFilename,factPrefix);
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
