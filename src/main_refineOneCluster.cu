#include "mergeClusteringGB.h"

int main(int argc,char *argv[])
{   
    if (argc != 6) 
    {
        cout << "usage: ./refineOneCluster [one cluster file] [signal floder path] [signal file root name] [output file] [threshold factor]" << endl;
        exit(-1);
    }
    /*define veribels*/
    vector<vector<int> > oneclusters;
    vector<int> allEle; 
    vector<vector<float>> sigs;
    const string oneClusterFile = argv[1];
    const string sigDirPath = argv[2];
    const string sigRootName = argv[3];
    float threshold = 85;
    const int maxLocalLength = 2048;
    const string outFile = argv[4];
    float devideIndex = atof(argv[5]);
    vector<vector<int>> gpuclusterResult;
    readGoodClusterFile(oneclusters, oneClusterFile);
    fromClusterGetEle(allEle, oneclusters);
    sort(allEle.begin(), allEle.end());
    getSigsOfList(sigs, allEle, sigDirPath, sigRootName);
    /*define veribels*/
    
    refineOneEleCluster(sigs, threshold, devideIndex, gpuclusterResult, maxLocalLength);
    // cout << gpuclusterResult.size() << endl;

    // cout << threshold << endl;
    ofstream disfile(outFile, ios::out);
    string temp;
    for (int i = 0; i < gpuclusterResult.size(); i++)
    {
      for (int j = 0; j < gpuclusterResult[i].size(); j++)
      {
        temp = std::to_string(allEle[gpuclusterResult[i][j]]);
        disfile << temp;
        disfile << " ";
      }
      disfile << endl;
    }

    return 0;
    
}