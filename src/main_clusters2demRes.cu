#include "mergeClusteringGB.h"

int main(int argc,char *argv[])
{   
    if(argc != 7) 
    {
        cout << "usage: ./mainCluster2DemRes [good cluster file] [signal floder path] [signal root name] [barcode signals path] [output file] [barcode number]" << endl;
        exit(-1);
    }
    /* define variables */
    const string goodClusterFile = argv[1]; // good cluster file.
    const string sigsPath = argv[2]; // signal path.
    const string sigRootName = argv[3]; // the name of signal file.
    const string trueBarcodeSigsPath = argv[4]; 
    const string outFile = argv[5]; // the output file.
    const int barcodeNum = atoi(argv[6]);
    
    /* define variables */

    /* get the index of clusters and define demuplexing result. */
    vector<vector<int> > goodClusterList;
    readGoodClusterFile(goodClusterList, goodClusterFile);
    vector<int> demulRes;
    /* get the index of clusters and define demuplexing result. */
    
    /* get the sampled signals from each cluster. */
    vector<int> allSampNumIndex;  // all sampling numbers.
    vector<int> allSampSigIndex;
    vector<vector<float>> sampSigs;
    for(int i = 0; i < goodClusterList.size(); i++)
    {
        int samNum = goodClusterList[i].size() > 5 ? 5 : goodClusterList[i].size();
        allSampNumIndex.push_back(samNum);
        for(int j = 0; j < samNum; j++) allSampSigIndex.push_back(goodClusterList[i][j]);
    }
    getSigsOfList(sampSigs, allSampSigIndex, sigsPath, sigRootName);
    /* get the sampled signals from each cluster. */

    /* get barcode signals. */
    vector<int> barcodeIdxList;
    vector<vector<float>> barcodeSigs;
    for(int i = 0; i < barcodeNum; i++) barcodeIdxList.push_back(i);
    getSigsOfList(barcodeSigs, barcodeIdxList, trueBarcodeSigsPath, sigRootName);
    /* get barcode signals. */

    /* get DTW matrix between barcode signals and nanopore signals. */
    vector<vector<float>> distMatrix;
    gpuCalculatemnDynamicTimeWarping_1024(barcodeSigs, sampSigs, distMatrix);
    distMatrix = transpose(distMatrix);
    /* get DTW matrix between barcode signals and nanopore signals. */

    int t = 0;
    for(int i = 0; i < allSampNumIndex.size(); i++)
    {
        vector<vector<float>> sampDistMatrix(distMatrix.begin() + t, distMatrix.begin() + t + allSampNumIndex[i]);
        vector<int> tempDemIndexList;
        for(int j = 0; j < allSampNumIndex[i]; j++)
        {
            int minPosition = min_element(sampDistMatrix[j].begin(),sampDistMatrix[j].end()) - sampDistMatrix[j].begin();
            tempDemIndexList.push_back(minPosition);
        }
        int thisDemRes = majorElemCandidate(tempDemIndexList);
        // cout << thisDemRes << endl;
        demulRes.push_back(thisDemRes);
        t += allSampNumIndex[i];
    }

    int readNum = 0;
    for(int i = 0; i < goodClusterList.size(); i++) readNum += goodClusterList[i].size();
    vector<int> demulResLabelList(readNum, 0);
    for(int i = 0; i < goodClusterList.size(); i++)
    {
        for(int j = 0; j < goodClusterList[i].size(); j++)
        {
            demulResLabelList[goodClusterList[i][j]] = demulRes[i];
        }
    }

    ofstream disfile(argv[5], ios::out);
    string temp;
    for (int i = 0; i < demulResLabelList.size(); i++)
    {
        temp = std::to_string(demulResLabelList[i]);
        disfile << temp;
        disfile << endl;
    }

    return 0;
}