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

    vector<int> allReadIdxList;
    for(int i = 0; i < goodClusterList.size(); i++)
    {
        for(int j = 0; j < goodClusterList[i].size(); j++)
        {
            allReadIdxList.push_back(goodClusterList[i][j]);
        }
    }

    sort(allReadIdxList.begin(), allReadIdxList.end());
    /* get the index of clusters and define demuplexing result. */

    /* get the threshold. */
    vector<int> IndexList;
    int count = 100 < goodClusterList[0].size() ? 100 : goodClusterList[0].size();
    for(int i = 0; i < count; i++)
    {
        IndexList.push_back(goodClusterList[0][i]);
    }

    vector<vector<float>> selSigs;
    getSigsOfList(selSigs, IndexList, sigsPath, sigRootName);

    vector<vector<float>> distList;
    gpuCalculatemnDynamicTimeWarping_1024(selSigs, selSigs, distList);
    
    vector<float> oneDdistList;
    oneT2D(oneDdistList, distList);
    double sum = std::accumulate(std::begin(oneDdistList), std::end(oneDdistList), 0.0);
    float meanDist = sum / oneDdistList.size();
    float threshold = meanDist * 1.6;
    cout << "threshold " << threshold << endl;
    // float minDist = *min_element(oneDdistList.begin(), oneDdistList.end());
    // float maxDist = *max_element(oneDdistList.begin(), oneDdistList.end());
    // float threshold = (minDist + maxDist) * 2.5;
    /* get the threshold. */
    
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
    
    if (barcodeIdxList.size() < 9999)
    {
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
            if (allSampNumIndex[i] == 1)
            {
                int minPosition = min_element(sampDistMatrix[0].begin(),sampDistMatrix[0].end()) - sampDistMatrix[0].begin();
                if(sampDistMatrix[0][minPosition] < threshold) tempDemIndexList.push_back(minPosition);
                else tempDemIndexList.push_back(-1);
            }
            else
            {
                for(int j = 0; j < allSampNumIndex[i]; j++)
                {
                    int minPosition = min_element(sampDistMatrix[j].begin(),sampDistMatrix[j].end()) - sampDistMatrix[j].begin();
                    tempDemIndexList.push_back(minPosition);
                }
            }
            
            int thisDemRes = majorElemCandidate(tempDemIndexList);
            // cout << thisDemRes << endl;
            demulRes.push_back(thisDemRes);
            t += allSampNumIndex[i];
        }
    }

    else
    {   
        int t = 0;
        for(int i = 0; i < allSampNumIndex.size(); i++)
        {   
            vector<vector<float>> sampDistMatrix;
            vector<vector<float>> subSampSigs(sampSigs.begin() + t, sampSigs.begin() + t + allSampNumIndex[i]);
            gpuCalculatemnDynamicTimeWarping_1024(subSampSigs, barcodeSigs, sampDistMatrix);
            // vector<vector<float>> sampDistMatrix(distMatrix.begin() + t, distMatrix.begin() + t + allSampNumIndex[i]);
            vector<int> tempDemIndexList;
            for(int j = 0; j < allSampNumIndex[i]; j++)
            {
                int minPosition = min_element(sampDistMatrix[j].begin(),sampDistMatrix[j].end()) - sampDistMatrix[j].begin();
                if(sampDistMatrix[j][minPosition] < threshold) tempDemIndexList.push_back(minPosition);
                else tempDemIndexList.push_back(-1);
            }
            int thisDemRes = majorElemCandidate(tempDemIndexList);
            // cout << thisDemRes << endl;
            demulRes.push_back(thisDemRes);
            t += allSampNumIndex[i];
        }
    }

    map<int, int> demRes;
    for(int i = 0; i < goodClusterList.size(); i++)
    {
        for(int j = 0; j < goodClusterList[i].size(); j++)
        {
            demRes[goodClusterList[i][j]] = demulRes[i];
        }
    }

    ofstream disfile(argv[5], ios::out);
    string barID;
    string readID;
    disfile << "read ID" << ": " << "barcode ID" << endl;
    map<int,int>::iterator it;
    for (it = demRes.begin(); it != demRes.end(); it++)
    {   
        barID = std::to_string(it->second);
        readID = std::to_string(it->first);
        disfile << readID << ": " << barID;
        disfile << endl;
    }

    return 0;
}