#include "mergeClusteringGB.h"

int main(int argc,char *argv[])
{   
    if(argc != 7) 
    {
        cout << "usage: ./mergeGoodCluster [good cluster file] [signal floder path] [signal root name] [sampling number] [output file] [threshold factor]" << endl;
        cout << "[sampling number] determines the merging conditions, the larger A is, the more accurate the merging is, and it is recommended to set it to 10." << endl;
        exit(-1);
    }
    cout << "Clusters are being merged..." << endl;
    /* define variables */
    const string goodClusterFile = argv[1]; // good cluster file.
    const string sigsPath = argv[2]; // signal path.
    const string sigRootName = argv[3]; // the name of signal file.
    const int sampLine = atoi(argv[4]); // the baselin of sampling number.
    const string outFile = argv[5]; // the output file.
    const float thresFactor = atof(argv[6]); // threshold factor.
    /* define variables */

    /* get the index of good cluster. */
    vector<vector<int> > goodClusterList;
    readGoodClusterFile(goodClusterList, goodClusterFile);

    

    int minSize;
    findMinSize(goodClusterList, minSize);
    const int randnum = sampLine < minSize ? sampLine : minSize; // sample number.
    // cout << randnum << endl;
    
    /* get the index of good cluster. */

    
    /* get the threshold. */
    vector<int> randomIndexList;
    int count = 200 < goodClusterList.size() ? 200 : goodClusterList.size();
    for(int i = 0; i < count; i++)
    {
        randomIndexList.push_back(goodClusterList[i][0]);
    }
    // cout << "add is end!" << endl;
    vector<vector<float>> randomSigs;
    getSigsOfList(randomSigs, randomIndexList, sigsPath, sigRootName);
    
    vector<vector<float>> distList;
    gpuCalculatemnDynamicTimeWarping_1024(randomSigs, randomSigs, distList);
    
    vector<float> oneDdistList;
    oneT2D(oneDdistList, distList);
    float minDist = *min_element(oneDdistList.begin(), oneDdistList.end());
    float maxDist = *max_element(oneDdistList.begin(), oneDdistList.end());
    float threshold = (minDist + maxDist) / thresFactor; // 2.85 is good.
    // cout << threshold << endl;
    /* get the threshold. */

    cout << "threshold is end!" << endl;

    /*tqdm*/
    // int i = 0;
	// char bar[102];
	// const char *lable = "|/-\\";
	// bar[0] = 0;
    /*tqdm*/

    /* try merge good cluster. */
    vector<vector<int> > MergedClusterList;
    while(goodClusterList.size() != 1)
    {   
        /*tqdm*/
        // printf("[%-100s][%d%%][%c]\r", bar, i, lable[i%4]);
		// fflush(stdout);
	    // bar[i] = '#';
		// i+= 100/goodClusterList.size();
		// bar[i] = 0;
        /*tqdm*/
        cout << "NO.of clusters: " << goodClusterList.size() << ". When it is equal to 1, the algorithm terminates." << endl;
        vector<vector<float>> startSigs; // the signals of top element.
        vector<int > topSigsRandIndex; // the randomly selected index of top element.
        randomSample(goodClusterList[0], randnum, topSigsRandIndex);
        // cout << goodClusterList[0][0] << endl;
        getSigsOfList(startSigs, topSigsRandIndex, sigsPath, sigRootName);

        vector<int> tempCluster(goodClusterList[0].begin(), goodClusterList[0].end()); // top cluster.
        // cout << tempCluster.size() << endl;
        auto iter = goodClusterList.erase(goodClusterList.begin()); // update goodClusterList.
        
        vector<vector<float>> sampledSigs;
        vector<int> samIndexList; // load the all index of random selected signals.
        for(int i = 0; i < goodClusterList.size(); i++) randomSample(goodClusterList[i], randnum, samIndexList);
        getSigsOfList(sampledSigs, samIndexList, sigsPath, sigRootName);
        // cout << goodClusterList.size() << endl;
        
        vector<vector<float>> distMatrix; // a distance matrix, size=3*(goodClusterList.size*3)
        gpuCalculatemnDynamicTimeWarping_1024(startSigs, sampledSigs, distMatrix);
        // cout << distMatrix.size() << distMatrix[0].size() << endl;
        
        vector<int> mergedIndexList;
        findMeredCluster(distMatrix, threshold, randnum, mergedIndexList);
        // cout << mergedIndexList.size() << endl;

        for(int i = 0; i < mergedIndexList.size(); i++)
        {   
            tempCluster.insert(tempCluster.end(), goodClusterList[mergedIndexList[i]].begin(), goodClusterList[mergedIndexList[i]].end());
            goodClusterList[mergedIndexList[i]].push_back(-1);
        }

        // for(int i = 0; i < mergedIndexList.size(); i++) cout << tempCluster[0] << " " << goodClusterList[mergedIndexList[i]][0] << endl;
        // cout << tempCluster.size() << endl;
        DeleteNoteOff(goodClusterList);
        // cout << "good size: " << goodClusterList.size() << endl;
        // cout << mergedIndexList.size() << "  " << goodClusterList[mergedIndexList[0]][0] << endl;
        MergedClusterList.push_back(tempCluster);
        if(goodClusterList.size() == 0) break;
        // break;
    }
    if(goodClusterList.size() == 1) MergedClusterList.push_back(goodClusterList[0]);
   
    /* try merge good cluster. */

    ofstream disfile(argv[5], ios::out);
    string temp;
    for (int i = 0; i < MergedClusterList.size(); i++)
    {
      for (int j = 0; j < MergedClusterList[i].size(); j++)
      {
        temp = std::to_string(MergedClusterList[i][j]);
        disfile << temp;
        disfile << " ";
      }
      disfile << endl;
    }

    return 0;
}