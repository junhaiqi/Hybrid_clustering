#include "mergeClusteringGB.h"

int main(int argc,char *argv[])
{   
    if(argc != 5) 
        {
        //  vector<int> goodClusterList; vector<int> badClusterList;
        //  const string clusterFile = "../tempoutput/newCleanedBarcode350ForGB_goupedIndex.txt";
        //  readClusterFile(goodClusterList, badClusterList, clusterFile);
        // //  cout << badClusterList.size() << goodClusterList.size() << endl;
        //     cout << badClusterList[0] << " " << badClusterList[badClusterList.size()-1] << endl;
        //     cout << goodClusterList[0] << " " << goodClusterList[goodClusterList.size()-1] << endl;
         cout << "usage ./mainMerge [cluster file path] [barcode signal floder path] [signal root name] [threshold factor]" << endl; 
         cout << "The first line of cluster file is the index of sequence in good cluster, and the second line is the index of sequence in bad cluster." << endl;
         exit(-1);
        }
    else
    {
        // initialization index list of good cluster and list of bad cluster.
        vector<int> goodClusterList; vector<int> badClusterList; 

        // get the index of good clusters and bad clusters.
        const string clusterFile = argv[1];
        readClusterFile(goodClusterList, badClusterList, clusterFile);
        
        // initialization signal list of good cluster and list of bad cluster.
        vector<vector<float> > goodClusterSigs; vector<vector<float>> badClusterSigs;

        const string sigDir = argv[2]; // path of signal folder.
        const string sigRootName = argv[3]; // root name of signal file.
        const float thresFactor = atof(argv[4]); // threshold factor.

        // get good cluster signals.
        cout << "Loading data..." << endl;
        for(int i = 0; i < goodClusterList.size(); i++)
        {   
            // cout << "good: " << i << " " << goodClusterList[i] << endl;
            vector<float> sig = GetSingleSignalData(goodClusterList[i], sigDir, sigRootName);
            goodClusterSigs.push_back(sig);
        }

        // get bad cluster signals.
        for(int i = 0; i < badClusterList.size(); i++)
        {   
            // cout << "bad: " << i << " " << badClusterList[i] << endl;
            vector<float> sig = GetSingleSignalData(badClusterList[i], sigDir, sigRootName);
            badClusterSigs.push_back(sig);
        }
        cout << "Data loading is complete!" << endl;
        

        // get DTW distance matrix.
        vector<vector<float>> gpudistResult;
        gpuCalculatemnDynamicTimeWarping_1024(goodClusterSigs, badClusterSigs, gpudistResult);

        // get threshold.
        gpudistResult = transpose(gpudistResult);
        // cout << "Data Size: " << gpudistResult.size() << "*" << gpudistResult[0].size() << endl;

        float minDist = *min_element(gpudistResult[0].begin(),gpudistResult[0].end());
        float maxDist = *max_element(gpudistResult[0].begin(),gpudistResult[0].end());
        float threshold = (minDist + maxDist) / thresFactor;
        const int failedMergeIndex = goodClusterList.size() + 1;
        // cout << "threshold: " << threshold << endl;

        // from DTW distance matrix to merge bad clusters into good clusters.
        vector<int> mergeIndexList;
        for(int i = 0; i<gpudistResult.size(); i++)
        {   
            float minValue = *min_element(gpudistResult[i].begin(),gpudistResult[i].end());
            if(minValue < threshold+1) // merge condition.
            {   
                int minPosition = min_element(gpudistResult[i].begin(),gpudistResult[i].end()) - gpudistResult[i].begin();
                mergeIndexList.push_back(minPosition);
            }
            else mergeIndexList.push_back(failedMergeIndex);
        }

        for(int i = 0; i < mergeIndexList.size(); i++)
        {
            cout << mergeIndexList[i] << endl;
        }
        return 0;
    }
}