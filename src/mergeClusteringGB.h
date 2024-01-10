#ifndef MERGECLUSTERINGGB
#define MERGECLUSTERINGGB

#include "cuda_def.cuh"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
//#include "cuda_kernels.h"
#include <omp.h>
#include <vector>
#include <set>
#include <map>
#define DBL_MAX 1.7976931348623158e+308 // max value
#define PI 3.1415926535898
#define MAX_NUM 2147483647
#define MAXLENGTH 1023
#define CPUTHREADS 36
// int CLUSTERNUM = 10;
// float THRESHOLD = 90;
#define ENDL printf("\n");
using namespace std;

// vector<vector<float>> GetAlignSignalList_1(vector<int> &SignalIdList, const string &argv_1);
void findMeredCluster(const vector<vector<float>> &distMatrix, const float &threshold, const int &sampNum, 
                      vector<int> &mergedIndexList);

int majorElemCandidate(const vector<int> A);

void findMinSize(vector<vector<int> > vec, int &minSize);

void randomSample(vector<int> &vec, const int &sampleNum, vector<int> &resList);

__global__ void cuDTW_1024(float *g_allColData, unsigned int *g_allColLength, float *g_allRowData,
                           unsigned int *g_allRowLength, float *g_odata);
__global__ void cuDTW_2048(float *g_allColData, unsigned int *g_allColLength, float *g_allRowData,
                           unsigned int *g_allRowLength, float *g_odata);
void TestDataGeneration(vector<vector<float>> &siga, vector<vector<float>> &sigb, int siga_n = 50,
                        int sigb_n = 50, bool aeb = false);
float gpuCalculatemnDynamicTimeWarping_1024(const vector<vector<float>> &siga,
                                            const vector<vector<float>> &sigb,
                                            vector<vector<float>> &gpudistResult);
float gpuCalculatemnDynamicTimeWarping_2048(const vector<vector<float>> &siga,
                                            const vector<vector<float>> &sigb,
                                            vector<vector<float>> &gpudistResult);
float cpuCalculatemnDynamicTimeWarping(const vector<vector<float>> &siga,
                                       const vector<vector<float>> &sigb,
                                       vector<vector<float>> &cpudistResult);
float cpuDynamicTimeWarping(const std::vector<float> &seq1, const std::vector<float> &seq2);
// void gpuCluster(const vector<vector<float>> &sigb, float threshold,
//                 vector<vector<int>> &gpuclusterResult, int maxLocalLength);

// vector<int> GetAlignList();
vector<vector<float>> transpose(vector<vector<float>>& A);

vector<vector<double>> GetAlignSignalList(vector<int> &SignalIdList, const string &argv_1);

// vector<string> GetargvList();
void readClusterFile(vector<int> &goodClusterList, vector<int> &badClusterList, const string &clusterFile);

void readGoodClusterFile(vector<vector<int> > &goodClusterList, const string &goodClusterFile);

vector<float> GetSingleSignalData(int &SignalId, const string &argv_1, const string &sigRootName);

void getSigsOfList(vector<vector<float> > &res, vector<int> &idxList, const string &sigDirPath, const string &sigRootName);

void DeleteNoteOff(vector<vector<int> > &eventStore);

// vector<vector<float>> GetCenterSignalList(const int &SignalScale);
void ZScoreNormalize(std::vector<float> &signals, float *avg = NULL, float *stdev = NULL);
template <class T>
void convertFromString(T &value, string &s);

void oneT2D(vector<float> &oneDList, vector<vector<float>> &twoDList);

void fromClusterGetEle(vector<int> &oneDList, vector<vector<int>> &twoDList);

void refineOneEleCluster(const vector<vector<float>> &sigb, float &threshold, float &devideIndex,
                        vector<vector<int>> &gpuclusterResult, int maxLocalLength);

#endif
