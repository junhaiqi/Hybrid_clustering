//#include "cuda_def.cuh"
//#include "cuda_proc.h"
#include "cuda_def.cuh"
#include <iostream>
#include <vector>
//#include "cuda_kernels.h"
#include <omp.h>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <algorithm>
using namespace std;
#define DBL_MAX 1.7976931348623158e+308 // max value
#define PI 3.1415926535898
#define MAX_NUM 2147483647
#define MAXLENGTH 1023
#define CPUTHREADS 5

using namespace std;
//////////////////////////////////
//函数声明
//////////////////////////////////

__global__ void cuDTW_2048(float *g_allColData, unsigned int *g_allColLength, float *g_allRowData,
                           unsigned int *g_allRowLength, float *g_odata);

float gpuCalculatemnDynamicTimeWarping_2048(const vector<vector<float>> &siga,
                                            const vector<vector<float>> &sigb,
                                            vector<vector<float>> &gpudistResult);
__global__ void cuDTW_ultimate(float *g_allColData,
                               unsigned int *g_allColLength,
                               float *g_allRowData,
                               unsigned int *g_allRowLength, float *g_odata);
void TestDataGeneration(vector<vector<float>> &siga,
                        vector<vector<float>> &sigb, int siga_n = 1000,
                        int sigb_n = 50);
float gpuCalculate1nDynamicTimeWarping(const vector<vector<float>> &siga,
                                       const vector<vector<float>> &sigb,
                                       vector<vector<float>> &gpudistResult);
float cpuCalculate1nDynamicTimeWarping(const vector<vector<float>> &siga,
                                       const vector<vector<float>> &sigb,
                                       vector<vector<float>> &cpudistResult);
float cpuDynamicTimeWarping(const std::vector<float> &seq1,
                            const std::vector<float> &seq2);

// using namespace std;

vector<double> GetSingleSignalData(int &SignalId, const string &argv_1);

template <class T>
void convertFromString(T &value, string &s);

vector<double> GetSingleSignalData1(int &SignalId, const string &argv_2);

vector<int> GetAlignList();

vector<string> GetargvList();

vector<vector<double>> GetCenterSignalList(const vector<int> &AlignList, const string &argv_2);

vector<vector<double>> GetAlignSignalList(vector<int> &SignalIdList, const string &argv_1);

vector<vector<float>> FromDoubleToFloat(vector<vector<double>> a);

void GetminAndpos(vector<vector<float>> &dislist, vector<vector<float>> MinAndPosresultList);

void ZScoreNormalize(vector<double> &signals);

int main()
{
  // int id = 10;
  // int scale = 6126;
  // vector<double > test = GetSingleSignalData1(id,scale);
  // for(int i = 0; i < test.size(); i++){
  //   cout << test[i] << endl;
  // }
  vector<vector<float>> gpudtwresult;
  vector<vector<float>> MinAndPosresultList;
  vector<int> Alignlist = GetAlignList();
  // cout << Alignlist.size() << endl;
  // for(int i = 0; i < Alignlist.size();i++){
  //   cout << Alignlist[i] << endl;
  // }
  // cout << Alignlist[2] << endl;
  int scale = Alignlist[1];
  vector<string> argvList = GetargvList();
  // for(int i = 0; i < test.size(); i++){
  //   cout << test[i] << endl;
  // }
  vector<vector<double>> AlignSignallist = GetAlignSignalList(Alignlist, argvList[0]);
  vector<vector<double>> CenterSignallist = GetCenterSignalList(Alignlist, argvList[1]);
  vector<vector<float>> AlignSignallist1 = FromDoubleToFloat(AlignSignallist);
  vector<vector<float>> CenterSignallist1 = FromDoubleToFloat(CenterSignallist);

  gpuCalculatemnDynamicTimeWarping_2048(AlignSignallist1, CenterSignallist1, gpudtwresult);
  // GetminAndpos(gpudtwresult,MinAndPosresultList);
  // for(int i = 0; i < gpudtwresult.size(); i++){
  //   for(int j = 0; j < gpudtwresult[i].size(); j++){
  //     cout << gpudtwresult[i][j] << endl;
  //   }
  // }
  ofstream disfile("OnetoNdisfile.txt", ios::out);
  string temp;
  for (int i = 0; i < gpudtwresult.size(); i++)
  {
    for (int j = 0; j < gpudtwresult[i].size(); j++)
    {
      temp = std::to_string(gpudtwresult[i][j]);
      disfile << temp;
      disfile << " ";
    }
    disfile << endl;
  }
  // cout << gpudtwresult.size() << endl;
  disfile.close();
  // cout << MinAndPosresultList[0][0] << " " << MinAndPosresultList[1][0] << endl;
  cout << "OnetoNdisfile created successfully!" << endl;
  return 0;
}

//////////////////////////////////
//function
//////////////////////////////////

void ZScoreNormalize(vector<double> &signals)
{
  double sum = accumulate(signals.begin(), signals.end(), 0.0);
  double mean = sum / signals.size();

  double acc = 0.0;
  for (size_t i = signals.size(); i--;)
  {
    signals[i] = signals[i] - mean;
    acc += signals[i] * signals[i];
  }

  double deviation = sqrt(acc / signals.size());

  for (size_t i = signals.size(); i--;)
  {
    signals[i] /= deviation;
  }
}

vector<vector<float>> FromDoubleToFloat(vector<vector<double>> a)
{
  vector<vector<float>> temp;
  for (int i = 0; i < a.size(); i++)
  {
    vector<float> templist;
    for (int j = 0; j < a[i].size(); j++)
    {
      float tempvalue = (float)a[i][j];
      templist.push_back(tempvalue);
    }
    temp.push_back(templist);
  }
  return temp;
}


// void GetminAndpos(const vector<vector<float > >& dislist,vector<vector<float > > MinAndPosresultList){
//    vector<float> minvaluelist;
//    vector<float> positionlist;
//    omp_set_num_threads(CPUTHREADS);
//    #pragma omp parallel for
//    for(int i = 0; i < dislist.size(); i++){
//       float minvalue = dislist[i][0];
//       float position = 0;
//       for(int j = 0; j < dislist[i].size(); j++){
//          if(minvalue > dislist[i][j]){
//            minvalue = dislist[i][j];
//            position = (float)j;
//          }
//       minvaluelist.push_back(minvalue);
//       positionlist.push_back(position);
//       }
//     MinAndPosresultList.push_back(minvaluelist);
//     MinAndPosresultList.push_back(positionlist);
//    }
// }

void GetminAndpos(vector<vector<float>> &dislist, vector<vector<float>> MinAndPosresultList)
{
  vector<float> minvaluelist;
  vector<float> positionlist;
  omp_set_num_threads(CPUTHREADS);
#pragma omp parallel for
  for (int i = 0; i < dislist.size(); i++)
  {
    // float minvalue = 0;
    // float position = 0;
    vector<float>::iterator minvalue = min_element(dislist[i].begin(), dislist[i].end());
    cout << (int)*minvalue << endl;
    minvaluelist.push_back((int)*minvalue);
    positionlist.push_back(distance(dislist[i].begin(), minvalue));
  }
  MinAndPosresultList.push_back(minvaluelist);
  MinAndPosresultList.push_back(positionlist);
}

vector<double> GetSingleSignalData(int &SignalId, const string &argv_1)
{
  // cout << SignalId << endl;
  ifstream SignalFile;
  string TempId = std::to_string(SignalId);
  // string TempId2 = std::to_string(SignalScale);
  string TempString = argv_1 + "/" + "signal_" + TempId + ".txt";
  const char *FileName = TempString.data();
  SignalFile.open(FileName, ios::in);
  if (!SignalFile.is_open())
  {
    cout << "Signal file open error!" << endl;
    cout << argv_1 << endl;
    cout << SignalId << endl;
  }
  string FileLine;
  vector<double> SignalData;
  while (getline(SignalFile, FileLine))
  {
    double SignalValue;
    convertFromString(SignalValue, FileLine);
    SignalData.push_back(SignalValue);
  }
  ZScoreNormalize(SignalData);
  SignalFile.close();
  return SignalData;
}


vector<double> GetSingleSignalData1(int &SignalId, const string &agrv_2)
{
  ifstream SignalFile;
  string TempId = std::to_string(SignalId);
  // string TempId2 = std::to_string(SignalScale);
  string TempString = agrv_2 + "/" + "consensus_sig_" + TempId + ".txt";
  const char *FileName = TempString.data();
  SignalFile.open(FileName, ios::in);
  if (!SignalFile.is_open())
  {
    cout << "CenterSignal file open error!" << endl;
  }
  string FileLine;
  vector<double> SignalData;
  while (getline(SignalFile, FileLine))
  {
    double SignalValue;
    convertFromString(SignalValue, FileLine);
    SignalData.push_back(SignalValue);
  }
  ZScoreNormalize(SignalData);
  SignalFile.close();
  return SignalData;
}

vector<int> GetAlignList()
{
  ifstream Alignfile;
  // const char* filename = "ReadyToSort.txt".data();
  Alignfile.open("ReadyToSortfile.txt", ios::in);
  if (!Alignfile.is_open())
  {
    cout << "Align file open error!" << endl;
  }
  string FileLine;
  vector<int> AlignList;
  while (getline(Alignfile, FileLine))
  {
    int SignalId;
    convertFromString(SignalId, FileLine);
    AlignList.push_back(SignalId);
  }
  return AlignList;
}

vector<string> GetargvList()
{
  ifstream argvfile;
  // const char* filename = "ReadyToSort.txt".data();
  argvfile.open("argv_file.txt", ios::in);
  if (!argvfile.is_open())
  {
    cout << "argv_file.txt open error!" << endl;
  }
  string FileLine;
  vector<string> argvList;
  while (getline(argvfile, FileLine))
  {
    argvList.push_back(FileLine);
  }
  return argvList;
}

vector<vector<double>> GetCenterSignalList(const vector<int> &AlignList, const string &argv_2)
{
  vector<vector<double>> SignalList;
  for (int i = 0; i < AlignList[0]; i++)
  {
    vector<double> TempSignal = GetSingleSignalData1(i, argv_2);
    SignalList.push_back(TempSignal);
  }
  return SignalList;
}

vector<vector<double>> GetAlignSignalList(vector<int> &SignalIdList, const string &argv_1)
{
  vector<vector<double>> CenterSignalList;
  for (int i = 2; i < SignalIdList.size(); i++)
  {
    vector<double> TempSignal = GetSingleSignalData(SignalIdList[i], argv_1);
    CenterSignalList.push_back(TempSignal);
  }
  return CenterSignalList;
}

void TestDataGeneration(vector<vector<float>> &siga,
                        vector<vector<float>> &sigb, int siga_n, int sigb_n)
{
  siga.resize(siga_n);
  sigb.resize(sigb_n);

  omp_set_num_threads(CPUTHREADS);
#pragma omp parallel for
  for (int i = 0; i < siga_n; i++)
  {
    int siga_length = rand() % 200 + 700; // siga_i的长度在700-900之间
    siga[i].resize(siga_length);
    for (int j = 0; j < siga_length; j++)
    {
      siga[i][j] = rand() % 400 + 400; // siga_i_j的范围在400-800之间
    }
  }

#pragma omp parallel for
  for (int i = 0; i < sigb_n; i++)
  {
    int sigb_length = rand() % 200 + 700; // sigb_i的长度在700-900之间
    sigb[i].resize(sigb_length);
    for (int j = 0; j < sigb_length; j++)
    {
      sigb[i][j] = rand() % 400 + 400; // sigb_i_j的范围在400-800之间
    }
  }
}

float gpuCalculate1nDynamicTimeWarping(const vector<vector<float>> &siga,
                                       const vector<vector<float>> &sigb,
                                       vector<vector<float>> &gpudistResult)
{
  // return 0;
  int siga_n = siga.size();
  int sigb_n = sigb.size();
  gpudistResult.resize(siga_n);
  for (int i = 0; i < siga_n; i++)
  {
    gpudistResult[i].resize(sigb_n);
  }
  // int siga_length = 0;
  float *d_distResult = NULL;
  float *d_allColData = NULL;
  float *d_allRowData = NULL;
  unsigned int *d_allRowLength;
  unsigned int *d_allColLength;
  // vector<float *> rowDataList(sigb_n);
  cudaMalloc((void **)&d_distResult, siga_n * sigb_n * sizeof(float));
  CUERR
  cudaMalloc((void **)&d_allColData, siga_n * 1024 * sizeof(float));
  CUERR
  cudaMalloc((void **)&d_allRowData, sigb_n * 1024 * sizeof(float));
  CUERR

  vector<unsigned int> h_allRowLength(sigb_n);
  for (int i = 0; i < sigb_n; i++)
  {
    h_allRowLength[i] = sigb[i].size();
    cudaMemcpy(&d_allRowData[1024 * i], &sigb[i][0],
               h_allRowLength[i] * sizeof(float), cudaMemcpyHostToDevice);
    CUERR
  }
  cudaMalloc((void **)&d_allRowLength, sigb_n * sizeof(unsigned int));
  CUERR
  cudaMemcpy(d_allRowLength, &h_allRowLength[0], sigb_n * sizeof(unsigned int),
             cudaMemcpyHostToDevice);
  CUERR

  vector<unsigned int> h_allColLength(siga_n);
  for (int i = 0; i < siga_n; i++)
  {
    h_allColLength[i] = siga[i].size();
    cudaMemcpy(&d_allColData[1024 * i], &siga[i][0],
               h_allColLength[i] * sizeof(float), cudaMemcpyHostToDevice);
    CUERR
  }
  cudaMalloc((void **)&d_allColLength, siga_n * sizeof(unsigned int));
  CUERR
  cudaMemcpy(d_allColLength, &h_allColLength[0], siga_n * sizeof(unsigned int),
             cudaMemcpyHostToDevice);
  CUERR

  float timesum = 0;
  dim3 threadsPerBlock(1024);
  dim3 blocksPerGrid(sigb_n, siga_n);
  cuDTW_ultimate<<<blocksPerGrid, threadsPerBlock>>>(
      d_allColData, d_allColLength, d_allRowData, d_allRowLength, d_distResult);
  CUERR

  for (int i = 0; i < siga_n; i++)
  {
    cudaMemcpy(&gpudistResult[i][0], &d_distResult[sigb_n * i],
               sigb_n * sizeof(float), cudaMemcpyDeviceToHost);
    CUERR
  }

  cudaFree(d_allColData);
  CUERR
  cudaFree(d_distResult);
  CUERR
  cudaFree(d_allRowData);
  CUERR
  cudaFree(d_allRowLength);
  CUERR

  return timesum;
}

float cpuCalculate1nDynamicTimeWarping(const vector<vector<float>> &siga,
                                       const vector<vector<float>> &sigb,
                                       vector<vector<float>> &cpudistResult)
{

  int siga_n = siga.size();
  int sigb_n = sigb.size();
  cpudistResult.resize(siga_n);
  for (int i = 0; i < siga_n; i++)
  {
    cpudistResult[i].resize(sigb_n);
  }
  omp_set_num_threads(CPUTHREADS);
#pragma omp parallel for
  for (int i = 0; i < siga_n; i++)
  {
    for (int j = 0; j < sigb_n; j++)
    {
      cpudistResult[i][j] = cpuDynamicTimeWarping(siga[i], sigb[j]);
    }
  }

  float timesum = 0;
  // printf("cpu Average time use of DTW= %f sec\n", timesum / sigb_n);
  return timesum;
}

float cpuDynamicTimeWarping(const std::vector<float> &seq1,
                            const std::vector<float> &seq2)
{
  vector<vector<float>> score(seq1.size());

  for (int i = 0; i < seq1.size(); i++)
  {
    score[i].resize(seq2.size());
  }

  for (int i = 0; i < seq1.size(); i++)
  {
    for (int j = 0; j < seq2.size(); j++)
    {
      score[i][j] = std::fabs(seq1[i] - seq2[j]);
    }
  }

  for (int i = 1; i < seq1.size(); i++)
  {
    score[i][0] += score[i - 1][0];
  }

  for (int j = 1; j < seq2.size(); j++)
  {
    score[0][j] += score[0][j - 1];
  }

  for (int i = 1; i < seq1.size(); i++)
  {
    for (int j = 1; j < seq2.size(); j++)
    {
      score[i][j] += std::min(std::min(score[i - 1][j], score[i][j - 1]),
                              score[i - 1][j - 1]);
    }
  }

  float diff = score[seq1.size() - 1][seq2.size() - 1];

  return diff;
}

template <class T>
void convertFromString(T &value, string &s)
{
  std::stringstream ss(s);
  ss >> value;
}

__global__ void cuDTW_ultimate(float *g_allColData,
                               unsigned int *g_allColLength,
                               float *g_allRowData,
                               unsigned int *g_allRowLength, float *g_odata)
{

  unsigned int inblockThreadIdx = threadIdx.x;
  unsigned int rowLength =
      g_allRowLength[blockIdx.x]; 
  unsigned int colLength = g_allColLength[blockIdx.y];
  float myNum = 0, myColNum;
  __shared__ unsigned int s_turn;
  __shared__ float preNum[1024], prepreNum[1024], rowData[1024];
  rowData[inblockThreadIdx] =
      g_allRowData[blockIdx.x * blockDim.x + threadIdx.x];
  if (inblockThreadIdx == 0)
  {
    s_turn = 0;
  }
  __syncthreads();
  if (inblockThreadIdx < colLength)
  {
    myColNum = g_allColData[blockIdx.y * 1024 + inblockThreadIdx];
    prepreNum[inblockThreadIdx] = preNum[inblockThreadIdx] = 0;
    int col;
    while (s_turn < colLength + rowLength)
    {
      col = s_turn - inblockThreadIdx;
      if (col >= 0 && col < rowLength)
      {
        if (inblockThreadIdx == 0)
        {
          myNum = preNum[inblockThreadIdx] + fabs(myColNum - rowData[col]);
        }
        else
        {
          if (col == 0)
          {
            myNum =
                preNum[inblockThreadIdx - 1] + fabs(myColNum - rowData[col]);
          }
          else
          {
            myNum = min(min(prepreNum[inblockThreadIdx - 1],
                            preNum[inblockThreadIdx - 1]),
                        preNum[inblockThreadIdx]) +
                    fabs(myColNum - rowData[col]);
          }
        }
      }
      __syncthreads();
      prepreNum[inblockThreadIdx] = preNum[inblockThreadIdx];
      preNum[inblockThreadIdx] = myNum;
      if (inblockThreadIdx == 0)
      {
        // printf("--\nI am first thread of %d block,myIdx=%d,turn=%d,col=%d\n",
        // blockIdx.x,
        //      inblockThreadIdx, s_turn, col);
        s_turn++;
      }
      __syncthreads();
    }
  }
  if (inblockThreadIdx == colLength - 1)
  {
    // if (myNum < 50000)
    //    printf("my blockIdx=%d,my result=%f\n", blockIdx.x,
    //    preNum[inblockThreadIdx]);
    g_odata[blockIdx.y * gridDim.x + blockIdx.x] = myNum;
  }
}

float gpuCalculatemnDynamicTimeWarping_2048(const vector<vector<float>> &siga,
                                            const vector<vector<float>> &sigb,
                                            vector<vector<float>> &gpudistResult)
{
  // return 0;
  int siga_n = siga.size();
  int sigb_n = sigb.size();
  gpudistResult.resize(siga_n);
  for (int i = 0; i < siga_n; i++)
  {
    gpudistResult[i].resize(sigb_n);
  }
  // int siga_length = 0;
  float *d_distResult = NULL;
  float *d_allColData = NULL;
  float *d_allRowData = NULL;
  unsigned int *d_allRowLength;
  unsigned int *d_allColLength;
  // vector<float *> rowDataList(sigb_n);
  cudaMalloc((void **)&d_distResult, siga_n * sigb_n * sizeof(float));
  CUERR
  cudaMalloc((void **)&d_allColData, siga_n * 2048 * sizeof(float));
  CUERR
  cudaMalloc((void **)&d_allRowData, sigb_n * 2048 * sizeof(float));
  CUERR

  
  vector<unsigned int> h_allRowLength(sigb_n);
  for (int i = 0; i < sigb_n; i++)
  {
    h_allRowLength[i] = min(int(sigb[i].size()), 2048);
    cudaMemcpy(&d_allRowData[2048 * i], &sigb[i][0], h_allRowLength[i] * sizeof(float),
               cudaMemcpyHostToDevice);
    CUERR
  }
  cudaMalloc((void **)&d_allRowLength, sigb_n * sizeof(unsigned int));
  CUERR
  cudaMemcpy(d_allRowLength, &h_allRowLength[0], sigb_n * sizeof(unsigned int),
             cudaMemcpyHostToDevice);
  CUERR


  vector<unsigned int> h_allColLength(siga_n);
  for (int i = 0; i < siga_n; i++)
  {
    h_allColLength[i] = min(int(siga[i].size()), 2048);
    cudaMemcpy(&d_allColData[2048 * i], &siga[i][0], h_allColLength[i] * sizeof(float),
               cudaMemcpyHostToDevice);
    CUERR
  }
  cudaMalloc((void **)&d_allColLength, siga_n * sizeof(unsigned int));
  CUERR
  cudaMemcpy(d_allColLength, &h_allColLength[0], siga_n * sizeof(unsigned int),
             cudaMemcpyHostToDevice);
  CUERR

  float timesum = 0;
  dim3 threadsPerBlock(1024);
  dim3 blocksPerGrid(sigb_n, siga_n); 
  cuDTW_2048<<<blocksPerGrid, threadsPerBlock>>>(d_allColData, d_allColLength, d_allRowData,
                                                 d_allRowLength, d_distResult);
  CUERR

  for (int i = 0; i < siga_n; i++)
  {
    cudaMemcpy(&gpudistResult[i][0], &d_distResult[sigb_n * i], sigb_n * sizeof(float),
               cudaMemcpyDeviceToHost);
    CUERR
  }

  cudaFree(d_allColData);
  CUERR
  cudaFree(d_distResult);
  CUERR
  cudaFree(d_allRowData);
  CUERR
  cudaFree(d_allRowLength);
  CUERR

  return timesum;
}

__global__ void cuDTW_2048(float *g_allColData, unsigned int *g_allColLength, float *g_allRowData,
                           unsigned int *g_allRowLength, float *g_odata)
{

  unsigned int rowLength = g_allRowLength[blockIdx.x]; 
  unsigned int colLength = g_allColLength[blockIdx.y];
  float myNum1 = 0, myNum2 = 0, myColNum1, myColNum2;
  __shared__ unsigned int s_turn;
  __shared__ float preNum1[1024], preNum2[1024], prepreNum2[1024], rowData[2048];
  
  rowData[threadIdx.x] = g_allRowData[blockIdx.x * 2048 + threadIdx.x];
  __syncthreads();
  rowData[threadIdx.x + 1024] = g_allRowData[blockIdx.x * 2048 + threadIdx.x + 1024];
  if (threadIdx.x == 0)
  {
    s_turn = 0;
  }
  __syncthreads();
  if (threadIdx.x < (colLength - 1) / 2 + 1)
  {
    
    myColNum1 = g_allColData[blockIdx.y * 2048 + (threadIdx.x) * 2];
    myColNum2 = g_allColData[blockIdx.y * 2048 + (threadIdx.x) * 2 + 1];
    prepreNum2[threadIdx.x] = preNum2[threadIdx.x] = preNum1[threadIdx.x] = 0; //初始化
    int col;
    while (s_turn < (colLength - 1) / 2 + 1 + rowLength)
    {                             
      col = s_turn - threadIdx.x; 
      if (col >= 0 && col < rowLength)
      {
        if (threadIdx.x == 0)
        {
          myNum1 = preNum1[0] + fabs(myColNum1 - rowData[col]);

          if (col == 0)
          { 
            myNum2 = myNum1 + fabs(myColNum2 - rowData[col]);
          }
          else
          {
            myNum2 = min(min(myNum1, preNum1[0]), preNum2[0]) +
                     fabs(myColNum2 - rowData[col]);
          }
        }
        else
        {
          if (col == 0)
          {
            myNum1 = preNum2[threadIdx.x - 1] + fabs(myColNum1 - rowData[col]);
            myNum2 = myNum1 + fabs(myColNum2 - rowData[col]);
          }
          else
          {
            myNum1 = min(min(prepreNum2[threadIdx.x - 1], preNum2[threadIdx.x - 1]),
                         preNum1[threadIdx.x]) +
                     fabs(myColNum1 - rowData[col]);
            myNum2 = min(min(myNum1, preNum1[threadIdx.x]), preNum2[threadIdx.x]) +
                     fabs(myColNum2 - rowData[col]);
          }
        }
      }
      __syncthreads();
      prepreNum2[threadIdx.x] = preNum2[threadIdx.x];
      preNum2[threadIdx.x] = myNum2;
      preNum1[threadIdx.x] = myNum1;
      if (threadIdx.x == 0)
      {
        // printf("--\nI am first thread of %d block,myIdx=%d,turn=%d,col=%d\n",
        // blockIdx.x,
        //      threadIdx.x, s_turn, col);
        s_turn++;
      }
      __syncthreads();
    }
  }
  if (threadIdx.x == (colLength - 1) / 2)
  {
    // if (myNum1 < 50000)
    //    printf("my blockIdx=%d,my result=%f\n", blockIdx.x,
    //    preNum[threadIdx.x]);
    if (colLength % 2 == 0)
      g_odata[blockIdx.y * gridDim.x + blockIdx.x] = myNum2;
    else
      g_odata[blockIdx.y * gridDim.x + blockIdx.x] = myNum1;
  }
}