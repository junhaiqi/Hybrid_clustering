//#include "cuda_def.cuh"
//#include "cuda_proc.h"
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
#define DBL_MAX 1.7976931348623158e+308 // max value
#define PI 3.1415926535898
#define MAX_NUM 2147483647
#define MAXLENGTH 1023
#define CPUTHREADS 36
int CLUSTERNUM = 10;
float THRESHOLD = 90;
#define ENDL printf("\n");

using namespace std;
//////////////////////////////////
//function
//////////////////////////////////

vector<vector<float>> GetAlignSignalList_1(vector<int> &SignalIdList, const string &argv_1);

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
void gpuCluster(const vector<vector<float>> &sigb, float threshold,
                vector<vector<int>> &gpuclusterResult, int maxLocalLength);

vector<int> GetAlignList();

vector<vector<double>> GetAlignSignalList(vector<int> &SignalIdList, const string &argv_1);

vector<string> GetargvList();

vector<float> GetSingleSignalData(int &SignalId, const string &argv_1);

vector<vector<float>> GetCenterSignalList(const int &SignalScale);
void ZScoreNormalize(std::vector<float> &signals, float *avg = NULL, float *stdev = NULL);
template <class T>
void convertFromString(T &value, string &s)
{
    std::stringstream ss(s);
    ss >> value;
}

// using namespace std;

int main(int argc, char **argv)
{
    vector<vector<float>> siga;
    vector<vector<float>> sigb;

    int siga_n = 200;
    int sigb_n = 10000;
    if (argc >= 3)
    {
        siga_n = atoi(argv[1]);
        sigb_n = atoi(argv[2]);
        if (argc == 4)
        {
            THRESHOLD = atoi(argv[3]);
        }
    }

    float gputimeuse, cputimeuse;

    // printf("reading data...\n");

    vector<int> Alignlist = GetAlignList();

    vector<string> argvList = GetargvList();

    sigb = GetAlignSignalList_1(Alignlist, argvList[0]);

    // printf("TestDataGeneration Finished!... \n");
    vector<vector<int>> realgpuclusterResult;
    vector<vector<int>> gpuclusterResult;
    {

        int mod = 2048;
        // printf("choose a mod: 1024 or 2048 ? ");
        // cin >> mod;
        // printf("\nstart cuda cluster... \n");
        CLOCKSTART
        gpuCluster(sigb, Alignlist[1], gpuclusterResult, mod);
        // printf("cuda cluster ");
        CLOCKSTOP
    }

    // cout << gpuclusterResult[0][1] << endl;
    // cout << Alignlist.size() << endl;

    vector<int> templist;
    for (int p = 0; p < gpuclusterResult.size(); p++)
    {
        vector<int> templist;
        for (int q = 0; q < gpuclusterResult[p].size(); q++)
        {
            // cout << i;
            // cout << gpuclusterResult[p][q] << endl;
            int temp = gpuclusterResult[p][q] + 2;
            // cout << gpuclusterResult[p][q] << endl;
            // cout << Alignlist[temp] << endl;
            templist.push_back(Alignlist[temp]);
            // cout << realgpuclusterResult[p][q] << endl;
            // realgpuclusterResult[p].push_back(Alignlist[temp]);
            // cout << "--------------------";
            // cout <<  realgpuclusterResult[i][j] << endl;
            // cout << 1 << endl;
        }
        realgpuclusterResult.push_back(templist);
    }

    vector<vector<float>> gpudistResult;
    vector<vector<float>> cpudistResult;

    ofstream disfile("gpuClusterResult.txt", ios::out);
    string temp;
    for (int i = 0; i < realgpuclusterResult.size(); i++)
    {
        for (int j = 0; j < realgpuclusterResult[i].size(); j++)
        {
            temp = std::to_string(realgpuclusterResult[i][j]);
            disfile << temp;
            disfile << " ";
        }
        disfile << endl;
    }
    // cout << gpudtwresult.size() << endl;
    disfile.close();

    // exit(0);

    return 0;
}

//////////////////////////////////
// function
//////////////////////////////////

vector<vector<float>> GetAlignSignalList_1(vector<int> &SignalIdList, const string &argv_1)
{

    vector<vector<float>> CenterSignalList;
    for (int i = 2; i < SignalIdList.size(); i++)
    {
        vector<float> TempSignal = GetSingleSignalData(SignalIdList[i], argv_1);
        CenterSignalList.push_back(TempSignal);
    }
    return CenterSignalList;
}

vector<int> GetAlignList()
{
    ifstream Alignfile;
    // const char* filename = "ReadyToSort.txt".data();
    Alignfile.open("ReadyToSortfile_1.txt", ios::in);
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

void TestDataGeneration(vector<vector<float>> &siga, vector<vector<float>> &sigb, int siga_n,
                        int sigb_n, bool aeb)
{
    siga.resize(siga_n);
    sigb.resize(sigb_n);
  
    omp_set_num_threads(CPUTHREADS);
#pragma omp parallel for
    for (int i = 0; i < siga_n; i++)
    {
        int siga_length = rand() % 200 + 700; 
        siga[i].resize(siga_length);
        for (int j = 0; j < siga_length; j++)
        {
            siga[i][j] = rand() % 400 + 400; 
        }
    }
    if (siga_n == sigb_n && aeb)
    {
        for (int i = 0; i < sigb_n; i++)
        {
            sigb[i] = siga[i];
        }
        return;
    }
#pragma omp parallel for
    for (int i = 0; i < sigb_n; i++)
    {
        int sigb_length = rand() % 200 + 700; 
        sigb[i].resize(sigb_length);
        for (int j = 0; j < sigb_length; j++)
        {
            sigb[i][j] = rand() % 400 + 400; 
        }
    }
}

float gpuCalculatemnDynamicTimeWarping_1024(const vector<vector<float>> &siga,
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
        cudaMemcpy(&d_allRowData[1024 * i], &sigb[i][0], h_allRowLength[i] * sizeof(float),
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
        h_allColLength[i] = min(int(siga[i].size()), 1024);
        cudaMemcpy(&d_allColData[1024 * i], &siga[i][0], h_allColLength[i] * sizeof(float),
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
    cuDTW_1024<<<blocksPerGrid, threadsPerBlock>>>(d_allColData, d_allColLength, d_allRowData,
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

float cpuCalculatemnDynamicTimeWarping(const vector<vector<float>> &siga,
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

float cpuDynamicTimeWarping(const std::vector<float> &seq1, const std::vector<float> &seq2)
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
            score[i][j] +=
                std::min(std::min(score[i - 1][j], score[i][j - 1]), score[i - 1][j - 1]);
        }
    }

    float diff = score[seq1.size() - 1][seq2.size() - 1];

    return diff;
}
vector<float> GetSingleSignalData(int &SignalId, const string &argv_1)
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
    vector<float> SignalData;
    while (getline(SignalFile, FileLine))
    {
        float SignalValue;
        convertFromString(SignalValue, FileLine);
        SignalData.push_back(SignalValue);
    }
    ZScoreNormalize(SignalData);
    SignalFile.close();
    return SignalData;
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

// vector<vector<float>> GetCenterSignalList(const int &SignalScale)
// {
//     vector<vector<float>> SignalList;
//     for (int i = 0; i < SignalScale; i++)
//     {
//         vector<float> TempSignal = GetSingleSignalData(i, SignalScale);
//         SignalList.push_back(TempSignal);
//     }
//     return SignalList;
// }

void gpuCluster(const vector<vector<float>> &sigb, float threshold,
                vector<vector<int>> &gpuclusterResult, int maxLocalLength)
{
    vector<vector<float>> siga, remainSigb(sigb);
    int sigb_n = sigb.size();
    int remainSigb_n = remainSigb.size();
    int siga_n, usedCount = 0;
    vector<int> label(remainSigb_n), aindex(remainSigb_n), bindex(remainSigb_n),
        used(remainSigb_n, 0);

    {
        for (int i = 0; i < remainSigb_n; i++)
        {
            label[i] = i;
            bindex[i] = i;
        }
    }

    for (int loop = 0; loop < 10; loop++)
    {
        // printf("loop=%d-----------------------\n", loop);
        if (loop > 0)
        {
            remainSigb.clear();
            for (int i = 0; i < sigb_n; i++)
            {
                if (used[i] == 0)
                {
                    remainSigb.push_back(sigb[i]);
                    bindex[remainSigb.size() - 1] = i;
                }
            }
        }
        remainSigb_n = remainSigb.size();
        // printf("sigb_length=%d\n", remainSigb_n);
        if (remainSigb.size() == 0)
        {
            break;
        }
        {
            //通过随机的方式来初始化这一次要用的数据
            int randnum;
            siga.clear();
            for (int i = 0; i < remainSigb_n; i++)
            {
                randnum = rand() % 1000;
                // printf("%d ", randnum);
                
                if (randnum < 100000.0 / remainSigb_n)
                {
                    siga.push_back(remainSigb[i]);
                    aindex[siga.size() - 1] = bindex[i];
                    // printf("picked i=%d\n", i);
                }
            }
            // printf("siga_length=%d\n", siga.size());
        }

        siga_n = siga.size();
        vector<vector<float>> gpudistResult;
        if (maxLocalLength == 1024)
        {
            gpuCalculatemnDynamicTimeWarping_1024(siga, remainSigb, gpudistResult);
        }
        else if (maxLocalLength == 2048)
        {
            gpuCalculatemnDynamicTimeWarping_2048(siga, remainSigb, gpudistResult);
        }
        else
        {
            printf("wrong max local length!\n");
            exit(-2);
        }
        if (loop == 0)
        {

            float maxVal = *max_element(gpudistResult[0].begin(), gpudistResult[0].end());
            float minVal = *min_element(gpudistResult[0].begin(), gpudistResult[0].end());
            THRESHOLD = (maxVal + minVal) / 2;
            // printf("maxDist,minDist,thredhold=%f %f %f\n", maxVal, minVal, THRESHOLD);
        }

        for (int i = 0; i < siga_n; i++)
        {
            if (used[aindex[i]] == 0)
            {
                for (int j = 0; j < remainSigb_n; j++)
                {
                    if (gpudistResult[i][j] < THRESHOLD)
                    {
                        label[bindex[j]] = label[aindex[i]];
                        used[bindex[j]] = 1;
                    }
                }
            }
        }
    }
    for (int i = 0; i < sigb_n; i++)
    {
        // if (i % CLUSTERNUM == 0) {
        //     printf("\n");
        // }
        // printf("%d ", label[i]);
        bool notFound = true;
        for (int j = 0; j < gpuclusterResult.size(); j++)
        {
            if (label[i] == gpuclusterResult[j][0])
            {
                if (label[i] != i)
                {
                    gpuclusterResult[j].push_back(i);
                }
                notFound = false;
                break;
            }
        }
        if (notFound)
        {
            vector<int> thisCluster;
            thisCluster.push_back(label[i]);
            if (label[i] != i)
            {
                thisCluster.push_back(i);
            }
            gpuclusterResult.push_back(thisCluster);
        }
    }
    // printf("cluster numbers=%d\n", gpuclusterResult.size());

    // // for (int i = 0; i < gpuclusterResult.size(); i++) {
    // //     for (int j = 0; j < gpuclusterResult[i].size(); j++) {
    // //         printf("%d ", gpuclusterResult[i][j]);
    // //     }
    //     ENDL;
    // // }

    // printf("\n");
}
void ZScoreNormalize(std::vector<float> &signals, float *avg, float *stdev)
{
    // CLOCKSTART
    float sum = std::accumulate(signals.begin(), signals.end(), 0.0);
    float mean = sum / signals.size();

    float acc = 0.0;
    for (size_t i = signals.size(); i--;)
    {
        signals[i] = signals[i] - mean;
        acc += signals[i] * signals[i];
    }

    float deviation = std::sqrt(acc / signals.size());

    for (size_t i = signals.size(); i--;)
    {
        signals[i] /= deviation;
    }

    if (avg)
    {
        *avg = mean;
    }
    if (stdev)
    {
        *stdev = deviation;
    }
    // printf("ZScoreNorm ");
    // CLOCKSTOP
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

__global__ void cuDTW_1024(float *g_allColData, unsigned int *g_allColLength, float *g_allRowData,
                           unsigned int *g_allRowLength, float *g_odata)
{

    unsigned int rowLength = g_allRowLength[blockIdx.x]; 
    unsigned int colLength = g_allColLength[blockIdx.y];
    float myNum = 0, myColNum;
    __shared__ unsigned int s_turn;
    __shared__ float preNum[1024], prepreNum[1024], rowData[1024];
    rowData[threadIdx.x] =
        g_allRowData[blockIdx.x * 1024 + threadIdx.x]; 
    if (threadIdx.x == 0)
    {
        s_turn = 0;
    }
    __syncthreads();
    if (threadIdx.x < colLength)
    {
        myColNum = g_allColData[blockIdx.y * 1024 + threadIdx.x];
        prepreNum[threadIdx.x] = preNum[threadIdx.x] = 0;
        int col;
        while (s_turn < colLength + rowLength)
        {
            col = s_turn - threadIdx.x;
            if (col >= 0 && col < rowLength)
            {
                if (threadIdx.x == 0)
                {
                    myNum = preNum[threadIdx.x] + fabs(myColNum - rowData[col]);
                }
                else
                {
                    if (col == 0)
                    {
                        myNum = preNum[threadIdx.x - 1] + fabs(myColNum - rowData[col]);
                    }
                    else
                    {
                        myNum = min(min(prepreNum[threadIdx.x - 1], preNum[threadIdx.x - 1]),
                                    preNum[threadIdx.x]) +
                                fabs(myColNum - rowData[col]);
                    }
                }
            }
            __syncthreads();
            prepreNum[threadIdx.x] = preNum[threadIdx.x];
            preNum[threadIdx.x] = myNum;
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
    if (threadIdx.x == colLength - 1)
    {
        // if (myNum < 50000)
        //    printf("my blockIdx=%d,my result=%f\n", blockIdx.x,
        //    preNum[threadIdx.x]);
        g_odata[blockIdx.y * gridDim.x + blockIdx.x] = myNum;
    }
}

__global__ void cuDTW_2048(float *g_allColData, unsigned int *g_allColLength, float *g_allRowData,
                           unsigned int *g_allRowLength, float *g_odata)
{

    unsigned int rowLength = g_allRowLength[blockIdx.x]; // block的x是sigb，block的y是siga
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
        prepreNum2[threadIdx.x] = preNum2[threadIdx.x] = preNum1[threadIdx.x] = 0;
        int col;
        while (s_turn < (colLength - 1) / 2 + rowLength)
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

// vector<vector<float>> GetAlignSignalList(vector<int> &SignalIdList, const int &SignalScale)
// {
//     vector<vector<float>> CenterSignalList;
//     for (int i = 2; i < SignalIdList.size(); i++)
//     {
//         vector<float> TempSignal = GetSingleSignalData(SignalIdList[i], SignalScale);
//         CenterSignalList.push_back(TempSignal);
//     }
//     return CenterSignalList;
// }