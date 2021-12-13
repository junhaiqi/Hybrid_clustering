#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <algorithm>
#include <omp.h>
#include <stdio.h>

#define CPUTHREADS 30
using namespace std;

template <class T> 
void convertFromString(T &value, string &s);

vector<double> GetSingleSignalData(int &SignalId, const string &argv_1);

vector<string> GetargvList();

double dtw(const vector<double> &s1, const vector<double> &s2);

void ZScoreNormalize(vector<double> &signals);

double GetZscoreDTW(const string &argv_1,int& SignalId1,int& SignalId2);

vector<vector<int> > GetInitClusterResult(int& DataScale);

vector<int > Get_clustering_info();

vector<vector<int> > MergedGoodCluster(vector<vector<int> >& InitClusteringResult,
                     const string &argv_1,const float& threshold,int& RandomCount,int& goodIndex);

int main(){
    vector<int > Th_Scale_List = Get_clustering_info();
    int scale = Th_Scale_List[0];
    float threshold = Th_Scale_List[1];
    int RandomCount = Th_Scale_List[2];
    int goodIndex = Th_Scale_List[3];
    vector<string > argvList = GetargvList();
    vector<vector<int> > InitResult = GetInitClusterResult(scale); 
    vector<vector<int> > goodcluster = MergedGoodCluster(InitResult,argvList[0],threshold,RandomCount,goodIndex);//将好的类合并,放进goodcluster.txt文件中

    for(int i = 0; i < goodcluster.size(); i++){
        for(int j = 0; j < goodcluster[i].size(); j++){
            if(goodcluster[i][j] > scale){
                cout << "算法异常值:" << goodcluster[i][j] << endl;
            }
        }
    }
    ofstream disfile("goodclusterfile.txt",ios::out);
    string temp;
    for(int i = 0; i < goodcluster.size(); i++){
        for(int j = 0; j < goodcluster[i].size(); j++){
            temp = std::to_string(goodcluster[i][j]);
            disfile << temp;
            disfile << " ";
        }
        disfile << endl;
    }
    disfile.close();
    return 0;
}

template <class T> 
void convertFromString(T &value, string &s) {
  std::stringstream ss(s);
  ss >> value;
}

vector<int > Get_clustering_info(){
    ifstream Th_Scale_file;
    Th_Scale_file.open("th_goodcount_file.txt",ios::in);
    if(!Th_Scale_file.is_open()){
       cout << "th_goodcount_file open error!" << endl;
    }
    string FileLine;
    vector<int > Th_Scale_List;
    while(getline(Th_Scale_file,FileLine)){
        int th_scale_temp;
        convertFromString(th_scale_temp,FileLine);
        Th_Scale_List.push_back(th_scale_temp);
    }

    Th_Scale_file.close();

    return Th_Scale_List;
}

vector<string> GetargvList(){
  ifstream argvfile;
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

double dtw(const vector<double> &s1, const vector<double> &s2)
{
   int len1 = s1.size();
   int len2 = s2.size();
   double dp[len2 + 1], dppre[len2 + 1];

   //printf("len1=%d;len2=%d\n", len1, len2);
   for (int j = 0; j < len2; j++)
   {
      dppre[j] = fabs(s1[0] - s2[j]);
   }
   for (int i = 1; i < len1; i++)
   {
      dp[0] = fabs(s1[i] - s2[0]);
      for (int j = 1; j < len2; j++)
      {
         //    printf("%.2f ",dp[j]);        //
         dp[j] = std::min(dp[j - 1], dppre[j - 1]);
         dp[j] = std::min(dp[j], dppre[j]) + fabs(s1[i] - s2[j]);
      }
      //printf("\n");         //
      for (int j = 0; j < len2; j++)
      {
         dppre[j] = dp[j];
      }
   }
   return dp[len2 - 1];
}

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

double GetZscoreDTW(const string &argv_1,int& SignalId1,int& SignalId2){
    vector<double> Signal1 = GetSingleSignalData(SignalId1,argv_1);
    vector<double> Signal2 = GetSingleSignalData(SignalId2,argv_1);
    ZScoreNormalize(Signal1);
    ZScoreNormalize(Signal2);
    double DTWDistance;
    DTWDistance = dtw(Signal1,Signal2);
    return DTWDistance;
}

vector<vector<int> > GetInitClusterResult(int& DataScale){
    ifstream InitRes;
    vector<int> SingleClstr;
    vector<vector<int> > Allclstr;
    string flag = ">";
    string flag1 = ".";
    char flag2 = *flag.data();
    char flag3 = *flag1.data();
    string temp = std::to_string(DataScale);
    temp = temp + ".clstr";
    InitRes.open(temp,ios::in);
    while (!InitRes.is_open())
    {
        cout << "File2 open error!" << endl;
    }
    string FileLine;
    int indicator = 0;
    while(getline(InitRes,FileLine)){
        int start = 0;
        int end = 0;
        int DNAIndicator = 0;
        if(FileLine[0] == flag2 && indicator != 0) {
            Allclstr.push_back(SingleClstr);
            SingleClstr.clear();
        }
        indicator += 1;
        if(FileLine[0] != flag2){
        for(int i = 0; i<FileLine.size(); i++){
            start = FileLine.find(">");
            end = FileLine.find(".") - start;
        }
        string NewLine = FileLine.substr(start + 1,end - 1);
        convertFromString(DNAIndicator,NewLine);
        SingleClstr.push_back(DNAIndicator-1);
        }
    }
    Allclstr.push_back(SingleClstr);
    return Allclstr;
}

vector<vector<int> > MergedGoodCluster(vector<vector<int> >& InitClusteringResult,const string &argv_1,const float& threshold,int& RandomCount,int& goodIndex){
      vector<vector<int> > InitClusteringResultTemp;
      // 将符合规模的类放进列表中
      for(int i = 0; i < InitClusteringResult.size(); i++){ 
          if(InitClusteringResult[i].size() > goodIndex){
            InitClusteringResultTemp.push_back(InitClusteringResult[i]);
          }
       }
       vector<vector<int > > goodcluster;
       vector<vector<int > > RealGoodcluster;
       
       for(int i = 0; i < InitClusteringResultTemp.size() - 1; i++){
         int count_1 = std::count(InitClusteringResultTemp[i].begin(), InitClusteringResultTemp[i].end(), -1);
          if(std::count(InitClusteringResultTemp[i].begin(), InitClusteringResultTemp[i].end(), -1)){
            //   cout << InitClusteringResultTemp[i][InitClusteringResultTemp[i].size() - 1] << endl;
              continue;
          }
          else{
              random_shuffle(InitClusteringResultTemp[i].begin(), InitClusteringResultTemp[i].end());
              vector<int > RandomCheck_1;
              for(int j  = 0; j < RandomCount; j++){
                  RandomCheck_1.push_back(InitClusteringResultTemp[i][j]);
              }
            
            
            for(int t = i + 1; t < InitClusteringResultTemp.size(); t++){
                int count_2 = std::count(InitClusteringResultTemp[t].begin(), InitClusteringResultTemp[t].end(), -1);
                if(std::count(InitClusteringResultTemp[t].begin(), InitClusteringResultTemp[t].end(), -1)){
                //    cout << InitClusteringResultTemp[t][InitClusteringResultTemp[t].size() - 1] << endl;
                   continue;
                }
                int bool_1 = 0; 
                random_shuffle(InitClusteringResultTemp[t].begin(),InitClusteringResultTemp[t].end());
                vector<int > RandomCheck_2;
                for(int j = 0; j < RandomCount; j++){
                  RandomCheck_2.push_back(InitClusteringResultTemp[t][j]);
                }
                
                omp_set_num_threads(30);
                #pragma omp parallel for
                for(int n = 0; n < RandomCount; n++){
                    // if(bool_1 == 1){
                    //     // break;
                    // }
                    for(int m = 0; m < RandomCount && bool_1 != 1; m++){
                        // RandomCheck_1[n] = RandomCheck_1[n] - 1;
                        // RandomCheck_1[m] = RandomCheck_1[m] - 1;
                        if(GetZscoreDTW(argv_1, RandomCheck_1[n],RandomCheck_2[m]) > threshold){
                           bool_1 = 1;
                           break; 
                        }
                    }
                }
                
                if(bool_1 == 0){
                    for(int w = 0; w < InitClusteringResultTemp[t].size();w++){
                        InitClusteringResultTemp[i].push_back(InitClusteringResultTemp[t][w]);
                        int check_count = std::count(InitClusteringResultTemp[i].begin(),InitClusteringResultTemp[i].end(),-1);
                        if(check_count){
                           cout << "出现问题的位置:" << i << endl;
                        }
                        
                    }
                    InitClusteringResultTemp[t].push_back(-1);

                }
            }
          }
       } 
       for(int p = 0; p < InitClusteringResultTemp.size(); p++){
           int temp = InitClusteringResultTemp[p].size() - 1;
           if(InitClusteringResultTemp[p][temp] != -1){
               goodcluster.push_back(InitClusteringResultTemp[p]);
           }
       }
       for(int i = 0; i < goodcluster.size(); i++){
           vector<int> list;
           for(int j = 0; j < goodcluster[i].size(); j++){
               if(goodcluster[i][j] != -1){
                   list.push_back(goodcluster[i][j]);
               }
           }
           RealGoodcluster.push_back(list);
       }
    return RealGoodcluster;      
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
