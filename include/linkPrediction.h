// 链路预测算法
#ifndef _LINKPREDICTION_
#define _LINKPREDICTION_

#include <vector>
#include <string>
#include <unordered_set>
#include <ctime>
#include <array>
#include "graph.h"
#include <set>
#include <map>
#include "matrix.h"
#include "bioAlgorith.h"
#include "dataDeal.h"
#include <cstdlib>

using namespace std;

using fun_ptr = vector<vector<double>>(*)(UDG &G);

// 基于邻接矩阵的链路预测算法
class LinkPrediction{
public:
    LinkPrediction();
    LinkPrediction(UDG &G);
    vector<vector<double>> mmain(const std::string &name,UDG &G,double parameter = 0,double parameter2 = 2,double parameter3 = 1,string datasetName = "");
    vector<vector<double>> PA(UDG &G);
    // 基于二阶路径（公共邻居）的算法
    vector<vector<double>> CN(UDG &G);
    vector<vector<double>> RA(UDG &G);
    vector<vector<double>> JC(UDG &G);
    vector<vector<double>> CCLP(UDG &G);
    // 基于三阶路径的算法
    vector<vector<double>> A3(UDG &G);
    vector<vector<double>> LP(UDG &G, double lam = 0.01);
    vector<vector<double>> L3(UDG &G);
    vector<vector<double>> CH2_L3(UDG &G);
    vector<vector<double>> SIM(UDG &G);
    vector<vector<double>> Sim(UDG &G, std::vector<std::vector<double>> &similar);
    vector<vector<double>> maxSIM(UDG &G);
    vector<vector<double>> SMS(UDG &G,const string &rw = "JC",double degree = 1);
    vector<vector<double>> SMS(UDG &G,std::vector<std::vector<double>> &similar);
    vector<vector<double>> SMS_KNCN(UDG &G,double degree = 0,double parameter = 1);
    vector<vector<double>> maxSMS(UDG &G,const string &rw = "CN",bool degree = false);
    vector<vector<double>> maxSMS(UDG &G,const std::vector<std::vector<double>> &similar);
    // 基于序列的算法
    vector<vector<double>> SMS_sequence(UDG &G,map<pair<int,int>,double> &data,map<int,string>&seq, ofstream &f);
    vector<vector<double>> DSMS_sequence(UDG &G,map<pair<int,int>,double> &data,map<int,string>&seq, ofstream &f);
    vector<vector<double>> maxSMS_sequence(UDG &G,map<pair<int,int>,double> &data,map<int,string>&seq, ofstream &f);
    // 基于三阶路径和序列的混合算法
    vector<vector<double>> SMS_mix(UDG &G,map<pair<int,int>,double> &data,map<int,string>&seq, ofstream &f,double lam = 0);
    vector<vector<double>> DSMS_mix(UDG &G,map<pair<int,int>,double> &data,map<int,string>&seq, ofstream &f,double lam = 0);
    vector<vector<double>> maxSMS_mix(UDG &G,map<pair<int,int>,double> &data,map<int,string>&seq, ofstream &f,double lam = 0);
private:
    map<string,function<vector<vector<double>>(UDG &)>> _similar;
    double _A3(UDG &G,int x,int y);
    double _L3(UDG &G,int x,int y);
    double _C3(UDG &G,int x,int y,vector<vector<double>> &similarMatrix);
    double _CH2_L3(UDG &G,int x,int y);
    double _SMS(UDG &G,vector<vector<double>> &similarMatrix,int x,int y,double degree = 1);
    vector<vector<double>> _CN(UDG &G);
    vector<vector<double>> _AA(UDG &G);
    vector<vector<double>> _RA(UDG &G);
    vector<vector<double>> _JC(UDG &G);
    vector<vector<pair<int,double>>> _JC_(UDG &G);
    vector<vector<double>> _Sorenson(UDG &G);
    vector<vector<double>> _HPI(UDG &G);
    vector<vector<double>> _HDI(UDG &G);
    vector<vector<double>> _LHNI(UDG &G);
    vector<vector<double>> _Salton(UDG &G);
    double acc_similarBySequence(int x,int y,map<pair<int,int>,double> &similar,map<int,string>&seq);
    double acc_similarBySequence2(int x,int y,map<pair<int,int>,double> &similar,map<int,string>&seq,ofstream &f);
};

#endif