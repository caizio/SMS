#ifndef _EVALUATOR_DEFINE
#define _EVALUATOR_DEFINE
#include <cmath>
#include <ctime>
#include <vector>
#include <utility>
#include <algorithm>
#include "baseStruct.hpp"

struct train_test_data{
    std::vector<std::vector<int>> train;
    std::vector<std::vector<int>> test;
    std::vector<std::pair<int,int>> &trainPair;
    std::vector<std::pair<int,int>> &testPair;
};

// 边的数据，xy：边，TF：1表示正样本（测试集P），0表示负样本（U-T-P），score为算法预测得分
struct edgeScore{
    int x;
    int y;
    bool TF;
    double score;
};

struct ROCandPR{
    std::vector<std::pair<double,double>> roc;
    double roc_area;
    std::vector<std::pair<double,double>> pr;
    double pr_area;
    std::vector<double> precision;
    std::vector<double> recall;  
    std::vector<double> fpr;  
};
class Evaluator{
public:
    Evaluator();
    ~Evaluator();
    double accAUC(const std::vector<std::vector<int>> &train,const std::vector<std::vector<int>> &test,const std::vector<std::vector<double>> &similarMatrix,int AUCNum = 672400);
    double _IDGC(int L,int base = -1);
    double _nDCG(const std::vector<int> &test,int base = -1);
    std::vector<std::pair<int,int>> _similiarToRankALL2(std::vector<std::vector<double>> &similarM,std::vector<std::vector<int>> &train);
    std::vector<edgeScore> _similiarToRankALL3(std::vector<std::vector<double>> &similarM,std::vector<std::vector<int>> &train,std::vector<std::vector<int>> &test);
    std::vector<int> _find_test_rank(const std::vector<std::pair<int,int>> &predictionEdge, const std::vector<std::pair<int,int>> &test);
    double area(std::vector<std::pair<double,double>> &xy);
    ROCandPR _ROC_PR(std::vector<edgeScore> &predictionEdge,int T,int F);
    ROCandPR acc_ROC_PR(std::vector<std::vector<double>> &similar,std::vector<std::vector<int>> &train,std::vector<std::vector<int>> &test,std::vector<std::pair<int,int>> &trainPair,std::vector<std::pair<int,int>> &testPair);
    double acc_nDCG(std::vector<std::vector<double>> &similar,std::vector<std::vector<int>> &train,std::vector<std::pair<int,int>> &testPair);
    double _rank_score(const std::vector<int> &testRank,int H);
    double acc_rank_score(std::vector<std::vector<double>> &similar,std::vector<std::vector<int>> &train,std::vector<std::pair<int,int>> &trainPair,std::vector<std::pair<int,int>> &testPair);
    baseStruct::LinkRall acc_all_evaluators(std::vector<std::vector<double>> &similar,std::vector<std::vector<int>> &train,std::vector<std::vector<int>> &test,std::vector<std::pair<int,int>> &trainPair,std::vector<std::pair<int,int>> &testPair);
};

#endif