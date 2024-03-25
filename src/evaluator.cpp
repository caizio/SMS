#include "evaluator.h"
using namespace std;

// 重载pair的等号判断，避免边的节点先大后小的情况
bool operator==(const pair<int,int> &pair1,const pair<int,int> &pair2){
    if((pair1.first == pair2.first && pair1.second == pair2.second) || (pair1.first == pair2.second && pair1.second == pair2.first)){
        return true;
    }
    return false;
}

Evaluator::Evaluator(){};
Evaluator::~Evaluator(){};
// 通过随机抽样（分别从预测边和要预测的边抽样，预测边分值大计1分，相等计0.5分）的方式近似计算auc
double Evaluator::accAUC(const vector<vector<int>> &train,const vector<vector<int>> &test,const vector<vector<double>> &similarMatrix,int AUCNum){
    int n = train.size();
    vector<pair<int,int>> test_edge;
    vector<pair<int,int>> non_edge;
    if(n == 0 || test.size() == 0) {
        return -1;
    };
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(test[i][j] != 0){
                test_edge.emplace_back(make_pair(i,j));
            }
            if(train[i][j] == 0 && test[i][j] == 0){
                non_edge.emplace_back(make_pair(i,j));
            }
        }
    }
    // cout << "测试集边数:"  << test_edge.size() << endl << "未观察边数:" << non_edge.size() << endl; 
    double sum = 0;
    srand(time(0));
    int test_size = test_edge.size(), non_size = non_edge.size();
    for(int i = 0; i < AUCNum; i++){
        int rand1 = rand() % test_size, rand2 = rand() % non_size;
        double test_score = similarMatrix[test_edge[rand1].first][test_edge[rand1].second];
        double non_score = similarMatrix[non_edge[rand2].first][non_edge[rand2].second];
        if(test_score > non_score){
            sum += 1;
        }else if(test_score == non_score){
            sum += 0.5;
        }
    }
    double auc = sum / AUCNum;
    return auc;
}

// 计算nDCG的分母
double Evaluator::_IDGC(int L,int base){
    double res = 0;
    for(int i = 1; i <= L; i++){
        res += pow(log2(i + 1),base);
    }
    return res;
};
// 计算nDCG指标，test为测试集即正样本在所有预测样本中的排序，base为-1
// 2022_11_28
double Evaluator::_nDCG(const std::vector<int> &test,int base){
    int L = test.size();
    double DCG = 0;
    double IDCG = this->_IDGC(L);
    for(auto &p:test){
       DCG +=  pow(log2(p + 1),base);
    }
    return DCG / IDCG;
}
// 返回测试集边的排序，predictionEdge：算法预测的边，按照得分从到到小的排序
// 2022_11_28
std::vector<int> Evaluator::_find_test_rank(const std::vector<std::pair<int,int>> &predictionEdge,const std::vector<std::pair<int,int>> &test){
    vector<int> res;
    int n = predictionEdge.size();
    for(auto &p:test){
        for(int i = 0; i < n; i++){
            if(predictionEdge[i] == p){
                res.push_back(i+1);
            }
        }
    }
    return res;
}
// 计算曲线的面积
double Evaluator::area(std::vector<std::pair<double,double>> &xy){
    int n = xy.size();
    double res = 0;
    for(int i = 0; i < n - 1; i++){
        res += (xy[i].second + xy[i+1].second) * (xy[i+1].first - xy[i].first) / 2;
    }
    return res;
}
// 计算nDCG,参数：相似度矩阵，训练集，测试集边
// 2022_11_30
double Evaluator::acc_nDCG(std::vector<std::vector<double>> &similar,std::vector<std::vector<int>> &train,std::vector<std::pair<int,int>> &testPair){
    auto pairRank = this->_similiarToRankALL2(similar,train);
    auto pairIndex = this->_find_test_rank(pairRank,testPair);
    auto nDCG = this->_nDCG(pairIndex);
    return nDCG;
}

// 计算排序分,tsetRank：测试集样本在所有要预测样本中的排序，H为要预测样本的数量,ep：测试集样本数量
double Evaluator::_rank_score(const std::vector<int> &testRank,int H){
    int ep = testRank.size();
    double res = 0;
    for(int i = 0; i < ep; i++){
        res += testRank[i];
    }
    return res / ep / H;
}

// 计算排序分
double Evaluator::acc_rank_score(std::vector<std::vector<double>> &similar,std::vector<std::vector<int>> &train,std::vector<std::pair<int,int>> &trainPair,std::vector<std::pair<int,int>> &testPair){
    auto pairRank = this->_similiarToRankALL2(similar,train);
    auto pairIndex = this->_find_test_rank(pairRank,testPair);
    int n = train.size();
    int H = n * (n - 1) / 2 - trainPair.size();
    return this->_rank_score(pairIndex,H);
}

// 从相似性矩阵中挑选出要预测的边（全集-训练集），用sort
// 2022_11_30
vector<pair<int,int>> Evaluator::_similiarToRankALL2(vector<vector<double>> &similarM,vector<vector<int>> &train){
    int n = similarM.size();
    vector<pair<int,int>> res;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(train[i][j] == 0){
                res.push_back(make_pair(i,j));
            }
        }
    }
    sort(res.begin(),res.end(),[&similarM](const pair<int,int> &edge1,const pair<int,int> &edge2){
        if(similarM[edge1.first][edge1.second] <= similarM[edge2.first][edge2.second]){
            return false;
        }
        return true;
    });
    return res;
}
// 从相似性矩阵中挑选出要预测的边（全集-训练集），用sort
// 2022_11_30
std::vector<edgeScore> Evaluator::_similiarToRankALL3(std::vector<std::vector<double>> &similarM,std::vector<std::vector<int>> &train,std::vector<std::vector<int>> &test){
    int n = similarM.size();
    vector<edgeScore> res;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(train[i][j] == 0){
                edgeScore es;
                es.x = i;
                es.y = j;
                es.score = similarM[i][j];
                if(test[i][j] != 0){
                    es.TF = 1;
                }else{
                    es.TF = 0;
                }
                res.push_back(es);
            }
        }
    }
    sort(res.begin(),res.end(),[](const edgeScore &edge1,const edgeScore &edge2){
        if(edge1.score <= edge2.score){
            return false;
        }
        return true;
    });
    return res;
}

// 计算ROC和PR，参数：prediction为预测边，按照得分排序,T为测试集数量（所有正样本），F为要U - 训练集-测试集数量（所有负样本）
// 2022_11_30
ROCandPR Evaluator::_ROC_PR(std::vector<edgeScore> &predictionEdge,int T,int F){
    int n = predictionEdge.size();
    // 最开始，所有数据判为反例
    int TP = 0,FP = 0,FN = T,TN = F;
    double TPR = 0;
    double FPR = 0;
    double P = 1;
    double R = 0;
    ROCandPR roc_pr;
    roc_pr.pr.push_back(make_pair(R,P));
    roc_pr.roc.push_back(make_pair(FPR,TPR));
    for(int i = 0; i < n; i++){
        if(predictionEdge[i].TF == true){
            TP +=1;
            FN -=1;
        }else{
            FP +=1;
            TN -=1;
        }
        TPR = double(TP) / T;
        FPR = double(FP) / F;
        P = double(TP) / (TP + FP);
        R = double(TP) / T;
        roc_pr.roc.push_back(make_pair(FPR,TPR));
        roc_pr.pr.push_back(make_pair(R,P));
        roc_pr.precision.push_back(P);
        roc_pr.recall.push_back(R);
        roc_pr.fpr.push_back(FPR);
    }
    roc_pr.roc_area = this->area(roc_pr.roc);
    roc_pr.pr_area = this->area(roc_pr.pr);
    return roc_pr;
}
// 计算ROC和PR
ROCandPR Evaluator::acc_ROC_PR(std::vector<std::vector<double>> &similar,std::vector<std::vector<int>> &train,std::vector<std::vector<int>> &test,std::vector<std::pair<int,int>> &trainPair,std::vector<std::pair<int,int>> &testPair){
    auto predictionEdge = this->_similiarToRankALL3(similar,train,test);
    int T = testPair.size();
    int F = similar.size() * (similar.size() - 1) / 2 - trainPair.size() - testPair.size();
    auto roc_pr = this->_ROC_PR(predictionEdge,T,F);
    return roc_pr;
}

baseStruct::LinkRall Evaluator::acc_all_evaluators(std::vector<std::vector<double>> &similar,std::vector<std::vector<int>> &train,std::vector<std::vector<int>> &test,std::vector<std::pair<int,int>> &trainPair,std::vector<std::pair<int,int>> &testPair){
    baseStruct::LinkRall temp;
    temp.auc = accAUC(train,test,similar);
    auto p =  acc_ROC_PR(similar,train,test,trainPair,testPair);
    temp.AUROC = p.roc_area;
    temp.AUPRC = p.pr_area;
    // temp.rankScore = acc_rank_score(similar,train,trainPair,testPair);
    temp.nDCG = acc_nDCG(similar,train,testPair);
    temp.precision = p.precision;
    temp.recall = p.recall;
    temp.FPR = p.fpr;
    temp.roc = p.roc;
    temp.pr = p.pr;
    return temp;
}

