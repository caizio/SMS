// 基础数据结构

#ifndef BASE_STRUCT_
#define BASE_STRUCT_
#include <vector>
#include <string>
#include <utility>
#include <cstring>
#include "matrix.h"
namespace baseStruct
{
    // ROC曲线
    struct ROC{
        std::string name;
        std::vector<double> FPR;
        std::vector<double> TPR;
    };
    // PR曲线
    struct PR{
        std::string name;
        std::vector<double> precision;
        std::vector<double> recall;
    };
    // 链路预测的结果
    struct LinkRall{
        std::string name;
        // 以抽样的方式计算的AUC
        double auc;
        // 通过ROC曲线绘制的曲线求的面积（理论上等同于auc，但是实际不相等，没有找到为什么？）
        double AUROC;
        // 通过PR曲线绘制的曲线求的面积
        double AUPRC;
        // 排序分
        double rankScore;
        // 归一化折损累计增益
        double nDCG;
        // 精度，查准率,ps：precision[49]表示前50的精度,因为数组从0开始
        std::vector<double> precision;
        // TPR/召回率，查全率,ps：recall[49]表示前50的召回率
        std::vector<double> recall;
        // 假正例，ROC曲线的横坐标
        std::vector<double> FPR;
        // roc、pr曲线的值
        std::vector<std::pair<double,double>> roc;
        std::vector<std::pair<double,double>> pr;
        public:
        LinkRall(){
            this->auc = 0;
            this->AUROC = 0;
            this->AUPRC = 0;
            this->rankScore = 0;
            this->nDCG = 0;
            this->precision.resize(500,0);
            this->recall.resize(500,0);
            this->FPR.resize(500,0);
        }
        void operator+(LinkRall &lr){
            this->auc += lr.auc;
            this->AUROC += lr.AUROC;
            this->AUPRC += lr.AUPRC;
            this->rankScore += lr.rankScore;
            this->nDCG += lr.nDCG;
            this->precision = this->precision + lr.precision;
            this->recall = this->recall + lr.recall;
            this->FPR = this->FPR + lr.FPR;
            int n = this->roc.size();
            if(n == 0){
                this->roc = lr.roc;
                this->pr = lr.pr;
                return;
            }
            if(lr.roc.size() != n || this->roc.size() != this->pr.size() || lr.roc.size() != lr.pr.size()){
                // This warning can be ignored and has no impact on the experimental results. 
                // 该错误可以忽视，对实验结果没有影响，roc和pr的数据仅用于画图
                // 造成的原因是十折交叉验证中，第十份测试集的边数量和其他九份不一样
                cout << "LinkRall::operator+::warning::size wrong" << endl;      
            }
            int m = lr.roc.size();
            n = min(n,m);
            for(int i = 0; i < n; i++){
                this->roc[i].first += lr.roc[i].first;
                this->roc[i].second += lr.roc[i].second;
                this->pr[i].first += lr.pr[i].first;
                this->pr[i].second += lr.pr[i].second;
            }
        }
        void getAverage(int n){
            this->auc /= n;
            this->AUROC /= n;
            this->AUPRC /= n;
            this->rankScore /= n;
            this->nDCG /= n;
            this->precision = this->precision * (1.0 / n);
            this->recall = this->recall * (1.0 / n);
            this->FPR = this->FPR * (1.0 / n);
            int m = this->roc.size();
            for(int i = 0; i < m; i++){
                this->roc[i].first = this->roc[i].first / n;
                this->roc[i].second = this->roc[i].second / n;
                this->pr[i].first = this->pr[i].first / n;
                this->pr[i].second = this->pr[i].second / n;
            }   
        }
        void show(){
            cout << name << endl;
            if(precision.size() < 500 || recall.size() < 500 || FPR.size() < 500){
                cout << "LinkRall::show()::error" << endl;
                return;
            }
            cout << "P:" << precision[99] << "   " << precision[199] << "   " << precision[499] << endl; 
            cout << "R\\TPR:" << recall[99] << "   " << recall[199] << "   " << recall[499] << endl; 
            cout << "FPR:" << FPR[99] << "   " << FPR[199] << "   " << FPR[499] << endl; 
            cout << "nDCG:" << nDCG << endl;  
            cout << "auc:" << auc << endl;
            cout << "AUROC:" << AUROC << endl;
            cout << "AUPRC:" << AUPRC << endl;
            cout << "rankScore:" << rankScore << endl;
        }
    };

} // namespace baseStruct


#endif