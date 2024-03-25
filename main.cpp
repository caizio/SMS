#include <iostream>
#include <vector>
#include <string>

#include "main_function.h"

using namespace std;

/*
    @brief 通过choose的值运行不同算法。Run different algorithms based on the value of choose
    @param [methods] The name of the algorithm being run
    @param [parameter] Parameters of different algorithms
    @param [parameter2] Parameters2 of different algorithms
*/
int main(int argc,char *argv[]){
    int choose = 0;
    DataDeal dd;
    string dir = "";
    string stime = getDate();
    string outFile = "../result/" + stime;
    if(argc >= 2){
        choose = stoi(argv[1]);
    }
    std::vector<std::string> methods;
    std::vector<double> parameter;
    std::vector<double> parameter2;
    if(choose == 0){
        // for(int i = 0; i < 1; i++){
        //     dd.cross_validation(LNBL3DataSet[i],10);
        // }
        // 交叉验证划分数据集
        // Cross validation partitioning dataset
        for(auto &p:LNBL3DataSet){
            dd.cross_validation(p,10);
        }
        // // 随机划分
        // for(auto &p:LNBL3DataSet){
        //     dd.random_divide(p,0.5,10);
        // }
    }else if(choose == 1){
        // 不同算法的对比实验
        // Comparative experiments of different algorithms
        methods = {"CN","A3","L3","CH2_L3","SIM","maxSIM","SMS","SMS","maxSMS"};
        parameter = {-1,  -1,  -1,      -1,   -1,      -1,   9,    9,       9};
        parameter2 ={-1,  -1,  -1,      -1,   -1,      -1,   0,    1,       1};
        for(int i = 0; i < LNBL3DataSet.size();i++){
            dir = LNBL3DataSet[i];
            m_LinkPrediction(dir,methods,parameter,outFile,10,parameter2);
        } 
    }else if(choose == 2){
        // 相似性指标的选择
        // parameter:0-9表示10中相似性指标
        // Selection of similarity indicators
        // Parameters: 0-9 represent 10 similarity indicators
        methods = {"SMS","SMS","SMS","SMS","SMS","SMS","SMS","SMS","SMS","SMS"};
        parameter =  {0,1,2,3,4,5,6,7,8,9};
        
        // parameter2: 0表示不除以度，1表示除以度，2表示除以CN
        // Parameter2: 0 represents not dividing by degrees, 1 represents dividing by degrees, and 2 represents dividing by CN
        parameter2 = {0,0,0,0,0,0,0,0,0,0};
        // parameter2 = {1,1,1,1,1,1,1,1,1,1};
        // parameter2 = {2,2,2,2,2,2,2,2,2,2};

        // methods = {"maxSMS","maxSMS","maxSMS","maxSMS","maxSMS","maxSMS","maxSMS","maxSMS","maxSMS","maxSMS"};
        // parameter =  {0,1,2,3,4,5,6,7,8,9};
        // parameter2 = {0,0,0,0,0,0,0,0,0,0}; 
        for(int i = 0; i < LNBL3DataSet.size(); i++){
            dir = LNBL3DataSet[i];
            m_LinkPrediction(dir,methods,parameter,outFile,10,parameter2);
        } 
    }else if(choose == 3){
        // 互补传递的验证 
        // Verification of complementary transfer
        // m_3();
        // m_3_2();
        m_3maxDCN_l();
    }else if(choose == 4){
        // 寻找不同数据集划分下，DSMS的最优参数
        // Finding the optimal parameters for DSMS under different dataset partitions
        parameter = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        for(int i = 0; i < 1;i++){
            dir = LNBL3DataSet[i];
            m_find_DSMS_parameter(dir,parameter,outFile,4,90);
        }
    }else if(choose == 5){
        // SMS+SeqS,DSMS+SeqS,maxSMS+SeqS
        // SMS algorithm based solely on sequences
        int dataset_repeat_count = 10;
        for(int i = 0; i < LNBL3DataSet.size(); i++){
            m_LinkPrediction_seq(LNBL3DataSet[i],outFile,dataset_repeat_count);
        } 
    }else if(choose == 6){
        // The SMS algorithm considering sequence similarity and topological similarity
        int dataset_repeat_count = 1;
        vector<vector<double>> lam = {
            {2}
        };
        for(int i = 0; i < LNBL3DataSet.size(); i++){
            m_LinkPrediction_SMS_mix(LNBL3DataSet[i],outFile,dataset_repeat_count,lam[0]);
        } 
    }else if(choose == 7){
        // The DSMS algorithm considering sequence similarity and topological similarity
        int dataset_repeat_count = 10;
        vector<vector<double>> lam = {
            {0.5}
        };
        for(int i = 0; i < LNBL3DataSet.size(); i++){
            m_LinkPrediction_DSMS_mix(LNBL3DataSet[i],outFile,dataset_repeat_count,lam[0]);
        } 
    }else if(choose == 8){
        // The maxSMS algorithm considering sequence similarity and topological similarity
        int dataset_repeat_count = 10;
        vector<vector<double>> lam = {
            {0.2}
        };
        for(int i = 0; i < LNBL3DataSet.size(); i++){
            m_LinkPrediction_maxSMS_mix(LNBL3DataSet[i],outFile,dataset_repeat_count,lam[0]);
        } 
    }else if(choose == 9){
        // Calculate the distribution of edge density
        // 计算相连边和不相连边的共同邻居数量和边密度
        for(int i = 0; i < 6; i++){
            auto dataset= LNBL3DataSet[i];
            auto adj = dd.readGraph(dataset);
            UDG G(adj);
            string datasetName = dd.getDataSetName(dataset);
            cout << datasetName << ",";
            auto res = G.acc_xy_quadrangle_graph_cn_and_edges_density();
            std::string outfile = "../temp/" + datasetName + "_cn_and_edges_density_and_ave_interdimate_nodes_degree.txt";
            ofstream f(outfile);
            for(auto &p:res){
                f << get<0>(p) << "," << get<1>(p) << "," << get<2>(p) << ","  << endl;
            }
            f.close();
        }
    }
    return 0;
}
