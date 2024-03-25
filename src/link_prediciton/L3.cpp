#include "linkPrediction.h"

using namespace std;


// A3算法，计算两节点三阶路径的数量
vector<vector<double>> LinkPrediction::A3(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> res(n,vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                res[i][j] = _A3(G,i,j);
                res[j][i] = res[i][j];
            }
        }
    }
    return res;
}

vector<vector<double>> LinkPrediction::LP(UDG &G,double lam){
    int n = G.mg.vertex_num;
    vector<vector<double>> res(n,vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                res[i][j] = G.gd.cn_number[i][j] + lam * _A3(G,i,j);
                res[j][i] = res[i][j];
            }
        }
    }
    return res;
}

// L3算法
vector<vector<double>> LinkPrediction::L3(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> res(n,vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                res[i][j] = _L3(G,i,j);
                res[j][i] = res[i][j];
            }
        }
    }
    return res;
};

// CH2_L3,三阶路径和LCP范式的结合
vector<vector<double>> LinkPrediction::CH2_L3(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,0));
    int similar = 0;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                similar = _CH2_L3(G,i,j);
                similarMatrix[i][j] = similar;
                similarMatrix[j][i] = similar;                
            }
        }
    }
    return similarMatrix;      
}

// SIM
vector<vector<double>> LinkPrediction::SIM(UDG &G){
    int n = G.mg.vertex_num;
    auto JC = this->_JC(G);
    // JC = this->_CN(G);
    vector<vector<double>> res(n,vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                double si = 0;
                for(auto &ii:G.mg.allNei[i]){
                    si += JC[ii][j];
                }
                for(auto &jj:G.mg.allNei[j]){
                    si += JC[jj][i];
                }
                res[i][j] = si;
                res[j][i] = si;
            }
        }
    }
    return res;
}

vector<vector<double>> LinkPrediction::maxSIM(UDG &G){
    int n = G.mg.vertex_num;
    auto JC = this->_JC(G);
    vector<vector<double>> res(n,vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                double si = 0;
                for(auto &ii:G.mg.allNei[i]){
                    si = max(JC[ii][j],si);
                }
                double di = 0;
                for(auto &jj:G.mg.allNei[j]){
                    di = max(JC[jj][i],di);
                }
                res[i][j] = si + di;
                res[j][i] = res[i][j];
            }
        }
    }
    return res;    
}