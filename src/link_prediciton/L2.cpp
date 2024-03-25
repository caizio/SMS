#include "linkPrediction.h"

using namespace std;

vector<vector<double>> LinkPrediction::PA(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,0));
    double similar = 0;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                similar = G.gd.degree[i] * G.gd.degree[j];
                similarMatrix[i][j] = similar;
                similarMatrix[j][i] = similar;                
            }
        }
    }
    return similarMatrix; 
}

// cn
vector<vector<double>> LinkPrediction::CN(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,0));
    int similar = 0;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            // 只考虑原本不存在的链路,对其计算CN值
            if(G.mg.edge[i][j] == 0){
                similar = G.gd.cn_number[i][j];
                similarMatrix[i][j] = similar;
                similarMatrix[j][i] = similar;                
            }
        }
    }
    return similarMatrix;    
}

// sum(1/kz)
vector<vector<double>> LinkPrediction::RA(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,0));
    auto degree = G.acc_alldegree();
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                auto v = G.gd.cn[i][j];
                int v_sum = v.size();
                double similar = 0;
                for(int z = 0; z < v_sum; z++){
                    similar += 1.0 / degree[v[z]];
                }
                similarMatrix[i][j] = similar;
                similarMatrix[j][i] = similar;
            }
        }
    }
    return similarMatrix;
}

// Jaccard
vector<vector<double>> LinkPrediction::JC(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                int similar1 = G.gd.cn_number[i][j];
                int similar2 = G.commonUnion_xy(i,j).size();
                if(similar2 == 0) continue;
                double similar = similar1 / double(similar2);
                similarMatrix[i][j] = similar;
                similarMatrix[j][i] = similar;                
            }
        }
    }
    return similarMatrix;    
}

// 共同邻居的聚类系数之和
vector<vector<double>> LinkPrediction::CCLP(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,0));
    auto CC = G.acc_allClustering();
    for(int x = 0; x < n; x++){
        for(int y = x + 1; y < n; y++){
            double similar_xy = 0;
            if(G.mg.edge[x][y] == 0){
                for(auto &w:G.gd.cn[x][y]){
                    similar_xy += CC[w];
                }
            }
            similarMatrix[x][y] = similar_xy;
            similarMatrix[y][x] = similar_xy;
        }
    }
    return similarMatrix;
}