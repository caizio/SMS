#include "linkPrediction.h"
using namespace std;

LinkPrediction::LinkPrediction(){

};
LinkPrediction::LinkPrediction(UDG &G){
    // _similar.insert({"JC",this->_JC});
    // _similar["CN"] = _CN;
    // _similar["JC"] = _JC;

}


std::vector<std::vector<double>> LinkPrediction::mmain(const std::string &name,UDG &G,double parameter,double parameter2,double parameter3,string dataset_name){
    int n = G.mg.vertex_num;
    std::vector<std::vector<double>> similar(n,std::vector<double>(n,0));
    if(name == "CN"){
        return CN(G);
    }else if(name == "PA"){
        return PA(G);
    }else if(name == "RA"){
        return RA(G);
    }else if(name == "JC"){
        return JC(G);
    }else if(name == "A3"){
        return A3(G);
    }else if(name == "L3"){
        return L3(G);
    }else if(name == "CH2_L3"){
        return CH2_L3(G);
    }else if(name == "SIM"){
        return SIM(G);
    }else if(name == "maxSIM"){
        return maxSIM(G);
    }else if(name == "SMS"){
        if(parameter == 0){
            return SMS(G,"CN",parameter2);
        }else if(parameter == 1){
            return SMS(G,"RA",parameter2);
        }else if(parameter == 2){
            return SMS(G,"JC",parameter2);
        }else if(parameter == 3){
            return SMS(G,"AA",parameter2);
        }else if(parameter == 4){
            return SMS(G,"Sorenson",parameter2);
        }else if(parameter == 5){
            return SMS(G,"HPI",parameter2);
        }else if(parameter == 6){
            return SMS(G,"HDI",parameter2);
        }else if(parameter == 7){
            return SMS(G,"LHNI",parameter2);
        }else if(parameter == 8){
            return SMS(G,"Salton",parameter2);
        }else if(parameter == 9){
            return SMS_KNCN(G,parameter2,parameter3);
        }
        cout << "waring: parameter is wrong" << endl;
        return SMS(G);
    }else if(name == "maxSMS"){
        if(parameter == 0){
            return maxSMS(G,"CN");
        }else if(parameter == 1){
            return maxSMS(G,"RA");
        }else if(parameter == 2){
            return maxSMS(G,"JC");
        }else if(parameter == 3){
            return maxSMS(G,"AA");
        }else if(parameter == 4){
            return maxSMS(G,"Sorenson");
        }else if(parameter == 5){
            return maxSMS(G,"HPI");
        }else if(parameter == 6){
            return maxSMS(G,"HDI");
        }else if(parameter == 7){
            return maxSMS(G,"LHNI");
        }else if(parameter == 8){
            return maxSMS(G,"Salton");
        }else if(parameter == 9){
            return maxSMS(G,"CN",true);
        };
        cout << "参数错数，已计算maxSMS_KNCN" << endl;
        return maxSMS(G,"CN",true);
    }else if(name == "LP"){
        return LP(G);
    }else if(name == "CCLP"){
        return CCLP(G);
    }
    cout << "LinkPrediction::mmain::方法参数:" << name << " 错误,已默认计算CN" << endl;
    return CN(G);
}

// 根据不同相似度的链路预测算法
vector<vector<double>> LinkPrediction::Sim(UDG &G, std::vector<std::vector<double>> &similar){
    int n = G.mg.vertex_num;
    vector<vector<double>> res(n,vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                double si = 0;
                for(auto &ii:G.mg.allNei[i]){
                    si += similar[ii][j];
                }
                for(auto &jj:G.mg.allNei[j]){
                    si += similar[jj][i];
                }
                res[i][j] = si;
                res[j][i] = si;
            }
        }
    }
    return res;
}

// SMS
vector<vector<double>> LinkPrediction::SMS(UDG &G,const string &rw,double degree){
    int n = G.mg.vertex_num;
    vector<vector<double>> similar(n,vector<double>(n,0));
    vector<vector<double>> similarMatrix;
    if(rw == "JC"){
        similarMatrix = _JC(G);
    }else if(rw == "CN"){
        similarMatrix = _CN(G);
    }else if(rw == "RA"){
        similarMatrix = _RA(G);
    }else if(rw == "AA"){
        similarMatrix = _AA(G);
    }else if(rw == "Sorenson"){
        similarMatrix = _Sorenson(G);
    }else if(rw == "HPI"){
        similarMatrix = _HPI(G);
    }else if(rw == "HDI"){
        similarMatrix = _HDI(G);
    }else if(rw == "LHNI"){
        similarMatrix = _LHNI(G);
    }else if(rw == "Salton"){
        similarMatrix = _Salton(G);
    }
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                similar[i][j] = _SMS(G,similarMatrix,i,j,degree);
                similar[j][i] = similar[i][j];
            }
        }
    }
    return similar;    
};

vector<vector<double>> LinkPrediction::SMS(UDG &G,std::vector<std::vector<double>> &similarMatrix){
    int n = G.mg.vertex_num;
    vector<vector<double>> similar(n,vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                similar[i][j] = _SMS(G,similarMatrix,i,j,0);
                similar[j][i] = similar[i][j];
            }
        }
    }
    return similar; 
}


vector<vector<double>> LinkPrediction::SMS_KNCN(UDG &G,double degree,double parameter){
    int n = G.mg.vertex_num;
    vector<vector<double>> similar(n,vector<double>(n,0));
    auto similarMatrix = _CN(G);
    for(int x = 0; x < n; x++){
        for(int y = x + 1; y < n; y++){
            if(G.mg.edge[x][y] == 0){
                double res = 0;
                for(auto &xx:G.mg.allNei[x]){
                    for(auto &yy:G.mg.allNei[y]){
                        if(G.mg.edge[xx][yy] != 0 && xx != yy){
                            double sim = similarMatrix[xx][y] * similarMatrix[yy][x] / (G.gd.degree[xx] * G.gd.degree[yy]);
                            double d = 1;
                            if(degree == 1){
                                d = 1.0 / (pow(G.gd.degree[xx],parameter) * pow(G.gd.degree[yy],parameter));                  
                            }else if(degree == 2){
                                d = 1.0 / (G.gd.cn_number[xx][y] * G.gd.cn_number[yy][x]);
                            }
                            res = res + d * sim;
                        }
                    }
                } 
                similar[x][y] = res;
                similar[y][x] = res;
            }
        }
    }
    return similar;       
}

vector<vector<double>> LinkPrediction::maxSMS(UDG &G,const string &rw,bool degree){
    int n = G.mg.vertex_num;
    vector<vector<double>> similar(n,vector<double>(n,0));
    vector<vector<double>> similarMatrix;
    if(rw == "JC"){
        similarMatrix = _JC(G);
    }else if(rw == "CN"){
        similarMatrix = _CN(G);
    }else if(rw == "RA"){
        similarMatrix = _RA(G);
    }else if(rw == "AA"){
        similarMatrix = _AA(G);
    }else if(rw == "Sorenson"){
        similarMatrix = _Sorenson(G);
    }else if(rw == "HPI"){
        similarMatrix = _HPI(G);
    }else if(rw == "HDI"){
        similarMatrix = _HDI(G);
    }else if(rw == "LHNI"){
        similarMatrix = _LHNI(G);
    }else if(rw == "Salton"){
        similarMatrix = _Salton(G);
    }
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                double d = 0;
                for(auto &xx:G.mg.allNei[i]){
                    for(auto &yy:G.mg.allNei[j]){
                        if(G.mg.edge[xx][yy] != 0 && xx != yy){
                            double deg = 1;
                            if(degree == true){
                                deg = 1.0 / (G.gd.degree[xx] * G.gd.degree[yy]);
                            }
                            d = max(deg * similarMatrix[xx][j] * similarMatrix[yy][i],d);
                        }
                    }
                }
                similar[i][j] = d;
                similar[j][i] = d;
            }
        }
    }
    return similar;      
}

vector<vector<double>> LinkPrediction::maxSMS(UDG &G,const std::vector<std::vector<double>> &similarMatrix){
    int n = G.mg.vertex_num;
    bool degree = true;
    vector<vector<double>> similar(n,vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                double d = 0;
                for(auto &xx:G.mg.allNei[i]){
                    for(auto &yy:G.mg.allNei[j]){
                        if(G.mg.edge[xx][yy] != 0 && xx != yy && similarMatrix[i][yy] != 0 && similarMatrix[j][xx] != 0){
                            double deg = 1;
                            // if(degree == true){
                            //     deg = 1.0 / (G.gd.degree[xx] * G.gd.degree[yy]);
                            // }
                            d = max(deg * similarMatrix[xx][j] * similarMatrix[yy][i],d);
                        }
                    }
                }
                similar[i][j] = d;
                similar[j][i] = d;
            }
        }
    }
    return similar;     
}

// 计算节点xy之间三阶路径的数量
double LinkPrediction::_A3(UDG &G,int x,int y){
    if(G.mg.allNei[x].size() == 0 || G.mg.allNei[y].size() == 0) return 0;
    int n = G.mg.vertex_num;
    double res = 0;
    for(auto &xx:G.mg.allNei[x]){
        for(auto &yy:G.mg.allNei[y]){
            if(G.mg.edge[xx][yy] != 0 && xx != yy){
                res += 1.0;
            }
        }
    }
    return res;
}

// 计算节点xy之间三阶路径的数量,且除以中间节点的度的0.5次方（根号）
// sum 1.0 / pow((degree[xx] * degree[yy]),0.5)
double LinkPrediction::_L3(UDG &G,int x,int y){
    if(G.mg.allNei[x].size() == 0 || G.mg.allNei[y].size() == 0) return 0;
    int n = G.mg.vertex_num;
    double res = 0;
    for(auto &xx:G.mg.allNei[x]){
        for(auto &yy:G.mg.allNei[y]){
            if(G.mg.edge[xx][yy] != 0 && xx != yy){
                res += 1.0 / pow((G.gd.degree[xx] * G.gd.degree[yy]),0.5);
            }
        }
    }
    return res;
}

// 计算x和y之间的L3社区数量
double LinkPrediction::_C3(UDG &G,int x,int y,vector<vector<double>> &similarMatrix){
    if(G.mg.allNei[x].size() == 0 || G.mg.allNei[y].size() == 0) return 0;
    int n = G.mg.edge_num;
    vector<set<int>> C;
    vector<double> similars;
    for(auto &xx:G.mg.allNei[x]){
        for(auto &yy:G.mg.allNei[y]){
            if(G.mg.edge[xx][yy] != 0 && xx != yy){
                // 寻找C中是否有重复的节点
                int i = 0;
                for(; i < C.size(); i++){
                    if(C[i].count(xx) == 1 || C[i].count(yy) == 1){
                        C[i].insert(xx);
                        C[i].insert(yy);
                        double similar = similarMatrix[x][yy] * similarMatrix[xx][y] / G.gd.degree[xx] / G.gd.degree[yy];
                        similars[i] = max(similar,similars[i]);
                        break;
                    }
                }
                // 说明当前的社区中不存在这条边
                if(i == C.size()){
                    set<int> temp;
                    temp.insert(xx);
                    temp.insert(yy);
                    C.push_back(temp);
                    double similar = similarMatrix[x][yy] * similarMatrix[xx][y] / G.gd.degree[xx] / G.gd.degree[yy];
                    similars.push_back(similar);
                }
            }
        }
    }
    double result = 0;
    int i = 0;
    for(auto &s:similars){
        if(C[i].size() > 2){
            result += s;
        }
        i++;
    }
    return result;
}

// 计算xy之间三阶路径的CH2_L3
double LinkPrediction::_CH2_L3(UDG &G,int x,int y){
    if(G.mg.allNei[x].size() == 0 || G.mg.allNei[y].size() == 0) return 0;
    int n = G.mg.vertex_num;
    double res = 0;
    set<int> lcl;
    for(auto &xx:G.mg.allNei[x]){
        for(auto &yy:G.mg.allNei[y]){
            if(G.mg.edge[xx][yy] != 0 && xx != yy){
                lcl.insert(xx);
                lcl.insert(yy);                
            }
        }
    }
    for(auto &xx:G.mg.allNei[x]){
        for(auto &yy:G.mg.allNei[y]){
            if(G.mg.edge[xx][yy] != 0 && xx != yy){
                double up1 = 1;
                double up2 = 1;
                for(auto &p:lcl){
                    if(G.mg.edge[xx][p] != 0){
                        up1++;
                    }
                    if(G.mg.edge[yy][p] != 0){
                        up2++;
                    }
                }
                double down1 = G.gd.degree[xx] + 2 - up1;
                double down2 = G.gd.degree[yy] + 2 - up2;
                double si = pow(up1 * up2 / down1 / down2,0.5);
                res += si;
            }
        }
    }
    return res;    
}

double LinkPrediction::_SMS(UDG &G,vector<vector<double>> &similarMatrix,int x,int y,double degree){
    int n = G.mg.vertex_num;
    double res = 0;
    for(auto &xx:G.mg.allNei[x]){
        for(auto &yy:G.mg.allNei[y]){
            if(G.mg.edge[xx][yy] != 0 && xx != yy){
                double d = 1;
                if(degree == 1){
                    d = 1.0 / (G.gd.degree[xx] * G.gd.degree[yy]);                  
                }else if(degree == 2){
                    d = 1.0 / (G.gd.cn_number[xx][y] * G.gd.cn_number[yy][x]);
                }else if(degree == 3){
                    d = 1.0 / pow(G.gd.degree[xx] * G.gd.degree[yy],0.5);
                }
                res = res + d * (similarMatrix[xx][y] * similarMatrix[yy][x]);
            }
        }
    }
    return res;      
}

vector<vector<double>> LinkPrediction::_CN(UDG &G){
    return Matrix_my::intToDouble(G.gd.cn_number);
}
vector<vector<double>> LinkPrediction::_AA(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,0));
    auto degree = G.acc_alldegree();
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            auto v = G.gd.cn[i][j];
            int v_sum = v.size();
            double similar = 0;
            for(int z = 0; z < v_sum; z++){
                similar += 1.0 / log(degree[v[z]]);
            }
            similarMatrix[i][j] = similar;
            similarMatrix[j][i] = similar;
        }
    }
    return similarMatrix;    
}


vector<vector<double>> LinkPrediction::_RA(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,0));
    auto degree = G.acc_alldegree();
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            auto v = G.gd.cn[i][j];
            // auto v = G.commonNeighbor_xy(i,j);
            int v_sum = v.size();
            double similar = 0;
            for(int z = 0; z < v_sum; z++){
                similar += 1.0 / degree[v[z]];
            }
            similarMatrix[i][j] = similar;
            similarMatrix[j][i] = similar;
        }
    }
    return similarMatrix;    
}
vector<vector<double>> LinkPrediction::_JC(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            int similar1 = G.gd.cn_number[i][j];
            int similar2 = G.commonUnion_xy(i,j).size();
            if(similar2 == 0) continue;
            double similar = similar1 / double(similar2);
            similarMatrix[i][j] = similar;
            similarMatrix[j][i] = similar;                
        }
    }
    return similarMatrix;   
}

vector<vector<pair<int,double>>> LinkPrediction::_JC_(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<pair<int,double>>> result(n);
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            int similar1 = G.gd.cn_number[i][j];
            int similar2 = G.commonUnion_xy(i,j).size();
            if(similar2 == 0) continue;
            double s = similar1 / double(similar2);
            if(s != 0){
                pair<int,double> sij = make_pair(j,s);
                pair<int,double> sji = make_pair(i,s);
                result[i].push_back(sij);
                result[j].push_back(sji);
            }             
        }
    }
    return result;  
}

vector<vector<double>> LinkPrediction::_Sorenson(UDG &G){
    int n = G.mg.vertex_num;
    auto degree = G.acc_alldegree();
    vector<vector<double>> similarMatrix(n,vector<double>(n,-1));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                int commonNeighboor = G.gd.cn_number[i][j];
                int ki = degree[i];
                int kj = degree[j];
                double similar = 0;
                if(ki + kj != 0){
                    similar = 2.0 * commonNeighboor / (ki + kj);
                }
                similarMatrix[i][j] = similar;
                similarMatrix[j][i] = similar;
            }
        }
    }
    return similarMatrix;    
}

vector<vector<double>> LinkPrediction::_HPI(UDG &G){
    int n = G.mg.vertex_num;
    auto degree = G.acc_alldegree();
    vector<vector<double>> similarMatrix(n,vector<double>(n,-1));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                int commonNeighboor = G.gd.cn_number[i][j];
                int ki = degree[i];
                int kj = degree[j];
                int mink = min(kj,kj);
                double similar = 0;
                if(mink != 0){
                    similar = commonNeighboor / double(mink);                    
                }
                similarMatrix[i][j] = similar;
                similarMatrix[j][i] = similar;
            }
        }
    }
    return similarMatrix;
}

// HDI
vector<vector<double>> LinkPrediction::_HDI(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,-1));
    auto degree = G.acc_alldegree();
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                int commonNeighboor = G.gd.cn_number[i][j];
                int ki = degree[i];
                int kj = degree[j];
                int maxk = max(kj,kj);
                double similar = 0;
                if(maxk != 0){
                    similar = commonNeighboor / double(maxk);
                }
                similarMatrix[i][j] = similar;
                similarMatrix[j][i] = similar;
            }
        }
    }
    return similarMatrix;
}
// LHNI
vector<vector<double>> LinkPrediction::_LHNI(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,-1));
    auto degree = G.acc_alldegree();
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                int commonNeighboor = G.gd.cn_number[i][j];
                int ki = degree[i];
                int kj = degree[j];
                double similar = 0;
                if(ki * kj != 0){
                    similar = commonNeighboor / double(ki * kj);
                }
                similarMatrix[i][j] = similar;
                similarMatrix[j][i] = similar;
            }
        }
    }
    return similarMatrix;
}

vector<vector<double>> LinkPrediction::_Salton(UDG &G){
    int n = G.mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,-1));
    auto degree = G.acc_alldegree();
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] == 0){
                int commonNeighboor = G.gd.cn_number[i][j];
                int ki = degree[i];
                int kj = degree[j];
                double similar = 0;
                if(ki * kj != 0){
                    similar = commonNeighboor / pow(ki * kj,0.5);
                }
                similarMatrix[i][j] = similar;
                similarMatrix[j][i] = similar;
            }
        }
    }
    return similarMatrix;
}