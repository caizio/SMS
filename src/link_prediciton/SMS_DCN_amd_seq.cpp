#include "linkPrediction.h"

using namespace std;

// 基于序列相似性的SMS算法
vector<vector<double>> LinkPrediction::SMS_sequence(UDG &G,map<pair<int,int>,double> &seq_similar, map<int,string> &seq,ofstream &file){
    int n = G.mg.vertex_num;
    vector<vector<double>> complement(n,vector<double>(n,0));  
    for(int x = 0; x < n; x++){
        // if(x % 10 == 0) cout << x << endl;
        for(int y = x + 1; y < n; y++){
            if(G.mg.edge[x][y] == 0){
                double c = 0;
                for(auto &xx:G.mg.allNei[x]){
                    for(auto &yy:G.mg.allNei[y]){
                        if(G.mg.edge[xx][yy] != 0 && xx != yy){ 
                            auto similar_x_yy = acc_similarBySequence2(x,yy,seq_similar,seq,file);
                            auto similar_y_xx = acc_similarBySequence2(y,xx,seq_similar,seq,file);
                            c += similar_x_yy * similar_y_xx;
                        }
                    }
                }
                complement[x][y] = c;
                complement[y][x] = c;
            }
        }
    }
    return complement;
}

vector<vector<double>> LinkPrediction::DSMS_sequence(UDG &G,map<pair<int,int>,double> &seq_similar, map<int,string> &seq,ofstream &file){
    int n = G.mg.vertex_num;
    vector<vector<double>> complement(n,vector<double>(n,0));
    double d = 0;
    for(int x = 0; x < n; x++){
        // if(x % 10 == 0) cout << x << endl;
        for(int y = x + 1; y < n; y++){
            if(G.mg.edge[x][y] == 0){
                double c = 0;
                for(auto &xx:G.mg.allNei[x]){
                    for(auto &yy:G.mg.allNei[y]){
                        if(G.mg.edge[xx][yy] != 0 && xx != yy){ 
                            d = 1.0 / (G.gd.degree[xx] * G.gd.degree[yy]);  
                            auto similar_x_yy = acc_similarBySequence2(x,yy,seq_similar,seq,file);
                            auto similar_y_xx = acc_similarBySequence2(y,xx,seq_similar,seq,file);
                            c += d * similar_x_yy * similar_y_xx;
                        }
                    }
                }
                complement[x][y] = c;
                complement[y][x] = c;
            }
        }
    }
    return complement;
}

vector<vector<double>> LinkPrediction::maxSMS_sequence(UDG &G,map<pair<int,int>,double> &seq_similar, map<int,string> &seq,ofstream &file){
    int n = G.mg.vertex_num;
    vector<vector<double>> complement(n,vector<double>(n,0));
    for(int x = 0; x < n; x++){
        // if(x % 10 == 0) cout << x << endl;
        for(int y = x + 1; y < n; y++){
            if(G.mg.edge[x][y] == 0){
                double c = 0;
                for(auto &xx:G.mg.allNei[x]){
                    for(auto &yy:G.mg.allNei[y]){
                        if(G.mg.edge[xx][yy] != 0 && xx != yy){  
                            auto similar_x_yy = acc_similarBySequence2(x,yy,seq_similar,seq,file);
                            auto similar_y_xx = acc_similarBySequence2(y,xx,seq_similar,seq,file);
                            c = max(similar_x_yy * similar_y_xx,c);
                        }
                    }
                }
                complement[x][y] = c;
                complement[y][x] = c;
            }
        }
    }
    return complement;
}


// S = DCN * (SeqS ^ lam)
vector<vector<double>> LinkPrediction::SMS_mix(UDG &G,map<pair<int,int>,double> &seq_similar, map<int,string> &seq,ofstream &file,double lam){
    int n = G.mg.vertex_num;
    vector<vector<double>> complement(n,vector<double>(n,0));
    auto similarMatrix = _CN(G);  
    // double d = 1;
    for(int x = 0; x < n; x++){
        // if(x % 10 == 0) cout << x << endl;
        for(int y = x + 1; y < n; y++){
            if(G.mg.edge[x][y] == 0){
                double c = 0;
                for(auto &xx:G.mg.allNei[x]){
                    for(auto &yy:G.mg.allNei[y]){
                        if(G.mg.edge[xx][yy] != 0 && xx != yy){ 
                            auto similar_x_yy = acc_similarBySequence2(x,yy,seq_similar,seq,file);
                            auto similar_y_xx = acc_similarBySequence2(y,xx,seq_similar,seq,file);
                            double s_x_yy = similarMatrix[x][yy] / (G.gd.degree[yy]) * pow(similar_x_yy,lam);
                            double s_y_xx = similarMatrix[y][xx] / (G.gd.degree[xx]) * pow(similar_y_xx,lam);
                            c += s_x_yy * s_y_xx;
                        }
                    }
                }
                complement[x][y] = c;
                complement[y][x] = c;
            }
        }
    }
    return complement;
}

// S = DCN * (SeqS ^ lam) + maxSMS
vector<vector<double>> LinkPrediction::maxSMS_mix(UDG &G,map<pair<int,int>,double> &seq_similar, map<int,string> &seq,ofstream &file,double lam){
    int n = G.mg.vertex_num;
    vector<vector<double>> complement(n,vector<double>(n,0));
    auto similarMatrix = _CN(G);  
    for(int x = 0; x < n; x++){
        // if(x % 10 == 0) cout << x << endl;
        for(int y = x + 1; y < n; y++){
            if(G.mg.edge[x][y] == 0){
                double c = 0;
                for(auto &xx:G.mg.allNei[x]){
                    for(auto &yy:G.mg.allNei[y]){
                        if(G.mg.edge[xx][yy] != 0 && xx != yy){ 
                            auto similar_x_yy = acc_similarBySequence2(x,yy,seq_similar,seq,file);
                            auto similar_y_xx = acc_similarBySequence2(y,xx,seq_similar,seq,file);
                            double s_x_yy = similarMatrix[x][yy] / (G.gd.degree[yy]) * pow(similar_x_yy,lam);
                            double s_y_xx = similarMatrix[y][xx] / (G.gd.degree[xx]) * pow(similar_y_xx,lam);
                            c = max(s_x_yy * s_y_xx,c);
                        }
                    }
                }
                complement[x][y] = c;
                complement[y][x] = c;
            }
        }
    }
    return complement;
}

// S = DCN / k * (SeqS ^ lam)
vector<vector<double>> LinkPrediction::DSMS_mix(UDG &G,map<pair<int,int>,double> &seq_similar, map<int,string> &seq,ofstream &file,double lam){
    int n = G.mg.vertex_num;
    vector<vector<double>> complement(n,vector<double>(n,0));
    auto similarMatrix = _CN(G);  
    for(int x = 0; x < n; x++){
        // if(x % 10 == 0) cout << x << endl;
        for(int y = x + 1; y < n; y++){
            if(G.mg.edge[x][y] == 0){
                double c = 0;
                for(auto &xx:G.mg.allNei[x]){
                    for(auto &yy:G.mg.allNei[y]){
                        if(G.mg.edge[xx][yy] != 0 && xx != yy){ 
                            auto similar_x_yy = acc_similarBySequence2(x,yy,seq_similar,seq,file);
                            auto similar_y_xx = acc_similarBySequence2(y,xx,seq_similar,seq,file);
                            double s_x_yy = similarMatrix[x][yy] / (G.gd.degree[yy] * G.gd.degree[yy]) * pow(similar_x_yy,lam);
                            double s_y_xx = similarMatrix[y][xx] / (G.gd.degree[xx] * G.gd.degree[xx]) * pow(similar_y_xx,lam);
                            c += s_x_yy * s_y_xx;
                        }
                    }
                }
                complement[x][y] = c;
                complement[y][x] = c;
            }
        }
    }
    return complement;
}

// 计算两个序列的相似性
double LinkPrediction::acc_similarBySequence(int a,int b, map<pair<int,int>,double> &similar,map<int,string>&seq){
    int x = min(a,b);
    int y = max(a,b);
    // 如果序列间的相似性已经计算，则直接寻找
    auto similar_end = similar.end();
    if(similar.find({x,y}) != similar_end){
        return similar[{x,y}];
    }
    //计算两个序列的相似性 
    bio::SequenceCompare sc;
    string s1 = seq[x];
    string s2 = seq[y];
    auto p = sc.blosub(seq[x],seq[y]);
    auto similar_x_y = sc.find_path(s1,s2,p,s1.size(),s2.size()).similar;   
    similar[{x,y}] = similar_x_y;
    return similar_x_y;
}

double LinkPrediction::acc_similarBySequence2(int a,int b,map<pair<int,int>,double> &similar,map<int,string>&seq,ofstream &file){
    int x = min(a,b);
    int y = max(a,b);
    // 如果序列间的相似性已经计算，则直接寻找
    auto similar_end = similar.end();
    if(similar.find({x,y}) != similar_end){
        return similar[{x,y}];
    }
    //计算两个序列的相似性 
    bio::SequenceCompare sc;
    string s1 = seq[x];
    string s2 = seq[y];
    auto p = sc.blosub(seq[x],seq[y]);
    auto similar_x_y = sc.find_path(s1,s2,p,s1.size(),s2.size()).similar;   
    similar[{x,y}] = similar_x_y;
    file << to_string(x) << "," << to_string(y) << "," << to_string(similar_x_y) << endl;
    return similar_x_y;
}
