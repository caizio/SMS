#ifndef _Graph_
#define _Graph_
#include <iostream>
#include <vector>
#include <numeric>
#include <unordered_set>
#include <algorithm>
#include <set>
#include <fstream>
#include <queue>
#include <cmath>
#include <random>
#include <fstream>
#include "myAlgorithm.h"
#include "matrix.h"

using namespace std;

// 图的邻接矩阵
struct MGraph{
    vector<vector<int>>edge;            //邻接矩阵
    vector<vector<int>>allNei;          //储存某节点的邻居(该储存方式方便了寻找邻居)
    vector<set<int>>allNei2;
    int vertex_num;                     //顶点数量
    int edge_num;                       //边数量
};
// 储存计算过程的数据，避免重复计算
// sign = 1表示该数据已经计算过，不用重复计算
struct GraphData{
    vector<int> degree;
    bool sign_degree;
    vector<vector<vector<int>>> cn;
    bool sign_cn;
    vector<vector<int>> cn_number;   
    GraphData();
};
// 相似节点及得分
struct nei_jc{
    int node;
    double similar;
};

struct sms_connection{
    // 左边是开区间，右边闭区间
    double range_left;
    double range_right;
    double connect;
    double non_connect;
};

struct connection_and_non{
    int number;
    double connect;
    double non_connect;    
};

// 无向图
class UDG{
public:
    MGraph mg;
    GraphData gd;
    UDG();
    UDG(const vector<vector<int>> &data);
    vector<int> commonNeighbor_xy(int x,int y);
    vector<int> commonUnion_xy(int x,int y);
    vector<vector<int>> commonNeighbors();
    int acc_degree(int x); 
    vector<int> acc_alldegree();
    double acc_average_degree();
    double acc_clustering(int x);
    vector<double> acc_allClustering();
    double acc_closure(int x);
    double acc_average_clustering(bool is_all = false);
    void reset_gd();
    vector<sms_connection> acc_sms_connection_JC(int divide = 10);
    vector<sms_connection> acc_sms_connection_DCN(int divide = 10);
    map<double,pair<int,int>> acc_sms_connection_DCN2();
    map<double,pair<int,int>> acc_sms_connection_CN();
    int acc_distance(int x,int y);
    int acc_distance2(int x,int y,int L = 2);
    int acc_distance3(int x,int y);
    int acc_distance4(int x,int y);
    set<pair<int,int>> find_k_hop_subgraph(int i, int j, int k = 2);
    double acc_average_cn_number();
    double acc_degree_proportion_less_int(int k);
    void sample_neg(int number);
    void acc_average_CN_DCN();
    void acc_max_CN_DCN();
    void acc_xy_quadrangle_graph_edges_density();
    vector<vector<double>> JC();
    int acc_similar_drived_edges_number();
private:
    vector<vector<nei_jc>> _CN();
    vector<vector<nei_jc>> _JC();
    
};

#endif