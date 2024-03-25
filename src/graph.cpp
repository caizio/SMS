#include "graph.h"

using namespace std;

GraphData::GraphData(){
    sign_degree = 0;
    sign_cn = 0;
}
// 不要使用默认的构造函数
UDG::UDG(){
    cout << "不要使用默认的构造函数" << endl;
}
UDG::UDG(const vector<vector<int>> &data){
    int m = data.size(), n = data[0].size();
    if (m != n) cout << "UDG::UDG::warning：数据的尺寸不匹配" << endl;
    m = m <= n ? m : n;
    this->mg.vertex_num = m;
    this->mg.edge.resize(m);
    this->mg.allNei.resize(m);
    this->mg.allNei2.resize(m);
    for(int i = 0; i < m; i++){
        this->mg.edge[i].resize(m,0);
    }
    // 添加图的边
    int edge_count = 0;
    for (int i = 0; i < m; i++){
        for (int j = 0; j < m; j++){
            // 不考虑自边
            if (data[i][j] != 0 && i != j){
                this->mg.edge[i][j] = data[i][j];
                this->mg.allNei[i].push_back(j);
                this->mg.allNei2[i].insert(j);
                edge_count++;
            }
        }
    }
    // 因为一条边添加了两次，所以要除以2
    edge_count /= 2;
    this->mg.edge_num = edge_count;
    this->acc_alldegree();
    this->commonNeighbors();
}
//计算xy两个节点的交集（共同邻居）
vector<int> UDG::commonNeighbor_xy(int x,int y){    
    // 通过哈希函数，优化求共同邻居的算法
    unordered_set<int> set1,set2;
    for (auto &xx:mg.allNei[x]) {
        set1.insert(xx);
    }
    for (auto &yy:mg.allNei[y]) {
        set2.insert(yy);
    }
    return myAlgorithm::getIntersection(set1,set2);
}
//计算xy两个节点的并集
vector<int> UDG::commonUnion_xy(int x,int y){
    vector<int> res;
    unordered_set<int> set1;
    for(auto &xx:mg.allNei[x]){
        set1.insert(xx);
    }
    for(auto &yy:mg.allNei[y]){
        set1.insert(yy);
    }
    for(auto &pp:set1){
        res.push_back(pp);
    }
    return res;
}
// 计算任意两节点的交集数量
vector<vector<int>> UDG::commonNeighbors(){
    int n = this->mg.vertex_num;
    vector<vector<int>> res(n,vector<int>(n,0));
    if(this->gd.sign_cn == 1){
        return gd.cn_number;
    }
    this->gd.cn.resize(n,vector<vector<int>>(n));
    this->gd.cn_number.resize(n,vector<int>(n,0));
    this->gd.sign_cn = 1;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            auto cn = this->commonNeighbor_xy(i,j);
            this->gd.cn[i][j] = cn;
            this->gd.cn[j][i] = cn;
            int similar = cn.size();
            gd.cn_number[i][j] = similar;
            gd.cn_number[j][i] = similar;                
        }
    }
    return gd.cn_number;
}
//计算节点x的度
int UDG::acc_degree(int x){
    return this->mg.allNei[x].size();
}
//计算所有节点的度
vector<int> UDG::acc_alldegree(){
    // 如果sign_degree==1,说明度已经计算过了
    if(this->gd.sign_degree == 1){
        return this->gd.degree;
    }
    int n = this->mg.vertex_num;
    vector<int> degree(n,0);
    for(int i = 0; i < n; i++){
        degree[i] = this->mg.allNei[i].size();
    }
    this->gd.degree = degree;
    this->gd.sign_degree = 1;
    return degree;
}
// 计算度的平均值
double UDG::acc_average_degree(){
    return accumulate(this->gd.degree.begin(),this->gd.degree.end(),0) / double(this->mg.vertex_num);
}
// 计算节点x的聚类系数
double UDG::acc_clustering(int x){
    auto degree_x = this->gd.degree[x];
    if(degree_x <= 1) return 0;
    int l = 0;
    int n = this->mg.allNei[x].size();
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            int xNei_i = this->mg.allNei[x][i];
            int xNei_j = this->mg.allNei[x][j];
            if(this->mg.edge[xNei_i][xNei_j] != 0) l++;
        }
    }
    return double(2 * l) / (degree_x * (degree_x - 1));
}
// 计算所有节的聚类系数
vector<double> UDG::acc_allClustering(){
    int n = this->mg.vertex_num;
    vector<double> res(n,0);
    for(int i = 0; i < n; i++){
        res[i] = this->acc_clustering(i);
        // cout << res[i] << " ";
    }
    // cout << endl;
    return res;
}
// 计算节点x的闭包系数
double UDG::acc_closure(int x){
    auto degree_x = this->gd.degree[x];
    if(degree_x <= 1) return 0;
    int l = 0;
    int n = this->mg.allNei[x].size();
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            int xNei_i = this->mg.allNei[x][i];
            int xNei_j = this->mg.allNei[x][j];
            if(this->mg.edge[xNei_i][xNei_j] != 0) l++;
        }
    }
    double dj = 0;
    for(auto &p:mg.allNei[x]){
        dj += (gd.degree[p] - 1);
    }
    return 2 * double(l) / dj;
}
// 计算所有节点的平均聚类系数,注意，如果is_all == false,除以的是所有度大于2的节点数量，因为度为1的节点不能计算聚类系数
double UDG::acc_average_clustering(bool is_all){
    int n = this->mg.vertex_num;
    auto degree = this->acc_alldegree();
    double aveClustering = 0;
    for(int i = 0; i < n; i++){
        aveClustering += this->acc_clustering(i);
    }
    int count = 0;
    if(is_all == true) return aveClustering / n;
    for(int &d:degree){
        if(d >= 2) count++;
    }
    return aveClustering / count;
}

// 重新计算和cn
void UDG::reset_gd(){
    gd.sign_cn = 0;
    gd.sign_degree = 0;
    this->acc_alldegree();
    this->commonNeighbors();
}

// 按照JC相似度统计的相连比例
vector<sms_connection> UDG::acc_sms_connection_JC(int divide){
    vector<sms_connection> group(divide);
    for(int i  = 0; i < divide; i++){
        group[i].range_left = 1.0 / divide * i;
        group[i].range_right =  1.0 / divide * (i + 1);
        group[i].connect = 0;
        group[i].non_connect = 0;
    }
    auto JC = _JC();
    int n = mg.vertex_num;
    // 遍历相连的边ij
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(mg.edge[i][j] != 0){
                // 遍历与ij相似的节点xy
                for(auto &x:JC[i]){
                    for(auto &y:JC[j]){
                        if(x.node != y.node && x.node != j && y.node != i){
                            // 按照相似性区间统计相连和不相连边的数量
                            double similar = x.similar * y.similar;
                            for(int z = 0; z < divide; z++){
                                if(similar > group[z].range_left && similar <= group[z].range_right){
                                    if(mg.edge[x.node][y.node] != 0){
                                        group[z].connect++;
                                    }else{
                                        group[z].non_connect++;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return group;
}
vector<sms_connection> UDG::acc_sms_connection_DCN(int divide){
    vector<sms_connection> group(divide);
    for(int i  = 0; i < divide; i++){
        group[i].range_left = 1.0 / divide * i;
        group[i].range_right =  1.0 / divide * (i + 1);
        group[i].connect = 0;
        group[i].non_connect = 0;
    }
    auto CN = _CN();
    int n = mg.vertex_num;
    // 遍历相连的边ij
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(mg.edge[i][j] != 0){
                // 遍历与ij相似的节点xy
                for(auto &x:CN[i]){
                    for(auto &y:CN[j]){
                        if(x.node != y.node && x.node != j && y.node != i){
                            // 按照相似性区间统计相连和不相连边的数量
                            double similar = x.similar * y.similar / (gd.degree[i] * gd.degree[j]);
                            for(int z = 0; z < divide; z++){
                                if(similar > group[z].range_left && similar <= group[z].range_right){
                                    if(mg.edge[x.node][y.node] != 0){
                                        group[z].connect++;
                                    }else{
                                        group[z].non_connect++;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return group;    
}
// max_DCN的传递，按照相似度统计
map<double,pair<int,int>> UDG::acc_sms_connection_DCN2(){
    map<double,pair<int,int>> count;
    auto DCN = _CN();
    int n = mg.vertex_num;
    // 遍历相连的边ij
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(mg.edge[i][j] != 0){
                // 遍历与ij相似的节点xy
                double similar_max = 0;
                double x_node_max = 0;
                double y_node_max = 0;
                for(auto &x:DCN[i]){
                    for(auto &y:DCN[j]){
                        if(x.node != y.node && x.node != j && y.node != i){
                            // 按照相似性区间统计相连和不相连边的数量
                            double similar = x.similar * y.similar / (gd.degree[i] * gd.degree[j]);
                            if(similar_max < similar){
                                similar_max = similar;
                                x_node_max = x.node;
                                y_node_max = y.node;
                            }
                        }
                    }
                }
                if(similar_max == 0) continue;
                if(mg.edge[x_node_max][y_node_max] != 0){
                    count[similar_max].first++;
                }else{
                    count[similar_max].second++;
                }
            }
        }
    }
    return count; 
}

// 按照CN统计连边比例
map<double,pair<int,int>> UDG::acc_sms_connection_CN(){
    map<double,pair<int,int>> count;
    auto CN = _CN();
    int n = mg.vertex_num;
    // 遍历相连的边ij
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(mg.edge[i][j] != 0){
                // 遍历与ij相似的节点xy
                double similar_max = 0;
                double x_node_max = 0;
                double y_node_max = 0;
                for(auto &x:CN[i]){
                    for(auto &y:CN[j]){
                        if(x.node != y.node && x.node != j && y.node != i){
                            // 按照相似性区间统计相连和不相连边的数量
                            double similar = x.similar * y.similar;
                            if(similar_max < similar){
                                similar_max = similar;
                                x_node_max = x.node;
                                y_node_max = y.node;
                            }
                        }
                    }
                }
                if(similar_max == 0) continue;
                if(mg.edge[x_node_max][y_node_max] != 0){
                    count[similar_max].first++;
                }else{
                    count[similar_max].second++;
                }
            }
        }
    }
    return count;
}

// 计算两节点之间的最短距离
int UDG::acc_distance(int source,int dest){
    int n = mg.vertex_num;
    queue<int> q;
    unordered_set<int> visited;
    q.push(source);
    visited.insert(source);
    int depth = 0;
    while (!q.empty()) {
        int size = q.size();
        for (int i = 0; i < size; ++i) {
            int u = q.front();
            q.pop();
            if (u == dest) {
                return depth;
            }
            for(int v = 0; v < n; v++){
                if(mg.edge[u][v] != 0){
                    if (visited.find(v) == visited.end()) {
                        q.push(v);
                        visited.insert(v);
                    }
                }
            }
        }
        depth++;
    }
    return -1;
}
// 计算两节点之间的最短距离,寻找距离小于L的两个，否则返回-1
int UDG::acc_distance2(int source,int dest,int L){
    int n = mg.vertex_num;
    queue<int> q;
    unordered_set<int> visited;
    q.push(source);
    visited.insert(source);
    int depth = 0;
    while (!q.empty()) {
        int size = q.size();
        for (int i = 0; i < size; ++i) {
            int u = q.front();
            q.pop();
            if (u == dest) {
                return depth;
            }
            for(int v = 0; v < n; v++){
                if(mg.edge[u][v] != 0){
                    if (visited.find(v) == visited.end()) {
                        q.push(v);
                        visited.insert(v);
                    }
                }
            }
        }
        if(depth == L){
            return -1;
        }
        depth++;
    }
    return -1;
}

// 计算两节点之间的最短距离,是1 2 3 还是不可达
int UDG::acc_distance3(int x,int y){
    int n = mg.vertex_num;
    if(x == y){
        return 0;
    }else if(mg.edge[x][y] != 0){
        return 1;
    }else if(gd.cn_number[x][y] != 0){
        return 2;
    }else{
        for(auto &xx:mg.allNei[x]){
            for(auto &yy:mg.allNei[y]){
                if(mg.edge[xx][yy] != 0 && xx != yy){
                    return 3;
                }
            }
        }
    }
    return -1;
}

// 计算路径=2的节点对，是否存在路径=3的节点对，
int UDG::acc_distance4(int source,int dest){
    int n = mg.vertex_num;
    queue<int> q;
    unordered_set<int> visited;
    q.push(source);
    visited.insert(source);
    int depth = 0;
    while (!q.empty()) {
        int size = q.size();
        for (int i = 0; i < size; ++i) {
            int u = q.front();
            q.pop();
            if (u == dest) {
                return depth;
            }
            for(int v = 0; v < n; v++){
                if(mg.edge[u][v] != 0){
                    if (visited.find(v) == visited.end()) {
                        q.push(v);
                        visited.insert(v);
                    }
                }
            }
        }
        depth++;
    }
    return -1;
}

// 计算CN相似度，并以邻接表的形式储存
vector<vector<nei_jc>> UDG::_CN(){
    int n = mg.vertex_num;
    vector<vector<nei_jc>> result(n);
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            int similar1 = gd.cn_number[i][j];
            if(similar1 == 0) continue;
            nei_jc njc;
            njc.node = j;
            njc.similar = similar1;
            result[i].push_back(njc);  
            nei_jc njc2;
            njc2.node = i;
            njc2.similar = similar1;
            result[j].push_back(njc2);       
        }
    }
    return result; 
}
// 计算jc相似度，并以邻接表的形式储存
vector<vector<nei_jc>> UDG::_JC(){
    int n = mg.vertex_num;
    vector<vector<nei_jc>> result(n);
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            int similar1 = gd.cn_number[i][j];
            int similar2 = commonUnion_xy(i,j).size();
            if(similar2 == 0 || similar1 == 0) continue;
            double similar = similar1 / double(similar2);
            nei_jc njc;
            njc.node = j;
            njc.similar = similar;
            result[i].push_back(njc);  
            nei_jc njc2;
            njc2.node = i;
            njc2.similar = similar;
            result[j].push_back(njc2);       
        }
    }
    return result; 
}

// 计算节点i和j的xy四边形图的平均共同邻居
double UDG::_acc_xy_quadrangle_graph_i_j_ave_cn(int i, int j){
    int n = this->mg.vertex_num;
    set<int> V1;
    set<int> V2;
    for(auto &xx:this->mg.allNei[i]){
        for(auto &yy:this->mg.allNei[j]){
            if(xx == j or yy == i) continue;
            if(this->mg.edge[xx][yy] != 0 && xx != yy){
                V1.insert(xx);
                V2.insert(yy);
            }
        }
    }
    if(V1.size() == 0 || V2.size() == 0) return 0;
    int common_nei_number = 0;
    for(auto &node:V1){
        auto com = this->gd.cn_number[j][node];
        common_nei_number += com;
    }
    for(auto &node:V2){
        auto com = this->gd.cn_number[i][node];
        common_nei_number += com;
    }
    double ave = double(common_nei_number) / (V1.size() + V2.size()); 
    return ave;
}


// 计算节点i和j的xy四边形图的边密度
double UDG::_acc_xy_quadrangle_graph_i_j_edges_density(int i, int j){
    int n = this->mg.vertex_num;
    set<int> V1;
    set<int> V2;
    set<pair<int,int>> Edges;
    V1.insert(j);
    V2.insert(i);
    if(mg.edge[i][j] != 0){
        Edges.insert({min(i,j),max(i,j)});
    }
    for(auto &xx:this->mg.allNei[i]){
        for(auto &yy:this->mg.allNei[j]){
            if(xx == j or yy == i) continue;
            if(this->mg.edge[xx][yy] != 0 && xx != yy){
                V1.insert(xx);
                V2.insert(yy);
                Edges.insert({min(i,xx),max(i,xx)});
                Edges.insert({min(j,yy),max(j,yy)});
                Edges.insert({min(xx,yy),max(xx,yy)});
            }
        }
    }
    double density = double(Edges.size()) / V1.size() / V2.size();
    return density;
}

// 返回节点i和j的k跳子图，k只能等于2或者3
set<pair<int,int>> UDG::find_k_hop_subgraph(int i, int j, int k){
    set<int> nodes;
    set<pair<int,int>> edges;
    if(k == 2 && gd.cn_number[i][j] != 0){
        nodes.insert(i);
        nodes.insert(j);
        for(auto &p:gd.cn[i][j]){
            nodes.insert(p);
        }
        for(int node:nodes){
            for(auto &p:mg.allNei[node]){
                if(nodes.count(p) == 1){
                    int min_n = min(node,p);
                    int max_n = max(node,p);
                    edges.insert(make_pair(min_n,max_n));
                }
            }
        }
    }
    if(k == 3 && gd.cn_number[i][j] != 0){
        nodes.insert(i);
        nodes.insert(j);
        bool flag = 0;
        for(auto u:mg.allNei[i]){
            nodes.insert(u);
            for(auto v:mg.allNei[j]){
                nodes.insert(v);
                if(mg.edge[u][v] != 0 && u != j && v != i){
                    flag = 1;
                }
            }
        }
        if(flag == 0){
            return set<pair<int,int>>();
        }
        edges.insert(make_pair(i,j));
        for(int node:nodes){
            for(auto &p:mg.allNei[node]){
                if(nodes.count(p) == 1){
                    int min_n = min(node,p);
                    int max_n = max(node,p);
                    edges.insert(make_pair(min_n,max_n));
                }
            }
        }
    }
    return edges;
}

// 计算平均的共同邻居数量
double UDG::acc_average_cn_number(){
    int n = mg.vertex_num;
    int count = 0;
    double sum = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(gd.cn_number[i][j] != 0){
                sum += gd.cn_number[i][j];
                count++;
            }
        }
    }
    return sum / count;
}

// 计算小于等于k的度的节点占比
double UDG::acc_degree_proportion_less_int(int k){
    int n = mg.vertex_num;
    int less_int_node_number = 0;
    int count = 0;
    for(int i = 0; i < n; i++){
        if(gd.degree[i] > 0){
            count ++;
            if(gd.degree[i] <= k){
                less_int_node_number += 1;
            }
        }
    }
    double result = double(less_int_node_number) / count;
    return result;
}

// 采样距离大于3的样本为负样本
void UDG::sample_neg(int number){
    // cout << 1 << endl;
    set<pair<int,int>> neg_edges;

    std::random_device rd;
    std::mt19937 gen(rd());
    // 定义随机数分布范围
    std::uniform_int_distribution<int> distribution(0, this->mg.vertex_num-1);
    while(neg_edges.size() != mg.edge_num){
        int rand1 = distribution(gen);
        int rand2 = distribution(gen);
        while(rand1 == rand2){
            rand2 = distribution(gen);
        }
        int min1 = min(rand1,rand2);
        int max1 = max(rand1,rand2);
        int d = acc_distance3(min1,max1);
        if(d == -1){
            neg_edges.insert({min1,max1});
        }
    }
    ofstream f;
    f.open("./neg.csv");
    for(auto &p:neg_edges){
        f << p.first << " " << p.second << endl;
    }
    f.close();
}

// 计算平均共同邻居数量和平均DCN    
void UDG::acc_average_CN_DCN(){
    int n = this->mg.vertex_num;
    int xy_number = 0;
    int xy_cn4 = 0;
    double all_CN = 0;
    double all_DCN = 0;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(this->mg.edge[i][j] != 0){
                set<int> V1;
                set<int> V2;
                for(auto &xx:this->mg.allNei[i]){
                    for(auto &yy:this->mg.allNei[j]){
                        if(xx == j or yy == i) continue;
                        if(this->mg.edge[xx][yy] != 0 && xx != yy){
                            V1.insert(xx);
                            V2.insert(yy);
                        }
                    }
                }
                if(V1.size() == 0 || V2.size() == 0) continue;
                xy_number++;
                int common_nei_number = 0;
                double xy_DCN = 0;
                for(auto &node:V1){
                    auto com = this->gd.cn_number[j][node];
                    auto com_DCN = double(this->gd.cn_number[j][node]) / this->gd.degree[node];
                    common_nei_number += com;
                    xy_DCN += com_DCN;
                }
                for(auto &node:V2){
                    auto com = this->gd.cn_number[i][node];
                    auto com_DCN = double(this->gd.cn_number[i][node]) / this->gd.degree[node];
                    common_nei_number += com;
                    xy_DCN += com_DCN;
                }
                double ave = double(common_nei_number) / (V1.size() + V2.size()); 
                if(ave >= 4) xy_cn4++;
                double ave_DCN = xy_DCN / (V1.size() + V2.size());
                all_CN += ave;
                all_DCN += ave_DCN;
            }
        }
    }
    // printf("%f,%f,",all_CN / xy_number,all_DCN / xy_number);
    printf("%f\n",double(xy_cn4) / xy_number);
    // ofstream f("../temp/average_CN_DCN.csv",ios::app);
    // f << all_CN / xy_number << "," << all_DCN / xy_number << endl;
    // f.close();
    return;
}

void UDG::acc_max_CN_DCN(){
    int n = this->mg.vertex_num;
    int xy_number = 0;
    double all_CN = 0;
    double all_DCN = 0;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(this->mg.edge[i][j] != 0){
                set<int> V1;
                set<int> V2;
                for(auto &xx:this->mg.allNei[i]){
                    for(auto &yy:this->mg.allNei[j]){
                        if(xx == j or yy == i) continue;
                        if(this->mg.edge[xx][yy] != 0 && xx != yy){
                            V1.insert(xx);
                            V2.insert(yy);
                        }
                    }
                }
                if(V1.size() == 0 || V2.size() == 0) continue;
                xy_number++;
                int max_xy_CN = 0;
                double max_xy_DCN = 0;
                for(auto &node:V1){
                    auto com = this->gd.cn_number[j][node];
                    auto com_DCN = double(this->gd.cn_number[j][node]) / this->gd.degree[node];
                    max_xy_CN  = max(com,max_xy_CN);
                    max_xy_DCN = max(com_DCN,max_xy_DCN);
                }
                for(auto &node:V2){
                    auto com = this->gd.cn_number[i][node];
                    auto com_DCN = double(this->gd.cn_number[i][node]) / this->gd.degree[node];
                    max_xy_CN  = max(com,max_xy_CN);
                    max_xy_DCN = max(com_DCN,max_xy_DCN);
                }
                // double ave = double(common_nei_number) / (V1.size() + V2.size()); 
                // double ave_DCN = xy_DCN / (V1.size() + V2.size());
                all_CN += max_xy_CN;
                all_DCN += max_xy_DCN;
            }
        }
    }
    printf("%f,%f\n",all_CN / xy_number,all_DCN / xy_number);
    ofstream f("../temp/max_CN_DCN.csv",ios::app);
    f << all_CN / xy_number << "," << all_DCN / xy_number << endl;
    f.close();
    return;
}
// 计算xy四边形图 边的稠密程度
void UDG::acc_xy_quadrangle_graph_edges_density(){
    int n = this->mg.vertex_num;
    double sum_density = 0;
    int numbers = 0;
    int numbers_dayu_12 = 0;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(this->mg.edge[i][j] != 0){
                set<int> V1;
                set<int> V2;
                set<pair<int,int>> Edges;
                V1.insert(j);
                V2.insert(i);
                Edges.insert({min(i,j),max(i,j)});
                for(auto &xx:this->mg.allNei[i]){
                    for(auto &yy:this->mg.allNei[j]){
                        if(xx == j or yy == i) continue;
                        if(this->mg.edge[xx][yy] != 0 && xx != yy){
                            V1.insert(xx);
                            V2.insert(yy);
                            Edges.insert({min(i,xx),max(i,xx)});
                            Edges.insert({min(j,yy),max(j,yy)});
                            Edges.insert({min(xx,yy),max(xx,yy)});
                        }
                    }
                } 
                if(Edges.size() == 1) continue;
                double density = double(Edges.size()) / V1.size() / V2.size();
                sum_density += density;
                if(density >= 0.5){
                    numbers_dayu_12++;
                }
                numbers++;
            }
        }
    }
    // cout << numbers;
    printf("%f,%f\n",sum_density / numbers,double(numbers_dayu_12) / numbers);
    // ofstream f("../temp/xy_quadrangle_graph_edges_density.csv",ios::app);
    // f << sum_density / numbers << double(numbers_dayu_12) / numbers << endl;
    // f.close();
    return;
}

vector<vector<double>> UDG::JC(){
    int n = this->mg.vertex_num;
    vector<vector<double>> similarMatrix(n,vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            int similar1 = this->gd.cn_number[i][j];
            int similar2 = this->commonUnion_xy(i,j).size();
            if(similar2 == 0) continue;
            double similar = similar1 / double(similar2);
            similarMatrix[i][j] = similar;
            similarMatrix[j][i] = similar;                
        }
    }
    return similarMatrix;  
}

// 计算相似性驱动形成边的比例
int UDG::acc_similar_drived_edges_number(){
    int n = this->mg.vertex_num;
    int sim_edges = 0;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(this->mg.edge[i][j] != 0){
                double sim1 = double(gd.cn_number[i][j]) / double(gd.degree[i]);
                double sim2 = double(gd.cn_number[i][j]) / double(gd.degree[j]);
                if(sim1 > 0.5 || sim2 > 0.5){
                    sim_edges++;
                }
            }
        }
    }
    return sim_edges;
}

std::vector<tuple<double,double,int>> UDG::acc_xy_quadrangle_graph_cn_and_edges_density(){
    std::vector<std::tuple<double,double,int>> res;
    int n = mg.vertex_num;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            auto ave_cn = _acc_xy_quadrangle_graph_i_j_ave_cn(i,j);
            auto edges_density = _acc_xy_quadrangle_graph_i_j_edges_density(i,j);
            if(ave_cn == 0 || edges_density == 0) continue;
            res.push_back(make_tuple(ave_cn,edges_density,mg.edge[i][j]));
        }
    }
    return res;
}