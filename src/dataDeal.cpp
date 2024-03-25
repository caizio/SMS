#include "dataDeal.h"
using namespace std;

// 采用正则匹配对字符串进行分割
vector<string> DataDeal::split(string s,string label){
    regex ws_re(label);
    vector<string> v(sregex_token_iterator(s.begin(),s.end(),ws_re,-1),sregex_token_iterator());
    return v;
}
// 读取dir路径下的第i个训练集和测试集，并以边和矩阵的形式储存
void DataDeal::readDatas(string dir,int i,vector<vector<int>> &train,vector<vector<int>> &test,vector<pair<int,int>>&trainPair,vector<pair<int,int>>&testPair){
    DataDeal dd;
    string tempPath = dir.substr(0,dir.size()-4); 
    string trainPath = tempPath + "Train_" + to_string(i) + ".txt";
    string testPath = tempPath + "Test_" + to_string(i) + ".txt";
    auto trainSize = dd.readData(trainPath,trainPair);
    auto testSize = dd.readData(testPath,testPair);
    int size = trainSize >= testSize ? trainSize : testSize;
    train = dd.pairToMatrix(trainPair,size);
    test = dd.pairToMatrix(testPair,size);
}

void DataDeal::readDatas_cs6(string dir,int i,vector<vector<int>> &train,vector<vector<int>> &test,vector<pair<int,int>>&trainPair,vector<pair<int,int>>&testPair,string temp_dir){
    DataDeal dd;
    string tempPath = dir.substr(0,dir.size()-4) + "_" + temp_dir + "_"; 
    string trainPath = tempPath + "Train_50%" + to_string(i) + ".txt";
    string testPath = tempPath + "Test_50%" + to_string(i) + ".txt";
    auto trainSize = dd.readData(trainPath,trainPair);
    auto testSize = dd.readData(testPath,testPair);
    int size = trainSize >= testSize ? trainSize : testSize;
    train = dd.pairToMatrix(trainPair,size);
    test = dd.pairToMatrix(testPair,size);
}

void DataDeal::readDatas(string training_set,string test_set, vector<vector<int>> &train,vector<vector<int>> &test,vector<pair<int,int>>&trainPair,vector<pair<int,int>>&testPair){
    DataDeal dd;
    auto trainSize = dd.readData(training_set,trainPair);
    auto testSize = dd.readData(test_set,testPair);
    int size = trainSize >= testSize ? trainSize : testSize;
    train = dd.pairToMatrix(trainPair,size);
    test = dd.pairToMatrix(testPair,size);
}

// 读取dir路径下的第i个训练集和测试集，并以边和矩阵的形式储存
void DataDeal::readDatas2(string dir,int i,vector<vector<int>> &train,vector<vector<int>> &test,vector<pair<int,int>>&trainPair,vector<pair<int,int>>&testPair,int ratio){
    DataDeal dd;
    string tempPath = dir.substr(0,dir.size()-4); 
    string trainPath = tempPath + "Train_" + to_string(i) + ".txt";
    string testPath = tempPath + "Test_" + to_string(i) + ".txt";
    auto trainSize = dd.readData(trainPath,trainPair);
    auto testSize = dd.readData(testPath,testPair);
    int size = trainSize >= testSize ? trainSize : testSize;
    train = dd.pairToMatrix(trainPair,size);
    test = dd.pairToMatrix(testPair,size);
}

// 读取dir路径下的第i个训练集和测试集，并以边和矩阵的形式储存，适用于按照不同比例的跑的算法的读取数据
void DataDeal::readDatas3(string dir,int i,vector<vector<int>> &train,vector<vector<int>> &test,vector<pair<int,int>>&trainPair,vector<pair<int,int>>&testPair,int ratio){
    DataDeal dd;
    string tempPath = dir.substr(0,dir.size()-4); 
    string trainPath = tempPath + "Train_" + to_string(ratio) + "%" +  to_string(i) + ".txt";
    string testPath = tempPath + "Test_" +  to_string(100 - ratio) + "%" +  to_string(i) + ".txt";
    auto trainSize = dd.readData(trainPath,trainPair);
    auto testSize = dd.readData(testPath,testPair);
    int size = trainSize >= testSize ? trainSize : testSize;
    train = dd.pairToMatrix(trainPair,size);
    test = dd.pairToMatrix(testPair,size);
}

// 将读取的数据以边的形式储存（data），一个pair是一条边,并返回节点的最大编号+1，即节点数量，因为编号是从0开始编号的，方便转矩阵
int DataDeal::readData(string file,vector<pair<int,int>> &data){
    ifstream infile;
    infile.open(file,ios::in);
    if(!infile.is_open()){
        cout << "readData::文件打开失败" << endl;
        return -1;
    }
    string s;
    int maxNode = 0;
    data.clear();
    while(getline(infile,s)){
        vector<string> text = DataDeal::split(s,",");
        pair<int,int> edge;
        edge.first = stoi(text[0]);
        edge.second = stoi(text[1]);
        data.push_back(edge);
        if(edge.first > maxNode) maxNode = edge.first;
        if(edge.second > maxNode) maxNode = edge.second;
    }
    infile.close();
    return maxNode + 1;
}
// 读取图数据，返回矩阵的形式
vector<vector<int>> DataDeal::readGraph(string file){
    vector<vector<int>> edge;
    vector<pair<int,int>> edgePair;
    auto edge_size = readData(file,edgePair);
    return pairToMatrix(edgePair,edge_size);
}

// 将读取的边转化成邻接矩阵，注意矩阵的大小是edge中最大的边+1，或者手动设置Nodenum1为最大节点
vector<vector<int>> DataDeal::pairToMatrix(vector<pair<int,int>> edge,int NodeNum1){
    int n = edge.size();
    // 寻找最大的节点
    int NodeNum2 = 0;
    if(NodeNum1 == 0){
        for(int i = 0; i < n; i++){
            if(edge[i].first > NodeNum2){NodeNum2 = edge[i].first;}
            if(edge[i].second > NodeNum2){NodeNum2 = edge[i].second;}
        }
        NodeNum2 += 1;
    }else{
        NodeNum2 = NodeNum1;
    }
    vector<vector<int>> result(NodeNum2,vector<int>(NodeNum2,0));
    for(int i = 0; i < n; i++){
        // cout << edge[i].first << "," << edge[i].second << endl;
        result[edge[i].first][edge[i].second] = 1;
        result[edge[i].second][edge[i].first] = 1;
    }
    return result;
}
// 获取数据集的名字
string DataDeal::getDataSetName(string dir){
    string lab = "/";
    auto substring = this->split(dir,lab);
    string res = substring[substring.size() - 1];
    return res.substr(0,res.size()-4);
}
// 写出一个string数据
int DataDeal::writeData(const string &dir,const string &data){
    ofstream f;
    f.open(dir, ios::app);
    if(!f.is_open()){
        cout << "error::DataDeal::writeData::写出数据失败！" << endl;
        return -1;
    }
    f << data;
    f.close();
    return 0;
}
// 写出一个向量数据
void DataDeal::writeData(string dir,const vector<double> &data){
    ofstream f;
    f.open(dir);
    if(!f.is_open()){
        cout << "error::DataDeal::writeData::写出数据失败！" << endl;
        return;
    }
    for(auto &p:data){
        f << p << ",";
    }
    f << endl;
    f.close();
}
// 将边数据写出到file文件中
void DataDeal::writeData(vector<pair<int,int>> data,string file, bool is_continue){
    int n = data.size();
    string sOut = "";
    for(int i = 0; i < n; i++){
        string temp = to_string(data[i].first) + "," + to_string(data[i].second) + "\n";
        sOut += temp;
    }
    ofstream f;
    if(is_continue == true){
        f.open(file,ios::app);
    }else{
        f.open(file);
    }
    f << sOut;
    f.close();
}

void DataDeal::writeData(map<int,int> data,string dir){
    ofstream f;
    f.open(dir,ios::app);
    if(!f.is_open()){
        cout << "error::DataDeal::writeData::写出数据失败！" << endl;
        return;
    }
    for(auto p:data){
        f << p.first << "," << p.second << endl;
    } 
    f.close();
}

void DataDeal::writeString(string data,string dir){
    ofstream f(dir,ios::app);
    if(!f.is_open()){
        cout << "error::DataDeal::writeData::写出数据失败！" << endl;
        return;
    }
    f << data;
    f.close();
}


// 写出一个数据集的P-R曲线结果,(R到0.1为止)
int DataDeal::writeLinkRall(const string &dir,const string &datasetName,const vector<baseStruct::LinkRall> &data){
    ofstream f;
    f.open(dir,ios::app);
    if(!f.is_open()){
        cout << "error::DataDeal::writeLinkRall::写出数据失败！" << endl;
        return -1;
    }
    int n = data.size();
    f << datasetName << endl;
    f << to_string(n) << endl;
    for(int i = 0; i < n; i++){
        f << data[i].name << endl;
        int m = data[i].pr.size();
        for(int j = 0; j < m - 1; j++){
            // 只打印R 0-0.1的部分
            f << data[i].pr[j].first << ",";
            if(data[i].pr[j].first > 0.1){
                m = j + 1;
                break;
            }
        }
        // f << data[i].pr[m-1].first;
        f << endl;
        for(int j = 0; j < m; j++){
            f << data[i].pr[j].second << ",";
        }
        f << endl;
    }
    f.close();
    return 0;
}
// 写出精度100、200、prauc、ndcg
int DataDeal::writeLinkRall2(const string &dir,const string &datasetName,const vector<baseStruct::LinkRall> &data){
    ofstream f;
    f.open(dir,ios::app);
    if(!f.is_open()){
        cout << "error::DataDeal::writeLinkRall2::写出数据失败！" << endl;
        return -1;
    }
    int n = data.size();
    f << datasetName << endl;
    f << ",P@100,P@200,P@500,AUPRC,AUROC,nDCG,AUC," << endl;
    for(int i = 0; i < n; i++){
        f << data[i].name << "," << data[i].precision[99] << "," << data[i].precision[199]  << "," << data[i].precision[499]  << ","
        << data[i].AUPRC << "," << data[i].AUROC << "," << data[i].nDCG << "," <<  data[i].auc << "," << endl;;
    }
    f.close();
    return 0;    
}
// 按照交叉验证的思路划分数据集，dir为原始数据的路径，div_size为划分的数量，2022_11_4
void DataDeal::cross_validation(string dir,int div_size){
    cout << "处理:" << endl;
    cout << dir << endl;  
    DataDeal dd;
    vector<pair<int,int>> edge;
    int nodeNum = dd.readData(dir,edge);
    if(nodeNum == -1){
        cout << dir << "error::cross_validation::数据为空" << endl;
        return;
    }
    vector<vector<pair<int,int>>> newEdges(div_size,vector<pair<int,int>>());
    int edgeNum = edge.size();
    int everyEdgeNums = edgeNum / div_size;
    srand((int)time(0));
    // 从edge中挑选edgeNum/div_size的边，放入newEdges[i]中,剩下的边是最后一份的
    for(int i = 0; i < div_size - 1; i++){
        for(int j = 0; j < everyEdgeNums; j++){
            int random = rand() % edgeNum;
            newEdges[i].push_back(edge[random]);
            swap(edge[random],edge[edge.size() - 1]);
            edge.pop_back();
            edgeNum--;
        }
    }
    newEdges[div_size-1] = edge;
    // 合并newEdges中的div_size-1份作为训练集，剩下的1份作为测试集，并写出到文件中
    for(int i = 0; i < div_size; i++){
        vector<pair<int,int>> train;
        vector<pair<int,int>> test;
        for(int j = 0; j < div_size; j++){
            // i==j则表示放入测试集，否则放入训练集
            if(i == j){
                for(auto &p:newEdges[j]){
                    test.push_back(p);
                }
            }else{
                for(auto &p:newEdges[j]){
                    train.push_back(p);
                }                
            }
        }
        auto path = dir.substr(0,dir.size() - 4);
        string trainPath =  path + "Train_" + to_string(i) + ".txt";
        string testPath =  path + "Test_" + to_string(i) + ".txt";
        dd.writeData(train,trainPath,false);
        dd.writeData(test,testPath,false);
        cout << "第" + to_string(i) + "轮划分：" << "train:" << train.size() << "  ";
        cout << "test:" << test.size() << endl;
    }
}
// 随机划分数据集（不是按照交叉验证）
// ratio为划分的比例，n为划分的次数
void DataDeal::random_divide(string dir,double ratio,int n,string outPath){
    srand((int)time(NULL)); 
    for(int j = 0; j < n; j++){
        cout << "处理:" << endl;
        cout << dir << endl;  
        DataDeal dd;
        vector<pair<int,int>> edge;
        int nodeNum = dd.readData(dir,edge);
        if(nodeNum == -1){
            cout << dir << "error::cross_validation::数据为空" << endl;
            return;
        }
        // 将edge中的边分别放到训练集trainEdge、测试集testEdge中
        int trainNumber = int(edge.size() * ratio);
        // cout << trainNumber << endl;
        int randSize = edge.size();
        vector<pair<int,int>> trainEdge;
        vector<pair<int,int>> testEdge;
        for(int i = 0; i < trainNumber; i++){
            int random = rand() % randSize;
            trainEdge.push_back(edge[random]);
            edge.erase(edge.begin() + random);
            randSize--;
        }
        testEdge = edge;
        string tempPath;
        if(outPath == "Null"){
            tempPath = dir.substr(0,dir.size()-4); 
        }else{
            tempPath = dir.substr(0,dir.size()-4) + "_" + outPath + "_"; 
        }
        string path1 = tempPath + "Train_" + to_string(int(ratio * 100)) + "%" + to_string(j) +".txt";
        string path2 = tempPath + "Test_" + to_string(100 - int(ratio * 100)) + "%" + to_string(j)+ ".txt";
        writeData(trainEdge,path1);
        writeData(testEdge,path2);
        int train_size = trainEdge.size();
        int test_size = testEdge.size();
        cout << "train:" << train_size << endl;
        cout << "test:" << test_size << endl;
    }

}

void DataDeal::TxtToNet(string file){
    ifstream infile;
    infile.open(file);
    if(!infile.is_open()){
        cout << file << "：读取失败！" <<endl;
        return;
    };
    vector<pair<int,int>> data;
    int nodeN = 0;
    string s;
    while(getline(infile,s)){
        string label = ",";
        auto temp = split(s,label);
        auto edge = make_pair(stoi(temp[0]),stoi(temp[1]));
        data.push_back(edge);
        nodeN = max(max(edge.first,edge.second),nodeN);
    }
    auto matrix = this->pairToMatrix(data);
    int nodeNum = nodeN;
    int edgeNum = data.size();
    string node = string("*Vertices") + string("\t") + to_string(nodeNum) + string("\n");
    for(int i = 1; i <= nodeNum; i++){
        node += to_string(i) + string("\t") + string("\"") + to_string(i) + string("\"") + string("\n");
    }
    string edge = string("*Edges") + string("\t") + to_string(edgeNum) + string("\n");
    for(int i = 0; i < data.size(); i++){
        edge += to_string(data[i].first) + string("\t") + to_string(data[i].second) + string("\n"); 
    }
    string res = node + edge;
    string outPath = file.substr(0,file.size() - 4) + string(".net");
    ofstream f;
    f.open(outPath);
    f << res;
    f.close();
    cout << outPath << ":已完成！" << endl;     
}

void DataDeal::writeSubGraph(int i,int j, string dir, const set<pair<int,int>> &edges){
    ofstream f;
    f.open(dir);
    if(!f.is_open()){
        cout << "error::DataDeal::writeData::写出数据失败！" << endl;
        return;
    }
    f << i << ":" << j << endl;
    for(auto &p:edges){
        f << p.first << "," << p.second << endl;
    }
    f.close();   
}
// 随机添加rate比例的边，并按照9：1的交叉验证的方式划分数据集
void DataDeal::randomAddEdges(const string &dir, double rate){
    // 正常交叉验证，分成10份
    // this->cross_validation(dir,10);
    // 将添加的边全部放入训练集中
    auto path = dir.substr(0,dir.size() - 4);
    for(int i = 0; i < 10; i++){
        vector<pair<int,int>> graphPair;
        auto size = this->readData(dir,graphPair);
        auto train = this->pairToMatrix(graphPair,size);
        int n = graphPair.size() * rate;
        vector<pair<int,int>> addEdges;
        for(int i = 0; i < n; i++){
            bool tag = false;
            while(tag == false){
                int random1 = rand() % size;
                int random2 = rand() % size;
                while(random1 == random2){
                    random2 = rand() % size;
                }
                if(random1 > random2){
                    auto tem = random1;
                    random1 = random2;
                    random2 = tem;
                }
                if(train[random1][random2] == 0){
                    addEdges.push_back(make_pair(random1,random2));
                    tag = true;
                }
            }
        }
        string trainPath =  path + "Train_" + to_string(i) + ".txt";
        this->writeData(addEdges,trainPath);
    }
}

// 读取序列的比对结果
map<pair<int,int>,double> DataDeal::read_seq_compare_result(const string &path){
    map<pair<int,int>,double> SeqCompareResult;
    ifstream f(path);
    if(!f){
        cout << "文件无法打开:" << path << endl;
        return SeqCompareResult;
    }
    string line;
    while(getline(f,line)){
        auto temp = split(line,",");
        if(temp.size() != 3){
            if(temp.size() != 3){
                cout << "seq compare size is wrong" << endl;
                for(auto p:temp){cout << p << endl;
            };
            continue;
        }
        };
        double similar = 0;
        try{
            similar = stod(temp[2]);
        }catch(const std::invalid_argument& e){
            std::cout << "转换失败，无效的参数: " << e.what() << std::endl;
            similar = 0;
        }
        pair<int,int> ids(stoi(temp[0]),stoi(temp[1]));
        // if(stod(temp[2]) == 0) cout << "0";
        SeqCompareResult.insert({ids,similar});
    }
    // for (const auto &pair : SeqCompareResult) {  
    //     std::cout << pair.first.first << "-"  << pair.first.second << ": " << pair.second << std::endl;  
    // }
    return SeqCompareResult;
}

// 保存序列比对的结果
void DataDeal::save_seq_compare_result(const map<pair<int,int>,double> &data, const string &sava_path){
    // 打开文件流，如果文件不存在则创建  
    ofstream file(sava_path, ios::out);  
  
    if (!file) {  
        // 如果文件打开失败，输出错误并退出  
        cerr << "Error opening file:" << sava_path << endl;  
        return;  
    }  
  
    // 遍历map，将数据写入文件  
    for (const auto &item : data) {  
        // 将pair的第一个元素（一个pair）转换为字符串  
        string key_str = to_string(item.first.first) + "," + to_string(item.first.second);  
        // 将double转换为字符串  
        string value_str = to_string(item.second);  
        // 将键值对写入文件，以逗号分隔  
        file << key_str << "," << value_str << endl;  
    }  
  
    // 关闭文件流  
    file.close(); 
}

map<int,string> DataDeal::readSeq(const string &path){
    ifstream f;
    f.open(path);
    if(!f){
        cout << "DataDeal::readSeq:文件读取失败:" << path << endl;
        exit(1);
    }
    map<int,string> res;
    string line;
    while(getline(f,line)){
        auto data = split(line,",");
        if(data.size() != 3){
            cout << "seq size is wrong" << endl;
            for(auto p:data){cout << p << endl;};
            continue;
        }
        data[2].erase(std::remove_if(data[2].begin(), data[2].end(), [](char c) { 
            if(c =='J' || c == 'O' || c == 'U'){
                return true;
            }else{
                return false;
            }
             }), data[2].end());
        res.insert({stoi(data[0]),data[2]});
    }
    return res;
}

void DataDeal::sava_one_ppi_similar(int x,int y,double similar,string path){
    // 打开文件流，如果文件不存在则创建  
    ofstream file(path, ios::app);  
    if (!file) {  
        // 如果文件打开失败，输出错误并退出  
        cerr << "Error opening file:" << path << endl;  
        return;  
    }  
    // int temp_min = min(x,y);
    // int temp_max = max(x,y);
    file << to_string(x) << "," << to_string(y) << "," << to_string(similar) << endl;  
    // 关闭文件流  
    file.close();     
}

// 采样出负样本
set<pair<int,int>> DataDeal::sample_neg(const vector<pair<int,int>>&edges,int n){
    set<pair<int,int>> neg_edges;
    int edge_number = edges.size();
    if(n != 0) edge_number = n;
    int node_number = 0;
    for(auto &edge:edges){
        node_number = max(node_number,edge.first);
        node_number = max(node_number,edge.second);
    }
    node_number++;
    std::random_device rd;
    std::mt19937 gen(rd());
    // 定义随机数分布范围
    std::uniform_int_distribution<int> distribution(0,node_number - 1);
    while(neg_edges.size() != edge_number){
        int rand1 = distribution(gen);
        int rand2 = distribution(gen);
        while(rand1 == rand2){
            rand2 = distribution(gen);
        }
        int min1 = min(rand1,rand2);
        int max1 = max(rand1,rand2);
        neg_edges.insert({min1,max1});
    }
    return neg_edges;
}