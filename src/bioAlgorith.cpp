#include "bioAlgorith.h"

using namespace std;
bio::SequenceCompare::SequenceCompare(){};
bio::SequenceCompare::~SequenceCompare(){};



// 计算两个字符的匹配得分
int  bio::SequenceCompare::getScore(char a, char b) {
    return scoreMatrix[acid[a]][acid[b]];
}
// 统计字符串a和字符串b共有的字符/最大字符串长度
double bio::SequenceCompare::JC(std::string &a,std::string &b){
    map<char,int> couta;
    map<char,int> coutb;
    // 统计a和b的字符频率
    for(int i = 0; i < a.size(); i++){
        if(couta.count(a[i]) == 0){
            couta[a[i]] = 1;
        }else{
            couta[a[i]] ++;
        }
    }
    for(int i = 0; i < b.size(); i++){
        if(coutb.count(b[i]) == 0){
            coutb[b[i]] = 1;
        }else{
            coutb[b[i]] ++;
        }
    }
    int common = 0;
    int uni = max(a.size(),b.size());
    for(auto item = couta.begin(); item != couta.end(); item++){
        common += min(item->second,coutb[item->first]);
    }
    return common / double(uni);
}
// 经典LCS，寻找最大子序列
double bio::SequenceCompare::LCS(std::string &text1,std::string &text2){
    int m = text1.length(), n = text2.length();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1));
    for (int i = 1; i <= m; i++) {
        char c1 = text1.at(i - 1);
        for (int j = 1; j <= n; j++) {
            char c2 = text2.at(j - 1);
            if (c1 == c2) {
                dp[i][j] = dp[i - 1][j - 1] + 1;
            } else {
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
            }
        }
    }
    return dp[m][n] / double(max(m,n));
}
// 基于blosub得分矩阵比对的蛋白质比对算法
bio::result_blosub bio::SequenceCompare::blosub(std::string &a,std::string &b){
    int m = a.size(), n = b.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1,0));
    vector<vector<int>> path(m + 1, vector<int>(n + 1,0));
    // 初始化
    for(int i = 1; i <= m; i++){
        char c = a[i-1];
        if(acid.count(c) == 0){
            cout << "非法字符:" << c << endl;
            return {};
        }
        dp[i][0] = dp[i-1][0] + scoreMatrix[acid[c]][acid['-']];
    }
    for(int j = 1; j <= n; j++){
        char c = b[j-1];
        if(acid.count(c) == 0){
            cout << "非法字符:" << c << endl;
            return {};
        }
        dp[0][j] = dp[0][j-1] + scoreMatrix[acid['-']][acid[c]];
    }
    for (int i = 1; i <= m; i++){
        char c1 = a[i-1];
        for (int j = 1; j <= n; j++) {
            char c2 = b[j-1];
            int t1 = dp[i-1][j] + scoreMatrix[acid[c1]][acid['-']];
            int t2 = dp[i][j-1] + scoreMatrix[acid['-']][acid[c2]];
            int t3 = dp[i-1][j-1] + scoreMatrix[acid[c1]][acid[c2]];
            dp[i][j] = Algorithm::maxInt(t1,t2,t3);
            if(dp[i][j] == t3){
                path[i][j] = 2;
            }else if(dp[i][j] == t2){
                // 1表示左走
                path[i][j] = 1;
            }else{
                // 0表示上走
                path[i][j] = 0;
            }
        }
    }
    bio::result_blosub rb;
    rb.score = dp[m][n];
    rb.path = path;
    return rb;
}
// 将比对的过程还原，输出比对的序列
bio::result_path bio::SequenceCompare::find_path(std::string &a,std::string &b,bio::result_blosub &rb,int i,int j){
    string aa,bb;
    double same = 0;
    while(i > 0 && j > 0){
        if(rb.path[i][j] == 2){
            aa.push_back(a[i-1]);
            bb.push_back(b[j-1]);
            i--;
            j--;
            // cout << a[i-1] << "," << b[j-1] << endl;
            if(a[i] == b[j]){
                same++;
            }
        }else if(rb.path[i][j] == 1){
            aa.push_back('-');
            bb.push_back(b[j-1]);
            j--;
        }else{
            aa.push_back(a[i-1]);
            bb.push_back('-');
            i--;
        }
    }
    if(i > 0){
        for(int z = 0; z < i; z++){
            aa.push_back(a[z]);
            bb.push_back('-');
        }
    }
    if(j > 0){
        for(int z = 0; z < j; z++){
            aa.push_back('-');
            bb.push_back(b[z]);
        }
    }
    reverse(aa.begin(),aa.end());
    reverse(bb.begin(),bb.end());
    result_path rp;
    // cout << aa << endl;
    // cout << bb << endl;
    rp.protein1 = aa;
    rp.protein2 = bb; 
    // cout << aa.size() << "," << same << endl;
    rp.similar = same / double(aa.size());
    return rp;   
}

// 经典NW序列比对算法
double bio::SequenceCompare::NW(std::string &seq_a, std::string &seq_b){
    auto p = this->blosub(seq_a,seq_b);
    auto pp = this->find_path(seq_a,seq_b,p,seq_a.size(),seq_b.size());
    return pp.similar;
}

// 根据比对的区间计算在整个序列中的占比
double bio::SequenceCompare::similarByCombine(std::vector<std::pair<int,int>> &merged,std::string &seq){
    int l = seq.size();
    int sum = 0;
    for(auto &p:merged){
        sum += p.second - p.first + 1;
    }
    return double(sum) / l;
}

double bio::SequenceCompare::similarByCombine(vector<pair<int,int>> &merged,string &seq_a,string &seq_b){
    int l = max(seq_a.size(),seq_b.size());
    int sum = 0;
    for(auto &p:merged){
        sum += p.second - p.first + 1;
    }
    return double(sum) / l;
}

int bio::SequenceCompare::SequenceCompare::maxInt(int a,int b,int c){
    int d = max(a,b);
    return max(d,c);
}
int bio::SequenceCompare::SequenceCompare::maxInt(int a,int b,int c,int d){
    int e = max(a,b);
    int f = max(c,d);
    return max(e,f);
}

// 输入储存序列的数据和编号ij，输出序列比对的得分，其中data为读入的编号、蛋白质、序列数据，ij是要比较的蛋白质编号
double bio::SequenceCompare::compare(std::vector<bio::protein_sequence> &data,int i,int j){
    bio::SequenceCompare seq;
    string a = data[i-1].sequence, b = data[j-1].sequence;
    auto p = seq.blosub(a,b);
    auto pp = seq.find_path(a,b,p,a.size(),b.size());
    return pp.similar;
}
// 输入储存序列的数据和编号ij，输出序列比对的得分，其中data为读入的编号、蛋白质、序列数据，ij是要比较的蛋白质编号
double bio::SequenceCompare::compare2(std::vector<std::string> &data,int i,int j){
    bio::SequenceCompare seq;
    auto p = seq.blosub(data[i],data[j]);
    auto pp = seq.find_path(data[i],data[j],p,data[i].size(),data[j].size());
    return pp.similar;
}
bio::ReadProteins::ReadProteins(){};
bio::ReadProteins::~ReadProteins(){};
// 字符串分割
std::vector<std::string> bio::ReadProteins::split(std::string s,std::string label){
    regex ws_re(label);
    vector<string> v(sregex_token_iterator(s.begin(),s.end(),ws_re,-1),sregex_token_iterator());
    return v;
}
// 处理非法字符X
std::string bio::ReadProteins::judge(const std::string &s){
    // string isos = "ACDEFGHIKLMNPQRSTVWY";
    string res;
    for(auto &ss:s){
        for(auto it = s.begin(); it != s.end(); it++){
            if(*it != 'X'){
                res.push_back(*it);
            }
        }
    }
    return res;
}

// 读取蛋白质、序列数据
vector<bio::protein_sequence> bio::ReadProteins::read_protein_sequence(std::string path){
    ifstream f;
    f.open(path,ios::in);
    if(!f.is_open()){
        cout << "bio::ReadProteins::readTxt:" << path << ":文件打开失败！" << endl;
        return {};
    }
    string s;
    vector<bio::protein_sequence> res;
    int line = 0;
    while(getline(f,s)){
        auto temp_data = this->split(s,",");
        if(temp_data.size() != 3){
            cout << "bio::ReadProteins::readTxt:第" << line << "行数据读取失败！已返回前面成功读取的数据。" << endl;
            return res;
        }
        bio::protein_sequence temp;
        temp.number = stoi(temp_data[0]);
        temp.protein = temp_data[1];
        temp.sequence = judge(temp_data[2]);
        res.push_back(temp);
        line++;
    }
    return res;
}
// 读取蛋白质、序列数据
vector<string> bio::ReadProteins::read_protein_sequence2(std::string path){
    ifstream f;
    f.open(path,ios::in);
    if(!f.is_open()){
        cout << "bio::ReadProteins::readTxt:" << path << ":文件打开失败！" << endl;
        return {};
    }
    string s;
    // 因为蛋白质编号是从1开始的，保持一致性，第一行为空
    vector<string> res = {"null"};
    int line = 0;
    while(getline(f,s)){
        auto temp_data = this->split(s,",");
        if(temp_data.size() != 3){
            cout << "bio::ReadProteins::readTxt:第" << line << "行数据读取失败！已返回前面成功读取的数据。" << endl;
            return res;
        }
        auto temp = judge(temp_data[2]);
        res.push_back(temp);
        line++;
    }
    return res;
}

bio::Algorithm::Algorithm(){};
bio::Algorithm::~Algorithm(){};
int bio::Algorithm::maxInt(int a,int b,int c){
    int d = max(a,b);
    return max(d,c);
}

