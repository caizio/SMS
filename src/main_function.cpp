#include "main_function.h"

using namespace std;
vector<baseStruct::LinkRall> m_LinkPrediction(string dir,std::vector<std::string>names,vector<double> parameter,string outfile,int n,vector<double> parameter2){
        if(parameter.size() < names.size() || parameter2.size() < names.size()){
            cout << "error::m_LinkPrediction::parameter.size() < names.size() or parameter2.size() < names.size()";
        }
        LinkPrediction lp;
        DataDeal dd;
        Evaluator ev;
        cout << getCurrentSystemTime() << endl;
        string datasetName = dd.getDataSetName(dir);
        cout << datasetName << endl;
        int methods_number = names.size();
        vector<baseStruct::LinkRall> lkrs(methods_number);
        for(int i = 0; i < methods_number; i++){
            lkrs[i].name = names[i];
        }
        for(int i = 0; i < n; i++){
            // 读取数据
            vector<vector<int>> train;vector<vector<int>> test;
            vector<pair<int,int>> trainPair;vector<pair<int,int>> testPair;
            dd.readDatas(dir,i,train,test,trainPair,testPair);
            UDG udg(train);
            // 计算相似性及评价指标
            vector<vector<double>> similar;
            baseStruct::LinkRall result;
            for(int j = 0; j < methods_number; j++){
                similar = lp.mmain(names[j],udg,parameter[j],parameter2[j]);
                result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
                lkrs[j] + result;
            }
        }
        // 交叉验证的结果取均值并打印
        for(int i = 0; i < methods_number; i++){
            lkrs[i].getAverage(n);   
            lkrs[i].show();         
        }
        // dd.writeLinkRall(outfile + "s.txt",datasetName,lkrs);
        // 写出各性能的结果
        dd.writeLinkRall2(outfile + ".csv",datasetName,lkrs);
        return lkrs; 
}

vector<baseStruct::LinkRall> m_LinkPrediction_test(string dir,vector<string> names, vector<string> parameter,string outFile,int n){
        LinkPrediction lp;
        DataDeal dd;
        Evaluator ev;
        cout << getCurrentSystemTime() << endl;
        string datasetName = dd.getDataSetName(dir);
        cout << datasetName << endl;
        int methods_number = names.size();
        vector<baseStruct::LinkRall> lkrs(methods_number);
        for(int i = 0; i < methods_number; i++){
            lkrs[i].name = names[i];
        }
        for(int i = 0; i < n; i++){
            // 读取数据
            vector<vector<int>> train;vector<vector<int>> test;
            vector<pair<int,int>> trainPair;vector<pair<int,int>> testPair;
            dd.readDatas(dir,i,train,test,trainPair,testPair);
            UDG udg(train);
            // 计算相似性及评价指标
            vector<vector<double>> similar;
            baseStruct::LinkRall result;
            for(int j = 0; j < methods_number; j++){
                similar = lp.mmain(names[j],udg);
                result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
                lkrs[j] + result;
            }
        }
        // 交叉验证的结果取均值并打印
        for(int i = 0; i < methods_number; i++){
            lkrs[i].getAverage(n);   
            lkrs[i].show();         
        }
        // 写出PR曲线的数据,用于绘图（见tool里的python代码）
        // dd.writeLinkRall(outfile + "s.txt",datasetName,lkrs);
        // 写出各性能的结果
        dd.writeLinkRall2(outFile + ".csv",datasetName,lkrs);
        return lkrs;     
}

// 计算每条边对应两个节点的相似性并写出
void m_LinkPrediction_edge_seq_sim(string dir){
    LinkPrediction lp;
    DataDeal dd;
    bio::SequenceCompare sc;

    cout << getCurrentSystemTime() << endl;
    string datasetName = dd.getDataSetName(dir);
    cout << datasetName << endl;

    auto adj = dd.readGraph(dir);
    UDG G(adj);
    auto seq = dd.readSeq("../dataset/seq/seq_" + datasetName + "Map2.txt");

    map<pair<int,int>,double> sim;
    int n = G.mg.vertex_num;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(G.mg.edge[i][j] != 0){
                string sa = seq[i];
                string sb = seq[j];
                double sim_ij = sc.NW(sa,sb);
                sim[{i,j}] = sim_ij;
                // cout << sim_ij << ",";
            }
        }
        if(i % 10 == 0){
            cout << i << ",";
        }
    }

    ofstream file("../temp/edge_SW_sim/" + datasetName + ".txt", ios::app);  
    if (!file) {  
            // 如果文件打开失败，输出错误并退出  
            cerr << "Error opening file:" << "../temp/edge_SW_sim/" + datasetName + ".txt" << endl;  
    } 
    for(auto &s:sim){
        file << to_string(s.first.first) << "," << to_string(s.first.second) << "," << to_string(s.second) << endl;
    }
    file.close();
}

// 计算基于序列相似的SMS，DSMS，maxSMS算法的性能
vector<baseStruct::LinkRall> m_LinkPrediction_seq(string dir,string outfile,int dataset_repeat_count){
    LinkPrediction lp;
    DataDeal dd;
    Evaluator ev;
    cout << getCurrentSystemTime() << endl;
    string datasetName = dd.getDataSetName(dir);
    cout << datasetName << endl;
    int methods_number = 3;
    vector<baseStruct::LinkRall> lkrs(methods_number);
    lkrs[0].name = "SMS_seq";
    lkrs[1].name = "DSMS_seq";
    lkrs[2].name = "maxSMS_Seq";
    // lkrs[3].name = "maxDSMS_Seq";
    for(int i = 0; i < dataset_repeat_count; i++){
        // 读取数据
        vector<vector<int>> train;vector<vector<int>> test;
        vector<pair<int,int>> trainPair;vector<pair<int,int>> testPair;
        dd.readDatas(dir,i,train,test,trainPair,testPair);
        UDG udg(train);
        // 计算相似性及评价指标
        vector<vector<double>> similar;
        baseStruct::LinkRall result;
        bio::SequenceCompare sc;
        auto similar_seq = dd.read_seq_compare_result("../temp/sc_" + datasetName + ".txt");
        auto seq = dd.readSeq("../dataset/seq/seq_" + datasetName + "Map2.txt");
        ofstream file("../temp/sc_" + datasetName + ".txt", ios::app);  
        if (!file) {  
            // 如果文件打开失败，输出错误并退出  
            cerr << "Error opening file:" << "../temp/sc_" + datasetName + ".txt" << endl;  
        }    
        similar = lp.SMS_sequence(udg,similar_seq,seq,file);
        // dd.save_seq_compare_result(similar_seq,"../temp/sc_" + datasetName + ".txt");
        result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
        lkrs[0] + result;
        similar = lp.DSMS_sequence(udg,similar_seq,seq,file);
        // dd.save_seq_compare_result(similar_seq,"../temp/sc_" + datasetName + ".txt");
        result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
        lkrs[1] + result;
        similar = lp.maxSMS_sequence(udg,similar_seq,seq,file);
        // dd.save_seq_compare_result(similar_seq,"../temp/sc_" + datasetName + ".txt");
        result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
        lkrs[2] + result;
    }
    // 交叉验证的结果取均值并打印
    for(int i = 0; i < methods_number; i++){
        lkrs[i].getAverage(dataset_repeat_count);   
        lkrs[i].show();         
    }
    // // 写出各性能的结果
    // dd.writeLinkRall(outfile + "_SMSseq62.txt",datasetName,lkrs);
    // 写出PR曲线的数据,用于绘图（见tool里的python代码）
    dd.writeLinkRall2(outfile + "_SMSseq62.csv",datasetName,lkrs);
    return lkrs; 
}

void m_LinkPrediction2(string dir,std::vector<std::string>names,vector<double> parameter,string outfile,int n,vector<double> parameter2,int ratio){
        if(parameter.size() < names.size() || parameter2.size() < names.size()){
            cout << "error::m_LinkPrediction::parameter.size() < names.size() or parameter2.size() < names.size()";
        }
        LinkPrediction lp;
        DataDeal dd;
        Evaluator ev;
        cout << getCurrentSystemTime() << endl;
        string datasetName = dd.getDataSetName(dir);
        cout << datasetName << endl;
        int methods_number = names.size();
        for(int i = 0; i < n; i++){
            vector<baseStruct::LinkRall> lkrs(methods_number);
            for(int j = 0; j < methods_number; j++){
                lkrs[j].name = names[j];
            }
            // 读取数据
            vector<vector<int>> train;vector<vector<int>> test;
            vector<pair<int,int>> trainPair;vector<pair<int,int>> testPair;
            dd.readDatas2(dir,i,train,test,trainPair,testPair,ratio);
            UDG udg(train);
            // 计算相似性及评价指标
            vector<vector<double>> similar;
            baseStruct::LinkRall result;
            for(int j = 0; j < methods_number; j++){
                similar = lp.mmain(names[j],udg,parameter[j],parameter2[j]);
                result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
                lkrs[j] + result;
            }
            // 写出PR曲线
            dd.writeLinkRall(outfile + ".txt",datasetName,lkrs);
            // 各算法的性能
            dd.writeLinkRall2(outfile + ".csv",datasetName,lkrs);
        }
}

// S = DCN * (SeqS ^ lam)
vector<baseStruct::LinkRall> m_LinkPrediction_SMS_mix(string dir,string outfile,int dataset_repeat_count,const vector<double> &lam){
    LinkPrediction lp;
    DataDeal dd;
    Evaluator ev;
    cout << getCurrentSystemTime() << endl;
    string datasetName = dd.getDataSetName(dir);
    cout << datasetName << endl;
    int methods_number = lam.size();
    vector<baseStruct::LinkRall> lkrs(methods_number);
    for(int i = 0; i < lam.size(); i++){
        lkrs[i].name = "SMS_mix_lam = " + to_string(lam[i]);
    }
    for(int i = 0; i < dataset_repeat_count; i++){
        // 读取数据
        vector<vector<int>> train;vector<vector<int>> test;
        vector<pair<int,int>> trainPair;vector<pair<int,int>> testPair;
        dd.readDatas(dir,i,train,test,trainPair,testPair);
        UDG udg(train);
        // 计算相似性及评价指标
        vector<vector<double>> similar;
        baseStruct::LinkRall result;
        bio::SequenceCompare sc;
        auto similar_seq = dd.read_seq_compare_result("../temp/sc_" + datasetName + ".txt");
        auto seq = dd.readSeq("../dataset/seq/seq_" + datasetName + "Map2.txt");
        ofstream file("../temp/sc_" + datasetName + ".txt", ios::app);  
        if (!file) {  
            // 如果文件打开失败，输出错误并退出  
            cerr << "Error opening file:" << "../temp/sc_" + datasetName + ".txt" << endl;  
        }    
        for(int j = 0; j < lam.size(); j++){
            // cout << lam[j];
            similar = lp.SMS_mix(udg,similar_seq,seq,file,lam[j]);
            result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
            lkrs[j] + result;
        }

    }
    // 交叉验证的结果取均值并打印
    for(int i = 0; i < methods_number; i++){
        lkrs[i].getAverage(dataset_repeat_count);   
        lkrs[i].show();         
    }
    // 写出各性能的结果
    dd.writeLinkRall(outfile + "_SMS_mix.txt",datasetName,lkrs);
    // 写出PR曲线的数据,用于绘图（见tool里的python代码）
    dd.writeLinkRall2(outfile + "_SMS_mix.csv",datasetName,lkrs);
    return lkrs;     
}

vector<baseStruct::LinkRall> m_LinkPrediction_DSMS_mix(string dir,string outfile,int dataset_repeat_count,const vector<double> &lam){
    LinkPrediction lp;
    DataDeal dd;
    Evaluator ev;
    cout << getCurrentSystemTime() << endl;
    string datasetName = dd.getDataSetName(dir);
    cout << datasetName << endl;
    vector<baseStruct::LinkRall> lkrs(lam.size());
    for(int i = 0; i < lam.size(); i++){
        lkrs[i].name = "DSMS_mix_lam = " + to_string(lam[i]);
    }

    for(int i = 0; i < dataset_repeat_count; i++){
        // 读取数据
        vector<vector<int>> train;vector<vector<int>> test;
        vector<pair<int,int>> trainPair;vector<pair<int,int>> testPair;
        dd.readDatas(dir,i,train,test,trainPair,testPair);
        UDG udg(train);
        // 计算相似性及评价指标
        vector<vector<double>> similar;
        baseStruct::LinkRall result;
        bio::SequenceCompare sc;
        auto similar_seq = dd.read_seq_compare_result("../temp/sc_" + datasetName + ".txt");
        auto seq = dd.readSeq("../dataset/seq/seq_" + datasetName + "Map2.txt");
        ofstream file("../temp/sc_" + datasetName + ".txt", ios::app);  
        if (!file) {  
            // 如果文件打开失败，输出错误并退出  
            cerr << "Error opening file:" << "../temp/sc_" + datasetName + ".txt" << endl;  
        }    
        for(int j = 0; j < lam.size(); j++){
            // cout << lam[j];
            similar = lp.DSMS_mix(udg,similar_seq,seq,file,lam[j]);
            result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
            lkrs[j] + result;
        }

    }
    // 交叉验证的结果取均值并打印
    for(int i = 0; i < lam.size(); i++){
        lkrs[i].getAverage(dataset_repeat_count);   
        lkrs[i].show();         
    }
    // 写出各性能的结果
    dd.writeLinkRall(outfile + "_DSMS_mix.txt",datasetName,lkrs);
    // 写出PR曲线的数据,用于绘图（见tool里的python代码）
    dd.writeLinkRall2(outfile + "_DSMS_mix.csv",datasetName,lkrs);
    return lkrs;     
}

vector<baseStruct::LinkRall> m_LinkPrediction_maxSMS_mix(string dir,string outfile,int dataset_repeat_count,const vector<double> &lam){
    LinkPrediction lp;
    DataDeal dd;
    Evaluator ev;
    cout << getCurrentSystemTime() << endl;
    string datasetName = dd.getDataSetName(dir);
    cout << datasetName << endl;
    vector<baseStruct::LinkRall> lkrs(lam.size());
    for(int i = 0; i < lam.size(); i++){
        lkrs[i].name = "maxSMS_mix_lam = " + to_string(lam[i]);
    }

    for(int i = 0; i < dataset_repeat_count; i++){
        // 读取数据
        vector<vector<int>> train;vector<vector<int>> test;
        vector<pair<int,int>> trainPair;vector<pair<int,int>> testPair;
        dd.readDatas(dir,i,train,test,trainPair,testPair);
        UDG udg(train);
        // 计算相似性及评价指标
        vector<vector<double>> similar;
        baseStruct::LinkRall result;
        bio::SequenceCompare sc;
        auto similar_seq = dd.read_seq_compare_result("../temp/sc_" + datasetName + ".txt");
        auto seq = dd.readSeq("../dataset/seq/seq_" + datasetName + "Map2.txt");
        ofstream file("../temp/sc_" + datasetName + ".txt", ios::app);  
        if (!file) {  
            // 如果文件打开失败，输出错误并退出  
            cerr << "Error opening file:" << "../temp/sc_" + datasetName + ".txt" << endl;  
        }    
        for(int j = 0; j < lam.size(); j++){
            similar = lp.maxSMS_mix(udg,similar_seq,seq,file,lam[j]);
            result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
            lkrs[j] + result;
        }

    }
    // 交叉验证的结果取均值并打印
    for(int i = 0; i < lam.size(); i++){
        lkrs[i].getAverage(dataset_repeat_count);   
        lkrs[i].show();         
    }
    // 写出各性能的结果
    dd.writeLinkRall(outfile + "_maxSMS_mix.txt",datasetName,lkrs);
    // 写出PR曲线的数据,用于绘图（见tool里的python代码）
    dd.writeLinkRall2(outfile + "_maxSMS_mix.csv",datasetName,lkrs);
    return lkrs;       
}

vector<baseStruct::LinkRall> m_LinkPrediction_SMS_parameter(string dir,string outfile,int dataset_repeat_count,const vector<double> &lam,string parameter){
    LinkPrediction lp;
    DataDeal dd;
    Evaluator ev;
    cout << getCurrentSystemTime() << endl;
    string datasetName = dd.getDataSetName(dir);
    cout << datasetName << endl;
    vector<baseStruct::LinkRall> lkrs(lam.size());
    for(int i = 0; i < lam.size(); i++){
        lkrs[i].name = parameter + "_lam = " + to_string(lam[i]);
    }

    for(int i = 0; i < dataset_repeat_count; i++){
        // 读取数据
        vector<vector<int>> train;vector<vector<int>> test;
        vector<pair<int,int>> trainPair;vector<pair<int,int>> testPair;
        dd.readDatas(dir,i,train,test,trainPair,testPair);
        UDG udg(train);
        // 计算相似性及评价指标
        vector<vector<double>> similar;
        baseStruct::LinkRall result;
        bio::SequenceCompare sc;
        auto similar_seq = dd.read_seq_compare_result("../temp/sc_" + datasetName + ".txt");
        auto seq = dd.readSeq("../dataset/seq/seq_" + datasetName + "Map2.txt");
        ofstream file("../temp/sc_" + datasetName + ".txt", ios::app);  
        if (!file) {  
            // 如果文件打开失败，输出错误并退出  
            cerr << "Error opening file:" << "../temp/sc_" + datasetName + ".txt" << endl;  
        }    
        for(int j = 0; j < lam.size(); j++){
            if(parameter == "SMS_mix"){
                similar = lp.SMS_mix(udg,similar_seq,seq,file,lam[j]);
                result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
                lkrs[j] + result;                
            }else if(parameter == "DSMS_mix"){
                similar = lp.DSMS_mix(udg,similar_seq,seq,file,lam[j]);
                result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
                lkrs[j] + result;               
            }else if(parameter == "maxSMS_mix"){
                similar = lp.maxSMS_mix(udg,similar_seq,seq,file,lam[j]);
                result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
                lkrs[j] + result;
            }else if(parameter == "SMS_DCN_mul_SeqS"){
                similar = lp.SMS_mix(udg,similar_seq,seq,file,lam[j]);
                result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
                lkrs[j] + result;  
            }

        }

    }
    // 交叉验证的结果取均值并打印
    for(int i = 0; i < lam.size(); i++){
        lkrs[i].getAverage(dataset_repeat_count);   
        lkrs[i].show();         
    }
    // 写出各性能的结果
    // dd.writeLinkRall(outfile + ".txt",datasetName,lkrs);
    // 写出PR曲线的数据,用于绘图（见tool里的python代码）
    dd.writeLinkRall2(outfile + "_p.csv",datasetName,lkrs);
    return lkrs;    
}

// 寻找DSMS_DCN的最优参数
void m_find_DSMS_parameter(string dataset, vector<double> parameter, string outfile, int repeat, int ratio){
        LinkPrediction lp;
        DataDeal dd;
        Evaluator ev;
        cout << getCurrentSystemTime() << endl;
        string datasetName = dd.getDataSetName(dataset);
        cout << datasetName << endl;
        int methods_number = parameter.size();
        for(int i = 0; i < repeat; i++){
            vector<baseStruct::LinkRall> lkrs(methods_number);
            for(int j = 0; j < methods_number; j++){
                lkrs[j].name = to_string(parameter[j]);
            }
            // 读取数据
            vector<vector<int>> train;vector<vector<int>> test;
            vector<pair<int,int>> trainPair;vector<pair<int,int>> testPair;
            dd.readDatas2(dataset,i,train,test,trainPair,testPair,ratio);
            UDG udg(train);
            // 计算相似性及评价指标
            vector<vector<double>> similar;
            baseStruct::LinkRall result;
            for(int j = 0; j < methods_number; j++){
                similar = lp.mmain("SMS",udg,9,1,parameter[j]);
                result = ev.acc_all_evaluators(similar,train,test,trainPair,testPair);
                lkrs[j] + result;
            }
            // 写出各性能的结果
            dd.writeLinkRall(outfile + ".txt",datasetName,lkrs);
            // 写出PR曲线的数据,用于绘图（见tool里的python代码）
            dd.writeLinkRall2(outfile + ".csv",datasetName,lkrs);
        }
}

// 统计JC或者DCN的连接比例
void m_3(string out_dir){
    DataDeal dd;
    for(int i = 0; i < LNBL3DataSet.size(); i++){
        string name = LNBL3DataSet[i];
        string datasetName = dd.getDataSetName(name);
        auto edge = dd.readGraph(name);
        UDG G(edge);
        // auto result = G.acc_sms_connection_JC(10);
        auto result = G.acc_sms_connection_DCN(10);
        cout << datasetName << endl;
        ofstream f;
        f.open(out_dir,ios::app);
        f << datasetName << endl;
        for(auto &p:result){
            double p_connection = p.connect / (p.non_connect + p.connect);
            cout << p_connection << endl;
            f << p.range_left << "," << p.range_right << "," << p.connect << "," << p.non_connect << "," << p_connection << endl;
        }
        f << endl;
        f.close();
    }    
}
// 统计maxSMS_CN的连接比例
void m_3_2(string out_dir){
    DataDeal dd;
    for(int i = 0; i < LNBL3DataSet.size(); i++){
        string name = LNBL3DataSet[i];
        string datasetName = dd.getDataSetName(name);
        auto edge = dd.readGraph(name);
        UDG G(edge);
        auto count = G.acc_sms_connection_CN();
        cout << datasetName << endl;
        ofstream f;
        f.open(out_dir,ios::app);
        f << datasetName << endl;
        for(auto &p:count){
            f << p.first << "," << p.second.first << "," << p.second.second << endl;
        }
        f << endl;
        f.close();
    }      
}
// 统计累积DCN的连接比例
void m_3maxDCN_l(string out_dir){
    DataDeal dd;
    for(int i = 0; i < LNBL3DataSet.size(); i++){
        string name = LNBL3DataSet[i];
        string datasetName = dd.getDataSetName(name);
        // string name = "D:\\code\\edit\\c++\\newGraph3_easy_SMS\\dataset\\test_graph2.txt";
        auto edge = dd.readGraph(name);
        UDG G(edge);
        auto count = G.acc_sms_connection_DCN2();
        cout << datasetName << endl;
        ofstream f;
        f.open(out_dir,ios::app);
        f << datasetName << endl;
        int connect = 0;
        int non_connect = 0;
        f << "maxDCN,connect edge,nonconnect edge,ratio," << endl;
        for(auto &p:count){
            connect += p.second.first;
            non_connect += p.second.second;
            f << p.first << "," << connect << "," << non_connect <<"," <<  endl;
        }
        f << endl;
        f.close();
    }       
}

void acc_average_cn_number(string dir){
    DataDeal dd;
    auto datasetName = dd.getDataSetName(dir);
    auto adj = dd.readGraph(dir);
    UDG G(adj);
    double average_cn_number = G.acc_average_cn_number();
    cout << datasetName << "," << average_cn_number << endl;
}

// 计算相连边的序列相似性
void acc_connected_edges_seq_similar(string dataset_name){
    string dir = "../dataset/LNBL3DataSet/" + dataset_name + "/" + dataset_name + ".txt";
    string seq_dir = "../dataset/seq/" + dataset_name + "/seq_" + dataset_name + "Map2.txt";
    string outfile = "../dataset/edge_seq_sim_" + dataset_name + ".txt";
    ofstream file(outfile);
    if(!file.is_open()) cout << dataset_name + " is not open";
    DataDeal dd;
    vector<pair<int,int>> edgePair;
    dd.readData(dir,edgePair);
    auto seq = dd.readSeq("../dataset/seq/seq_" + dataset_name + "Map2.txt");
    bio::SequenceCompare sc;
    for(auto &edge:edgePair){
        string s1 = seq[edge.first];
        string s2 = seq[edge.second];
        auto p = sc.blosub(s1,s2);
        auto similar_x_y = sc.find_path(s1,s2,p,s1.size(),s2.size()).similar;   
        file << to_string(edge.first) << "," << to_string(edge.second) << "," << to_string(similar_x_y) << endl;
    }
    file.close();
}