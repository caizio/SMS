#ifndef _DATADEAL_
#define _DATADEAL_

#include "baseStruct.hpp"
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <regex>
#include <set>
#include <random>


using namespace std;
class DataDeal{
public:
    void readDatas(string dir,int i,vector<vector<int>> &train,vector<vector<int>> &test,vector<pair<int,int>>&trainPair,vector<pair<int,int>>&testPair);
    void readDatas_cs6(string dir,int i,vector<vector<int>> &train,vector<vector<int>> &test,vector<pair<int,int>>&trainPair,vector<pair<int,int>>&testPair,string temp_dir = "cs6");
    void readDatas2(string dir,int i,vector<vector<int>> &train,vector<vector<int>> &test,vector<pair<int,int>>&trainPair,vector<pair<int,int>>&testPair,int ratio = 90);
    void readDatas3(string dir,int i,vector<vector<int>> &train,vector<vector<int>> &test,vector<pair<int,int>>&trainPair,vector<pair<int,int>>&testPair,int ratio = 90);
    void readDatas(string training_set,string test_set, vector<vector<int>> &train,vector<vector<int>> &test,vector<pair<int,int>>&trainPair,vector<pair<int,int>>&testPair);
    vector<string> split(string s,string label);
    int readData(string file,vector<pair<int,int>> &data);
    vector<vector<int>> readGraph(string file);
    vector<vector<int>> pairToMatrix(vector<pair<int,int>> edge,int nodeNum = 0);
    string getDataSetName(string dir);
    int writeData(const string &dir,const string &data);
    void writeData(vector<pair<int,int>> data,string file,bool is_continue = true);
    void writeData(map<int,int>,string file);
    void writeData(string path,const vector<double> &data);
    void writeString(string data,string dir);
    int writeLinkRall(const string &dir,const string &datasetName,const vector<baseStruct::LinkRall> &data);
    int writeLinkRall2(const string &dir,const string &datasetName,const vector<baseStruct::LinkRall> &data);
    void cross_validation(string dir,int div_size);
    void random_divide(string dir,double rate,int n = 1,string outPath = "Null");
    void TxtToNet(string file);
    void writeSubGraph(int i, int j, string file, const set<pair<int,int>> &edges);
    void randomAddEdges(const string &dir, double rate = 0.1);
    map<pair<int,int>,double> read_seq_compare_result(const string &path);
    void save_seq_compare_result(const map<pair<int,int>,double> &data, const string &sava_path);
    map<int,string> readSeq(const string &path);
    void sava_one_ppi_similar(int x,int y,double similar,string path);
    set<pair<int,int>> sample_neg(const vector<pair<int,int>> &edges ,int n = 0);
};

#endif