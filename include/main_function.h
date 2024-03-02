#ifndef _MAIN_FUNCTION_
#define _MAIN_FUNCTION_

#include <iostream>
#include <vector>
#include <string>

#include "matrix.h"
#include "graph.h"
#include "dataDeal.h"
#include "linkPrediction.h"
#include "evaluator.h"
#include "bioAlgorith.h"

#include "baseStruct.hpp"
#include "date.hpp"
#include "dataDir.h"

// 实验1、2,不同算法的性能比较
/*
    @brief Main function of comparative experiment, need to partition the dataset first
    @param [dir] The file path of the dataset
    @param [names] The name of the algorithm being run
    @param [parameter] Parameters of different algorithms
    @param [outFile] The path of experimental results
    @param [n] The number of times the cross validation algorithm runs
    @param [parameter2] Parameters2 of different algorithms
    @result Container of experimental results for different algorithms
*/
vector<baseStruct::LinkRall> m_LinkPrediction(string dir,vector<string> names, vector<double> parameter,string outFile = "./result/0",int n = 10,vector<double> parameter2 = {});
// 不同算法的比较，非交叉验证
void m_LinkPrediction2(string dir,vector<string> names, vector<double> parameter,string outFile = "./result/0",int n = 10,vector<double> parameter2 = {},int ratio = 90);
// 实验3 统计相连uv的JC或DCN相似性节点xy相连的概率
void m_3(string out_dir = "../result/m_3_jc.csv");
// 统计累积的DCN
void m_3maxDCN_l(string out_dir = "../result/maxCN.csv");
// 统计maxCN
void m_3_2(string out_dir = "../result/m_3_cn.csv");
// 寻找DSMS的最优参数
void m_find_DSMS_parameter(string dataset, vector<double> parameter, string oufFile = "./result/0", int repeat = 1, int raito = 90);

// 计算数据集的序列相似性
void m_LinkPrediction_edge_seq_sim(string dir);

vector<baseStruct::LinkRall> m_LinkPrediction_seq(string dir,string outfile,int dataset_repeat_count);
vector<baseStruct::LinkRall> m_LinkPrediction_SMS_mix(string dir,string outfile,int dataset_repeat_count,const vector<double> &lam);
vector<baseStruct::LinkRall> m_LinkPrediction_DSMS_mix(string dir,string outfile,int dataset_repeat_count,const vector<double> &lam);
vector<baseStruct::LinkRall> m_LinkPrediction_maxSMS_mix(string dir,string outfile,int dataset_repeat_count,const vector<double> &lam);

vector<baseStruct::LinkRall> m_LinkPrediction_SMS_parameter(string dir,string outfile,int dataset_repeat_count,const vector<double> &lam,string parameter);

vector<baseStruct::LinkRall> m_LinkPrediction_test(string dir,vector<string> names, vector<string> parameter = {},string outFile = "./result/0.csv",int n = 10);

void acc_average_cn_number(string dir);

void acc_connected_edges_seq_similar(string dataset);
#endif