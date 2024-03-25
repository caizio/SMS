#pragma once
#include<vector>
#include<string>
#include<map>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<regex>
#include<utility>

namespace bio{
    // 储存蛋白质序列数据的数据结构
    struct protein_sequence{
        int number;
        std::string protein;
        std::string sequence;        
    };
    struct result_blosub{
        int score;
        std::vector<std::vector<int>> path;
    };
    struct result_path{
        std::string protein1;
        std::string protein2;
        double similar;
    };

    // 序列比对算法
    class SequenceCompare{
    private:
        // BLOSUM62矩阵
        std::vector<std::vector<int>> scoreMatrix = {
           // _, A, B, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, X, Y, Z
            { 1,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,}, //_
            {-4, 4,-2, 0,-2,-1,-2, 0,-2,-1,-1,-1,-1,-2,-1,-1,-1, 1, 0, 0,-3, 0,-2,-1,}, //A
            {-4,-2, 4,-3, 4, 1,-3,-1, 0,-3, 0,-4,-3, 3,-2, 0,-1, 0,-1,-3,-4,-1,-3, 1,}, //B
            {-4, 0,-3, 9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2,-2,-3,}, //C
            {-4,-2, 4,-3, 6, 2,-3,-1,-1,-3,-1,-4,-3, 1,-1, 0,-2, 0,-1,-3,-4,-1,-3, 1,}, //D
            {-4,-1, 1,-4, 2, 5,-3,-2, 0,-3, 1,-3,-2, 0,-1, 2, 0, 0,-1,-2,-3,-1,-2, 4,}, //E
            {-4,-2,-3,-2,-3,-3, 6,-3,-1, 0,-3, 0, 0,-3,-4,-3,-3,-2,-2,-1, 1,-1, 3,-3,}, //F
            {-4, 0,-1,-3,-1,-2,-3, 6,-2,-4,-2,-4,-3, 0,-2,-2,-2, 0,-2,-3,-2,-1,-3,-2,}, //G
            {-4,-2, 0,-3,-1, 0,-1,-2, 8,-3,-1,-3,-2, 1,-2, 0, 0,-1,-2,-3,-2,-1, 2, 0,}, //H
            {-4,-1,-3,-1,-3,-3, 0,-4,-3, 4,-3, 2, 1,-3,-3,-3,-3,-2,-1, 3,-3,-1,-1,-3,}, //I
            {-4,-1, 0,-3,-1, 1,-3,-2,-1,-3, 5,-2,-1, 0,-1, 1, 2, 0,-1,-2,-3,-1,-2, 1,}, //K
            {-4,-1,-4,-1,-4,-3, 0,-4,-3, 2,-2, 4, 2,-3,-3,-2,-2,-2,-1, 1,-2,-1,-1,-3,}, //L
            {-4,-1,-3,-1,-3,-2, 0,-3,-2, 1,-1, 2, 5,-2,-2, 0,-1,-1,-1, 1,-1,-1,-1,-1,}, //M
            {-4,-2, 3,-3, 1, 0,-3, 0, 1,-3, 0,-3,-2, 6,-2, 0, 0, 1, 0,-3,-4,-1,-2, 0,}, //N
            {-4,-1,-2,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2, 7,-1,-2,-1,-1,-2,-4,-2,-3,-1,}, //P
            {-4,-1, 0,-3, 0, 2,-3,-2, 0,-3, 1,-2, 0, 0,-1, 5, 1, 0,-1,-2,-2,-1,-1, 3,}, //Q
            {-4,-1,-1,-3,-2, 0,-3,-2, 0,-3, 2,-2,-1, 0,-2, 1, 5,-1,-1,-3,-3,-1,-2, 0,}, //R
            {-4, 1, 0,-1, 0, 0,-2, 0,-1,-2, 0,-2,-1, 1,-1, 0,-1, 4, 1,-2,-3, 0,-2, 0,}, //S
            {-4, 0,-1,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1, 0,-1,-1,-1, 1, 5, 0,-2, 0,-2,-1,}, //T
            {-4, 0,-3,-1,-3,-2,-1,-3,-3, 3,-2, 1, 1,-3,-2,-2,-3,-2, 0, 4,-3,-1,-1,-2,}, //V
            {-4,-3,-4,-2,-4,-3, 1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,-2, 2,-3,}, //W
            {-4, 0,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,-1,-1, 0, 0,-1,-2,-1,-1,-1,}, //X
            {-4,-2,-3,-2,-3,-2, 3,-3, 2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1, 2,-1, 7,-2,}, //Y
            {-4,-1, 1,-3, 1, 4,-3,-2, 0,-3, 1,-3,-1, 0,-1, 3, 0, 0,-1,-2,-3,-1,-2, 4,}, //Z
        };
        std::map<char,int> acid = {
            {'A',1},{'B',2},{'C',3},{'D',4},{'E',5},{'F',6},{'G',7},{'H',8},{'I',9},{'K',10},{'L',11},{'M',12},
            {'N',13},{'P',14},{'Q',15},{'R',16},{'S',17},{'T',18},{'V',19},{'W',20},{'X',21},{'Y',22},{'Z',23},{'-',0}
        };
        int maxInt(int a,int b,int c);
        int maxInt(int a,int b,int c,int d);
    public:
        SequenceCompare();
        ~SequenceCompare();
        int getScore(char a, char b);
        double JC(std::string &a,std::string &b);
        double LCS(std::string &a,std::string &b);
        result_blosub blosub(std::string &a,std::string &b);
        result_path find_path(std::string &a,std::string &b,result_blosub &rb,int i,int j);
        double NW(std::string &seq_a, std::string &seq_b);
        double similarByCombine(std::vector<std::pair<int,int>> &merged,std::string &seq);
        double similarByCombine(std::vector<std::pair<int,int>> &merged,std::string &seq_a,std::string &seq_b);
        double compare(std::vector<protein_sequence> &data,int i,int j);
        double compare2(std::vector<std::string> &data,int i,int j);
    };

    // 读取数据
    class ReadProteins{
    public:
        ReadProteins();
        ~ReadProteins();
        std::vector<protein_sequence> read_protein_sequence(std::string path);
        std::vector<std::string> read_protein_sequence2(std::string path);
        std::vector<std::string> split(std::string s,std::string label);
        std::string judge(const std::string &s);
    };
    // 一些常规算法
    class Algorithm{
    public:
        Algorithm();
        ~Algorithm();
        static int maxInt(int a,int b,int c);
    };
}


