// 自己写的矩阵库
#ifndef _Matrix_
#define _Matrix_

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <climits>
#include <map>

using namespace std;
// 对向量相加进行重载
template<typename T>
std::vector<T> operator+(const std::vector<T> &data1, const std::vector<T> &data2){
    int n = data1.size(), m = data2.size();
    if(n == 0) return data2;
    if(m == 0) return data1;
    n = n >= m ? m : n;
    std::vector<T> res(n,0);
    for(int i = 0; i < n; i++){
        res[i] = data1[i] + data2[i];
    }
    return res;
}
// 数乘重载
template<typename T>
std::vector<T> operator*(const std::vector<T> &vec, double a){
    int n = vec.size();
    vector<T> res(n);
    for(int i = 0; i < n; i++){
        res[i] = vec[i] * a;
    }
    return res;
}

class Matrix_my{
public:
    Matrix_my();
    ~Matrix_my();
    template <typename T> static void showVector(const std::vector<T> &vec1);
    template <typename T> static void showMatrix(const std::vector<std::vector<T>> &data,int way = 0);
    template <typename T> static double dot(std::vector<T> vec1,std::vector<T> vec2);
    template <typename T> static vector<vector<T>> zeroToMax(const vector<vector<T>> &data); 
    template <typename T> static bool isInVector(T d,const vector<T> &data);
    template <typename T> static double cos(const vector<int> &vec1,const vector<int> &vec2);
    template <typename T> static double jaccard(const vector<int> &vec1,const vector<int> &vec2);
    static std::vector<std::vector<int>> doubleToInt(const std::vector<std::vector<double>> &data);
    static std::vector<std::vector<double>> intToDouble(const std::vector<std::vector<int>> &data);
    template <typename T>
    static vector<vector<double>> MatrixNormalization(const vector<vector<T>> &data,int catagory = 1);
    template <typename T>
    static vector<vector<double>> typetoDouble(const vector<vector<T>> &data);
    static vector<vector<double>> transpose(const vector<vector<double>> &vec1);
    static vector<vector<int>> eye(int x);
    static vector<vector<double>> eyeD(int x);
    static vector<vector<double>> gauss(vector<vector<double>> data);
};

// 打印向量
template <typename T>
void Matrix_my::showVector(const std::vector<T> &vec1){
    int n = vec1.size();
    std::cout << "[";
    for(int i = 0; i < n - 1; i++){
        std::cout << vec1[i] << ",";
    }
    std::cout << vec1[n-1] << "]" << std::endl;
}
// 打印矩阵，参数：way==0或1 2种打印方式
template <typename T> 
void Matrix_my::showMatrix(const std::vector<std::vector<T>> &data,int way){
    if(way == 0){
        std::cout << "[";
        for(auto &vec:data){
            int n = vec.size();
            std::cout << "[";
            for(int i = 0; i < n - 1; i++){
                std::cout << vec[i] << ",";
            }
            std::cout << vec[n-1] << "]";
        }
        std::cout << "]" << std::endl;
    }
    if(way == 1){
        std::cout << "---------------------------------" << std::endl;
        int n = data.size();
        for(int i = 0; i < n; i++){
            showVector(data[i]);
        }
        std::cout << "---------------------------------" << std::endl;
    }    
}
// 向量点乘
template <typename T>
double Matrix_my::dot(std::vector<T> vec1,std::vector<T> vec2){   
    int n = vec1.size();
    int m = vec2.size();
    if(m == 0 || n == 0 || m != n){
        std::cout << "Matrix::dot::error:" << std::endl;
        std::cout << "m.size() != n.size()" << std::endl;
        return 0;
    }else{
        double sum = 0;
        for(int i = 0; i < m; i++){
            double temp = vec1[i] * vec2[i];
            sum += temp;
        }
        return sum;
    }
}

// 将矩阵中的0变成最大,对角线不变
template <typename T>
vector<vector<T>> Matrix_my::zeroToMax(const vector<vector<T>> &data){
    int n = data.size();
    vector<vector<T>> res = data;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(res[i][j] == 0 && i != j) res[i][j] = 99999999;
        }
    }
    return res;
}

// 判断d是否在data中
template <typename T>
bool Matrix_my::isInVector(T d,const vector<T> &data){
    for(auto &dd:data){
        if(dd == d){
            return true;
        }
    }
    return false;
}

// 对矩阵归一化处理，catagroy = 1,按行归一化，catagory=2,按列归一化(未写)
template <typename T>
vector<vector<double>> Matrix_my::MatrixNormalization(const vector<vector<T>> &data,int catagory){
    int n = data.size(); int m = data[0].size();
    vector<vector<double>> res(n,vector<double>(m,0));
    for(int i = 0; i < n; i++){
        double sum = 0;
        for(int j = 0; j < m; j++){
            sum += data[i][j];
        }
        if(sum != 0){
            for(int j = 0; j < m; j++){
                res[i][j] = data[i][j] / sum;
            }
        }
    }
    return res;
}

// 将矩阵的数据类型转为double
template <typename T>
vector<vector<double>> Matrix_my::typetoDouble(const vector<vector<T>> &data){
    int n = data.size();int m = data[0].size();
    vector<vector<double>> res(n,vector<double>(m,0));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            res[i][j] = double(data[i][j]);
        }
    }
    return res;
}



#endif