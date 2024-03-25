#include "matrix.h"

std::vector<std::vector<int>> Matrix_my::doubleToInt(const std::vector<std::vector<double>> &data){
    int n = data.size();
    std::vector<std::vector<int>> res(n,std::vector<int>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n;j++){
            res[i][j] = data[i][j];
        }
    }
    return res;
}

std::vector<std::vector<double>> Matrix_my::intToDouble(const std::vector<std::vector<int>> &data){
    int n = data.size();
    std::vector<std::vector<double>> res(n,std::vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n;j++){
            res[i][j] = data[i][j];
        }
    }
    return res;
}
// 矩阵转置
vector<vector<double>> Matrix_my::transpose(const vector<vector<double>> &vec1){
    int n = vec1.size(), m = vec1[0].size();
    vector<vector<double>> ants(m,vector<double>(n,0));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            ants[j][i] = vec1[i][j];
        }
    }
    return ants;
}
// 生成一个x阶的单位阵
vector<vector<int>> Matrix_my::eye(int x){
    vector<vector<int>> ants(x,vector<int>(x,0));
    for(int i = 0; i < x; i++){
        ants[i][i] = 1;
    } 
    return ants;
}
// 生成一个x阶的单位阵,double类型
vector<vector<double>> Matrix_my::eyeD(int x){
    vector<vector<double>> ants(x,vector<double>(x,0));
    for(int i = 0; i < x; i++){
        ants[i][i] = 1.0;
    } 
    return ants;
}

// 高斯若当消元法求矩阵的逆
vector<vector<double>> Matrix_my::gauss(vector<vector<double>> data){
    int n = data.size();
    int k = 0;
    double max = 0, temp = 0;
    auto data2 = Matrix_my::typetoDouble(Matrix_my::eye(n));
    for(int i = 0; i < n; i++){
        max = data[i][i];
        k = i;
        // 寻找首元素最大的一行
        for(int j = i + 1; j < n; j++){
            if(abs(data[j][i] > abs(max))){
                max = data[j][i];
                k = j;
            }
        }
        // 将首元素最大的一行交换到i行
        if(k != i){
            for(int j = 0; j < n; j++){
                swap(data[i][j],data[k][j]);
                swap(data2[i][j],data2[k][j]);
            }
        }
        // 如果所有元素为0，说明无法转化为单位阵，无逆矩阵
        if(data[i][i] == 0){
            cout << "gauss-error:无逆矩阵！" << endl;
            return {};
        }
        // 对第i行单位化 
        temp = data[i][i];
        if(temp != 1){
            for(int j = 0; j < n; j++){
                data[i][j] = data[i][j] / temp;
                data2[i][j] = data2[i][j] / temp; 
            }
        }

        // 其他行减去第i行的倍数
        for(int j = 0; j < n; j++){
            if(j != i){
                temp = data[j][i];
                for(int z = 0; z < n; z++){
                    data[j][z] = data[j][z] - data[i][z] * temp;
                    data2[j][z] = data2[j][z] - data2[i][z] * temp;
                }
            }
        }
    }
    return data2;
}