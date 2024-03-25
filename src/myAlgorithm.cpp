#include "myAlgorithm.h"

using namespace std;
// 求两个集合的交集并返回向量
std::vector<int> myAlgorithm::getIntersection(std::unordered_set<int>& set1, std::unordered_set<int>& set2){
    if (set1.size() > set2.size()) {
        return getIntersection(set2, set1);
    }
    std::vector<int> intersection;
    for (auto& num : set1) {
        if (set2.count(num)) {
            intersection.push_back(num);
        }
    }
    return intersection;
}
// 求两个集合交集的数量
int myAlgorithm::getIntersection2(std::unordered_set<int>& set1, std::unordered_set<int>& set2){
    if (set1.size() > set2.size()) {
        return getIntersection2(set2, set1);
    }
    int res = 0;
    std::vector<int> intersection;
    for (auto& num : set1) {
        if (set2.count(num)) {
            res++;
        }
    }
    return res;
}

double myAlgorithm::sigmoid(double x){
    return 1 / (1 + std::pow(M_E,-x));
}

// 自定义比较函数，用于降序排序
bool myAlgorithm::compare(const std::pair<double, int>& a, const std::pair<double, int>& b) {
    return a.first > b.first;
}

std::vector<int> myAlgorithm::sortVectorDescending(const std::vector<double>& inputVector) {
    std::vector<std::pair<double, int>> indexedVector;
    std::vector<int> sortedIndices;

    // 构建带有下标的 vector
    for (int i = 0; i < inputVector.size(); i++) {
        indexedVector.push_back(std::make_pair(inputVector[i], i));
    }

    // 使用自定义比较函数对 vector 进行排序
    std::sort(indexedVector.begin(), indexedVector.end(), compare);

    // 提取排序后的下标
    for (const auto& pair : indexedVector) {
        sortedIndices.push_back(pair.second);
    }

    return sortedIndices;
}

// 计算两个交集的数量
int myAlgorithm::interSetNumber(const set<int> &set1,const set<int> &set2){
    std::set<int> intersection;
    // 使用 set_intersection 计算交集
    std::set_intersection(set1.begin(), set1.end(),
                          set2.begin(), set2.end(),
                          std::inserter(intersection, intersection.begin()));

    // 计算交集数量
    int intersectionCount = std::distance(intersection.begin(), intersection.end());
    return intersectionCount;
}