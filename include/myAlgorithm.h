#pragma once

#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <set>

namespace myAlgorithm{
    std::vector<int> getIntersection(std::unordered_set<int>& set1, std::unordered_set<int>& set2); 
    int getIntersection2(std::unordered_set<int>& set1, std::unordered_set<int>& set2);
    double sigmoid(double x);
    bool compare(const std::pair<double, int>& a, const std::pair<double, int>& b);
    std::vector<int> sortVectorDescending(const std::vector<double>& inputVector);
    int interSetNumber(const std::set<int> &a,const std::set<int> &b);
}