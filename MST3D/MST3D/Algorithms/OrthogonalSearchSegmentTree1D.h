#pragma once

#include <vector>
#include <algorithm>
//#include <iomanip>
#include <map>
#include <utility>
#include <string>
#include <complex>
#include <unordered_set>
#include <set>
#include "BaseClasses.h"

//using namespace std;


//struct Point {
//};

struct SegmentTreeCoreProperties {
    PointWithIdx queryPoint;
    PointWithIdx left;
    PointWithIdx right;
    std::vector<int> answer;
    bool need_clear;
    long long add_nums = 0;
    long long erase_nums = 0;
    long long get_nums = 0;
};


inline constexpr int VAL = 1000 * 1000 * 1000;


class SegmentTree1D {
public:
    SegmentTreeCoreProperties* coreProperties;
    const int n;
    const int std_length;
    std::vector<std::pair<int, int>> coord;
    std::vector<std::unique_ptr<std::unordered_set<int>>> indices;
    //std::vector<std::unique_ptr<std::set<int>>> indices;
    std::vector<bool> is_end_point;

    SegmentTree1D(SegmentTreeCoreProperties* coreProperties, std::vector<PointWithIdx>& points);


    void build(std::vector<int>& points, int v, int l, int r);

    void add_point(int v);

    void remove_point(int v);

    inline bool is_intersect(int v);

    inline bool is_inner(int v);

    inline void get_all_points_query(int v);

    void get_all_points(int v);
};


class SegmentTree1D_V2 {
public:
    SegmentTreeCoreProperties* coreProperties;
    const int n;
    const int std_length;
    std::vector<std::pair<int, int>> coord;
    std::vector<std::unique_ptr<std::unordered_set<int>>> indices;
    std::vector<bool> vbb;
    //std::vector<std::unique_ptr<std::set<int>>> indices;
    std::vector<bool> is_end_point;

    SegmentTree1D_V2(SegmentTreeCoreProperties* coreProperties, std::vector<PointWithIdx>& points);


    void build(std::vector<int>& points, int v, int l, int r);

    void add_point(int v);

    void remove_point(int v);

    inline bool is_intersect(int v);

    inline bool is_inner(int v);

    inline void get_all_points_query(int v);

    void get_all_points(int v);
};