#pragma once

#include "OrthogonalSearchSegmentTree1D.h"

class SegmentTree2D {
public:
    SegmentTreeCoreProperties* coreProperties;
    int n;
    int std_length;
    std::vector<std::pair<int, int>> coord;
    //std::vector<std::unique_ptr<SegmentTree1D>> segtree1Ds;
    std::vector<std::unique_ptr<SegmentTree1D_V2>> segtree1Ds;
    std::vector<bool> is_end_point;
    std::vector<bool> points_exists;

    SegmentTree2D(SegmentTreeCoreProperties* coreProperties, std::vector<PointWithIdx>& points);

    void init_segtree1D(
        std::vector<PointWithIdx>& points,
        std::vector<std::pair<int, std::pair<int, int>>>& vpp,
        int v, int l, int r);

    void build(
        std::vector<PointWithIdx>& points,
        std::vector<std::pair<int, std::pair<int, int>>>& vpp,
        int v, int l, int r);

    void add_point(int v);

    void remove_point(int v);

    inline bool is_intersect(int v);

    inline bool is_inner(int v);

    inline void get_all_points_query(int v);

    void get_all_points(int v);
};

class SegmentTree2D_V2 {
public:
    SegmentTreeCoreProperties* coreProperties;
    int n;
    int std_length;
    std::vector<std::pair<int, int>> coord;
    //std::vector<std::unique_ptr<SegmentTree1D>> segtree1Ds;
    std::vector<std::unique_ptr<SegmentTree1D_V2>> segtree1Ds;
    std::vector<bool> is_end_point;
    std::vector<bool> points_exists;

    SegmentTree2D_V2(SegmentTreeCoreProperties* coreProperties, std::vector<PointWithIdx>& points);

    void init_segtree1D(
        std::vector<PointWithIdx>& points,
        std::vector<std::pair<int, std::pair<int, int>>>& vpp,
        int v, int l, int r);

    void build(
        std::vector<PointWithIdx>& points,
        std::vector<std::pair<int, std::pair<int, int>>>& vpp,
        int v, int l, int r);

    void add_point(int v);

    void remove_point(int v);

    inline bool is_intersect(int v);

    inline bool is_inner(int v);

    inline void get_all_points_query(int v);

    void get_all_points(int v);
};
