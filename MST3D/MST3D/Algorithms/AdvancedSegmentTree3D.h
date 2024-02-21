#pragma once

#include "AdvancedSegmentTree2D.h"

class AdvancedSegmentTree3D {
public:
    SegmentTreeCoreProperties* coreProperties;
    int n;
    int std_length;
    std::vector<std::pair<int, int>> coord;
    std::vector<std::unique_ptr<AdvancedSegmentTree2D>> segtree2Ds;
    std::vector<bool> is_end_point;

    AdvancedSegmentTree3D(SegmentTreeCoreProperties* coreProperties, std::vector<PointWithIdx>& points);

    void init_segtree2D(
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