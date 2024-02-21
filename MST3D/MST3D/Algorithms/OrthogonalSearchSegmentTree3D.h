#pragma once

#include "OrthogonalSearchSegmentTree2D.h"

class SegmentTree3D {
public:
    SegmentTreeCoreProperties* coreProperties;
    int n;
    int std_length;
    std::vector<std::pair<int, int>> coord;
    std::vector<std::unique_ptr<SegmentTree2D>> segtree2Ds;
    std::vector<bool> is_end_point;
    std::vector<bool> points_exists;

    SegmentTree3D(SegmentTreeCoreProperties* coreProperties, std::vector<PointWithIdx>& points);

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

class OrthogonalSearch3DDummy {
public:
    SegmentTreeCoreProperties* coreProperties;
    std::vector<PointWithIdx> points;

    OrthogonalSearch3DDummy(SegmentTreeCoreProperties* coreProperties);

    void add_point(int v);
    void remove_point(int v);
    bool check_point(PointWithIdx& point);
    void get_all_points(int v);
};

class SegmentTree3D_V2 {
public:
    SegmentTreeCoreProperties* coreProperties;
    int n;
    int std_length;
    std::vector<std::pair<int, int>> coord;
    std::vector<std::unique_ptr<SegmentTree2D_V2>> segtree2Ds;
    std::vector<bool> is_end_point;
    std::vector<bool> points_exists;

    SegmentTree3D_V2(SegmentTreeCoreProperties* coreProperties, std::vector<PointWithIdx>& points);

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

class SegmentTree3D_V3 {
public:
    SegmentTreeCoreProperties* coreProperties;
    int n;
    int std_length;
    std::vector<std::pair<int, int>> coord;
    std::vector<std::unique_ptr<SegmentTree2D_V2>> segtree2Ds;
    std::vector<bool> is_end_point;
    std::vector<bool> points_exists;

    SegmentTree3D_V3(SegmentTreeCoreProperties* coreProperties, std::vector<PointWithIdx>& points);

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

class SegmentTree3D_Dummy {
public:
    SegmentTreeCoreProperties* coreProperties;
    int n;
    int std_length;
    std::vector<std::pair<int, int>> coord;
    std::vector<std::unique_ptr<OrthogonalSearch3DDummy>> segtree2Ds;
    std::vector<bool> is_end_point;
    std::vector<bool> points_exists;

    SegmentTree3D_Dummy(SegmentTreeCoreProperties* coreProperties, std::vector<PointWithIdx>& points);

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