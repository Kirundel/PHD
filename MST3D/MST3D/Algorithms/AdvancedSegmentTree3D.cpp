#include "BaseHeader.h"
#include "AdvancedSegmentTree3D.h"
#include <cassert>

AdvancedSegmentTree3D::AdvancedSegmentTree3D(SegmentTreeCoreProperties* coreProperties, vector<PointWithIdx>& points) :
    coreProperties(coreProperties)
{
    sort(
        points.begin(),
        points.end(),
        [](const auto& f, const auto& s) {return f.x < s.x; }
    );
    vector<pair<int, pair<int, int>>>vpp;
    vpp.push_back({ points[0].x, {0, 0} });
    for (int i = 1; i < points.size(); i++) {
        if (vpp.back().first == points[i].x) {
            vpp.back().second.second = i;
        }
        else {
            vpp.push_back({ points[i].x, {i, i} });
        }
    }

    n = vpp.size();
    std_length = 4 * n;
    coord.resize(std_length, { -INF, INF });
    segtree2Ds.resize(std_length);
    is_end_point.resize(std_length, false);


    build(points, vpp, 1, 0, n - 1);
}

void AdvancedSegmentTree3D::init_segtree2D(
    vector<PointWithIdx>& points,
    vector<pair<int, pair<int, int>>>& vpp,
    int v, int l, int r)
{
    vector<PointWithIdx> new_points(points.begin() + vpp[l].second.first, points.begin() + vpp[r].second.second + 1);
    coord[v] = { new_points.front().x, new_points.back().x };
    segtree2Ds[v] = make_unique<AdvancedSegmentTree2D>(coreProperties, new_points);
}

void AdvancedSegmentTree3D::build(
    vector<PointWithIdx>& points,
    vector<pair<int, pair<int, int>>>& vpp,
    int v, int l, int r)
{
    if (l >= r) {
        is_end_point[v] = true;
    }
    else {
        int mid = (l + r) / 2;
        build(points, vpp, v * 2, l, mid);
        build(points, vpp, v * 2 + 1, mid + 1, r);
    }

    init_segtree2D(points, vpp, v, l, r);
}

void AdvancedSegmentTree3D::add_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.x) {
            add_point(v * 2);
        }
        else {
            add_point(v * 2 + 1);
        }
    }

    segtree2Ds[v]->add_point(1);
    return;
}

void AdvancedSegmentTree3D::remove_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.x) {
            remove_point(v * 2);
        }
        else {
            remove_point(v * 2 + 1);
        }
    }

    segtree2Ds[v]->remove_point(1);
    return;
}

inline bool AdvancedSegmentTree3D::is_intersect(int v) {
    return coreProperties->left.x <= coord[v].second && coreProperties->right.x >= coord[v].first;
}

inline bool AdvancedSegmentTree3D::is_inner(int v) {
    return coreProperties->left.x <= coord[v].first && coreProperties->right.x >= coord[v].second;
}

inline void AdvancedSegmentTree3D::get_all_points_query(int v) {
    if (is_intersect(v) && segtree2Ds[v]->max_point[1] >= coreProperties->left.z) {
        get_all_points(v);
    }
}

void AdvancedSegmentTree3D::get_all_points(int v) {
    if (is_inner(v)) {
        segtree2Ds[v]->get_all_points(1);
        return;
    } else {
        if (!is_end_point[v]) {
            get_all_points_query(v * 2);
            get_all_points_query(v * 2 + 1);
        }
    }
}