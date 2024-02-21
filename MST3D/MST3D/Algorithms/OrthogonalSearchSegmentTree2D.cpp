#include "BaseHeader.h"
#include "OrthogonalSearchSegmentTree2D.h"


SegmentTree2D::SegmentTree2D(SegmentTreeCoreProperties* coreProperties, vector<PointWithIdx>& points) :
    coreProperties(coreProperties)
{
    sort(
        points.begin(),
        points.end(),
        [](const auto& f, const auto& s) {return f.y < s.y; }
    );
    vector<pair<int, pair<int, int>>>vpp;
    vpp.push_back({ points[0].y, {0, 0} });
    for (int i = 1; i < points.size(); i++) {
        if (vpp.back().first == points[i].y) {
            vpp.back().second.second = i;
        }
        else {
            vpp.push_back({ points[i].y, {i, i} });
        }
    }

    n = vpp.size();
    std_length = 4 * n;
    coord.resize(std_length, { -INF, INF });
    segtree1Ds.resize(std_length);
    is_end_point.resize(std_length, false);
    points_exists.resize(std_length, false);


    build(points, vpp, 1, 0, n - 1);
}

void SegmentTree2D::init_segtree1D(
    vector<PointWithIdx>& points,
    vector<pair<int, pair<int, int>>>& vpp,
    int v, int l, int r)
{
    vector<PointWithIdx> new_points(points.begin() + vpp[l].second.first, points.begin() + vpp[r].second.second + 1);
    coord[v] = { new_points.front().y, new_points.back().y };
    segtree1Ds[v] = make_unique<SegmentTree1D_V2>(coreProperties, new_points);
}

void SegmentTree2D::build(
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

    init_segtree1D(points, vpp, v, l, r);
}

void SegmentTree2D::add_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.y) {
            add_point(v * 2);
        }
        else {
            add_point(v * 2 + 1);
        }
    }

    segtree1Ds[v]->add_point(1);
    points_exists[v] = true;
    return;
}

void SegmentTree2D::remove_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.y) {
            remove_point(v * 2);
        }
        else {
            remove_point(v * 2 + 1);
        }
    }

    segtree1Ds[v]->remove_point(1);
    points_exists[v] = segtree1Ds[v]->vbb[1];
    //points_exists[v] = !segtree1Ds[v]->indices[1]->empty();
    return;
}

inline bool SegmentTree2D::is_intersect(int v) {
    return coreProperties->left.y <= coord[v].second && coreProperties->right.y >= coord[v].first;
}

inline bool SegmentTree2D::is_inner(int v) {
    return coreProperties->left.y <= coord[v].first && coreProperties->right.y >= coord[v].second;
}

inline void SegmentTree2D::get_all_points_query(int v) {
    if (is_intersect(v) && points_exists[v]) {
        get_all_points(v);
    }
}

void SegmentTree2D::get_all_points(int v) {
    if (is_inner(v)) {
        segtree1Ds[v]->get_all_points(1);
        return;
    }
    if (!is_end_point[v]) {
        get_all_points_query(v * 2);
        get_all_points_query(v * 2 + 1);
    }
}

//--------------------------------------------------------------------------------

SegmentTree2D_V2::SegmentTree2D_V2(SegmentTreeCoreProperties* coreProperties, vector<PointWithIdx>& points) :
    coreProperties(coreProperties)
{
    sort(
        points.begin(),
        points.end(),
        [](const auto& f, const auto& s) {return f.y < s.y; }
    );
    vector<pair<int, pair<int, int>>>vpp;
    vpp.push_back({ points[0].y, {0, 0} });
    for (int i = 1; i < points.size(); i++) {
        if (vpp.back().first == points[i].y) {
            vpp.back().second.second = i;
        }
        else {
            vpp.push_back({ points[i].y, {i, i} });
        }
    }

    n = vpp.size();
    std_length = 4 * n;
    coord.resize(std_length, { -INF, INF });
    segtree1Ds.resize(std_length);
    is_end_point.resize(std_length, false);
    points_exists.resize(std_length, false);


    build(points, vpp, 1, 0, n - 1);
}

void SegmentTree2D_V2::init_segtree1D(
    vector<PointWithIdx>& points,
    vector<pair<int, pair<int, int>>>& vpp,
    int v, int l, int r)
{
    vector<PointWithIdx> new_points(points.begin() + vpp[l].second.first, points.begin() + vpp[r].second.second + 1);
    segtree1Ds[v] = make_unique<SegmentTree1D_V2>(coreProperties, new_points);
}

void SegmentTree2D_V2::build(
    vector<PointWithIdx>& points,
    vector<pair<int, pair<int, int>>>& vpp,
    int v, int l, int r)
{
    coord[v] = { (points.begin() + vpp[l].second.first)->y, (points.begin() + vpp[r].second.second)->y };
    if (l >= r) {
        is_end_point[v] = true;
        init_segtree1D(points, vpp, v, l, r);
    }
    else {
        int mid = (l + r) / 2;
        build(points, vpp, v * 2, l, mid);
        build(points, vpp, v * 2 + 1, mid + 1, r);
    }
}

void SegmentTree2D_V2::add_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.y) {
            add_point(v * 2);
        }
        else {
            add_point(v * 2 + 1);
        }
    } else {
        segtree1Ds[v]->add_point(1);
    }
    points_exists[v] = true;
    return;
}

void SegmentTree2D_V2::remove_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.y) {
            remove_point(v * 2);
        }
        else {
            remove_point(v * 2 + 1);
        }
        points_exists[v] = points_exists[v * 2] || points_exists[v * 2 + 1];
    } else {
        segtree1Ds[v]->remove_point(1);
        points_exists[v] = segtree1Ds[v]->vbb[1];
    }
    return;
}

inline bool SegmentTree2D_V2::is_intersect(int v) {
    return coreProperties->left.y <= coord[v].second && coreProperties->right.y >= coord[v].first;
}

inline bool SegmentTree2D_V2::is_inner(int v) {
    return coreProperties->left.y <= coord[v].first && coreProperties->right.y >= coord[v].second;
}

inline void SegmentTree2D_V2::get_all_points_query(int v) {
    if (is_intersect(v) && points_exists[v]) {
        get_all_points(v);
    }
}

void SegmentTree2D_V2::get_all_points(int v) {
    if (is_end_point[v] && is_inner(v)) {
        segtree1Ds[v]->get_all_points(1);
        return;
    }
    if (!is_end_point[v]) {
        get_all_points_query(v * 2);
        get_all_points_query(v * 2 + 1);
    }
}