#include "BaseHeader.h"
#include "OrthogonalSearchSegmentTree3D.h"
#include <cassert>


SegmentTree3D::SegmentTree3D(SegmentTreeCoreProperties* coreProperties, vector<PointWithIdx>& points) :
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
    points_exists.resize(std_length, false);


    build(points, vpp, 1, 0, n - 1);
}

void SegmentTree3D::init_segtree2D(
    vector<PointWithIdx>& points,
    vector<pair<int, pair<int, int>>>& vpp,
    int v, int l, int r)
{
    vector<PointWithIdx> new_points(points.begin() + vpp[l].second.first, points.begin() + vpp[r].second.second + 1);
    coord[v] = { new_points.front().x, new_points.back().x };
    segtree2Ds[v] = make_unique<SegmentTree2D>(coreProperties, new_points);
}

void SegmentTree3D::build(
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

void SegmentTree3D::add_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.x) {
            add_point(v * 2);
        }
        else {
            add_point(v * 2 + 1);
        }
    }

    segtree2Ds[v]->add_point(1);
    points_exists[v] = true;
    return;
}

void SegmentTree3D::remove_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.x) {
            remove_point(v * 2);
        }
        else {
            remove_point(v * 2 + 1);
        }
    }

    segtree2Ds[v]->remove_point(1);
    points_exists[v] = segtree2Ds[v]->points_exists[1];
    return;
}

inline bool SegmentTree3D::is_intersect(int v) {
    return coreProperties->left.x <= coord[v].second && coreProperties->right.x >= coord[v].first;
}

inline bool SegmentTree3D::is_inner(int v) {
    return coreProperties->left.x <= coord[v].first && coreProperties->right.x >= coord[v].second;
}

inline void SegmentTree3D::get_all_points_query(int v) {
    if (is_intersect(v) && points_exists[v]) {
        get_all_points(v);
    }
}

void SegmentTree3D::get_all_points(int v) {
    if (is_inner(v)) {
        segtree2Ds[v]->get_all_points(1);
        return;
    }
    if (!is_end_point[v]) {
        get_all_points_query(v * 2);
        get_all_points_query(v * 2 + 1);
    }
}

//--------------------------------------------------------------------------------

OrthogonalSearch3DDummy::OrthogonalSearch3DDummy(SegmentTreeCoreProperties* coreProperties) :
    coreProperties(coreProperties)
{}

void OrthogonalSearch3DDummy::add_point(int v) {
    points.push_back(coreProperties->queryPoint);
}

void OrthogonalSearch3DDummy::remove_point(int v) {
    int idx = -1;
    for (int i = 0; i < points.size(); i++) {
        if (points[i].idx == coreProperties->queryPoint.idx) {
            idx = i;
            break;
        }
    }
    assert(idx != -1);

    points.erase(points.begin() + idx);
}

bool OrthogonalSearch3DDummy::check_point(PointWithIdx& point) {
    return coreProperties->left.x <= point.x && point.x <= coreProperties->right.x
        && coreProperties->left.y <= point.y && point.y <= coreProperties->right.y
        && coreProperties->left.z <= point.z && point.z <= coreProperties->right.z;
}


void OrthogonalSearch3DDummy::get_all_points(int v) {
    coreProperties->answer.clear();

    for (auto& point : points) {
        if (check_point(point)) {
            coreProperties->answer.push_back(point.idx);
        }
    }
}

//--------------------------------------------------------------------------------

SegmentTree3D_V2::SegmentTree3D_V2(SegmentTreeCoreProperties* coreProperties, vector<PointWithIdx>& points) :
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
    points_exists.resize(std_length, false);


    build(points, vpp, 1, 0, n - 1);
}

void SegmentTree3D_V2::init_segtree2D(
    vector<PointWithIdx>& points,
    vector<pair<int, pair<int, int>>>& vpp,
    int v, int l, int r)
{
    vector<PointWithIdx> new_points(points.begin() + vpp[l].second.first, points.begin() + vpp[r].second.second + 1);
    segtree2Ds[v] = make_unique<SegmentTree2D_V2>(coreProperties, new_points);
}

void SegmentTree3D_V2::build(
    vector<PointWithIdx>& points,
    vector<pair<int, pair<int, int>>>& vpp,
    int v, int l, int r)
{
    coord[v] = { (points.begin() + vpp[l].second.first)->x, (points.begin() + vpp[r].second.second)->x };
    if (l >= r) {
        is_end_point[v] = true;
        init_segtree2D(points, vpp, v, l, r);
    }
    else {
        int mid = (l + r) / 2;
        build(points, vpp, v * 2, l, mid);
        build(points, vpp, v * 2 + 1, mid + 1, r);
    }
}

void SegmentTree3D_V2::add_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.x) {
            add_point(v * 2);
        }
        else {
            add_point(v * 2 + 1);
        }
    }
    else {
        segtree2Ds[v]->add_point(1);
    }
    points_exists[v] = true;
    return;
}

void SegmentTree3D_V2::remove_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.x) {
            remove_point(v * 2);
        }
        else {
            remove_point(v * 2 + 1);
        }
        points_exists[v] = points_exists[v * 2] || points_exists[v * 2 + 1];
    }
    else {
        segtree2Ds[v]->remove_point(1);
        points_exists[v] = segtree2Ds[v]->points_exists[1];
    }
    return;
}

inline bool SegmentTree3D_V2::is_intersect(int v) {
    return coreProperties->left.x <= coord[v].second && coreProperties->right.x >= coord[v].first;
}

inline bool SegmentTree3D_V2::is_inner(int v) {
    return coreProperties->left.x <= coord[v].first && coreProperties->right.x >= coord[v].second;
}

inline void SegmentTree3D_V2::get_all_points_query(int v) {
    if (is_intersect(v) && points_exists[v]) {
        get_all_points(v);
    }
}

void SegmentTree3D_V2::get_all_points(int v) {
    if (is_end_point[v] && is_inner(v)) {
        segtree2Ds[v]->get_all_points(1);
        return;
    }
    if (!is_end_point[v]) {
        get_all_points_query(v * 2);
        get_all_points_query(v * 2 + 1);
    }
}

//--------------------------------------------------------------------------------

SegmentTree3D_V3::SegmentTree3D_V3(SegmentTreeCoreProperties* coreProperties, vector<PointWithIdx>& points) :
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
    points_exists.resize(std_length, false);


    build(points, vpp, 1, 0, n - 1);
}

void SegmentTree3D_V3::init_segtree2D(
    vector<PointWithIdx>& points,
    vector<pair<int, pair<int, int>>>& vpp,
    int v, int l, int r)
{
    vector<PointWithIdx> new_points(points.begin() + vpp[l].second.first, points.begin() + vpp[r].second.second + 1);
    coord[v] = { new_points.front().x, new_points.back().x };
    segtree2Ds[v] = make_unique<SegmentTree2D_V2>(coreProperties, new_points);
}

void SegmentTree3D_V3::build(
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

void SegmentTree3D_V3::add_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.x) {
            add_point(v * 2);
        }
        else {
            add_point(v * 2 + 1);
        }
    }

    segtree2Ds[v]->add_point(1);
    points_exists[v] = true;
    return;
}

void SegmentTree3D_V3::remove_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.x) {
            remove_point(v * 2);
        }
        else {
            remove_point(v * 2 + 1);
        }
    }

    segtree2Ds[v]->remove_point(1);
    points_exists[v] = segtree2Ds[v]->points_exists[1];
    return;
}

inline bool SegmentTree3D_V3::is_intersect(int v) {
    return coreProperties->left.x <= coord[v].second && coreProperties->right.x >= coord[v].first;
}

inline bool SegmentTree3D_V3::is_inner(int v) {
    return coreProperties->left.x <= coord[v].first && coreProperties->right.x >= coord[v].second;
}

inline void SegmentTree3D_V3::get_all_points_query(int v) {
    if (is_intersect(v) && points_exists[v]) {
        get_all_points(v);
    }
}

void SegmentTree3D_V3::get_all_points(int v) {
    if (is_inner(v)) {
        segtree2Ds[v]->get_all_points(1);
        return;
    }
    if (!is_end_point[v]) {
        get_all_points_query(v * 2);
        get_all_points_query(v * 2 + 1);
    }
}

//--------------------------------------------------------------------------------

SegmentTree3D_Dummy::SegmentTree3D_Dummy(SegmentTreeCoreProperties* coreProperties, vector<PointWithIdx>& points) :
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
    points_exists.resize(std_length, false);


    build(points, vpp, 1, 0, n - 1);
}

void SegmentTree3D_Dummy::init_segtree2D(
    vector<PointWithIdx>& points,
    vector<pair<int, pair<int, int>>>& vpp,
    int v, int l, int r)
{
    segtree2Ds[v] = make_unique<OrthogonalSearch3DDummy>(coreProperties);
}

void SegmentTree3D_Dummy::build(
    vector<PointWithIdx>& points,
    vector<pair<int, pair<int, int>>>& vpp,
    int v, int l, int r)
{
    coord[v] = { (points.begin() + vpp[l].second.first)->x, (points.begin() + vpp[r].second.second)->x };
    if (l >= r) {
        is_end_point[v] = true;
        init_segtree2D(points, vpp, v, l, r);
    }
    else {
        int mid = (l + r) / 2;
        build(points, vpp, v * 2, l, mid);
        build(points, vpp, v * 2 + 1, mid + 1, r);
    }
}

void SegmentTree3D_Dummy::add_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.x) {
            add_point(v * 2);
        }
        else {
            add_point(v * 2 + 1);
        }
    }
    else {
        segtree2Ds[v]->add_point(1);
    }
    points_exists[v] = true;
    return;
}

void SegmentTree3D_Dummy::remove_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.x) {
            remove_point(v * 2);
        }
        else {
            remove_point(v * 2 + 1);
        }
        points_exists[v] = points_exists[v * 2] || points_exists[v * 2 + 1];
    }
    else {
        segtree2Ds[v]->remove_point(1);
        points_exists[v] = !segtree2Ds[v]->points.empty();
    }
    return;
}

inline bool SegmentTree3D_Dummy::is_intersect(int v) {
    return coreProperties->left.x <= coord[v].second && coreProperties->right.x >= coord[v].first;
}

inline bool SegmentTree3D_Dummy::is_inner(int v) {
    return coreProperties->left.x <= coord[v].first && coreProperties->right.x >= coord[v].second;
}

inline void SegmentTree3D_Dummy::get_all_points_query(int v) {
    if (is_intersect(v) && points_exists[v]) {
        get_all_points(v);
    }
}

void SegmentTree3D_Dummy::get_all_points(int v) {
    if (is_end_point[v] && is_inner(v)) {
        segtree2Ds[v]->get_all_points(1);
        return;
    }
    if (!is_end_point[v]) {
        get_all_points_query(v * 2);
        get_all_points_query(v * 2 + 1);
    }
}