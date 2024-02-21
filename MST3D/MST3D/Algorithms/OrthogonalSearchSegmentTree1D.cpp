#include "BaseHeader.h"
#include "OrthogonalSearchSegmentTree1D.h"
#include <chrono>

SegmentTree1D::SegmentTree1D(SegmentTreeCoreProperties* coreProperties, std::vector<PointWithIdx>& points) :
    coreProperties(coreProperties),
    n(points.size()),
    std_length(n * 4),
    coord(std_length, { -INF, INF }),
    indices(std_length),
    is_end_point(std_length, false)
{
    //        sort(points.begin(), points.end(), [](auto& first, auto& second){return first.z < second.z;});
    std::vector<int> coordinates_split;
    for (const auto& point : points) {
        coordinates_split.push_back(point.z);
    }
    std::sort(coordinates_split.begin(), coordinates_split.end());
    std::unique(coordinates_split.begin(), coordinates_split.end());

    build(coordinates_split, 1, 0, n - 1);
}

void SegmentTree1D::build(std::vector<int>& points, int v, int l, int r) {
    indices[v] = make_unique<unordered_set<int>>();
    //indices[v] = make_unique<set<int>>();
    if (l >= r) {
        coord[v] = { points[l], points[l] };
        is_end_point[v] = true;
        return;
    }
    int mid = (l + r) / 2;
    build(points, v * 2, l, mid);
    build(points, v * 2 + 1, mid + 1, r);
    coord[v] = { coord[v * 2].first, coord[v * 2 + 1].second };
}

void SegmentTree1D::add_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.z) {
            add_point(v * 2);
        }
        else {
            add_point(v * 2 + 1);
        }
    }

    indices[v]->insert(coreProperties->queryPoint.idx);
    //if (is_end_point[v]) 
    {
        coreProperties->add_nums++;
        coreProperties->erase_nums += indices[v]->size();
    }
    return;
}

void SegmentTree1D::remove_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.z) {
            remove_point(v * 2);
        }
        else {
            remove_point(v * 2 + 1);
        }
    }

    indices[v]->erase(coreProperties->queryPoint.idx);
    return;
}

inline bool SegmentTree1D::is_intersect(int v) {
    return coreProperties->left.z <= coord[v].second && coreProperties->right.z >= coord[v].first;
}

inline bool SegmentTree1D::is_inner(int v) {
    return coreProperties->left.z <= coord[v].first && coreProperties->right.z >= coord[v].second;
}

inline void SegmentTree1D::get_all_points_query(int v) {
    if (is_intersect(v) && indices[v]->size() > 0) {
        get_all_points(v);
    }
}

void SegmentTree1D::get_all_points(int v) {
    if (is_inner(v)) {
        for (int idx : *indices[v]) {
            coreProperties->answer.push_back(idx);
        }
        return;
    }
    if (!is_end_point[v])
    {
        get_all_points_query(v * 2);
        get_all_points_query(v * 2 + 1);
    }
}

//--------------------------------------------------------------------------------


SegmentTree1D_V2::SegmentTree1D_V2(SegmentTreeCoreProperties* coreProperties, std::vector<PointWithIdx>& points) :
    coreProperties(coreProperties),
    n(points.size()),
    std_length(n * 4),
    coord(std_length, { -INF, INF }),
    indices(std_length),
    is_end_point(std_length, false),
    vbb(std_length, false)
{
    //        sort(points.begin(), points.end(), [](auto& first, auto& second){return first.z < second.z;});
    std::vector<int> coordinates_split;
    for (const auto& point : points) {
        coordinates_split.push_back(point.z);
    }
    std::sort(coordinates_split.begin(), coordinates_split.end());
    std::unique(coordinates_split.begin(), coordinates_split.end());

    build(coordinates_split, 1, 0, n - 1);
}

void SegmentTree1D_V2::build(std::vector<int>& points, int v, int l, int r) {
    if (l >= r) {
        coord[v] = { points[l], points[l] };
        is_end_point[v] = true;
        indices[v] = make_unique<unordered_set<int>>();
        return;
    }
    int mid = (l + r) / 2;
    build(points, v * 2, l, mid);
    build(points, v * 2 + 1, mid + 1, r);
    coord[v] = { coord[v * 2].first, coord[v * 2 + 1].second };
}

void SegmentTree1D_V2::add_point(int v) {
    //coreProperties->sum_tt++;
    vbb[v] = true;
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.z) {
            add_point(v * 2);
        }
        else {
            add_point(v * 2 + 1);
        }
    } else {
        indices[v]->insert(coreProperties->queryPoint.idx);
        coreProperties->add_nums++;
    }
    return;
}

void SegmentTree1D_V2::remove_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.z) {
            remove_point(v * 2);
        }
        else {
            remove_point(v * 2 + 1);
        }
        vbb[v] = vbb[v * 2] || vbb[v * 2 + 1];
    } else {
        indices[v]->erase(coreProperties->queryPoint.idx);
        vbb[v] = !indices[v]->empty();
        coreProperties->erase_nums++;
    }
    return;
}

inline bool SegmentTree1D_V2::is_intersect(int v) {
    return coreProperties->left.z <= coord[v].second && coreProperties->right.z >= coord[v].first;
}

inline bool SegmentTree1D_V2::is_inner(int v) {
    return coreProperties->left.z <= coord[v].first && coreProperties->right.z >= coord[v].second;
}

inline void SegmentTree1D_V2::get_all_points_query(int v) {
    if (is_intersect(v) && vbb[v]) {
        get_all_points(v);
    }
}

void SegmentTree1D_V2::get_all_points(int v) {
    if (is_end_point[v]) {
        if (is_inner(v))
        {
            coreProperties->get_nums++;
            for (int idx : *indices[v]) {
                coreProperties->answer.push_back(idx);
            }
        }
        return;
    } else {
        get_all_points_query(v * 2);
        get_all_points_query(v * 2 + 1);
    }
}