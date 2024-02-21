#include "BaseHeader.h"
#include "AdvancedSegmentTree2D.h"


AdvancedSegmentTree2D::AdvancedSegmentTree2D(SegmentTreeCoreProperties* coreProperties, vector<PointWithIdx>& points) :
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
    sets.resize(std_length);
    is_end_point.resize(std_length, false);
    max_point.resize(std_length, -INF);


    build(points, vpp, 1, 0, n - 1);
}

void AdvancedSegmentTree2D::build(
    vector<PointWithIdx>& points,
    vector<pair<int, pair<int, int>>>& vpp,
    int v, int l, int r)
{
    coord[v] = { (points.begin() + vpp[l].second.first)->y, (points.begin() + vpp[r].second.second)->y };
    if (l >= r) {
        is_end_point[v] = true;
        sets[v] = make_unique<set<pair<int, int>>>();
    }
    else {
        int mid = (l + r) / 2;
        build(points, vpp, v * 2, l, mid);
        build(points, vpp, v * 2 + 1, mid + 1, r);
    }
}

void AdvancedSegmentTree2D::add_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.y) {
            add_point(v * 2);
        }
        else {
            add_point(v * 2 + 1);
        }
        max_point[v] = max(max_point[v * 2], max_point[v * 2 + 1]);
    }
    else {
        sets[v]->insert({ coreProperties->queryPoint.z, coreProperties->queryPoint.idx });
        max_point[v] = max(max_point[v], coreProperties->queryPoint.z);
    }
    return;
}

void AdvancedSegmentTree2D::remove_point(int v) {
    if (!is_end_point[v]) {
        if (coord[v * 2].second >= coreProperties->queryPoint.y) {
            remove_point(v * 2);
        }
        else {
            remove_point(v * 2 + 1);
        }
        max_point[v] = max(max_point[v * 2], max_point[v * 2 + 1]);
    }
    else {
        sets[v]->erase({ coreProperties->queryPoint.z, coreProperties->queryPoint.idx });
        if (sets[v]->size() > 0) {
            max_point[v] = sets[v]->rbegin()->first;
        } else {
            max_point[v] = -INF;
        }
    }
    return;
}

inline bool AdvancedSegmentTree2D::is_intersect(int v) {
    return coreProperties->left.y <= coord[v].second && coreProperties->right.y >= coord[v].first;
}

inline bool AdvancedSegmentTree2D::is_inner(int v) {
    return coreProperties->left.y <= coord[v].first && coreProperties->right.y >= coord[v].second;
}

inline void AdvancedSegmentTree2D::get_all_points_query(int v) {
    if (is_intersect(v) && max_point[v] >= coreProperties->left.z) {
        get_all_points(v);
    }
}

void AdvancedSegmentTree2D::get_all_points(int v) {
    if (is_end_point[v] && is_inner(v)) {
        for (auto iter = sets[v]->rbegin(); iter != sets[v]->rend(); iter++) {
            if (iter->first >= coreProperties->left.z) {
                coreProperties->answer.push_back(iter->second);
            } else {
                break;
            }
        }
        return;
    }
    if (!is_end_point[v]) {
        get_all_points_query(v * 2);
        get_all_points_query(v * 2 + 1);
        max_point[v] = max(max_point[v * 2], max_point[v * 2 + 1]);
    }
}