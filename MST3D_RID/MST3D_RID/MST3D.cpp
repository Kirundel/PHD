#define _CRT_SECURE_NO_WARNINGS

using namespace std;

#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <utility>
#include <set>
#include <functional>
#include <chrono>

ifstream cin("input.txt");
ofstream cout("output.txt");


const int INF = 1000 * 1000 * 1000;

struct PointWithIdx {
    int x, y, z;
    int idx;
};

inline int l1_distance(const PointWithIdx& first, const PointWithIdx& second) {
    return abs(first.x - second.x) + abs(first.y - second.y) + abs(first.z - second.z);
}

struct Edge {
    int distance;
    int first, second;

    Edge()
    {}

    Edge(const PointWithIdx& first, const PointWithIdx& second) :
        distance(l1_distance(first, second)),
        first(first.idx),
        second(second.idx)
    {}

    bool operator ==(const Edge& other) const {
        return min(first, second) == min(other.first, other.second)
            && max(first, second) == max(other.first, other.second);
    }
};

struct Answer {
    long long sum = 0;
    vector<Edge> edges;
};

class MSTAlgorithm {
public:
    virtual Answer MST(int n, vector<Edge>& edges) = 0;
};

class MST3DAlgorithm {
public:
    virtual Answer MST3D(int n, vector<PointWithIdx>& points) = 0;
};


struct SegmentTreeCoreProperties {
    PointWithIdx queryPoint;
    PointWithIdx left;
    PointWithIdx right;
    vector<int> answer;
    bool need_clear;
    long long add_nums = 0;
    long long erase_nums = 0;
    long long get_nums = 0;
};


class AdvancedSegmentTree2D {
public:
    SegmentTreeCoreProperties* coreProperties;
    int n;
    int std_length;
    vector<pair<int, int>> coord;
    vector<unique_ptr<set<pair<int, int>>>> sets;
    vector<bool> is_end_point;
    vector<int> max_point;

    AdvancedSegmentTree2D(SegmentTreeCoreProperties* coreProperties, vector<PointWithIdx>& points);

    void build(
        vector<PointWithIdx>& points,
        vector<pair<int, pair<int, int>>>& vpp,
        int v, int l, int r);

    void add_point(int v);

    void remove_point(int v);

    inline bool is_intersect(int v);

    inline bool is_inner(int v);

    inline void get_all_points_query(int v);

    void get_all_points(int v);
};


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
        }
        else {
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
            }
            else {
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

class AdvancedSegmentTree3D {
public:
    SegmentTreeCoreProperties* coreProperties;
    int n;
    int std_length;
    vector<pair<int, int>> coord;
    vector<unique_ptr<AdvancedSegmentTree2D>> segtree2Ds;
    vector<bool> is_end_point;

    AdvancedSegmentTree3D(SegmentTreeCoreProperties* coreProperties, vector<PointWithIdx>& points);

    void init_segtree2D(
        vector<PointWithIdx>& points,
        vector<pair<int, pair<int, int>>>& vpp,
        int v, int l, int r);

    void build(
        vector<PointWithIdx>& points,
        vector<pair<int, pair<int, int>>>& vpp,
        int v, int l, int r);

    void add_point(int v);

    void remove_point(int v);

    inline bool is_intersect(int v);

    inline bool is_inner(int v);

    inline void get_all_points_query(int v);

    void get_all_points(int v);
};

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
    }
    else {
        if (!is_end_point[v]) {
            get_all_points_query(v * 2);
            get_all_points_query(v * 2 + 1);
        }
    }
}

class DSU {
    vector<int> parents;
    vector<int> sizes;

public:
    DSU(int n);

    int get_parent(int v);

    bool merge(int x, int y);
};

DSU::DSU(int n) : parents(n), sizes(n, 1) {
    for (int i = 0; i < n; i++) {
        parents[i] = i;
    }
}

int DSU::get_parent(int v) {
    return parents[v] == v ? v : (parents[v] = get_parent(parents[v]));
}

bool DSU::merge(int x, int y) {
    x = get_parent(x);
    y = get_parent(y);
    if (x == y) return false;
    if (sizes[x] < sizes[y]) swap(x, y);

    sizes[x] += sizes[y];
    parents[y] = x;
    return true;
}

class KruskalAlgorithm : public MSTAlgorithm {
public:
    Answer MST(int n, vector<Edge>& edges) override;
};

Answer KruskalAlgorithm::MST(int n, vector<Edge>& edges) {
    Answer answer;
    answer.edges.reserve(n - 1);

    sort(edges.begin(), edges.end(), [](const auto& f, const auto& s) { return f.distance < s.distance; });
    edges.resize(unique(edges.begin(), edges.end()) - edges.begin());

    DSU dsu(n);
    for (auto& edge : edges) {
        if (dsu.merge(edge.first, edge.second)) {
            answer.edges.push_back(edge);
            answer.sum += edge.distance;
        }
    }

    return answer;
}

class Orthogonal3DMSTPipeline : public MST3DAlgorithm {
private:
    MSTAlgorithm& mstAlgo;

public:
    Orthogonal3DMSTPipeline(MSTAlgorithm& algorithm);

    Answer MST3D(int n, vector<PointWithIdx>& points) override;
};

Orthogonal3DMSTPipeline::Orthogonal3DMSTPipeline(MSTAlgorithm& algorithm) : mstAlgo(algorithm)
{}


vector<vector<function<PointWithIdx(PointWithIdx&)>>> new_cones = {
    {
        [](PointWithIdx& p) {return PointWithIdx{p.x, p.y, p.z - p.x - p.y, p.idx}; } ,
        [](PointWithIdx& p) {return PointWithIdx{p.y, p.z, p.x - p.y - p.z, p.idx}; } ,
        [](PointWithIdx& p) {return PointWithIdx{p.x, p.z, p.y - p.x - p.z, p.idx}; } ,
        [](PointWithIdx& p) {return PointWithIdx{p.y + p.z - p.x, p.x + p.z - p.y, p.x + p.y - p.z, p.idx}; } ,
    },
    {
        [](PointWithIdx& p) {return PointWithIdx{p.x, p.y, -p.z - p.x - p.y, p.idx}; } ,
        [](PointWithIdx& p) {return PointWithIdx{p.y, -p.z, p.x - p.y + p.z, p.idx}; } ,
        [](PointWithIdx& p) {return PointWithIdx{p.x, -p.z, p.y - p.x + p.z, p.idx}; } ,
        [](PointWithIdx& p) {return PointWithIdx{p.x + p.y + p.z, p.y - p.z - p.x, p.x + p.z - p.y, p.idx}; } ,
    },
    {
        [](PointWithIdx& p) {return PointWithIdx{p.x, -p.y, p.z - p.x + p.y, p.idx}; } ,
        [](PointWithIdx& p) {return PointWithIdx{-p.y, p.z, p.x + p.y - p.z, p.idx}; } ,
        [](PointWithIdx& p) {return PointWithIdx{p.x, p.z, -p.y - p.x - p.z, p.idx}; } ,
        [](PointWithIdx& p) {return PointWithIdx{p.x - p.y - p.z, -p.y + p.z - p.x, p.x + p.z + p.y, p.idx}; } ,
    },
    {
        [](PointWithIdx& p) {return PointWithIdx{-p.x, p.y, p.z + p.x - p.y, p.idx}; } ,
        [](PointWithIdx& p) {return PointWithIdx{-p.x, p.z, p.y + p.x - p.z, p.idx}; } ,
        [](PointWithIdx& p) {return PointWithIdx{p.y, p.z, -p.x - p.y - p.z, p.idx}; } ,
        [](PointWithIdx& p) {return PointWithIdx{p.y + p.z + p.x, -p.x + p.z - p.y, -p.x + p.y - p.z, p.idx}; } ,
    },
};

vector<function<int(PointWithIdx&)>> cones_sort = {
    [](PointWithIdx& p) {return -p.x - p.y - p.z; },
    [](PointWithIdx& p) {return (p.z - p.x - p.y); } , //reverse z
    [](PointWithIdx& p) {return (p.y - p.x - p.z); } , //reverse y
    [](PointWithIdx& p) {return (p.x - p.y - p.z); } , //reverse x
};

vector<PointWithIdx> generate_translation(vector<PointWithIdx>& points, int x, int y, int z) {
    vector<PointWithIdx> result = points;
    for (int i = 0; i < points.size(); i++) {
        result[i].x *= x;
        result[i].y *= y;
        result[i].z *= z;
    }
    return result;
}

PointWithIdx ApplyTransitions(PointWithIdx point, int x, int y, int z, function<PointWithIdx(PointWithIdx&)>& cone) {
    point.x *= x;
    point.y *= y;
    point.z *= z;
    point = cone(point);
    return point;
}

vector<PointWithIdx> translate_via_cone(vector< PointWithIdx>& points, function<PointWithIdx(PointWithIdx&)>& cone) {
    vector<PointWithIdx> result = points;
    for (int i = 0; i < points.size(); i++) {
        result[i] = cone(result[i]);
    }
    return result;
}


Answer Orthogonal3DMSTPipeline::MST3D(int n, vector<PointWithIdx>& points)
{
    vector<Edge> edges;

    int cone_idx = 0;

    auto pipeline_begin = chrono::steady_clock::now();

    for (int m1 = -1; m1 < 2; m1 += 2) {
        for (int m2 = -1; m2 < 2; m2 += 2) {
            for (int m3 = -1; m3 < 2; m3 += 2) {
                vector<PointWithIdx> cur_points_prev = generate_translation(points, m1, m2, m3);

                for (int sc = 0; sc < 1; sc++)
                {
                    auto& sort_criterium = cones_sort[sc];
                    auto cur_points = cur_points_prev;
                    sort(cur_points.begin(), cur_points.end(),
                        [&](auto& first, auto& second)
                        { return sort_criterium(first) < sort_criterium(second); }
                    );


                    for (int cur_cone = 0; cur_cone < new_cones[sc].size(); cur_cone++) {
                        auto& cone = new_cones[sc][cur_cone];

                        int number_of_points = 0;
                        auto points_via_cone = cur_points;

                        points_via_cone = translate_via_cone(points_via_cone, cone);

                        SegmentTreeCoreProperties* coreProperties = NULL;
                        auto beg = chrono::steady_clock::now();
                        AdvancedSegmentTree3D* tree = NULL;
                        {
                            auto pointsCopy = points_via_cone;
                            coreProperties = new SegmentTreeCoreProperties();
                            tree = new AdvancedSegmentTree3D(coreProperties, pointsCopy);
                        }
                        long long tree_init = chrono::duration_cast<chrono::nanoseconds>(chrono::steady_clock::now() - beg).count();

                        long long add_points = 0;
                        long long get_points = 0;
                        long long remove_points = 0;

                        coreProperties->right = PointWithIdx{ .x = INF, .y = INF, .z = INF, .idx = -1 };

                        for (int i = 0; i < n; i++) {
                            auto& point = points_via_cone[i];

                            coreProperties->left = point;
                            coreProperties->answer.clear();

                            coreProperties->queryPoint = point;
                            tree->get_all_points(1);

                            number_of_points += coreProperties->answer.size();

                            for (auto pointIdx : coreProperties->answer) {
                                edges.push_back(Edge(points[point.idx], points[pointIdx]));
                                coreProperties->queryPoint = ApplyTransitions(points[pointIdx], m1, m2, m3, cone);

                                tree->remove_point(1);
                            }
                            coreProperties->queryPoint = point;
                            tree->add_point(1);
                        }

                        cone_idx++;

                        delete coreProperties;
                        delete tree;

                    }
                }
            }
        }
    }

    return mstAlgo.MST(n, edges);
}

void solve() {
    int n;
    cin >> n;

    vector<PointWithIdx> arr(n);

    for (int i = 0; i < n; i++) {
        auto& point = arr[i];
        cin >> point.x >> point.y >> point.z;
        point.idx = i;
    }

    KruskalAlgorithm kruskal;
    Orthogonal3DMSTPipeline orthogonal(kruskal);
    auto result = orthogonal.MST3D(n, arr);
    cout << result.sum << "\n";
    for (auto& edge : result.edges) {
        cout << edge.first + 1 << " " << edge.second + 1 << "\n";
    }
}

int main() {
    solve();
}