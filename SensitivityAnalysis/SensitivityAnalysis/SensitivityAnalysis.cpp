using namespace std;

#define _CRT_SECURE_NO_WARNINGS
#include <set>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <cassert>  

#include <fstream>
ifstream cin("input.txt");
ofstream cout("output.txt");
//#include <iostream>

using pii = pair<int, int>;
using ll = long long;

const int INF = 2009000999;

template <typename T>
inline void sort(vector<T>& x)
{
	sort(x.begin(), x.end());
}

template <typename T>
inline void rsort(vector<T>& x)
{
	sort(x.rbegin(), x.rend());
}

template <typename T>
inline void set_max(T& x, T y)
{
	x = max(x, y);
}

template<typename T>
void print_vector(vector<T>& a) {
	for (T& element : a) {
		cout << element << " ";
	}
	cout << endl;
}

class DSU {
	vector<int> parents;
	vector<int> sizes;

public:
	DSU(int n) : parents(n), sizes(n, 1) {
		for (int i = 0; i < n; i++) {
			parents[i] = i;
		}
	}

	int get_parent(int v) {
		return parents[v] == v ? v : (parents[v] = get_parent(parents[v]));
	}

	int merge(int x, int y) {
		x = get_parent(x);
		y = get_parent(y);
		if (x == y) return -1;
		if (sizes[x] < sizes[y]) swap(x, y);

		sizes[x] += sizes[y];
		parents[y] = x;
		return x;
	}

	pii merge_with_return(int x, int y) {
		x = get_parent(x);
		y = get_parent(y);
		if (x == y) return { -1, -1 };
		if (sizes[x] < sizes[y]) swap(x, y);

		sizes[x] += sizes[y];
		parents[y] = x;
		return { x, y };
	}
};

struct edge {
	int from;
	int to;
	int weight;
	int number;
	int lower_tolerance;
	int upper_tolerance;
	bool in_mst;
	int lower_tolerance_idx;
};

struct edge_2 {
	int to;
	int weight;
	int number;
};

int find_n(int n) {
	for (int i = 0; i < 31; i++) {
		if ((1 << i) >= n) {
			return i;
		}
	}
	throw "VALUE IS NOT FOUND";
}

int n, m;
vector<edge> all_edges;
vector<vector<edge_2>> graph;
vector<vector<pii>> mst; // to, weight
vector<pii> parent; // to, edge_num
ll mst_weight = 0;

void read_graph(){
	cin >> n >> m;

	all_edges.reserve(m);
	
	for (int i = 0; i < m; i++) {
		int from, to, weight;
		cin >> from >> to >> weight;
		from--; to--;
		all_edges.push_back({ min(from, to), max(from, to), weight, i, INF, -INF, false, -1 });
	}

	graph.resize(n);
	for (int i = 0; i < m; i++) {
		const auto& e = all_edges[i];
		graph[e.from].push_back({e.to, e.weight, e.number});
		graph[e.to].push_back({ e.from, e.weight, e.number });
	}
}

void calc_mst() {
	vector<pii> sorted_edges;
	sorted_edges.reserve(m);
	for (const auto& e : all_edges) {
		sorted_edges.push_back({ e.weight, e.number });
	}

	rsort(sorted_edges);
	DSU dsu(n);
	parent.resize(n, { -1, -1 });
	mst.resize(n);

	for (const pii& p : sorted_edges) {
		const auto& e = all_edges[p.second];
		auto f = dsu.merge(e.from, e.to);
		if (f != -1) {
			all_edges[e.number].in_mst = true;
			mst[e.from].push_back({ e.to, e.weight });
			mst[e.to].push_back({ e.from, e.weight });
			mst_weight += e.weight;
		}
	}
}

vector<int> depthes;
vector<int> traverse;
vector<int> tin;
vector<int> tout;
vector<int> traverse_indexing;
int t = 0;

void dfs(int v, int prev, int cur_depth) {
	int cur_idx = t;
	traverse_indexing[cur_idx] = v;
	t++;
	depthes[v] = cur_depth;
	tin[v] = traverse.size();
	tout[v] = traverse.size();
	traverse.push_back(cur_idx);
	for (auto& edge : graph[v]) {
		if (all_edges[edge.number].in_mst && edge.to != prev) {
			parent[edge.to] = { v, edge.number };
			dfs(edge.to, v, cur_depth + 1);
			tout[v] = traverse.size();
			traverse.push_back(cur_idx);
		}
	}
}

void calc_depthes() {
	depthes.resize(n, 0);
	traverse.reserve(2 * n);
	traverse_indexing.resize(n);
	tin.resize(n);
	tout.resize(n);
	dfs(0, -1, 0);
}

vector<vector<int>> sparse_table;
vector<int> st_sz;
vector<int> delta;

void build_sparse_table() {
	st_sz.resize(traverse.size());
	for (int i = 2; i < st_sz.size(); i++) {
		st_sz[i] = st_sz[i / 2] + 1;
	}
	int st_high = st_sz.back() + 1;
	sparse_table.resize(st_high, vector<int>(traverse.size()));
	sparse_table[0] = traverse;
	delta.resize(sparse_table.size());
	for (int i = 1; i < sparse_table.size(); i++) {
		int cur_d = 1 << (i - 1);
		delta[i] = cur_d * 2 - 1;
		for (int j = cur_d; j < traverse.size(); j++) {
			sparse_table[i][j - cur_d] = min(sparse_table[i - 1][j - cur_d], sparse_table[i - 1][j]);
		}
	}

}

int lca(int f, int s) {
	int l = min(tin[f], tin[s]);
	int r = max(tout[f], tout[s]);
	int d = r - l;
	int cur_i = st_sz[d];

	return traverse_indexing[min(sparse_table[cur_i][l], sparse_table[cur_i][r - delta[cur_i]])];
}

void precalc_lower_tolerances() {
	vector<pii> sorted_edges;
	sorted_edges.reserve(m - n + 1);

	for (const auto& e : all_edges) {
		if (!e.in_mst) {
			sorted_edges.push_back({ e.weight, e.number });
		}
	}

	rsort(sorted_edges);

	DSU dsu(n);

	vector<pair<int, pii>> edges_lt;
	edges_lt.reserve(m * 2);

	vector<int> cur_min(n);
	iota(cur_min.begin(), cur_min.end(), 0);

	for (const auto& p : sorted_edges) {
		auto& e = all_edges[p.second];

		int left = dsu.get_parent(e.from);
		int right = dsu.get_parent(e.to);

		while (left != right) {
			if (depthes[cur_min[left]] < depthes[cur_min[right]]) {
				swap(left, right);
			}

			auto& par = parent[cur_min[left]];
			int cur_p = dsu.get_parent(par.first);

			left = dsu.merge(left, cur_p);

			if (left != -1) {
				all_edges[par.second].lower_tolerance = e.weight;
				all_edges[par.second].lower_tolerance_idx = e.number;
				cur_min[left] = cur_min[cur_p];
			}

			right = dsu.get_parent(right);
		}
	}
}

void precalc_upper_tolerances() {
	vector<unordered_set<int>*> elements(n);
	for (int i = 0; i < n; i++) {
		elements[i] = new unordered_set<int>();
	}

	vector<pii> sorted_edges;

	for (auto &edge: all_edges) {
		if (edge.from != edge.to && !edge.in_mst) {
			elements[edge.from]->insert(edge.number);
			elements[edge.to]->insert(edge.number);
		} else {
			sorted_edges.push_back({ edge.weight, edge.number });
		}
	}

	rsort(sorted_edges);

	DSU dsu(n);

	for (auto& [cost, edge_idx] : sorted_edges) {
		auto& edge = all_edges[edge_idx];
		auto [f, s] = dsu.merge_with_return(edge.from, edge.to);
		if (f != -1) {
			if (elements[f]->size() < elements[s]->size()) {
				swap(elements[f], elements[s]);
			}

			for (int x : *(elements[s])) {
				auto temp_iter = elements[f]->find(x);
				if (temp_iter != elements[f]->end()) {
					elements[f]->erase(temp_iter);
					all_edges[x].upper_tolerance = cost;
				}
				else {
					elements[f]->insert(x);
				}
			}
		}
	}

	for (int i = 0; i < n; i++) {
		delete elements[i];
	}
}

void answer_queries() {
	int q;
	cin >> q;
	for (int i = 0; i < q; i++) {
		int l, r, e;
		cin >> l >> r >> e;
		l--; r--; e--;

		auto& edge = all_edges[e];
		int opt_vertex = (depthes[edge.from] > depthes[edge.to]) ? edge.from : edge.to;

		bool l_l = (lca(l, opt_vertex) == opt_vertex);
		bool r_l = (lca(r, opt_vertex) == opt_vertex);

		if (l_l != r_l) {
			if (edge.in_mst) {
				if (edge.lower_tolerance == INF) {
					cout << "INF ";
				}
				else {
					int additional_component = all_edges[edge.lower_tolerance_idx].upper_tolerance == -INF ? INF : all_edges[edge.lower_tolerance_idx].upper_tolerance;
					cout << edge.weight - min(edge.lower_tolerance, additional_component) << " ";
				}
				cout << "INF ";
			}
			else {
				cout << "INF ";
				if (edge.upper_tolerance == -INF) {
					cout << "INF";
				}
				else {
					cout << edge.upper_tolerance - edge.weight;
				}
			}

			cout << "\n";
		} else {
			cout << "INF INF\n";
		}
	}
}

void print_tree() {
	for (auto& edge : all_edges) {
		if (edge.in_mst) {
			if (edge.lower_tolerance == INF) {
				cout << "-1 ";
			} else {
				cout << mst_weight - edge.weight + edge.lower_tolerance << " ";
			}
			cout << mst_weight;
		} else {
			cout << mst_weight << " ";
			if (edge.upper_tolerance == -INF) {
				cout << "-1";
			} else {
				cout << mst_weight + edge.weight - edge.upper_tolerance;
			}
		}

		cout << "\n";
	}
}

void solve() {
	read_graph();
	calc_mst();
	calc_depthes();
	precalc_lower_tolerances();
	precalc_upper_tolerances();
	build_sparse_table();
	answer_queries();
}

int main(){
	ios_base::sync_with_stdio(0);
	cin.tie(0);
	cout.tie(0);

	solve();
}