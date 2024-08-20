using namespace std;

#define _CRT_SECURE_NO_WARNINGS
#include <set>
#include <vector>
#include <algorithm>
#include <math.h>
#include <cassert>  

#include <fstream>
//#include <iostream>
ifstream cin("input.txt");
ofstream cout("output.txt");


using pii = pair<int, int>;

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
		if (element == 0) {
			cout << 0 << " ";
		} else {
			cout << INF - element << " ";
		}
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

int find_n(int n) {
	for (int i = 0; i < 31; i++) {
		if ((1 << i) >= n) {
			return i;
		}
	}
	throw "VALUE IS NOT FOUND";
}

class LCA {
private:
	const int N;
	const int root;
	vector<vector<pii>>& tree;
	vector<vector<int>> parent;
	vector<vector<int>> result;
	vector<int>h;
	int answer = -1;

	void build(int v, int prev, int edge_cost) {
		h[v] = h[prev] + 1;
		parent[v][0] = prev;
		result[v][0] = edge_cost;
		for (int idx = 1; idx < N; idx++) {
			int prev_parent = parent[v][idx - 1];
			parent[v][idx] = parent[prev_parent][idx - 1];
			result[v][idx] = max(result[v][idx - 1], result[prev_parent][idx - 1]);
		}

		for (auto& [edge, cur_cost]: tree[v]) {
			if (edge != prev) {
				build(edge, v, cur_cost);
			}
		}
	}

public:
	LCA(vector<vector<pii>>& tree, int root) :
		tree(tree),
		N(find_n(tree.size()) + 2),
		parent(tree.size(), vector<int>(N, root)),
		result(tree.size(), vector<int>(N, 0)),
		h(tree.size(), 0),
		root(root)
	{
		build(root, root, 0);
	}

	int get_answer(int x, int y) {
		answer = 0;
		if (h[x] != h[y]) {
			if (h[x] > h[y]) swap(x, y);

			int cur_c = N - 1;
			while (cur_c >= 0 && h[y] > h[x]) {
				int edge = parent[y][cur_c];
				if (h[edge] >= h[x]) {
					set_max(answer, result[y][cur_c]);
					y = edge;
				}
				cur_c--;
			}
		}

		if (x == y) return answer;

		{
			int cur_c = N - 1;
			while (cur_c >= 0) {
				if (parent[x][cur_c] != parent[y][cur_c]) {
					set_max(answer, result[x][cur_c]);
					set_max(answer, result[y][cur_c]);
					x = parent[x][cur_c];
					y = parent[y][cur_c];
				}
				cur_c--;
			}
			return max(answer, max(result[x][0], result[y][0]));
		}
	}
};

int algo_type;
int n, m;
vector<vector<pii>> graph;
vector<pair<int, pii>> edges;
int q;
vector<pii> queries;

void read() {
	cin >> algo_type;
	cin >> n >> m;
	graph.resize(n);
	for (int i = 0; i < m; i++) {
		int a, b, c;
		cin >> a >> b >> c;
		c = INF - c;
		a--; b--;
		graph[a].push_back({ b, c });
		graph[b].push_back({ a, c });
		edges.push_back({ c, { a, b } });
	}
	
	cin >> q;
	for (int i = 0; i < q; i++) {
		int a, b;
		cin >> a >> b;
		a--; b--;
		queries.push_back({a, b});
	}
}

int dijkstra_from_point_to_point(int x, int y) {
	vector<int> distances(n, INF);
	set<pii> s;
	distances[x] = 0;
	s.insert({ 0, x });

	while (!s.empty()) {
		auto [dist, from] = *s.begin();
		s.erase(s.begin());

		for (auto [to, cost] : graph[from]) {
			int cur_cost = max(cost, dist);
			if (distances[to] > cur_cost) {
				s.erase({ distances[to], to });
				s.insert({ cur_cost, to });
				distances[to] = cur_cost;
			}
		}
	}

	return distances[y];
}

vector<int> solve_dijkstra() {
	vector<int> answer(q);

	for (int i = 0; i < q; i++) {
		answer[i] = dijkstra_from_point_to_point(queries[i].first, queries[i].second);
	}
	return answer;
}

vector<int> solve_merge_trees() {
	vector<set<int>*> elements(n);
	for (int i = 0; i < n; i++) {
		elements[i] = new set<int>();
	}

	vector<int> answers(q, -1);

	for (int idx = 0; idx < q; idx++) {
		auto& edge = queries[idx];
		if (edge.first != edge.second)
		{
			elements[edge.first]->insert(idx);
			elements[edge.second]->insert(idx);
		} else {
			answers[idx] = 0;
		}
	}

	auto sorted_edges = edges;
	sort(sorted_edges);

	DSU dsu(n);

	for (auto& [cost, edge] : sorted_edges) {
		auto [f, s] = dsu.merge_with_return(edge.first, edge.second);
		if (f != -1) {
			if (elements[f]->size() < elements[s]->size()) {
				swap(elements[f], elements[s]);
			}

			for (int x : *(elements[s])) {
				auto temp_iter = elements[f]->find(x);
				if (temp_iter != elements[f]->end()) {
					elements[f]->erase(temp_iter);
					answers[x] = cost;
				} else {
					elements[f]->insert(x);
				}
			}
		}
	}

	for (int i = 0; i < n; i++) {
		delete elements[i];
	}

	return answers;
}

vector<int> solve_lca() {
	auto sorted_edges = edges;
	sort(sorted_edges);

	vector<vector<pii>> mst(n);

	DSU dsu(n);

	for (auto& [cost, edge] : sorted_edges) {
		auto f = dsu.merge(edge.first, edge.second);
		if (f != -1) {
			mst[edge.first].push_back({ edge.second, cost });
			mst[edge.second].push_back({ edge.first, cost });
		}
	}

	LCA lca(mst, 0);

	vector<int> answers(q, -1);

	for (int idx = 0; idx < q; idx++) {
		auto& edge = queries[idx];
		answers[idx] = lca.get_answer(edge.first, edge.second);
	}

	return answers;
}

void solve() {
	read();

	vector<int> answer;
	switch (algo_type)
	{
	case 1:
		answer = solve_merge_trees();
		break;
	case 2:
		answer = solve_lca();
		break;
	case 3:
		answer = solve_dijkstra();
		break;
	default:
		throw "Error, Incorrect algo type";
		break;
	}
	
	print_vector(answer);
}


int main()
{
	solve();

	return 0;
}
