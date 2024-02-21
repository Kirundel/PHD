#include "BaseHeader.h"
#include "Prim3DAlgorithm.h"

Answer Prim3DAlgortihm::MST3D(int n, std::vector<PointWithIdx>& points) {
	vector<bool> used(n);
	//vector<int> min_e(n, 1000 * 1000 * 1000), sel_e(n, -1);
	vector<int> min_e(n, INF), sel_e(n, -1);
	min_e[0] = 0;

	Answer answer;

	for (int i = 0; i < n; ++i) {
		int v = -1;
		for (int j = 0; j < n; ++j) {
			if (!used[j] && (v == -1 || min_e[j] < min_e[v])) {
				v = j;
			}
		}

		used[v] = true;
		if (sel_e[v] != -1) {
			answer.edges.push_back(Edge(points[v], points[sel_e[v]]));
		}

		for (int to = 0; to < n; ++to)
		{
			const int edge_val = l1_distance(points[v], points[to]);
			if (edge_val < min_e[to]) {
				min_e[to] = edge_val;
				sel_e[to] = v;
			}
		}
	}

	for (auto& edge : answer.edges) {
		answer.sum += edge.distance;
	}

	return answer;
}