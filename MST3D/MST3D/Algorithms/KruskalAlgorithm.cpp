#include "BaseHeader.h"
#include "KruskalAlgorithm.h"
#include "DSU.h"


Answer KruskalAlgorithm::MST(int n, std::vector<Edge>& edges) {
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