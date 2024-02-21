#pragma once

#include <vector>
#include <utility>

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
		return std::min(first, second) == std::min(other.first, other.second)
			&& std::max(first, second) == std::max(other.first, other.second);
	}
};

struct Answer {
	long long sum = 0;
	std::vector<Edge> edges;
};

class MSTAlgorithm {
public:
	virtual Answer MST(int n, std::vector<Edge>& edges) = 0;
};

class MST3DAlgorithm {
public:
	virtual Answer MST3D(int n, std::vector<PointWithIdx>& points) = 0;
};