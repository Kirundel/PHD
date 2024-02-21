#pragma once

#include <vector>

class DSU {
	std::vector<int> parents;
	std::vector<int> sizes;

public:
	DSU(int n);

	int get_parent(int v);

	bool merge(int x, int y);
};