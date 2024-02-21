#include "DSU.h"

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
	if (sizes[x] < sizes[y]) std::swap(x, y);

	sizes[x] += sizes[y];
	parents[y] = x;
	return true;
}