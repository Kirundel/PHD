#include "Standard3DMSTPipeline.h"
#include "BaseHeader.h"

Standard3DMSTPipeline::Standard3DMSTPipeline(MSTAlgorithm& algorithm) : mstAlgo(algorithm)
{}

Answer Standard3DMSTPipeline::MST3D(int n, std::vector<PointWithIdx>& points)
{
	vector<Edge> edges;

	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			edges.push_back(Edge(points[i], points[j]));
		}
	}

	return mstAlgo.MST(n, edges);
}