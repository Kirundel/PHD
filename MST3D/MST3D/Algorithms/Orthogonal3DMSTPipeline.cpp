#include "Orthogonal3DMSTPipeline.h"
#include "OrthogonalSearchSegmentTree3D.h"
#include "AdvancedSegmentTree3D.h"
#include "BaseHeader.h"
#include <functional>
#include <iostream>
#include <chrono>

Orthogonal3DMSTPipeline::Orthogonal3DMSTPipeline(MSTAlgorithm& algorithm) : mstAlgo(algorithm)
{}


//vector<function<PointWithIdx(PointWithIdx&)>> cones = {
//	[](PointWithIdx& p) {return PointWithIdx{p.x, p.y, p.z - p.x - p.y, p.idx}; } ,
//	[](PointWithIdx& p) {return PointWithIdx{p.x, p.z, p.y - p.x - p.z, p.idx}; } ,
//	[](PointWithIdx& p) {return PointWithIdx{p.y, p.z, p.x - p.y - p.z, p.idx}; } ,
//	[](PointWithIdx& p) {return PointWithIdx{p.y + p.z - p.x, p.x + p.z - p.y, p.x + p.y - p.z, p.idx}; } ,
//};

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

Answer Orthogonal3DMSTPipeline::MST3D_main(int n, vector<PointWithIdx>& points)
{
	std::cout << "\n";
	vector<Edge> edges;

	int cone_idx = 0;

	//for (int m2 = -1; m2 < 2; m2 += 2) {
	//	for (int m3 = -1; m3 < 2; m3 += 2) {
	{
		{
			//vector<PointWithIdx> cur_points_prev = generate_translation(points, m2, m3);
			vector<PointWithIdx> cur_points_prev = points;

			for (int sc = 0; sc < cones_sort.size(); sc++)
			{
				auto& sort_criterium = cones_sort[sc];
				auto cur_points = cur_points_prev;
				sort(cur_points.begin(), cur_points.end(),
					[&](auto& first, auto& second)
					{ 
						const auto f_r = sort_criterium(first);
						const auto f_s = sort_criterium(second);
						if (f_r != f_s) return f_r < f_s;
						if (first.x != second.x) return first.x < second.x;
						if (first.y != second.y) return first.y < second.y;
						return first.z < second.z;
					}
				);


				for (int cur_cone = 0; cur_cone < new_cones[sc].size(); cur_cone++) {
					auto& cone = new_cones[sc][cur_cone];

					int number_of_points = 0;
					auto points_via_cone = cur_points;

					/*sort(points_via_cone.begin(), points_via_cone.end(),
						[&](auto& f, auto& s) { return cone_sort(f) < cone_sort(s); }
					);*/

					points_via_cone = translate_via_cone(points_via_cone, cone);

					auto beg = chrono::steady_clock::now();
					SegmentTreeCoreProperties* coreProperties = NULL;
					//OrthogonalSearch3DDummy* tree = NULL;
					SegmentTree3D_V3* tree = NULL;
					{
						auto pointsCopy = points_via_cone;
						coreProperties = new SegmentTreeCoreProperties();
						tree = new SegmentTree3D_V3(coreProperties, pointsCopy);
						//tree = new OrthogonalSearch3DDummy(coreProperties);
					}
					long long tree_init = std::chrono::duration_cast<std::chrono::nanoseconds>(chrono::steady_clock::now() - beg).count();


					coreProperties->right = PointWithIdx{ .x = INF, .y = INF, .z = INF, .idx = -1 };

					long long add_points = 0;
					long long get_points = 0;
					long long remove_points = 0;

					for (int i = 0; i < n; i++) {
						auto& point = points_via_cone[i];
						//std::cout << cur_cone << ": idx = " << point.idx << " x = " << point.x << " y = " << point.y << " z = " << point.z << "\n";

						coreProperties->queryPoint = point;
						coreProperties->answer.clear();

						beg = chrono::steady_clock::now();
						coreProperties->left = point;
						tree->get_all_points(1);
						get_points += std::chrono::duration_cast<std::chrono::nanoseconds>(chrono::steady_clock::now() - beg).count();

						number_of_points += coreProperties->answer.size();

						beg = chrono::steady_clock::now();
						for (auto pointIdx : coreProperties->answer) {
							edges.push_back(Edge(points[point.idx], points[pointIdx]));
							//coreProperties->queryPoint = ApplyTransitions(points[pointIdx], m2, m3, cone);
							coreProperties->queryPoint = ApplyTransitions(points[pointIdx], 1, 1, 1, cone);

							tree->remove_point(1);
						}
						remove_points += std::chrono::duration_cast<std::chrono::nanoseconds>(chrono::steady_clock::now() - beg).count();

						beg = chrono::steady_clock::now();
						coreProperties->queryPoint = point;
						tree->add_point(1);
						add_points += std::chrono::duration_cast<std::chrono::nanoseconds>(chrono::steady_clock::now() - beg).count();
					}

					std::cout << "CONE " << cone_idx << " ENDED NUM_POINTS: " << number_of_points
						<< " ADD POINTS: " << (add_points + 500000) / 1000000 
						<< " GET POINTS: " << (get_points + 500000) / 1000000
						<< " DEL POINTS: " << (remove_points + 500000) / 1000000
						<< " TREE INIT : " << (tree_init + 500000) / 1000000
						<< " INSERT    : " << coreProperties->add_nums
						<< " DELETE    : " << coreProperties->erase_nums
						<< " GET NUMS  : " << coreProperties->get_nums
						<< std::endl;					
					cone_idx++;

					delete coreProperties;
					delete tree;

				}
			}
		}
	}

	return mstAlgo.MST(n, edges);
}

Answer Orthogonal3DMSTPipeline::MST3D(int n, vector<PointWithIdx>& points)
{
	vector<Edge> edges;

	int cone_idx = 0;

	auto pipeline_begin = chrono::steady_clock::now();

	for (int m1 = -1; m1 < 2; m1+= 2) {
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
						//OrthogonalSearch3DDummy* tree = NULL;
						auto beg = chrono::steady_clock::now();
						AdvancedSegmentTree3D* tree = NULL;
						{
							auto pointsCopy = points_via_cone;
							coreProperties = new SegmentTreeCoreProperties();
							tree = new AdvancedSegmentTree3D(coreProperties, pointsCopy);
							//tree = new OrthogonalSearch3DDummy(coreProperties);
						}
						long long tree_init = std::chrono::duration_cast<std::chrono::nanoseconds>(chrono::steady_clock::now() - beg).count();

						long long add_points = 0;
						long long get_points = 0;
						long long remove_points = 0;

						coreProperties->right = PointWithIdx{ .x = INF, .y = INF, .z = INF, .idx = -1 };

						for (int i = 0; i < n; i++) {
							auto& point = points_via_cone[i];
							//std::cout << cur_cone << ": idx = " << point.idx << " x = " << point.x << " y = " << point.y << " z = " << point.z << "\n";

							coreProperties->left = point;
							coreProperties->answer.clear();

							//beg = chrono::steady_clock::now();
							coreProperties->queryPoint = point;
							tree->get_all_points(1);
							//get_points += std::chrono::duration_cast<std::chrono::nanoseconds>(chrono::steady_clock::now() - beg).count();

							number_of_points += coreProperties->answer.size();

							//beg = chrono::steady_clock::now();
							for (auto pointIdx : coreProperties->answer) {
								edges.push_back(Edge(points[point.idx], points[pointIdx]));
								coreProperties->queryPoint = ApplyTransitions(points[pointIdx], m1, m2, m3, cone);

								tree->remove_point(1);
							}
							//remove_points += std::chrono::duration_cast<std::chrono::nanoseconds>(chrono::steady_clock::now() - beg).count();

							//beg = chrono::steady_clock::now();
							coreProperties->queryPoint = point;
							tree->add_point(1);
							//add_points += std::chrono::duration_cast<std::chrono::nanoseconds>(chrono::steady_clock::now() - beg).count();
						}

						//std::cout << "CONE " << cone_idx << " ENDED NUM_POINTS: " << number_of_points
							//<< " ADD POINTS: " << (add_points + 500000) / 1000000 
							//<< " GET POINTS: " << (get_points + 500000) / 1000000
							//<< " DEL POINTS: " << (remove_points + 500000) / 1000000
							//<< " TREE INIT : " << (tree_init + 500000) / 1000000
							//<< " INSERT    : " << coreProperties->sum_tt 
							//<< std::endl;
						cone_idx++;

						delete coreProperties;
						delete tree;

					}
				}
			}
		}
	}

	//cout << "\nPIPELINE : " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono::steady_clock::now() - pipeline_begin).count() << std::endl;

	return mstAlgo.MST(n, edges);
}