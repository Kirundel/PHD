#include "BaseHeader.h"
#include "OrthogonalSearchChecker.h"
#include "OrthogonalSearchSegmentTree3D.h"
#include <random>
#include <cassert>
#include <iostream>

OrthogonalSearchChecker::OrthogonalSearchChecker() :
	cinput("input.txt"),
	coutput("output.txt")
{

}

void OrthogonalSearchChecker::print_point(PointWithIdx& point) {
	coutput << "Point: idx=" << point.idx << " x=" << point.x << " y=" << point.y << " z=" << point.z << "\n";
}

void OrthogonalSearchChecker::print_distances() {
	for (int i = 0; i < arr.size(); i++) {
		for (int j = i + 1; j < arr.size(); j++) {
			std::cout << arr[i].idx << " "  << arr[j].idx << " " << l1_distance(arr[i], arr[j]) << std::endl;
		}
	}
}

void OrthogonalSearchChecker::read_arr() {
	arr.clear();
	cinput >> n;

	arr.resize(n);

	for (int i = 0; i < n; i++) {
		auto& point = arr[i];
		cinput >> point.x >> point.y >> point.z;
		point.idx = i;
	}
}

const long long generator_seed = 543232244234;

void OrthogonalSearchChecker::generate_arr_full_random(int n) {
	arr.clear();
	this->n = n;
	
	arr.resize(n);

	mt19937 gen(generator_seed);

	for (int i = 0; i < n; i++) {
		auto& point = arr[i];
		point.idx = i;
		point.x = gen() % n;
		point.y = gen() % n;
		point.z = gen() % n;
	}
}

int OrthogonalSearchChecker::get_cubic_root(int n) {
	int l = 0;
	int r = n;

	while (r - l > 1) {
		int m = (l + r) / 2;
		if (m * m * m >= n) {
			r = m;
		}
		else {
			l = m;
		}
	}
	return r;
}

void OrthogonalSearchChecker::generate_arr_3D_cell(int n) {
	arr.clear();
	this->n = n;

	arr.reserve(n);

	mt19937 gen(generator_seed);

	long long number_of_coords = get_cubic_root(n);

	vector<int> ycoords(number_of_coords);
	vector<int> zcoords(number_of_coords);
	vector<int> xcoords;

	for (int i = 0; i < number_of_coords; i++) {
		ycoords[i] = gen() % n;
		zcoords[i] = gen() % n;
	}

	while (arr.size() < n)
	{
		int x = gen() % n;
		xcoords.push_back(x);
		for (auto y : ycoords) {
			for (auto z : zcoords) {
				auto point = PointWithIdx{ .x = x, .y = y, .z = z, .idx = (int)arr.size() };
				arr.push_back(
					point
				);

				if (arr.size() >= n)
					break;
			}
			if (arr.size() >= n)
				break;
		}
	}

	coutput << "X: ";
	sort(xcoords.begin(), xcoords.end());
	for (int i :xcoords) {
		coutput << i << " ";
	}
	coutput << "\nY: ";
	sort(ycoords.begin(), ycoords.end());
	for (int i : ycoords) {
		coutput << i << " ";
	}
	coutput << "\nZ: ";
	sort(zcoords.begin(), zcoords.end());
	for (int i : zcoords) {
		coutput << i << " ";
	}
	coutput << endl;
	

	std::shuffle(arr.begin(), arr.end(), gen);

	for (int i = 0; i < n; i++) {
		arr[i].idx = i;
	}
}

void OrthogonalSearchChecker::generate_arr_3D_cell_rotated(int n) {
	generate_arr_3D_cell(n);

	for (auto& point : arr) {
		int x = point.x + point.y + point.z;
		int y = point.x + point.y - point.z;
		int z = point.x - point.y + point.z;
	}
}

long long OrthogonalSearchChecker::get_working_time_mst_3d(MST3DAlgorithm& algorithm) {
	auto begin = std::chrono::steady_clock::now();

	Answer answer = algorithm.MST3D(this->n, this->arr);

	auto end = std::chrono::steady_clock::now();
	auto delta = end - begin;
	auto result_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(delta).count();

	return result_milliseconds;
}

pair<int, pair<int, int>> cmp_edge(Edge edge) {
	return { edge.distance, {min(edge.first, edge.second), max(edge.first, edge.second)} };
}

bool OrthogonalSearchChecker::compare_mst_3d(MST3DAlgorithm& first, MST3DAlgorithm& second) {
	auto begin = std::chrono::steady_clock::now();

	Answer first_answer = first.MST3D(this->n, this->arr);
	Answer second_answer = second.MST3D(this->n, this->arr);

	coutput << "F: " << first_answer.sum << "; S: " << second_answer.sum << "\n";

	sort(first_answer.edges.begin(), first_answer.edges.end(), [](auto& f, auto& s) { return cmp_edge(f) < cmp_edge(s); });
	/*coutput << "FIRST OUTPUT:" << " " << first_answer.edges.size() << "\n";
	for (auto& edge : first_answer.edges) {
		coutput << edge.distance << " " << edge.first << " " << edge.second << "\n";
	}*/

	sort(second_answer.edges.begin(), second_answer.edges.end(), [](auto& f, auto& s) { return cmp_edge(f) < cmp_edge(s); });
	//coutput << "SECOND OUTPUT:" << " " << second_answer.edges.size() << "\n";
	//for (auto& edge : second_answer.edges) {
	//	coutput << edge.distance << " " << edge.first << " " << edge.second << "\n";
	//}
	//coutput << "\n\n";

	int balance = 0;

	for (int idx = 0; idx < first_answer.edges.size(); idx++) {
		balance += (first_answer.edges[idx].distance - second_answer.edges[idx].distance);
		coutput << balance << "\t | "
			<< first_answer.edges[idx].distance << " " << second_answer.edges[idx].distance << " | \t"
			<< first_answer.edges[idx].first << " " << first_answer.edges[idx].second << " | "
			<< second_answer.edges[idx].first << " " << second_answer.edges[idx].second;
		coutput << "\n";
		if (first_answer.edges[idx].distance != second_answer.edges[idx].distance) {
			auto f_e = first_answer.edges[idx].first;
			auto s_e = first_answer.edges[idx].second;
			coutput << "(" << arr[f_e].x << " " << arr[f_e].y << " " << arr[f_e].z << ") - (" 
				<< arr[s_e].x << " " << arr[s_e].y << " " << arr[s_e].z << ") | ";
			f_e = second_answer.edges[idx].first;
			s_e = second_answer.edges[idx].second;
			coutput << "(" << arr[f_e].x << " " << arr[f_e].y << " " << arr[f_e].z << ") - ("
				<< arr[s_e].x << " " << arr[s_e].y << " " << arr[s_e].z << ")";
			coutput << "\n";
		}
	}
	

	return first_answer.sum == second_answer.sum;
}

bool OrthogonalSearchChecker::check_mst_pipeline() {

	auto copy_arr = arr;
	sort(copy_arr.begin(), copy_arr.end(),
		[](const auto& first, const auto& second)
		{ return -(first.x + first.y + first.z) < -(second.x + second.y + second.z); }
	);

	SegmentTreeCoreProperties* coreProperties = NULL;
	SegmentTreeCoreProperties* coreProperties_base = NULL;
	SegmentTree3D_V2* tree = NULL;
	OrthogonalSearch3DDummy* buddy = NULL;
	{
		auto pointsCopy = copy_arr;
		coreProperties = new SegmentTreeCoreProperties();
		tree = new SegmentTree3D_V2(coreProperties, pointsCopy);
		coreProperties_base = new SegmentTreeCoreProperties();
		buddy = new OrthogonalSearch3DDummy(coreProperties_base);
	}

	std::cout << "TREE INITIALIZED" << std::endl;

	long long number_of_points = 0;

	coreProperties->right = PointWithIdx{ .x = INF, .y = INF, .z = INF, .idx = -1 };
	coreProperties_base->right = PointWithIdx{ .x = INF, .y = INF, .z = INF, .idx = -1 };

	for (int i = 0; i < n; i++) {
		auto& point = copy_arr[i];

		coreProperties->queryPoint = point;
		coreProperties->left = point;
		coreProperties->answer.clear();
		coreProperties_base->queryPoint = point;
		coreProperties_base->left = point;
		coreProperties_base->answer.clear();

		tree->get_all_points(1);
		number_of_points += coreProperties->answer.size();
		buddy->get_all_points(1);

		assert(coreProperties->answer.size() == coreProperties_base->answer.size());
		
		//std::cout << "ADD: idx = " << point.idx << " x = " << point.x << " y = " << point.y << " z = " << point.z << "\n";
		tree->add_point(1);
		buddy->add_point(1);

		sort(coreProperties->answer.begin(), coreProperties->answer.end());
		sort(coreProperties_base->answer.begin(), coreProperties_base->answer.end());

		for (int idx = 0; idx < coreProperties->answer.size(); idx++) {
			assert(coreProperties->answer[idx] == coreProperties_base->answer[idx]);

			coreProperties->queryPoint = arr[coreProperties->answer[idx]];
			tree->remove_point(1);
			coreProperties_base->queryPoint = arr[coreProperties->answer[idx]];
			buddy->remove_point(1);

			//std::cout << "REMOVE: idx = " << coreProperties->queryPoint.idx << " x = " << coreProperties->queryPoint.x << " y = " << coreProperties->queryPoint.y << " z = " << coreProperties->queryPoint.z << "\n";
		}
	}

	delete tree;
	delete coreProperties;
	delete buddy;
	delete coreProperties_base;

	auto end = std::chrono::steady_clock::now();

	coutput << "NUMBER_OF_QUERIES_POINTS: " << number_of_points << "\n";
	return true;
}