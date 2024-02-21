#pragma once

#include <vector>
#include "OrthogonalSearchSegmentTree1D.h"
//#include "OrthogonalSearchSegmentTree2D.h"
#include <fstream>
#include <chrono>
#include <iostream>
#include "BaseHeader.h"

class OrthogonalSearchChecker {
private:
    std::ifstream cinput;
    std::ofstream coutput;

	int get_cubic_root(int n);

public:
    OrthogonalSearchChecker();

    void print_point(PointWithIdx& point);

    int n;
    std::vector<PointWithIdx> arr;

    void read_arr();

	void generate_arr_full_random(int n);

	void generate_arr_3D_cell(int n);
	void generate_arr_3D_cell_rotated(int n);

    template<typename T>
    void check_segment_tree();

	template<typename T>
	long long get_working_time_check_pipeline();

	template<typename T>
	long long get_working_time_mst_pipeline();

	bool check_mst_pipeline();

	void print_distances();

	long long get_working_time_mst_3d(MST3DAlgorithm& algorithm);	

	bool compare_mst_3d(MST3DAlgorithm& first, MST3DAlgorithm& second);
};

template<typename T>
void OrthogonalSearchChecker::check_segment_tree() {
	SegmentTreeCoreProperties* coreProperties = NULL;
	T* tree = NULL;
	{
		auto pointsCopy = arr;
		coreProperties = new SegmentTreeCoreProperties();
		tree = new T(coreProperties, pointsCopy);
	}

	coutput << "ADD:\n\n";
	for (int i = 0; i < n; i++) {
		auto& point = arr[i];
		print_point(point);

		coreProperties->queryPoint = point;
		coreProperties->left = point;
		coreProperties->right = PointWithIdx{ .x = INF, .y = INF, .z = INF, .idx = -1 };
		coreProperties->answer.clear();

		tree->get_all_points(1);

		coutput << "Result BEFORE: ";
		sort(coreProperties->answer.begin(), coreProperties->answer.end());
		for (int& idx : coreProperties->answer) {
			coutput << idx << " ";
		}
		coutput << "\n";

		tree->add_point(1);

		coreProperties->answer.clear();
		tree->get_all_points(1);

		coutput << "Result AFTER: ";
		sort(coreProperties->answer.begin(), coreProperties->answer.end());
		for (int& idx : coreProperties->answer) {
			coutput << idx << " ";
		}
		coutput << "\n\n";
	}

	coutput << "QUERY:\n\n";

	for (int i = 0; i < n; i++) {
		auto& point = arr[i];
		print_point(point);

		coreProperties->left = point;
		coreProperties->right = PointWithIdx{ .x = INF, .y = INF, .z = INF, .idx = -1 };
		coreProperties->answer.clear();

		tree->get_all_points(1);

		coutput << "Result: ";
		sort(coreProperties->answer.begin(), coreProperties->answer.end());
		for (int& idx : coreProperties->answer) {
			coutput << idx << " ";
		}
		coutput << "\n\n";
	}

	coutput << "REMOVE:\n\n";
	for (int i = 0; i < n; i++) {
		auto& point = arr[i];
		print_point(point);

		coreProperties->left = point;
		coreProperties->right = PointWithIdx{ .x = INF, .y = INF, .z = INF, .idx = -1 };
		coreProperties->queryPoint = point;
		coreProperties->answer.clear();
		tree->get_all_points(1);

		coutput << "Result BEFORE: ";
		sort(coreProperties->answer.begin(), coreProperties->answer.end());
		for (int& idx : coreProperties->answer) {
			coutput << idx << " ";
		}
		coutput << "\n";

		tree->remove_point(1);

		coreProperties->answer.clear();
		tree->get_all_points(1);

		coutput << "Result AFTER: ";
		sort(coreProperties->answer.begin(), coreProperties->answer.end());
		for (int& idx : coreProperties->answer) {
			coutput << idx << " ";
		}

		coutput << "\n\n";
	}

	delete tree;
	delete coreProperties;
}

template<typename T>
long long OrthogonalSearchChecker::get_working_time_check_pipeline() {
	auto begin = std::chrono::steady_clock::now();

	SegmentTreeCoreProperties* coreProperties = NULL;
	T* tree = NULL;
	{
		auto pointsCopy = arr;
		coreProperties = new SegmentTreeCoreProperties();
		tree = new T(coreProperties, pointsCopy);
	}
	
	long long number_of_points = 0;

	for (int i = 0; i < n; i++) {
		auto& point = arr[i];

		coreProperties->queryPoint = point;
		coreProperties->left = point;
		coreProperties->right = PointWithIdx{ .x = INF, .y = INF, .z = INF, .idx = -1 };
		coreProperties->answer.clear();

		tree->get_all_points(1);
		number_of_points += coreProperties->answer.size();

		tree->add_point(1);

		coreProperties->answer.clear();
		tree->get_all_points(1);
		number_of_points += coreProperties->answer.size();
	}

	for (int i = 0; i < n; i++) {
		auto& point = arr[i];

		coreProperties->left = point;
		coreProperties->right = PointWithIdx{ .x = INF, .y = INF, .z = INF, .idx = -1 };
		coreProperties->answer.clear();

		tree->get_all_points(1);
		number_of_points += coreProperties->answer.size();
	}

	for (int i = 0; i < n; i++) {
		auto& point = arr[i];

		coreProperties->left = point;
		coreProperties->right = PointWithIdx{ .x = INF, .y = INF, .z = INF, .idx = -1 };
		coreProperties->queryPoint = point;
		coreProperties->answer.clear();
		tree->get_all_points(1);
		number_of_points += coreProperties->answer.size();

		tree->remove_point(1);

		coreProperties->answer.clear();
		tree->get_all_points(1);
		number_of_points += coreProperties->answer.size();
	}

	delete tree;
	delete coreProperties;

	auto end = std::chrono::steady_clock::now();
	auto delta = end - begin;
	auto result_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(delta).count();

	coutput << "NUMBER_OF_QUERIES_POINTS: " << number_of_points << "\n";
	return result_milliseconds;
}

template<typename T>
long long OrthogonalSearchChecker::get_working_time_mst_pipeline() {
	auto copy_arr = arr;
	sort(copy_arr.begin(), copy_arr.end(),
		[](const auto& first, const auto& second) 
		{ return -(first.x + first.y + first.z) < -(second.x + second.y + second.z); }
	);

	auto begin = std::chrono::steady_clock::now();

	SegmentTreeCoreProperties* coreProperties = NULL;
	T* tree = NULL;
	{
		auto pointsCopy = copy_arr;
		coreProperties = new SegmentTreeCoreProperties();
		tree = new T(coreProperties, pointsCopy);
	}

	std::cout << "TREE INITIALIZED" << std::endl;

	long long number_of_points = 0;

	coreProperties->right = PointWithIdx{ .x = INF, .y = INF, .z = INF, .idx = -1 };

	for (int i = 0; i < n; i++) {
		auto& point = copy_arr[i];

		coreProperties->queryPoint = point;
		coreProperties->left = point;
		coreProperties->answer.clear();

		tree->get_all_points(1);
		number_of_points += coreProperties->answer.size();

		tree->add_point(1);

		for (auto pointIdx : coreProperties->answer) {
			coreProperties->queryPoint = arr[pointIdx];
			tree->remove_point(1);
		}

		/*if (i % 10000 == 0) {
			std::cout << "RESULTED " << i << std::endl;
		}*/
	}

	delete tree;
	delete coreProperties;

	auto end = std::chrono::steady_clock::now();
	auto delta = end - begin;
	auto result_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(delta).count();

	coutput << "NUMBER_OF_QUERIES_POINTS: " << number_of_points << "\n";
	return result_milliseconds;
}