#include "Algorithms/OrthogonalSearchChecker.h"
#include "Algorithms/OrthogonalSearchSegmentTree1D.h"
#include "Algorithms/OrthogonalSearchSegmentTree2D.h"
#include "Algorithms/OrthogonalSearchSegmentTree3D.h"
#include "Algorithms/KruskalAlgorithm.h"
#include "Algorithms/Standard3DMSTPipeline.h"
#include "Algorithms/Prim3DAlgorithm.h"
#include "Algorithms/Orthogonal3DMSTPipeline.h"
#include <iostream>


int main() {
    //checker.read_arr();
    //checker.check_segment_tree<SegmentTree1D>();
    //checker.check_segment_tree<SegmentTree2D>();
    //checker.check_segment_tree<SegmentTree3D>();
    //checker.generate_arr_full_random(320000);
    //checker.generate_arr_3D_cell(1000000);
    //auto result_1 = checker.get_working_time_check_pipeline<SegmentTree1D>();
    //std::cout << "1D: " << result_1 << std::endl;
    //auto result_2 = checker.get_working_time_check_pipeline<SegmentTree2D>();
    //std::cout << "2D: " << result_2 << std::endl;
    //auto result_3 = checker.get_working_time_check_pipeline<SegmentTree3D>();
    //std::cout << "3D: " << result_3 << std::endl;

    //auto result_1 = checker.get_working_time_mst_pipeline<SegmentTree1D>();
    //std::cout << "1D: " << result_1 << std::endl;
    //auto result_2 = checker.get_working_time_mst_pipeline<SegmentTree2D>();
    //std::cout << "2D: " << result_2 << std::endl;
    //auto result_3 = checker.get_working_time_mst_pipeline<SegmentTree3D>();
    //std::cout << "3D: " << result_3 << std::endl;

    //checker.generate_arr_3D_cell(10000);
    //checker.generate_arr_3D_cell_rotated(16);

    
    //checker.print_distances();

    for (int i = 9; i <= 18; i++) {
        OrthogonalSearchChecker checker;

        int num = 1 << i;

        checker.generate_arr_full_random(num);

        std::cout << num << std::endl;

        KruskalAlgorithm kruskal;
        //Standard3DMSTPipeline standard(kruskal);
        //std::cout << "Kruskal: " << checker.get_working_time_mst_3d(standard) << std::endl;

        //Prim3DAlgortihm prim;
        //std::cout << "Prim: " << checker.get_working_time_mst_3d(prim) << std::endl;

        Orthogonal3DMSTPipeline orthogonal(kruskal);
        std::cout << "Orthogonal: " << checker.get_working_time_mst_3d(orthogonal) << std::endl;
    }

    //std::cout << "COMPARE: " << checker.compare_mst_3d(standard, prim) << " " << checker.compare_mst_3d(prim, orthogonal) << std::endl;
    //std::cout << "COMPARE: " << checker.compare_mst_3d(prim, orthogonal) << std::endl;

    //std::cout << "FIND ERROR: " << checker.check_mst_pipeline() << std::endl;

    return 0;
}