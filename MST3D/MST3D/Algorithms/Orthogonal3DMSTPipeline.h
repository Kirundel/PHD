#pragma once

#include "BaseClasses.h"


class Orthogonal3DMSTPipeline : public MST3DAlgorithm {
private:
	MSTAlgorithm& mstAlgo;

public:
	Orthogonal3DMSTPipeline(MSTAlgorithm& algorithm);

	Answer MST3D(int n, std::vector<PointWithIdx>& points) override;
	Answer MST3D_main(int n, std::vector<PointWithIdx>& points);
};