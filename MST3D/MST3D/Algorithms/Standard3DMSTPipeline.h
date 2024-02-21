#pragma once

#include "BaseClasses.h"


class Standard3DMSTPipeline : public MST3DAlgorithm {
private:
	MSTAlgorithm& mstAlgo;

public:
	Standard3DMSTPipeline(MSTAlgorithm& algorithm);

	Answer MST3D(int n, std::vector<PointWithIdx>& points) override;
};