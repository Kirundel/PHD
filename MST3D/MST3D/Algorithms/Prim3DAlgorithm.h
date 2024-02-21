#pragma once

#include "BaseClasses.h"

class Prim3DAlgortihm : public MST3DAlgorithm {
public:
	Answer MST3D(int n, std::vector<PointWithIdx>& points) override;
};