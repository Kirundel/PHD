#pragma once

#include "BaseClasses.h"

class KruskalAlgorithm : public MSTAlgorithm {
public:
	Answer MST(int n, std::vector<Edge>& edges) override;
};
