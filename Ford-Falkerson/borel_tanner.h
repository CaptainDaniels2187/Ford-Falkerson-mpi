//
// Created by mhq on 13/03/23.
//

#pragma once
#include <random>
#include <tuple>

#define M_PI 3.1415926535
#define M_E 2.7182818284

// start snippet polya_1
template <class IntType = int, class DoubleType = double>
std::discrete_distribution<IntType>
borel_tanner(IntType R, DoubleType alpha, IntType size) {
	std::vector<DoubleType> probabilities;
	probabilities.reserve(size + 1);
	DoubleType p;
	for (int i = 0; i < R; ++i)
	{
		p = 0;
		probabilities.push_back(p);
	}
	p = pow(M_E, (-1 * alpha * R));
	probabilities.push_back(p);
	for (int i = R + 1; i <= size; ++i) {
		p = (R / sqrt(2.0 * M_PI * (i - R))) * (1.0 / i) * pow(((alpha * i) / (i - R)), (i - R)) * pow(M_E, ((1 - alpha) * i - R));
		probabilities.push_back(p);
	}
	return { probabilities.begin(), probabilities.end() };
}

void test_borel();
void borel_histograms();