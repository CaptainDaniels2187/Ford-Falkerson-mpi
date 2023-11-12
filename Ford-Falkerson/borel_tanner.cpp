//
// Created by mhq on 22/03/23.
//
#include "borel_tanner.h"
#include <map>
#include <iostream>

void test_borel() {
    auto dis = borel_tanner<>(2, 0.5, 5);
    std::mt19937 gen{std::random_device{}()};
    std::map<int, int> map;

    for (int n = 0; n < 1e5; ++n)
        ++map[dis(gen)];

    for(const auto& [num, count] : map)
        std::cout << num << "\t" << count << "\n";
}

// start snippet borel_histograms
void borel_histograms() {
    std::mt19937 gen{std::random_device{}()};
    constexpr auto N = 10000; // итераций
    constexpr auto n = 49;
    for (auto b : {7}) {
        auto distribution = borel_tanner(b, 0.2, n);
        std::vector<size_t> data(n + 1);
        std::cout << "b = " << b << '\n';
        for (int i = 0; i < N; ++i) {
            data[distribution(gen)]++;
        }
        std::cout << 'x' << '\t' << "p(x)" << '\n';
        for (int i = 0; i <= n; ++i) {
            std::cout << i << '\t' << static_cast<double>(data[i]) / N << '\n';
            //std::cout << data[i] * 1.0 / N << std::endl;
        }
    }
}
// end snippet borel_histograms