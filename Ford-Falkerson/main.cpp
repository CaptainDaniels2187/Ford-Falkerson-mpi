#include "borel_tanner.h"
#include "ford_falkerson.h"

using namespace graphs;

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>

std::mt19937 gen{ std::random_device{}() };

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<T>>& g) {
    for (const auto& row : g) {
        for (const auto& elem : row) {
            os << elem << ' ';
        }
        os << '\n';
    }
    return os;
}

int main() {

    return 0;
}