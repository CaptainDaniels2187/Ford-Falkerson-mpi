#pragma once

#include "borel_tanner.h"
#include <vector>
#include <set>
#include <cstddef>
#include <deque>

namespace graphs {
    // representation of a graph in a computer
    template<typename T = int>
    using adjacency_matrix = std::vector<std::vector<T>>;

    using Vertex = unsigned long;

    using path_t = std::vector<Vertex>;

    adjacency_matrix<> generate(size_t nVertices, std::mt19937& gen);

    std::vector<size_t> out_degrees(size_t nVertices, std::mt19937& gen);

    adjacency_matrix<> from_degrees(std::vector<size_t> vertex_degrees, std::mt19937& gen);

    constexpr auto INF = INT32_MAX;
    constexpr Vertex NO_VERTEX = -1u;

    // start snippet flow_graph_t
    struct flow_graph_t {
        // матрица пропускных способностей
        adjacency_matrix<> capacity;
        // вершина источник
        Vertex source;
        // вершина сток
        Vertex sink;
    };
    // end snippet flow_graph_t

    // adds to the graph supersource and supersink with ∞ edge capacities if needed,
    // returns new matrix and supersource & supersink vertexes

    // добавить к графу фиктивный исток и сток с бесконечной пропускной способностью
    // если стоков (истоков) несколько
    // возвращает новую матрицу и номера вершин истока и стока в новой матрице.
    flow_graph_t add_supersource_supersink(const adjacency_matrix<>& capacity);

    // start snippet flow_result_t
    struct flow_result_t {
        int max_flow; // величина потока
        adjacency_matrix<> flow; // матрица потока
    };
    // end snippet flow_result_t

    // find maximim flow of given capacities matrix. Ignores g.cost.
    // найти максимальный поток по данной матрице пропускных способностей (игнорируя стоимости).
    flow_result_t mpi_max_flow_ford_fulkerson(const flow_graph_t& g, int rank, int world_size);
}

