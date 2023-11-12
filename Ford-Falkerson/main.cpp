#include "borel_tanner.h"
#include "ford_falkerson.h"

using namespace graphs;

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>

#include "mpi.h"

std::mt19937 gen{ std::random_device{}() };

#define MASTER_RANK 0

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

template<typename T>
std::istream& operator>>(std::istream& is, const std::vector<std::vector<T>>& g) {
    for (const auto& row : g) {
        for (const auto& elem : row) {
            is >> elem;
        }
    }
    return is;
}

int main() {
    MPI::Init();
    int rank = MPI::COMM_WORLD.Get_rank();
    int world_size = MPI::COMM_WORLD.Get_size();

    size_t nVertexes;
    if (rank == MASTER_RANK) {
        std::cout << "Enter num of the vertexes: ";
        std::cin >> nVertexes;
    }
    MPI::COMM_WORLD.Bcast(&nVertexes, 1, MPI::UNSIGNED_LONG, MASTER_RANK);

    adjacency_matrix<> graph;

    // Random input(safe)
    if (rank == MASTER_RANK) {
        graph = generate(nVertexes, gen);
        std::cout << "Adjaceny matrix:\n" << graph;
    }

    // User input
    /*if (rank == MASTER_RANK) {
        graph.resize(nVertexes, std::vector<int>(nVertexes));
        std::cout << "Eneter the adjacency matrix:\n";
        std::cin >> graph;
    }*/

    // Transmit adjacency matrix to other nodes
    for (int i = 0; i < nVertexes; ++i)
    {
        std::vector<int> row(nVertexes);
        if (rank == MASTER_RANK) {
            row = graph[i];
        }
        MPI::COMM_WORLD.Bcast(row.data(), nVertexes, MPI::INT, MASTER_RANK);
        if (rank != MASTER_RANK) {
            graph.push_back(row);
        }
    }

    if (rank == MASTER_RANK) {
        std::cout << "Stage 1" << std::endl;
    }

    // Add suoersourse and supersink
    flow_graph_t flowgraph = add_supersource_supersink(graph);

    // Mpi ford-falkerson
    flow_result_t result = mpi_max_flow_ford_fulkerson(flowgraph, rank, world_size);

    // Program output
    if (rank == MASTER_RANK) {
        std::cout << "Max flow: " << result.max_flow << std::endl;
        std::cout << "Flow matrix:\n" << result.flow;
    }

    MPI::Finalize();
    return 0;
}