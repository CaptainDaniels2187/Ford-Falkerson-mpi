#include "ford_falkerson.h"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <queue>
#include <utility>
#include <map>
#include <set>

#include "mpi.h"
using namespace graphs;

#define MASTER_RANK 0

// start snippet mpi-generate
//adjacency_matrix<> graphs::generate(size_t nVertices, int rank, int world_size, std::mt19937& gen) {
//    int subsize;
//    if (rank < world_size - 1) {
//        subsize = nVertices / world_size;
//    }
//    else {
//        if (nVertices > 2) {
//            if (nVertices % world_size == 0) {
//                subsize = nVertices / world_size;
//            }
//            else {
//                subsize = nVertices % world_size;
//            }
//        }
//        else {
//            subsize = 0;
//        }
//    }
//    std::vector<size_t> vertex_degrees;
//    vertex_degrees = out_degrees(nVertices, subsize, rank, world_size, gen);
//    MPI::COMM_WORLD.Bcast(vertex_degrees.data(), vertex_degrees.size(), MPI::UNSIGNED_LONG, MASTER_RANK);
//    adjacency_matrix<> result;
//    if (vertex_degrees.size() > 1) {
//        adjacency_matrix<> subresult;
//        subresult = from_degrees(vertex_degrees, subsize, rank, world_size, gen);
//        
//    }
//    else {
//        if (rank == MASTER_RANK) {
//            if (vertex_degrees.size() == 1) {
//                result.resize(1);
//                result[0][0] = 1;
//            }
//        }
//        return result;
//    }
//
//}
//
//std::vector<size_t> graphs::out_degrees(size_t nVertices, int subsize, int rank, int world_size, std::mt19937& gen) {
//    if (nVertices == 0) {
//        return std::vector<size_t>{};
//    }
//    if (nVertices == 1) {
//        return std::vector<size_t>{0};
//    }
//    // степень вершины может быть в диапазоне [1 .. n-1]
//    // распределение выдает числа в диапазоне [0 .. n-2]
//    auto dis = borel_tanner<int, double>(1, 0.7, static_cast<int>(nVertices));
//     
//    std::vector<size_t> subresult(subsize);
//    for (size_t& x : subresult) {
//        x = dis(gen);
//    }
//
//    std::vector<int> sizes;
//    if (rank == MASTER_RANK) {
//        sizes.resize(world_size);
//    }
//    MPI::COMM_WORLD.Gather(&subsize, 1, MPI::INT, sizes.data(), 1, MPI::INT, MASTER_RANK);
//    std::vector<int> displacements;
//    if (rank == MASTER_RANK) {
//        displacements.resize(world_size, 0);
//        for (int i = 1; i < world_size; i++) {
//            if (sizes[i] > 0) {
//                displacements[i] = (displacements[i - 1] + sizes[i - 1]);
//            }
//        }
//    }
//
//    std::vector<size_t> result;
//    if (rank == MASTER_RANK) {
//        result.resize(nVertices);
//    }
//    MPI::COMM_WORLD.Gatherv(subresult.data(), subresult.size(), MPI::UNSIGNED_INT, result.data(), sizes.data(), displacements.data(), MPI::UNSIGNED_INT, MASTER_RANK);
//
//    if (rank == MASTER_RANK) {
//        std::sort(result.rbegin(), result.rend());
//
//        // чтобы обеспечить отсутствие циклов
//        for (size_t i = 0; i < nVertices; ++i) {
//            result[i] = std::min(nVertices - i - 1, result[i]);
//        }
//    }
//    return result;
//}
//
//adjacency_matrix<> graphs::from_degrees(std::vector<size_t> vertex_degrees, int subsize, int rank, int world_size, std::mt19937& gen) {
//    auto nVertices = vertex_degrees.size();
//    adjacency_matrix<> result(subsize, std::vector<int>(nVertices));
//    auto dis = borel_tanner<size_t, double>(7, 0.2, 50 - 1);
//
//    for (int i = 0; i < subsize; ++i) {
//        int displacement = rank * (nVertices / world_size) + i;
//        for (int j = displacement + 1; j < vertex_degrees[displacement] + displacement + 1; ++j) {
//            result[i][j] = dis(gen);
//            assert(result[i][j] > 0);
//        }
//        if (vertex_degrees[displacement] == nVertices - displacement - 1) {
//            continue;
//        }
//        // degrees[i] of ones, other are zeros.
//        // nVertices - i - 1 total
//        std::shuffle(result[i].begin() + displacement + 1, result[i].end(), gen);
//    }
//    return result;
//}
//// end snippet mpi-generate

// start snippet generate
adjacency_matrix<> graphs::generate(size_t nVertices, std::mt19937& gen) {
    return from_degrees(out_degrees(nVertices, gen), gen);
}

std::vector<size_t> graphs::out_degrees(size_t nVertices, std::mt19937& gen) {
    if (nVertices == 0) {
        return std::vector<size_t>{};
    }
    if (nVertices == 1) {
        return std::vector<size_t>{0};
    }
    // степень вершины может быть в диапазоне [1 .. n-1]
    // распределение выдает числа в диапазоне [0 .. n-2]
    auto dis = borel_tanner<int, double>(1, 0.7, static_cast<int>(nVertices));
    std::vector<size_t> result(nVertices);
    for (size_t& x : result) {
        x = dis(gen);
    }
    std::sort(result.rbegin(), result.rend());

    // чтобы обеспечить отсутствие циклов
    for (size_t i = 0; i < nVertices; ++i) {
        result[i] = std::min(nVertices - i - 1, result[i]);
    }
    return result;
}

adjacency_matrix<> graphs::from_degrees(std::vector<size_t> vertex_degrees, std::mt19937& gen) {
    auto nVertices = vertex_degrees.size();
    adjacency_matrix<> result(nVertices, std::vector<int>(nVertices));
    auto dis = borel_tanner<size_t, double>(7, 0.2, 50 - 1);
    for (int i = 0; i < nVertices; ++i) {
        for (int j = i + 1; j < vertex_degrees[i] + i + 1; ++j) {
            result[i][j] = dis(gen);
            assert(result[i][j] > 0);
        }
        if (vertex_degrees[i] == nVertices - i - 1) {
            continue;
        }
        // degrees[i] of ones, other are zeros.
        // nVertices - i - 1 total
        std::shuffle(result[i].begin() + i + 1, result[i].end(), gen);
    }
    return result;
}
// end snippet generate

// start snippet add_supersource_supersink
flow_graph_t graphs::add_supersource_supersink(const adjacency_matrix<>& capacity) {
    // identify sources and sinks
    std::vector<bool> is_source(capacity.size(), true);
    std::vector<bool> is_sink(capacity.size(), true);
    for (Vertex i{}; i < capacity.size(); ++i) {
        for (Vertex j{}; j < capacity.size(); ++j) {
            if (capacity[i][j] != 0) {
                is_source[j] = false;
                is_sink[i] = false;
            }
        }
    }
    std::vector<Vertex> sources;
    std::vector<Vertex> sinks;
    for (Vertex i{}; i < capacity.size(); ++i) {
        if (is_sink[i]) {
            sinks.push_back(i);
            continue;
        }
        if (is_source[i]) {
            sources.push_back(i);
        }
    }
    // fix matrices if needed
    flow_graph_t result;
    result.capacity = capacity;

    assert(!sources.empty());
    if (sources.size() == 1) {
        result.source = sources.front();
    }
    else {
        result.source = capacity.size();
        for (auto& row : result.capacity) {
            row.emplace_back();
        }
        result.capacity.emplace_back(result.source + 1, 0u);
        for (auto v : sources) {
            result.capacity[result.source][v] = INF;
        }
    }

    assert(!sinks.empty());
    if (sinks.size() == 1) {
        result.sink = sinks.front();
    }
    else {
        result.sink = capacity.size() + (sources.size() > 1);
        for (auto& row : result.capacity) {
            row.emplace_back();
        }
        result.capacity.emplace_back(result.sink + 1, 0u);
        for (auto v : sinks) {
            result.capacity[v][result.sink] = INF;
        }
    }
    return result;
}
// end snippet add_supersource_supersink

// start snippet mpi-max_flow_ford_fulkerson
flow_result_t graphs::mpi_max_flow_ford_fulkerson(const flow_graph_t& g, int rank, int world_size) {
    const auto s = g.source;
    const auto t = g.sink;
    if (rank == MASTER_RANK) {
        std::cout << s << ' ' << t << std::endl;
    }

    auto capacity = g.capacity;
    const auto n = g.capacity.size();
    auto parent = std::vector<Vertex>(n, NO_VERTEX);

    // Returns true if there is a path from
    // source `s` to sink `t` in residual graph.
    // Also fills parent[] to store the path.
    const auto bfs = [rank, world_size, n, s, t, &parent, &capacity = std::as_const(capacity)]() -> bool {
        // Mark all the vertices as not visited
        std::vector<bool> visited(n, false);
        // Create a queue for BFS
        std::deque<Vertex> queue{ s };
        // Mark the source node as visited and enqueue it
        visited[s] = true;
        parent[s] = s;

        // Parallel choice a part of all layer 1 vertexes
        Vertex u = queue.front();
        queue.pop_front();
        size_t count = 0;
        for (Vertex v{}; v < n; v++) {
            if (visited[v] || capacity[u][v] <= 0) {
                continue;
            }
            if (count % world_size == rank) {
                queue.push_back(v);
                parent[v] = u;
                visited[v] = true;
            }
            count++;
        }


        // Standard BFS loop
        while (!(queue.empty())) {
            Vertex u = queue.front();
            queue.pop_front();
            // Get all adjacent vertices of the dequeued vertex u
            // If an adjacent has not been visited, then mark it
            // visited and enqueue it
            for (Vertex v{}; v < n; v++) {
                if (visited[v] || capacity[u][v] <= 0) {
                    continue;
                }
                queue.push_back(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
        // If we reached sink in BFS starting from source, then return
        // true, else false
        return visited[t];
        };

    flow_result_t result;
    if (rank == MASTER_RANK) {
        result.max_flow = 0;
        result.flow = adjacency_matrix<>(g.capacity.size(), std::vector<int>(g.capacity.size()));
    }
    bool loop = true;

    while (loop) {
        // Try to find path to the sink and calculate the flow
        loop = bfs();
        if (rank == MASTER_RANK) {
            for (int j = 0; j < n; ++j) {
                std::cout << parent[j] << ' ';
            }
            std::cout << std::endl;
        }

        if (rank == MASTER_RANK) {
            std::cout << "Stage 2a" << std::endl;
        }
        int path_flow_temp = 0;
        if (loop) {
            path_flow_temp = INF;
            for (Vertex v = t; v != s; v = parent[v]) {
                path_flow_temp = std::min(path_flow_temp, capacity[parent[v]][v]);
            }
        }

        // Transmit the results to MASTER node
        std::vector<int> path_flows;
        if (rank == MASTER_RANK) {
            path_flows.resize(world_size);
        }
        MPI::COMM_WORLD.Gather(&path_flow_temp, 1, MPI::INT, path_flows.data(), 1, MPI::INT, MASTER_RANK);
        if (rank == MASTER_RANK) {
            std::cout << "Stage 2b" << std::endl;
        }

        std::vector<int> res;
        if (rank == MASTER_RANK) {
            res.resize(world_size);
        }

        // Calculate the exist path to the stick statement
        MPI::COMM_WORLD.Gather(&loop, 1, MPI::INT, res.data(), 1, MPI::INT, MASTER_RANK);
        if (rank == MASTER_RANK) {
            loop = false;
            for (auto p : res) {
                loop |= p;
            }
        }
        if (rank == MASTER_RANK) {
            std::cout << "Stage 2c" << std::endl;
        }
        // Transmit the exist path to the stick statement to other nodes
        MPI::COMM_WORLD.Bcast(&loop, 1, MPI::INT, MASTER_RANK);

        // Transmit the finded paths to the MASTER node
        path_t respath;
        if (rank == MASTER_RANK) {
            respath.resize(world_size * n);
        }
        MPI::COMM_WORLD.Gather(parent.data(), n, MPI::UNSIGNED_LONG, respath.data(), n, MPI::UNSIGNED_LONG, MASTER_RANK);

        if (rank == MASTER_RANK) {
            std::cout << "Stage 2d" << std::endl;
            for (int j = 0; j < n * world_size; ++j) {
                std::cout << respath[j] << ' ';
            }
            std::cout << std::endl;
        }
        if (rank == MASTER_RANK) {
            int path_flow = 0;
            // Parsing the recive paths
            std::vector<path_t> parents(world_size);
            for (int i = 0; i < world_size; ++i) {
                path_t temp(n);
                std::copy(respath.begin() + i * n, respath.begin() + (i + 1) * n, temp.begin());
                for (int j = 0; j < n; ++j) {
                    std::cout << temp[j] << ' ';
                }
                std::cout << std::endl;
                parents.push_back(temp);
            }
            if (rank == MASTER_RANK) {
                std::cout << "Stage 3a" << std::endl;
                std::cout << path_flows.data() << std::endl;
            }
            // Calculate max flow in the recived paths
            int maxp = 0;
            for (int i = 0; i < world_size; ++i) {           
                if (path_flows[i] > path_flow) {
                    path_flow = path_flow_temp;
                    maxp = i;
                }
            }
            if (rank == MASTER_RANK) {
                std::cout << "Stage 3aa" << std::endl;
            }
            // Modify the adjacency matrix - add the best flow
            for (Vertex v = t; v != s; v = parents[maxp][v]) {
                Vertex u = parents[maxp][v];
                capacity[u][v] -= path_flow;
                capacity[v][u] += path_flow;
                result.flow[u][v] += path_flow;
            }
            result.max_flow += path_flow;
            if (rank == MASTER_RANK) {
                std::cout << "Stage 3b" << std::endl;
            }
        }
        
        // Transmit the modified adjacency matrix to the other nodes
        for (int i = 0; i < n; ++i)
        {
            std::vector<int> row(n);
            if (rank == MASTER_RANK) {
                row = capacity[i];
            }
            MPI::COMM_WORLD.Bcast(row.data(), n, MPI::INT, MASTER_RANK);
            if (rank != MASTER_RANK) {
                std::copy(row.begin(), row.end(), capacity[i].begin());
            }
        }
        if (rank == MASTER_RANK) {
            std::cout << "Stage 3c" << std::endl;
        }
    }
    return result;
}
// end snippet mpi-max_flow_ford_fulkerson