#include "ford_falkerson.h"
#include <algorithm>
#include <cassert>
#include <queue>
#include <utility>
#include <map>
#include <set>


using namespace graphs;

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

// start snippet max_flow_ford_fulkerson
flow_result_t graphs::max_flow_ford_fulkerson(const flow_graph_t& g) {
    const auto s = g.source;
    const auto t = g.sink;

    auto capacity = g.capacity;
    const auto n = g.capacity.size();
    auto parent = std::vector<Vertex>(g.capacity.size(), -1u);

    // Returns true if there is a path from
    // source `s` to sink `t` in residual graph.
    // Also fills parent[] to store the path.
    const auto bfs = [n, s, t, &parent, &capacity = std::as_const(capacity)]() -> bool {
        // Mark all the vertices as not visited
        std::vector<bool> visited(n, false);
        // Create a queue for BFS
        std::deque<Vertex> queue{ s };
        // Mark the source node as visited and enqueue it
        visited[s] = true;
        parent[s] = s;
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
    result.max_flow = 0;
    result.flow = adjacency_matrix<>(g.capacity.size(), std::vector<int>(g.capacity.size()));
    while (bfs()) {
        int path_flow = INF;
        for (Vertex v = t; v != s; v = parent[v]) {
            path_flow = std::min(path_flow, capacity[parent[v]][v]);
        }
        for (Vertex v = t; v != s; v = parent[v]) {
            Vertex u = parent[v];
            capacity[u][v] -= path_flow;
            capacity[v][u] += path_flow;
            result.flow[u][v] += path_flow;
        }
        result.max_flow += path_flow;
    }
    return result;
}
// end snippet max_flow_ford_fulkerson