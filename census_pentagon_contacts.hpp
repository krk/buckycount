// Copyright 2026 Kerem KAT
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include "genfullerene.hpp"
#include <span>

namespace gf {

// The input is the decoder's oriented cubic rotation system. Pair counts at
// the two endpoints of an edge agree. At a cubic vertex every pair of
// incident faces shares exactly one incident edge, so summing over vertices
// counts every edge shared by two pentagons exactly twice.
inline unsigned census_pentagon_contacts(std::span<const Face> graph, std::vector<std::uint8_t> &state) {
    need(graph.size() >= 20 && graph.size() <= 65535 && graph.size() % 2 == 0,
         "invalid fullerene contact graph size");
    const int n = static_cast<int>(graph.size());
    // Low three bits mark traversed outgoing darts. The upper bits count
    // already encountered pentagons at this vertex.
    state.assign(n, 0);
    unsigned pentagons = 0, endpoints = 0;
    for (int start = 0; start < n; ++start) {
        for (int first = 0; first < 3; ++first) {
            if (state[start] & (1u << first)) continue;
            std::array<int, 6> vertices;
            int length = 0, v = start, slot = first;
            do {
                const auto flag = 1u << slot;
                need(length < 6 && !(state[v] & flag), "invalid fullerene face cycle");
                vertices[length++] = v;
                state[v] |= flag;
                const int w = graph[v][slot];
                need(w >= 0 && w < n && w != v, "invalid edge in contact graph");
                const auto& back = graph[w];
                if (back[0] == v) slot = 1;
                else if (back[1] == v) slot = 2;
                else {
                    need(back[2] == v, "nonreciprocal edge in contact graph");
                    slot = 0;
                }
                v = w;
            } while (v != start || slot != first);
            need(length == 5 || length == 6, "face is neither a pentagon nor a hexagon");
            if (length == 5) {
                ++pentagons;
                for (int j = 0; j < length; ++j) {
                    auto& occurrence = state[vertices[j]];
                    endpoints += occurrence >> 3;
                    occurrence += 8;
                }
            }
        }
    }
    need(pentagons == 12 && endpoints % 2 == 0 && endpoints <= 60,
         "inconsistent fullerene contact incidences");
    return endpoints / 2;
}

inline unsigned census_pentagon_contacts(std::span<const Face> graph) {
    std::vector<std::uint8_t> state;
    return census_pentagon_contacts(graph, state);
}

} // namespace gf
