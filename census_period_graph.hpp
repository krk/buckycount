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
#include "census_rank_selector.hpp"
#include "star_dense_decoder.hpp"

namespace gf::fullerene_rank {
template<class Visit> void visit_family_periods(I m, Point beta, Visit visit) {
    CENSUS_HEAP_SCOPE(owner);
    if (!small_height_detail::admissible(m, beta)) return;
    const CensusCentralStrip strip(beta);
    census_point_owners(m, beta, [&](const auto &v, const auto &x, const auto &y) {
        visit(strip.center(v.u, x.nearest, y.nearest, x.to_nearest, y.to_nearest));
    });
}

const std::vector<Face> &decode_period(I n, const IdentifierCusp &period, star_dense::Workspace &work) {
    CENSUS_HEAP_SCOPE(decode);
    auto tips = CENSUS_HEAP_EVAL(tips, census_reduced_tips(period));
    // Reduction certifies the nearest-tip inequalities. An owned, open
    // period is a triangulated cone sphere. The dense gluing still checks
    // its area, seam matches, and reciprocal simple edges as it builds them.
    auto star = CENSUS_HEAP_EVAL(star, star_boundary<true>(tips, 100000000, std::move(work.boundary)));
    need(star.area == n, "progressive period-to-star area mismatch");
    if (CENSUS_HEAP_EVAL(dense, decode_dense_star<false>(star, work))) {
        work.boundary = std::move(star.boundary);
        return work.graph;
    }
    CENSUS_HEAP_SCOPE(fallback);
    auto &graph = work.graph;
    auto triangulated = triangulate_star(std::move(star));
    graph = linear_star_graph(triangulated.atlas);
    work.boundary = std::move(triangulated.boundary);
    need(graph.size() == static_cast<std::size_t>(n), "progressive graph has incorrect size");
    return graph;
}

} // namespace gf::fullerene_rank
