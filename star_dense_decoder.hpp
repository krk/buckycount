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
#include "census_heap_scopes.hpp"
#include "star_center_decoder.hpp"
#include "star_linear_decoder.hpp"
#include <memory>
#include <span>
#include <cstring>

namespace gf {
namespace star_dense {
constexpr int shape_limit = 16;
constexpr std::array<Point, 3> directions{{{1, -2}, {1, 1}, {-2, 1}}};

// The geometry of a single directed lattice edge and the six center darts
// is translation invariant. Record scanline toggles and exact directed
// crossing order once for each displacement in a fixed square.
struct Mark {
    std::int8_t x, y, side;
    bool boundary;
};
struct Crossing {
    std::int8_t x, y, side, dart;
    std::uint16_t time;
    std::array<std::int8_t, 2> image_x, image_y, image_dart;
};
struct Stencil {
    std::vector<Mark> marks;
    std::vector<Crossing> crossings;
};
inline Stencil edge_stencil(Point d) {
    Stencil out;
    Point b{3 * d.x, 3 * d.y};
    for (I y = std::min<I>(0, d.y) - 1; y <= std::max<I>(0, d.y); ++y) {
        for (int side = 0; side < 2; ++side) {
            int r = side + 1;
            const I cy = 3 * y + r;
            if (cy > std::min<I>(0, b.y) && cy < std::max<I>(0, b.y)) {
                I num = b.x * cy - r * b.y, den = 3 * b.y;
                if (den < 0) {
                    num = -num;
                    den = -den;
                }
                const auto x = star_center_floor(num, den);
                out.marks.push_back({static_cast<std::int8_t>(x + 1), static_cast<std::int8_t>(y),
                                     static_cast<std::int8_t>(side), x * den == num});
            }
            for (I x = std::min<I>(0, d.x) - 1; x <= std::max<I>(0, d.x); ++x) {
                Point a{3 * x + r, cy};
                for (int k = 0; k < 3; ++k) {
                    auto v = side ? -directions[k] : directions[k];
                    const I den = cross(v, b), num = cross(b, a), where = cross(v, -a);
                    if (den <= 0 || num < 0 || num >= den || where > 0 || where < -den) continue;
                    // den <= 9*shape_limit=144. Distinct crossing fractions
                    // differ by at least 1/144^2, so this integer preserves
                    // their exact order, including ties and a crossing at zero.
                    out.crossings.push_back({static_cast<std::int8_t>(x), static_cast<std::int8_t>(y),
                                             static_cast<std::int8_t>(side), static_cast<std::int8_t>(k), static_cast<std::uint16_t>(65535 * num / den),
                                             {static_cast<std::int8_t>(x + y + side), static_cast<std::int8_t>(-y - 1)},
                                             {static_cast<std::int8_t>(-x - 1), static_cast<std::int8_t>(x + y + side)},
                                             {static_cast<std::int8_t>((k + 1) % 3), static_cast<std::int8_t>((k + 2) % 3)}});
                }
            }
        }
    }
    return out;
}
inline const Stencil& stencil(Point d) {
    constexpr int width = 2 * shape_limit + 1;
    thread_local std::array<std::optional<Stencil>, width * width> shapes;
    auto& found = shapes[(d.y + shape_limit) * width + d.x + shape_limit];
    if (!found) found.emplace(edge_stencil(d));
    return *found;
}
struct Position {
    int x, y, side;
};
struct Seam {
    int x, y;
    bool positive;
    Position image(Position p) const {
        // Scale-three rotation by +1 or -1 turn changes r to 3-r.
        if (positive) return {x - p.y - 1, y + p.x + p.y + p.side, 1 - p.side};
        return {x + p.x + p.y + p.side, y - p.x - 1, 1 - p.side};
    }
};
// A closed 22-edge walk with each coordinate increment at most L has
// coordinate span at most 11L. Every candidate source of a stencil hit
// comes from one of (L+2)^2 cells, two sides, and three darts.
static_assert(2 * (11 * shape_limit + 4) * (11 * shape_limit + 4) < 65535);
static_assert(22 * (shape_limit + 2) * (shape_limit + 2) * 6 < 65535);
struct Slot {
    // Zero denotes an absent hit, center or seam. Hit indices, vertex labels
    // and seam numbers are stored one greater than their outside representation.
    std::array<std::uint16_t, 3> first{};
    std::uint16_t label = 0;
    signed char boundary = 0, turn = 0;
    bool toggle = false;
};
struct Hit {
    std::uint16_t next, image, time;
    std::uint8_t dart;
};

struct Workspace {
    Polygon boundary;
    std::vector<Slot> cells;
    std::vector<Hit> hits;
    std::vector<std::uint16_t> reps;
    std::vector<Face> graph;
};

// Every constructed undirected edge installs its two reciprocal darts
// in matching slots. Their face successors form a permutation. Mark its
// cycles in place, then reuse the storage to check connectedness.
inline void check_sphere(const std::vector<Face>& graph, std::span<int> successor) {
    const int n = static_cast<int>(graph.size());
    int pentagons = 0, faces = 0;
    for (int v = 0; v < n; ++v) {
        const auto& row = graph[v];
        need(row[0] != row[1] && row[1] != row[2] && row[0] != row[2], "parallel center edge");
        for (int k = 0; k < 3; ++k) {
            const int start = 3 * v + k;
            if (successor[start] < 0) continue;
            int e = start, length = 0;
            do {
                int next = successor[e];
                need(next >= 0 && ++length <= 6, "invalid center face permutation");
                successor[e] = ~next;
                e = next;
            } while (e != start);
            need(length >= 5, "center face shorter than five");
            ++faces;
            pentagons += length == 5;
        }
    }
    need(pentagons == 12 && n - 3 * n / 2 + faces == 2, "center graph has incorrect Euler characteristic");
    auto seen = successor.subspan(n, n);
    successor[0] = 0;
    seen[0] = 0;
    int tail = 1;
    for (int head = 0; head < tail; ++head)
        for (int w : graph[successor[head]]) {
            if (seen[w] >= 0) continue;
            seen[w] = 0;
            successor[tail++] = w;
        }
    need(tail == n, "center graph is disconnected");
}
} // namespace star_dense

template<bool verify_sphere = true>
inline bool decode_dense_star(const Star& star, star_dense::Workspace &work) {
    using namespace star_dense;
    need(star.boundary.size() == 22 && star.area > 0, "invalid dense star chart");
    I xl = star.boundary[0].x, xh = xl, yl = star.boundary[0].y, yh = yl;
    std::size_t hit_count = 0;
    std::array<Seam, 22> seams;
    std::array<const Stencil*, 22> shapes;
    for (int i = 0; i < 22; ++i) {
        auto a = star.boundary[i], b = star.boundary[(i + 1) % 22];
        auto delta = b - a;
        if (std::abs(delta.x) > shape_limit || std::abs(delta.y) > shape_limit) return false;
        auto target = star.boundary[((i ^ 1) + 1) % 22];
        bool positive = !(i & 1);
        int phase = positive ? 1 : 5;
        need(rotate(delta, phase) == star.boundary[i ^ 1] - target, "invalid rotated star seam");
        auto offset = target - rotate(a, phase);
        // The first star tip is the origin. Translation gauges accepted by
        // the reference center decoder go to the ordinary linear path if
        // this 32-bit local geometry cannot represent them.
        if (std::abs(a.x) > 1000000 || std::abs(a.y) > 1000000 || std::abs(offset.x) > 4000000 || std::abs(offset.y) > 4000000) return false;
        seams[i] = {static_cast<int>(offset.x), static_cast<int>(offset.y), positive};
        shapes[i] = &stencil(delta);
        hit_count += shapes[i]->crossings.size();
        xl = std::min(xl, a.x);
        xh = std::max(xh, a.x);
        yl = std::min(yl, a.y);
        yh = std::max(yh, a.y);
    }
    const int xmin = static_cast<int>(xl) - 2, ymin = static_cast<int>(yl) - 2;
    const int width = static_cast<int>(xh - xl) + 4, height = static_cast<int>(yh - yl) + 4;
    const int size = 2 * width * height; // 22 bounded sides
    // The one-cell seam band decides when the dense layout pays off;
    // the outer padding additionally permits a checked neighbor step.
    if (2 * (width - 2) * (height - 2) > 8 * Wide(star.area) + 64) return false;
    auto &cells = work.cells;
    // Only the active rectangle is read. All-zero bytes represent the
    // empty cell; bulk clearing avoids per-field stores on every graph.
    if (cells.size() < static_cast<std::size_t>(size))
        CENSUS_HEAP_EVAL(dense_cells, cells.resize(size));
    std::fill_n(cells.begin(), size, Slot{});
    auto column = [&](int x, int y, int side) { return ((y - ymin) * width + (x - xmin)) * 2 + side; };
    auto &hits = work.hits;
    if (hits.size() < hit_count)
        CENSUS_HEAP_EVAL(dense_hits, hits.resize(hit_count));
    std::size_t written_hits = 0;
    for (int i = 0; i < 22; ++i) {
        auto a = star.boundary[i];
        const int offset = column(static_cast<int>(a.x), static_cast<int>(a.y), 0);
        for (const auto m : shapes[i]->marks) {
            int at = offset + 2 * (m.y * width + m.x) + m.side;
            cells[at].toggle = !cells[at].toggle;
            if (m.boundary) {
                auto& prev = cells[at - 2];
                need(prev.boundary == 0, "center lies on two star sides");
                prev.boundary = static_cast<std::uint8_t>(i + 1);
            }
        }
        const auto target = star.boundary[((i ^ 1) + 1) % 22];
        const int target_offset = column(static_cast<int>(target.x), static_cast<int>(target.y), 0);
        const int positive = !(i & 1);
        for (const auto& c : shapes[i]->crossings) {
            int at = offset + 2 * (c.y * width + c.x) + c.side;
            auto& first = cells[at].first[c.dart];
            // The image dart touches the twin segment, so its source lies at
            // most one cell outside the polygon box. The two-cell padding
            // covers this source and its neighbor without clipping.
            const int image = target_offset + 2 * (c.image_y[positive] * width + c.image_x[positive]) + (1 - c.side);
            need(image >= 2 * width + 1 && image < size - 2 * width - 1, "seam crossing outside dense box");
            hits[written_hits] = {first, static_cast<std::uint16_t>(image), c.time,
                                  static_cast<std::uint8_t>(c.image_dart[positive])};
            first = static_cast<std::uint16_t>(++written_hits);
        }
    }
    auto &reps = work.reps;
    if (reps.size() < static_cast<std::size_t>(star.area))
        CENSUS_HEAP_EVAL(dense_reps, reps.resize(static_cast<std::size_t>(star.area)));
    int n = 0;
    for (int y = static_cast<int>(yl); y < yh; ++y)
        for (int side = 0; side < 2; ++side) {
            bool inside = false;
            for (int x = xmin; x < xmin + width; ++x) {
                auto& slot = cells[column(x, y, side)];
                inside ^= slot.toggle;
                if (!inside && !slot.boundary) continue;
                if (slot.boundary) {
                    const auto& seam = seams[slot.boundary - 1];
                    auto p = seam.image({x, y, side});
                    need(p.x >= xl - 1 && p.x <= xh && p.y >= yl - 1 && p.y <= yh,
                         "paired center outside dense box");
                    auto& other = cells[column(p.x, p.y, p.side)];
                    need(other.boundary == (((slot.boundary - 1) ^ 1) + 1), "unpaired boundary center");
                    if (other.label) {
                        slot.label = other.label;
                        slot.turn = seam.positive ? 2 : 1;
                        continue;
                    }
                }
                need(n < star.area, "too many centers for dense star area");
                reps[n] = static_cast<std::uint16_t>(column(x, y, side));
                slot.label = static_cast<std::uint16_t>(++n);
            }
            need(!inside, "odd center scanline");
        }
    need(n == star.area, "dense star center quotient changed area");
    const std::array<int, 3> stride{1 - 2 * width, 1, -1};
    auto &graph = work.graph;
    CENSUS_HEAP_EVAL(adjacency, graph.assign(n, Face{-1, -1, -1}));
    std::unique_ptr<int[]> successor;
    if constexpr (verify_sphere) successor = std::make_unique_for_overwrite<int[]>(3 * n);
    for (int id = 0; id < n; ++id) {
        for (int k = 0; k < 3; ++k) {
            if (graph[id][k] >= 0) continue;
            int current = reps[id];
            int slot = k, t = 0;
            for (int crossings = 0;; ++crossings) {
                int best = -1;
                for (int h = cells[current].first[slot] - 1; h >= 0; h = hits[h].next - 1)
                    if (hits[h].time >= t && (best < 0 || hits[h].time < hits[best].time)) best = h;
                if (best < 0) break;
                if (crossings == 8) return false;
                t = hits[best].time;
                current = hits[best].image;
                slot = hits[best].dart;
            }
            const auto& end = cells[current + (current & 1 ? -stride[slot] : stride[slot])];
            int other = end.label - 1, back = slot + end.turn;
            if (back >= 3) back -= 3;
            need(other >= 0 && other != id && graph[other][back] < 0, "inconsistent traced star dart");
            graph[id][k] = other;
            graph[other][back] = id;
            if constexpr (verify_sphere) {
                successor[3 * id + k] = 3 * other + (back + 1) % 3;
                successor[3 * other + back] = 3 * id + (k + 1) % 3;
            }
        }
        if constexpr (!verify_sphere) {
            const auto &row = graph[id];
            need(row[0] != row[1] && row[1] != row[2] && row[0] != row[2], "parallel center edge");
        }
    }
    if constexpr (verify_sphere) check_sphere(graph, std::span(successor.get(), 3 * n));
    return true;
}
template<bool verify_sphere = true>
inline std::optional<std::vector<Face>> try_dense_star_graph(const Star& star) {
    star_dense::Workspace work;
    if (!decode_dense_star<verify_sphere>(star, work)) return {};
    return std::move(work.graph);
}
inline std::vector<Face> graph_from_reduced_star(const std::array<Point, 11>& tips, I limit = 1000000) {
    auto boundary = star_boundary(tips, limit);
    if (auto fast = try_dense_star_graph(boundary)) return std::move(*fast);
    return linear_star_graph(triangulate_star(std::move(boundary)).atlas);
}
} // namespace gf
