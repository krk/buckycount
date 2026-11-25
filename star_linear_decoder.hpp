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
#include "star_fragment_decoder.hpp"

namespace gf {
namespace star_linear {

// Allocate only the consecutive x slots meeting each unit-height strip
// of the convex chart. Unit padding keeps cells touching a strip boundary.
struct Row {
    I first_x;
    std::vector<std::array<int, 2>> ids;
};

struct Chart {
    I first_y{};
    std::vector<Row> rows;
    int begin{}, end{};

    int lookup(const StarFragmentKey &key) const {
        auto [x, y, orientation] = key;
        if (y < first_y || U(y - first_y) >= rows.size())
            return -1;
        const auto &row = rows[y - first_y];
        if (x < row.first_x || U(x - row.first_x) >= row.ids.size())
            return -1;
        return row.ids[x - row.first_x][orientation];
    }
};

Chart make_chart(const Polygon &p, std::vector<Patch> &triangles) {
    auto y_order = [](Point a, Point b) { return a.y < b.y; };
    auto [bottom, top] = std::minmax_element(p.begin(), p.end(), y_order);
    Chart chart;
    chart.first_y = bottom->y;
    chart.begin = static_cast<int>(triangles.size());
    for (I y = bottom->y; y < top->y; ++y) {
        I lo = std::numeric_limits<I>::max(), hi = std::numeric_limits<I>::min();
        // Only floor and ceiling are needed. Dividing unreduced numerators
        // and denominators avoids one Euclidean algorithm per intersection.
        auto bound = [&](I num, I den) {
            if (den < 0) {
                den = -den;
                num = -num;
            }
            I q = num / den, r = num % den;
            lo = std::min(lo, q - (r < 0));
            hi = std::max(hi, q + (r > 0));
        };
        for (std::size_t i = 0; i < p.size(); ++i) {
            Point a = p[i], b = p[(i + 1) % p.size()];
            if (a.y == y || a.y == y + 1)
                bound(a.x, 1);
            if (a.y != b.y)
                for (I level : {y, y + 1})
                    if (std::min(a.y, b.y) <= level && level <= std::max(a.y, b.y))
                        bound(a.x * (b.y - a.y) + (level - a.y) * (b.x - a.x), b.y - a.y);
        }
        need(lo <= hi, "empty strip inside convex chart");
        chart.rows.push_back({lo - 1, {}});
        auto &row = chart.rows.back();
        row.ids.resize(static_cast<std::size_t>(hi - row.first_x), {-1, -1});
        for (I x = row.first_x; x < hi; ++x)
            for (int side = 0; side < 2; ++side) {
                Patch triangle = star_fragment_triangle({x, y, side});
                if (overlaps(p, triangle)) {
                    need(triangles.size() < std::numeric_limits<int>::max(), "too many fragments");
                    row.ids[x - row.first_x][side] = static_cast<int>(triangles.size());
                    triangles.push_back(triangle);
                }
            }
    }
    chart.end = static_cast<int>(triangles.size());
    return chart;
}

// Strict overlap of a segment with an open triangle. Separating axes are
// the three triangle normals and the segment normal. A collinear seam
// has no interior overlap and is handled by the edge-link test instead.
bool interior_crossing(Point a, Point b, const Patch &t) {
    bool left = false, right = false;
    for (int i = 0; i < 3; ++i) {
        I side = cross(b - a, t[i] - a);
        left |= side > 0;
        right |= side < 0;
        Point edge = t[(i + 1) % 3] - t[i];
        if (std::max(cross(edge, a - t[i]), cross(edge, b - t[i])) <= 0)
            return false;
    }
    return left && right;
}

} // namespace star_linear

// Called only on the twenty triangular patches of one accepted star.
inline std::vector<Face> linear_star_graph(const Atlas &atlas) {
    need(atlas.rotation.empty() && atlas.patches.size() == 20 && atlas.seams.size() == 30,
         "expected a triangulated 22-sided star");
    std::vector<Patch> triangles;
    std::vector<star_linear::Chart> charts;
    int expected = 0;
    for (const auto &patch : atlas.patches) {
        need(patch.size() == 3, "nontriangular star patch");
        I area = polygon_area(patch);
        need(area > 0 && area <= std::numeric_limits<int>::max() / 3 - expected,
             "star mesh exceeds native index range");
        expected += static_cast<int>(area);
        charts.push_back(star_linear::make_chart(patch, triangles));
    }
    std::vector<StarFragmentLink> links;
    for (const auto &chart : charts)
        for (int i = chart.begin; i < chart.end; ++i) {
            const auto &t = triangles[i];
            for (int s = 0; s < 3; ++s) {
                Point a = t[s], b = t[(s + 1) % 3];
                auto key = star_fragment_key(Patch{a, b, a + b - t[(s + 2) % 3]});
                int j = chart.lookup(key);
                if (i < j)
                    links.push_back({i, s, j, star_fragment_edge(triangles[j], b, a)});
            }
        }
    // An arc (i,j,h) says slot k of i equals slot k+h of j.
    struct Arc { int to, shift; };
    std::vector<std::vector<Arc>> same(triangles.size());
    std::array<std::array<bool, 3>, 20> paired{};
    for (auto [f, i, g, j] : atlas.seams) {
        auto claim = [&](int face, int side) {
            need(face >= 0 && face < 20 && side >= 0 && side < 3 && !paired[face][side],
                 "invalid or repeated star seam");
            paired[face][side] = true;
        };
        claim(f, i);
        claim(g, j);
        Point a = atlas.patches[f][i], b = atlas.patches[f][(i + 1) % 3];
        Point c = atlas.patches[g][j], d = atlas.patches[g][(j + 1) % 3];
        int phase = unit_between(b - a, c - d);
        auto move = [&](Point p) { return d + rotate(p - a, phase); };
        for (int k = charts[f].begin; k < charts[f].end; ++k) {
            const auto &t = triangles[k];
            if (star_linear::interior_crossing(a, b, t)) {
                Patch image{move(t[0]), move(t[1]), move(t[2])};
                int h = charts[g].lookup(star_fragment_key(image));
                need(h >= 0, "missing transported star fragment");
                int shift = star_fragment_edge(triangles[h], image[0], image[1]);
                same[k].push_back({h, shift});
                same[h].push_back({k, -shift});
            } else {
                for (int s = 0; s < 3; ++s) {
                    Point p = t[s], q = t[(s + 1) % 3], r = t[(s + 2) % 3];
                    if (star_fragment_edge_overlap(a, b, p, q)) {
                        int h = charts[g].lookup(star_fragment_key(Patch{move(p), move(q), move(p + q - r)}));
                        need(h >= 0, "missing aligned star fragment");
                        links.push_back({k, s, h, star_fragment_edge(triangles[h], move(q), move(p))});
                    }
                }
            }
        }
    }
    // Resolve all identifications in one traversal. Local slot s at k
    // becomes slot s+offset[k] of its connected component.
    std::vector<int> label(triangles.size(), -1), offset(triangles.size()), stack;
    int n = 0;
    for (int root = 0; root < static_cast<int>(triangles.size()); ++root) {
        if (label[root] >= 0)
            continue;
        label[root] = n++;
        stack.push_back(root);
        while (!stack.empty()) {
            int k = stack.back();
            stack.pop_back();
            for (auto [h, shift] : same[k]) {
                int target = star_fragment_mod3(offset[k] - shift);
                if (label[h] < 0) {
                    label[h] = label[k];
                    offset[h] = target;
                    stack.push_back(h);
                } else {
                    need(label[h] == label[k] && offset[h] == target, "inconsistent fragment offsets");
                }
            }
        }
    }
    need(n == expected, "star fragment quotient has wrong area");
    std::vector<Face> graph(n, {-1, -1, -1}), mate = graph;
    for (auto [i, si, j, sj] : links) {
        int a = label[i], b = label[j];
        int s = star_fragment_mod3(si + offset[i]), t = star_fragment_mod3(sj + offset[j]);
        need(a != b && (graph[a][s] < 0 || (graph[a][s] == b && mate[a][s] == t)) &&
             (graph[b][t] < 0 || (graph[b][t] == a && mate[b][t] == s)),
             "inconsistent star fragment adjacency");
        graph[a][s] = b;
        graph[b][t] = a;
        mate[a][s] = t;
        mate[b][t] = s;
    }
    validate_star_fragment_graph(graph, mate);
    return graph;
}

} // namespace gf
