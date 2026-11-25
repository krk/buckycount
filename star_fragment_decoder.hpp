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

// Direct oriented unit-triangle quotient of a valid integral convex atlas.
#pragma once
#include "genfullerene.hpp"

namespace gf {
// Mathematical interface. The patches are occurrences, not geometric sets.
// Every patch corner is a lattice point; hence a unit-triangle interior contains
// no chart vertex. Its positive-area patch fragments are connected precisely
// through seams crossing that open interior. These transitions identify whole
// unit faces, with a cyclic corner offset. An edge-aligned seam instead joins
// opposite unit faces along the common edge; isolated point contacts do neither.
// Thus on a valid oriented integral convex cone atlas this quotient has exactly
// the unit faces and darts of assemble(atlas), without first constructing its
// geometric vertices and then rediscovering the dual edge incidence.
// Complexity after template construction is linear in fragment/link occurrences
// up to inverse-Ackermann factors, NOT a claimed all-input O(output-size) bound.

struct StarFragmentDecoderStats {
    U graphs = 0, vertices = 0, fragment_occurrences = 0;
    U patch_template_hits = 0, patch_templates = 0;
    U seam_template_hits = 0, seam_templates = 0;
    U interior_unions = 0, aligned_links = 0, local_links = 0;
    U self_incident_seams = 0;
};

using StarFragmentKey = std::array<I, 3>;
struct StarFragmentKeyHash {
    std::size_t operator()(const StarFragmentKey &key) const {
        U h = 1469598103934665603ULL;
        for (I x : key)
            h = (h ^ static_cast<U>(x)) * 1099511628211ULL;
        return static_cast<std::size_t>(h);
    }
};
Patch star_fragment_triangle(const StarFragmentKey &key) {
    I x = key[0], y = key[1];
    if (key[2] == 0)
        return {{{x, y}, {x + 1, y}, {x, y + 1}}};
    return {{{x + 1, y + 1}, {x, y + 1}, {x + 1, y}}};
}
StarFragmentKey star_fragment_key(const Patch &triangle) {
    I x = std::min({triangle[0].x, triangle[1].x, triangle[2].x});
    I y = std::min({triangle[0].y, triangle[1].y, triangle[2].y});
    bool lower = std::find(triangle.begin(), triangle.end(), Point{x, y}) != triangle.end();
    StarFragmentKey key{x, y, lower ? 0 : 1};
    auto canonical = star_fragment_triangle(key);
    need(triangle[0] != triangle[1] && triangle[1] != triangle[2] && triangle[0] != triangle[2],
         "collapsed fragment triangle");
    for (Point p : triangle)
        need(std::find(canonical.begin(), canonical.end(), p) != canonical.end(),
             "fragment is not an elementary lattice triangle");
    return key;
}
int star_fragment_edge(const Patch &triangle, Point a, Point b) {
    for (int j = 0; j < 3; ++j)
        if (triangle[j] == a && triangle[(j + 1) % 3] == b)
            return j;
    throw std::runtime_error("fragment seam reverses the unit orientation");
}
int star_fragment_mod3(int x) {
    x %= 3;
    return x < 0 ? x + 3 : x;
}

struct StarFragmentUnionFind {
    std::vector<int> parent, size, offset;
    explicit StarFragmentUnionFind(int count) : parent(count), size(count, 1), offset(count) {
        std::iota(parent.begin(), parent.end(), 0);
    }
    int find(int i) {
        if (parent[i] != i) {
            int old = parent[i];
            parent[i] = find(old);
            offset[i] = star_fragment_mod3(offset[i] + offset[old]);
        }
        return parent[i];
    }
    // Local slot k of i equals local slot k+shift of j.
    void merge(int i, int j, int shift) {
        int a = find(i), b = find(j);
        int delta = star_fragment_mod3(offset[i] - shift - offset[j]);
        if (a == b) {
            need(delta == 0, "nontrivial holonomy inside a unit face");
            return;
        }
        if (size[a] < size[b]) {
            std::swap(a, b);
            delta = star_fragment_mod3(-delta);
        }
        parent[b] = a;
        offset[b] = delta;
        size[a] += size[b];
    }
};

bool star_fragment_open_crossing(Point a, Point b, const Patch &triangle) {
    // interval() includes the triangle boundary. Remove the only positive-length
    // noninterior possibility: the entire segment lies on a triangle-side line.
    for (int j = 0; j < 3; ++j) {
        Point edge = triangle[(j + 1) % 3] - triangle[j];
        if (cross(edge, a - triangle[j]) == 0 && cross(edge, b - a) == 0)
            return false;
    }
    return interval(a, b, triangle).has_value();
}
bool star_fragment_edge_overlap(Point a, Point b, Point c, Point d) {
    if (cross(b - a, c - a) != 0 || cross(b - a, d - a) != 0)
        return false;
    I a0 = a.x, b0 = b.x, c0 = c.x, d0 = d.x;
    if (a.x == b.x) {
        a0 = a.y;
        b0 = b.y;
        c0 = c.y;
        d0 = d.y;
    }
    return std::max(std::min(a0, b0), std::min(c0, d0)) <
           std::min(std::max(a0, b0), std::max(c0, d0));
}

struct StarFragmentLink { int a, sa, b, sb; };
struct StarFragmentUnion { int a, b, shift; };
struct StarFragmentPatchTemplate {
    Polygon polygon;
    std::vector<Patch> triangles;
    std::unordered_map<StarFragmentKey, int, StarFragmentKeyHash> indices;
    std::vector<StarFragmentLink> links;
};
struct StarFragmentSeamTemplate {
    std::vector<StarFragmentUnion> unions;
    std::vector<StarFragmentLink> links;
};
struct StarFragmentPatchPosition { int type, start, arity; };

struct StarFragmentDecoderCache {
    std::map<Polygon, int> patch_ids;
    std::vector<StarFragmentPatchTemplate> patches;
    std::map<std::array<int, 4>, StarFragmentSeamTemplate> seams;

    StarFragmentPatchPosition patch(const Polygon &polygon, StarFragmentDecoderStats &stats) {
        // Orientation-preserving lattice isometries only. The returned start
        // transports original seam slots into the canonical cyclic ordering.
        int start = 0, best_phase = 0, n = static_cast<int>(polygon.size());
        for (int j = 0; j < n; ++j)
            for (int phase = 0; phase < 6; ++phase) {
                for (int k = 1; k < n; ++k) {
                    Point candidate = rotate(polygon[(j + k) % n] - polygon[j], phase);
                    Point old = rotate(polygon[(start + k) % n] - polygon[start], best_phase);
                    if (candidate != old) {
                        if (candidate < old) {
                            start = j;
                            best_phase = phase;
                        }
                        break;
                    }
                }
            }
        Polygon best;
        best.reserve(n);
        for (int k = 0; k < n; ++k)
            best.push_back(rotate(polygon[(start + k) % n] - polygon[start], best_phase));
        if (auto found = patch_ids.find(best); found != patch_ids.end()) {
            ++stats.patch_template_hits;
            return {found->second, start, n};
        }
        StarFragmentPatchTemplate entry;
        entry.polygon = best;
        entry.triangles = fragments(best);
        need(entry.triangles.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
             "too many patch fragments");
        for (int i = 0; i < static_cast<int>(entry.triangles.size()); ++i)
            need(entry.indices.emplace(star_fragment_key(entry.triangles[i]), i).second,
                 "duplicate patch fragment");
        for (int i = 0; i < static_cast<int>(entry.triangles.size()); ++i) {
            const auto &triangle = entry.triangles[i];
            for (int s = 0; s < 3; ++s) {
                Point a = triangle[s], b = triangle[(s + 1) % 3], c = triangle[(s + 2) % 3];
                auto key = star_fragment_key(Patch{a, b, a + b - c});
                auto neighbor = entry.indices.find(key);
                if (neighbor != entry.indices.end() && i < neighbor->second)
                    entry.links.push_back({i, s, neighbor->second,
                                           star_fragment_edge(star_fragment_triangle(key), b, a)});
            }
        }
        need(patches.size() < static_cast<std::size_t>(std::numeric_limits<int>::max()),
             "too many patch templates");
        int id = static_cast<int>(patches.size());
        patch_ids.emplace(std::move(best), id);
        patches.push_back(std::move(entry));
        ++stats.patch_templates;
        return {id, start, n};
    }

    const StarFragmentSeamTemplate &seam(int fi, int i, int gi, int j,
                                         StarFragmentDecoderStats &stats) {
        std::array<int, 4> key{fi, i, gi, j};
        if (auto found = seams.find(key); found != seams.end()) {
            ++stats.seam_template_hits;
            return found->second;
        }
        const auto &f = patches[fi];
        const auto &g = patches[gi];
        Point a = f.polygon[i], b = f.polygon[(i + 1) % f.polygon.size()];
        Point c = g.polygon[j], d = g.polygon[(j + 1) % g.polygon.size()];
        int phase = unit_between(b - a, c - d);
        Point offset = d - rotate(a, phase);
        auto transport = [&](Point p) { return rotate(p, phase) + offset; };
        StarFragmentSeamTemplate entry;
        for (int t = 0; t < static_cast<int>(f.triangles.size()); ++t) {
            const auto &triangle = f.triangles[t];
            if (star_fragment_open_crossing(a, b, triangle)) {
                Patch image{transport(triangle[0]), transport(triangle[1]), transport(triangle[2])};
                auto target_key = star_fragment_key(image);
                auto target = g.indices.find(target_key);
                need(target != g.indices.end(), "missing seam-interior triangle fragment");
                auto canonical = star_fragment_triangle(target_key);
                int shift = 0;
                while (shift < 3 && image[0] != canonical[shift])
                    ++shift;
                need(shift < 3, "invalid cyclic fragment offset");
                for (int s = 0; s < 3; ++s)
                    need(image[s] == canonical[(s + shift) % 3], "fragment orientation changed");
                entry.unions.push_back({t, target->second, shift});
            } else {
                for (int s = 0; s < 3; ++s) {
                    Point p = triangle[s], q = triangle[(s + 1) % 3], r = triangle[(s + 2) % 3];
                    if (!star_fragment_edge_overlap(a, b, p, q))
                        continue;
                    auto target_key = star_fragment_key(
                        Patch{transport(p), transport(q), transport(p + q - r)});
                    auto target = g.indices.find(target_key);
                    need(target != g.indices.end(), "missing edge-aligned seam neighbor");
                    entry.links.push_back({t, s, target->second,
                        star_fragment_edge(star_fragment_triangle(target_key), transport(q), transport(p))});
                }
            }
        }
        ++stats.seam_templates;
        return seams.emplace(key, std::move(entry)).first->second;
    }
};

void validate_star_fragment_graph(const std::vector<Face> &graph, const std::vector<Face> &mate) {
    int n = static_cast<int>(graph.size());
    need(n > 0, "empty fragment graph");
    for (const auto &row : graph)
        need(*std::min_element(row.begin(), row.end()) >= 0 && row[0] != row[1] &&
             row[1] != row[2] && row[0] != row[2], "fragment graph is not simple cubic");
    std::vector<bool> darts(3 * static_cast<std::size_t>(n));
    int pentagons = 0, faces = 0;
    for (int a = 0; a < n; ++a)
        for (int s = 0; s < 3; ++s)
            if (!darts[3 * a + s]) {
                int b = a, t = s, length = 0;
                do {
                    need(!darts[3 * b + t], "fragment graph has an invalid face permutation");
                    darts[3 * b + t] = true;
                    int next_vertex = graph[b][t];
                    t = (mate[b][t] + 1) % 3;
                    b = next_vertex;
                    need(++length <= 6, "fragment graph has a face longer than six");
                } while (b != a || t != s);
                need(length == 5 || length == 6, "fragment graph has a face shorter than five");
                pentagons += length == 5;
                ++faces;
            }
    need(pentagons == 12 && I(n) - 3 * I(n) / 2 + faces == 2,
         "fragment graph is not a fullerene sphere");
    std::vector<bool> seen(n);
    std::vector<int> queue{0};
    seen[0] = true;
    for (std::size_t i = 0; i < queue.size(); ++i)
        for (int b : graph[queue[i]])
            if (!seen[b]) {
                seen[b] = true;
                queue.push_back(b);
            }
    need(static_cast<int>(queue.size()) == n, "fragment graph is disconnected");
}

} // namespace gf
