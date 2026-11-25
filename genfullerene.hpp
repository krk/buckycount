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

// Exact arithmetic fullerene library: chart geometry, star periods, cap
// inventories and canonical cusp identifiers. C++20, no external dependency.
// Single-translation-unit header library: definitions are not inline, so a
// program includes it from exactly one .cpp file. Private research headers
// include this file and extend namespace gf; genfullerene_cli.hpp then adds
// the census drivers and the command-line interface (genfullerene_main).
#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace gf {
using I = std::int64_t;
using U = std::uint64_t;
__extension__ using Wide = __int128;
void need(bool yes, const char *message) {
    if (!yes)
        throw std::runtime_error(message);
}
int mod6(int k) {
    k %= 6;
    return k < 0 ? k + 6 : k;
}
struct Point {
    I x = 0, y = 0;
    auto operator<=>(const Point &) const = default;
};
Point operator+(Point a, Point b) { return {a.x + b.x, a.y + b.y}; }
Point operator-(Point a, Point b) { return {a.x - b.x, a.y - b.y}; }
Point operator-(Point a) { return {-a.x, -a.y}; }
I cross(Point a, Point b) { return a.x * b.y - a.y * b.x; }
I norm(Point a) { return a.x * a.x + a.x * a.y + a.y * a.y; }
I dot2(Point a, Point b) { return 2 * a.x * b.x + a.x * b.y + a.y * b.x + 2 * a.y * b.y; }
Point conjugate(Point a) { return {a.x + a.y, -a.y}; }
Point rotate(Point a, int k) {
    switch (mod6(k)) {
    case 0:
        return a;
    case 1:
        return {-a.y, a.x + a.y};
    case 2:
        return {-a.x - a.y, a.x};
    case 3:
        return -a;
    case 4:
        return {a.y, -a.x - a.y};
    default:
        return {a.x + a.y, -a.x};
    }
}
int unit_between(Point a, Point b) {
    need(norm(a) > 0, "zero seam vector");
    for (int k = 0; k < 6; ++k)
        if (rotate(a, k) == b)
            return k;
    throw std::runtime_error("seam vectors are not related by a lattice unit");
}
int next(int d) { return 3 * (d / 3) + (d + 1) % 3; }
int previous(int d) { return 3 * (d / 3) + (d + 2) % 3; }
struct Chart {
    std::array<int, 60> paired{}, phase{}, vertex{};
    std::array<std::vector<int>, 12> corners;
    int score = 0;
    explicit Chart(std::array<int, 60> p) : paired(p) {
        for (int d = 0; d < 60; ++d)
            need(p[d] >= 0 && p[d] < 60 && p[d] != d && p[p[d]] == d, "invalid dart pairing");
        vertex.fill(-1);
        int count = 0;
        for (int s = 0; s < 60; ++s)
            if (vertex[s] < 0) {
                need(count < 12, "too many cone vertices");
                int d = s;
                do {
                    need(vertex[d] < 0, "invalid vertex cycle");
                    vertex[d] = count;
                    corners[count].push_back(d);
                    d = next(p[d]);
                } while (d != s);
                score += std::abs(static_cast<int>(corners[count].size()) - 5);
                ++count;
            }
        need(count == 12, "cone map is not a sphere with twelve vertices");
        std::array<bool, 20> reached{};
        std::array<bool, 60> remaining{}, tree{};
        std::array<int, 20> queue{};
        int head = 0, tail = 1;
        reached[0] = true;
        while (head < tail) {
            int f = queue[head++];
            for (int d = 3 * f; d < 3 * f + 3; ++d)
                if (!reached[p[d] / 3]) {
                    reached[p[d] / 3] = true;
                    queue[tail++] = p[d] / 3;
                    tree[std::min(d, p[d])] = true;
                }
        }
        need(tail == 20, "disconnected cone map");
        for (int d = 0; d < 60; ++d)
            remaining[d] = (d < p[d] && !tree[d]);
        std::array<int, 12> residual;
        residual.fill(5);
        for (int step = 0; step < 11; ++step) {
            int edge = -1, sign = 0, leaf = -1;
            for (int v = 0; v < 12 && leaf < 0; ++v) {
                int degree = 0;
                for (int d : corners[v])
                    if (remaining[std::min(d, p[d])]) {
                        ++degree;
                        edge = std::min(d, p[d]);
                        sign = d < p[d] ? 1 : -1;
                    }
                if (degree == 1)
                    leaf = v;
            }
            need(leaf >= 0, "primal cotree is not a tree");
            int value = mod6(sign * residual[leaf]);
            phase[edge] = value;
            phase[p[edge]] = mod6(-value);
            for (int v = 0; v < 12; ++v)
                for (int d : corners[v])
                    if (std::min(d, p[d]) == edge)
                        residual[v] = mod6(residual[v] - (d < p[d] ? 1 : -1) * value);
            remaining[edge] = false;
        }
        for (int v = 0; v < 12; ++v) {
            int sum = 0;
            for (int d : corners[v])
                sum += phase[d];
            need(residual[v] == 0 && mod6(sum) == 5, "wrong cone holonomy");
        }
    }
};
Chart from_rotation(const std::vector<std::vector<int>> &adjacency) {
    need(adjacency.size() == 12, "cone map needs twelve vertices");
    int location[12][12];
    for (auto &row : location)
        std::fill(std::begin(row), std::end(row), -1);
    bool seen[12][12]{};
    std::vector<std::pair<int, int>> darts;
    for (int v = 0; v < 12; ++v)
        for (int w : adjacency[v]) {
            need(w >= 0 && w < 12 && w != v, "invalid cone edge");
            need(std::count(adjacency[v].begin(), adjacency[v].end(), w) == 1,
                 "parallel planar-code edge");
            need(std::count(adjacency[w].begin(), adjacency[w].end(), v) == 1,
                 "nonreciprocal cone edge");
            if (seen[v][w])
                continue;
            int a = v, b = w, face_size = 0;
            do {
                need(!seen[a][b], "repeated face dart");
                seen[a][b] = true;
                location[a][b] = static_cast<int>(darts.size());
                darts.emplace_back(a, b);
                ++face_size;
                auto it = std::find(adjacency[b].begin(), adjacency[b].end(), a);
                int c = adjacency[b][(static_cast<int>(it - adjacency[b].begin()) + 1) %
                                     adjacency[b].size()];
                a = b;
                b = c;
            } while (a != v || b != w);
            need(face_size == 3, "cone faces must be triangular");
        }
    need(darts.size() == 60, "cone map needs sixty darts");
    std::array<int, 60> p{};
    for (int d = 0; d < 60; ++d)
        p[d] = location[darts[d].second][darts[d].first];
    return Chart(p);
}
struct Triangle {
    std::array<Point, 3> side{};
    I area = 0;
    std::array<I, 3> cot{};
    std::array<int, 3> angle_lo{}, angle_hi{};
};
std::pair<int, int> angle_interval(Point current, Point previous_side) {
    Point a = -previous_side, b = conjugate(current);
    Point relative{a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x + a.y * b.y};
    constexpr std::array<Point, 7> rays{
        {{1, 0}, {1, 1}, {0, 1}, {-1, 2}, {-1, 1}, {-2, 1}, {-1, 0}}};
    need(relative.y > 0, "corner angle is not in (0,pi)");
    for (int k = 0; k <= 6; ++k) {
        I sign = cross(rays[k], relative);
        if (sign == 0 && dot2(rays[k], relative) > 0)
            return {k, k};
        if (k < 6 && sign > 0 && cross(rays[k + 1], relative) < 0)
            return {k, k + 1};
    }
    throw std::runtime_error("corner interval was not found");
}
Triangle triangle(std::array<Point, 3> sides) {
    need(sides[0] + sides[1] + sides[2] == Point{}, "triangle does not close");
    Triangle t;
    t.side = sides;
    t.area = cross(sides[0], sides[1]);
    need(t.area > 0, "nonpositive triangle");
    for (int j = 0; j < 3; ++j) {
        t.cot[j] = -dot2(sides[(j + 1) % 3], sides[(j + 2) % 3]);
        auto interval = angle_interval(sides[j], sides[(j + 2) % 3]);
        t.angle_lo[j] = interval.first;
        t.angle_hi[j] = interval.second;
    }
    return t;
}
std::vector<Triangle> triangles(int bound, int target) {
    need(bound > 0 && bound <= 1000000000, "edge norm must be in 1..1000000000");
    int radius = 2 * static_cast<int>(std::sqrt(bound)) + 2;
    std::vector<Point> vectors;
    for (int x = -radius; x <= radius; ++x)
        for (int y = -radius; y <= radius; ++y)
            if (norm({x, y}) > 0 && norm({x, y}) <= bound)
                vectors.push_back({x, y});
    std::vector<Triangle> result;
    for (Point a : vectors)
        for (Point b : vectors) {
            Point c = -(a + b);
            I area = cross(a, b);
            if (area > 0 && area <= target - 19 && norm(c) > 0 && norm(c) <= bound)
                result.push_back(triangle({a, b, c}));
        }
    return result;
}
struct VectorHash {
    std::size_t operator()(const std::vector<I> &a) const {
        U h = 0xcbf29ce484222325ULL;
        for (I x : a) {
            h ^= static_cast<U>(x);
            h *= 0x100000001b3ULL;
        }
        return static_cast<std::size_t>(h);
    }
};
// Automorphisms of the unmarked cone map, lifted to its normalized lattice
// connection. The six common units are handled by normalizing the image of
// the first side. Reversing a face sends its side vector to -conjugate(z).
struct ChartAction {
    std::array<int, 20> source{}, shift{}, rotation{};
    bool reflection = false;
};
std::vector<ChartAction> chart_actions(const Chart &c) {
    std::vector<ChartAction> result;
    for (int reverse = 0; reverse < 2; ++reverse)
        for (int root = 0; root < 60; ++root) {
            if (reverse == 0 && root == 0)
                continue;
            std::array<int, 60> image, inverse, queue{};
            image.fill(-1);
            inverse.fill(-1);
            image[0] = root;
            inverse[root] = 0;
            int head = 0, tail = 1;
            bool valid = true;
            while (head < tail && valid) {
                int d = queue[head++], e = image[d];
                for (auto [a, b] : {std::pair(next(d), reverse ? previous(e) : next(e)),
                                    std::pair(c.paired[d], c.paired[e])}) {
                    if (image[a] >= 0) {
                        valid &= image[a] == b;
                    } else if (inverse[b] >= 0) {
                        valid = false;
                    } else {
                        image[a] = b;
                        inverse[b] = a;
                        queue[tail++] = a;
                    }
                }
            }
            if (!valid)
                continue;
            need(tail == 60, "incomplete chart automorphism");
            std::array<int, 20> gauge, faces{};
            gauge.fill(-1);
            gauge[0] = 0;
            head = 0;
            tail = 1;
            while (head < tail) {
                int f = faces[head++];
                for (int d = 3 * f; d < 3 * f + 3; ++d) {
                    int g = c.paired[d] / 3;
                    int wanted =
                        mod6(gauge[f] + c.phase[image[d]] - (reverse ? -c.phase[d] : c.phase[d]));
                    if (gauge[g] < 0) {
                        gauge[g] = wanted;
                        faces[tail++] = g;
                    } else
                        need(gauge[g] == wanted, "automorphism does not lift to the connection");
                }
            }
            ChartAction action;
            action.reflection = reverse;
            for (int f = 0; f < 20; ++f) {
                int target = image[3 * f] / 3;
                action.source[target] = f;
                action.shift[target] = image[3 * f] % 3;
                action.rotation[target] = gauge[f];
            }
            result.push_back(action);
        }
    return result;
}
using Keys = std::unordered_set<std::vector<I>, VectorHash>;
struct Canonical {
    std::vector<I> key;
    int order = 0;
    bool duplicate = false;
};

struct Rational {
    I n = 0, d = 1;
    Rational(I a = 0, I b = 1) : n(a), d(b) {
        need(d != 0, "zero rational denominator");
        if (d < 0) {
            n = -n;
            d = -d;
        }
        I g = std::gcd(n, d);
        n /= g;
        d /= g;
    }
    auto operator<=>(const Rational &b) const {
        Wide a1 = Wide(n) * b.d, b1 = Wide(b.n) * d;
        return a1 <=> b1;
    }
    bool operator==(const Rational &b) const { return n == b.n && d == b.d; }
    I floor() const {
        I q = n / d;
        return q - (n % d < 0);
    }
    I ceil() const {
        I q = n / d;
        return q + (n % d > 0);
    }
};
using Patch = std::array<Point, 3>;
template <class Polygon> bool overlaps(const Polygon &a, const Patch &b) {
    auto separated = [](const auto &p, const auto &q) {
        for (std::size_t j = 0; j < p.size(); ++j) {
            I maximum = std::numeric_limits<I>::min();
            Point e = p[(j + 1) % p.size()] - p[j];
            for (Point v : q)
                maximum = std::max(maximum, cross(e, v - p[j]));
            if (maximum <= 0)
                return true;
        }
        return false;
    };
    return !separated(a, b) && !separated(b, a);
}
template <class Polygon> std::vector<Patch> fragments(const Polygon &p) {
    std::vector<Patch> result;
    I min_y = p[0].y, max_y = p[0].y;
    for (Point v : p) {
        min_y = std::min(min_y, v.y);
        max_y = std::max(max_y, v.y);
    }
    for (I y = min_y; y < max_y; ++y) {
        std::vector<Rational> xs;
        for (std::size_t j = 0; j < p.size(); ++j) {
            Point a = p[j], b = p[(j + 1) % p.size()];
            if (y <= a.y && a.y <= y + 1)
                xs.emplace_back(a.x);
            if (a.y != b.y)
                for (I level : {y, y + 1})
                    if (std::min(a.y, b.y) <= level && level <= std::max(a.y, b.y))
                        xs.emplace_back(a.x * (b.y - a.y) + (level - a.y) * (b.x - a.x), b.y - a.y);
        }
        if (xs.empty())
            continue;
        auto limits = std::minmax_element(xs.begin(), xs.end());
        for (I x = limits.first->floor() - 1; x < limits.second->ceil(); ++x) {
            Patch a{{{x, y}, {x + 1, y}, {x, y + 1}}}, b{{{x + 1, y + 1}, {x, y + 1}, {x + 1, y}}};
            if (overlaps(p, a))
                result.push_back(a);
            if (overlaps(p, b))
                result.push_back(b);
        }
    }
    return result;
}
std::optional<std::pair<Rational, Rational>> interval(Point a, Point b, const Patch &t) {
    Rational lo(0), hi(1);
    Point direction = b - a;
    for (int j = 0; j < 3; ++j) {
        Point e = t[(j + 1) % 3] - t[j];
        I constant = cross(e, a - t[j]), slope = cross(e, direction);
        if (slope > 0)
            lo = std::max(lo, Rational(-constant, slope));
        else if (slope < 0)
            hi = std::min(hi, Rational(-constant, slope));
        else if (constant < 0)
            return {};
    }
    if (lo < hi)
        return std::pair(lo, hi);
    return {};
}
struct UnionFind {
    std::vector<int> parent, size;
    explicit UnionFind(int n) : parent(n), size(n, 1) {
        std::iota(parent.begin(), parent.end(), 0);
    }
    int find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }
    void merge(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b)
            return;
        if (std::pair(size[a], -a) < std::pair(size[b], -b))
            std::swap(a, b);
        parent[b] = a;
        size[a] += size[b];
    }
};
using Face = std::array<int, 3>;
using Polygon = std::vector<Point>;
struct Atlas {
    std::vector<Polygon> patches;
    std::vector<std::array<int, 4>> seams;
    std::vector<int> rotation;
};
struct Mesh {
    std::vector<Face> faces;
    std::vector<int> rotation;
};
I polygon_area(const Polygon &p) {
    need(p.size() >= 3, "a patch needs at least three corners");
    need(std::set<Point>(p.begin(), p.end()).size() == p.size(), "repeated patch corner");
    I area = 0;
    for (std::size_t j = 0; j < p.size(); ++j) {
        area += cross(p[j], p[(j + 1) % p.size()]);
        need(cross(p[(j + 1) % p.size()] - p[j], p[(j + 2) % p.size()] - p[(j + 1) % p.size()]) >=
                 0,
             "patches must be convex and positively oriented");
        need(p[j] != p[(j + 1) % p.size()], "zero patch side");
        for (Point v : p)
            need(cross(p[(j + 1) % p.size()] - p[j], v - p[j]) >= 0,
                 "patches must be convex and positively oriented");
    }
    need(area > 0, "a patch must have positive area");
    return area;
}
Mesh assemble(const Atlas &atlas) {
    const auto &patches = atlas.patches;
    int nf = static_cast<int>(patches.size());
    need(nf > 0, "empty atlas");
    if (!atlas.rotation.empty()) {
        auto permutation = atlas.rotation;
        std::sort(permutation.begin(), permutation.end());
        need(static_cast<int>(permutation.size()) == nf, "patch rotation is not a permutation");
        for (int i = 0; i < nf; ++i)
            need(permutation[i] == i, "patch rotation is not a permutation");
    }
    I expected = 0;
    std::set<std::pair<int, int>> boundaries;
    for (int f = 0; f < nf; ++f) {
        I area = polygon_area(patches[f]);
        need(area <= std::numeric_limits<int>::max() - expected, "atlas is too large");
        expected += area;
        for (int j = 0; j < static_cast<int>(patches[f].size()); ++j)
            boundaries.emplace(f, j);
    }
    need(expected <= std::numeric_limits<int>::max(), "atlas is too large");
    for (auto seam : atlas.seams)
        for (int j : {0, 2})
            need(boundaries.erase({seam[j], seam[j + 1]}) == 1,
                 "every patch edge must be paired exactly once");
    need(boundaries.empty(), "every patch edge must be paired exactly once");
    std::vector<std::vector<Patch>> pieces(nf);
    std::vector<std::map<Point, int>> indices(nf);
    int tokens = 0;
    for (int f = 0; f < nf; ++f) {
        pieces[f] = fragments(patches[f]);
        for (const auto &t : pieces[f])
            for (Point p : t)
                indices[f].emplace(p, 0);
        for (auto &[p, index] : indices[f]) {
            (void)p;
            index = tokens++;
        }
    }
    UnionFind uf(tokens);
    struct Part {
        Rational lo, hi;
        int index;
    };
    for (auto [f, i, g, j] : atlas.seams) {
        Point a = patches[f][i], b = patches[f][(i + 1) % patches[f].size()], c = patches[g][j],
              e = patches[g][(j + 1) % patches[g].size()];
        int phase = unit_between(e - c, a - b);
        Point offset = b - rotate(c, phase);
        std::vector<Part> left, right;
        for (int k = 0; k < static_cast<int>(pieces[f].size()); ++k)
            if (auto v = interval(a, b, pieces[f][k]))
                left.push_back({v->first, v->second, k});
        for (int k = 0; k < static_cast<int>(pieces[g].size()); ++k)
            if (auto v = interval(c, e, pieces[g][k]))
                right.push_back({Rational(v->second.d - v->second.n, v->second.d),
                                 Rational(v->first.d - v->first.n, v->first.d), k});
        for (auto *parts : {&left, &right}) {
            std::sort(parts->begin(), parts->end(), [](const Part &x, const Part &y) {
                return std::tie(x.lo, x.hi, x.index) < std::tie(y.lo, y.hi, y.index);
            });
            Rational cursor;
            for (const auto &v : *parts) {
                need(v.lo == cursor, "seam intervals overlap or have a gap");
                cursor = v.hi;
            }
            need(cursor == Rational(1), "incomplete seam");
        }
        std::size_t l = 0, r = 0;
        while (l < left.size() && r < right.size()) {
            const auto &x = left[l];
            const auto &y = right[r];
            need(std::max(x.lo, y.lo) < std::min(x.hi, y.hi), "seam sweep lost overlap");
            int common = 0;
            for (Point p : pieces[f][x.index])
                for (Point v : pieces[g][y.index])
                    if (p == rotate(v, phase) + offset) {
                        ++common;
                        uf.merge(indices[f].at(p), indices[g].at(v));
                    }
            need(common >= 2, "incompatible grid fragments");
            bool advance_left = x.hi <= y.hi, advance_right = y.hi <= x.hi;
            l += advance_left;
            r += advance_right;
        }
    }
    std::map<int, int> labels;
    for (int i = 0; i < tokens; ++i)
        labels.emplace(uf.find(i), 0);
    int nv = 0;
    for (auto &[root, label] : labels) {
        (void)root;
        label = nv++;
    }
    std::set<Face> unique;
    for (int f = 0; f < nf; ++f)
        for (const auto &tri : pieces[f]) {
            Face face;
            for (int j = 0; j < 3; ++j)
                face[j] = labels.at(uf.find(indices[f].at(tri[j])));
            std::sort(face.begin(), face.end());
            need(face[0] != face[1] && face[1] != face[2], "collapsed unit triangle");
            unique.insert(face);
        }
    need(static_cast<I>(unique.size()) == expected, "mesh has wrong total area");
    Mesh mesh{{unique.begin(), unique.end()}, {}};
    if (!atlas.rotation.empty()) {
        mesh.rotation.assign(nv, -1);
        for (int f = 0; f < nf; ++f)
            for (const auto &[p, token] : indices[f]) {
                auto other = indices[atlas.rotation[f]].find(p);
                need(other != indices[atlas.rotation[f]].end(),
                     "rotation changes patch coordinates");
                int vertex = labels.at(uf.find(token)), image = labels.at(uf.find(other->second));
                need(mesh.rotation[vertex] < 0 || mesh.rotation[vertex] == image,
                     "patch rotation does not descend");
                mesh.rotation[vertex] = image;
            }
    }
    return mesh;
}
std::vector<Face> reconstruct(const Chart &chart, const std::array<Triangle, 20> &ts,
                              int expected) {
    Atlas atlas;
    for (const auto &t : ts)
        atlas.patches.push_back({Point{}, t.side[0], -t.side[2]});
    for (int d = 0; d < 60; ++d)
        if (d < chart.paired[d]) {
            int q = chart.paired[d];
            atlas.seams.push_back({d / 3, d % 3, q / 3, q % 3});
        }
    auto mesh = assemble(atlas);
    need(static_cast<int>(mesh.faces.size()) == expected, "mesh has wrong total area");
    return mesh.faces;
}
std::vector<Face> graph_from_faces(const std::vector<Face> &faces) {
    int n = static_cast<int>(faces.size()), nv = 0;
    for (auto f : faces)
        for (int v : f)
            nv = std::max(nv, v + 1);
    struct Side {
        int face, slot, sign;
    };
    std::map<std::pair<int, int>, std::vector<Side>> edges;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) {
            int a = faces[i][k], b = faces[i][(k + 1) % 3];
            edges[std::minmax(a, b)].push_back({i, k, a < b ? 1 : -1});
        }
    std::vector<int> degree(nv);
    std::vector<Face> adjacent(n), factors(n);
    std::set<std::pair<int, int>> graph_edges;
    for (const auto &[edge, incident] : edges) {
        need(incident.size() == 2, "nonclosed mesh");
        ++degree[edge.first];
        ++degree[edge.second];
        auto a = incident[0], b = incident[1];
        need(a.face != b.face && graph_edges.insert(std::minmax(a.face, b.face)).second,
             "nonsimple dual graph");
        adjacent[a.face][a.slot] = b.face;
        adjacent[b.face][b.slot] = a.face;
        factors[a.face][a.slot] = factors[b.face][b.slot] = -a.sign * b.sign;
    }
    need(nv - static_cast<int>(edges.size()) + n == 2, "nonspherical mesh");
    int pentagons = 0;
    for (int d : degree) {
        need(d == 5 || d == 6, "wrong mesh degree");
        pentagons += d == 5;
    }
    need(pentagons == 12, "wrong pentagon count");
    std::vector<int> sign(n), queue{0};
    sign[0] = 1;
    for (std::size_t i = 0; i < queue.size(); ++i)
        for (int k = 0; k < 3; ++k) {
            int a = queue[i], b = adjacent[a][k], wanted = sign[a] * factors[a][k];
            if (sign[b] == 0) {
                sign[b] = wanted;
                queue.push_back(b);
            } else
                need(sign[b] == wanted, "nonorientable mesh");
        }
    need(static_cast<int>(queue.size()) == n, "disconnected mesh");
    for (int i = 0; i < n; ++i)
        if (sign[i] < 0)
            std::reverse(adjacent[i].begin(), adjacent[i].end());
    return adjacent;
}
void write_graph(std::ostream &out, const std::vector<Face> &graph) {
    need(graph.size() <= 65535, "planar code supports at most 65535 vertices");
    const int n = static_cast<int>(graph.size());
    if (n <= 255) {
        std::array<char, 1 + 4 * 255> bytes;
        bytes[0] = static_cast<char>(n);
        for (int i = 0; i < n; ++i) {
            const auto &row = graph[i];
            auto *dst = bytes.data() + 1 + 4*i;
            dst[0] = static_cast<char>(row[0] + 1);
            dst[1] = static_cast<char>(row[1] + 1);
            dst[2] = static_cast<char>(row[2] + 1);
            dst[3] = 0;
        }
        out.write(bytes.data(), 1 + 4*n);
    } else {
        // Planar-code words are little-endian. Bound scratch space even for
        // the largest representable graph, and issue one write per block.
        constexpr int rows_per_block = 512;
        std::array<char, 3 + 8 * rows_per_block> bytes;
        bytes[0] = 0;
        bytes[1] = static_cast<char>(n & 255);
        bytes[2] = static_cast<char>(n >> 8);
        for (int start = 0; start < n; start += rows_per_block) {
            const int rows = std::min(n - start, rows_per_block);
            const int prefix = start == 0 ? 3 : 0;
            for (int i = 0; i < rows; ++i) {
                auto *dst = bytes.data() + prefix + 8*i;
                for (int j = 0; j < 3; ++j) {
                    const auto vertex = static_cast<unsigned>(graph[start + i][j] + 1);
                    dst[2*j] = static_cast<char>(vertex & 255);
                    dst[2*j + 1] = static_cast<char>(vertex >> 8);
                }
                dst[6] = dst[7] = 0;
            }
            out.write(bytes.data(), prefix + 8*rows);
        }
    }
    need(bool(out), "graph output failed");
}

// Direct high-order families. The ordered factor pair selects a lattice cap;
// this path does not visit the cone-chart search or an external graph census.
struct Family {
    int sector, a, b;
    std::string name, group;
    std::vector<int> parameters;
    int vertices() const { return (sector == 5 ? 10 : 2 * sector) * a * b; }
};
Family family(int sector, int a, int b) {
    need((sector == 5 || sector == 10 || sector == 12) && a > 0 && b > 0,
         "positive factor pairs and sector 5, 10, or 12 are required");
    need(a <= 1000000 && b <= 1000000 &&
             I(sector == 5 ? 10 : 2 * sector) * a * b <= std::numeric_limits<int>::max(),
         "family size is too large");
    Family f{sector, a, b, "", sector == 12 ? "D6d" : "D5d", {}};
    if (sector == 5) {
        need((a - b) % 2 == 0 && b != a && b != 3 * a,
             "sector 5 excludes opposite parity and the two walls");
        f.group = "D5h";
        if (a > b) {
            f.name = "H2";
            f.parameters = {(a - b) / 2, b};
        } else if (b < 3 * a) {
            f.name = "G1";
            f.parameters = {(b - a) / 2, (3 * a - b) / 2};
        } else {
            f.name = "H1";
            f.parameters = {a, (b - 3 * a) / 2};
        }
    } else if (b < a) {
        f.name = sector == 10 ? "A4" : "F4";
        f.parameters = {a - b, b};
    } else if (b == a) {
        f.name = sector == 10 ? "A1" : "F1";
        f.parameters = {a};
    } else if (b < 2 * a) {
        f.name = sector == 10 ? "I2" : "D3";
        f.parameters = {2 * a - b, b - a};
    } else if (b == 2 * a) {
        f.name = sector == 10 ? "I3" : "D4";
        f.parameters = {a};
    } else if (b < 3 * a) {
        f.name = sector == 10 ? "I1" : "D2";
        f.parameters = {b - 2 * a, 3 * a - b};
    } else if (b == 3 * a) {
        f.name = sector == 10 ? "A2" : "F2";
        f.parameters = {a};
    } else {
        f.name = sector == 10 ? "A5" : "F5";
        f.parameters = {a, b - 3 * a};
    }
    if (f.name == "A1" || f.name == "A2")
        f.group = "Ih";
    return f;
}
std::vector<Family> families(int vertices) {
    need(vertices >= 0, "negative family-list size");
    std::vector<Family> result;
    for (auto [sector, scale] : {std::pair(5, 10), std::pair(10, 20), std::pair(12, 24)}) {
        if (vertices == 0 || vertices % scale)
            continue;
        int value = vertices / scale;
        std::vector<int> divisors;
        for (int d = 1; d <= value / d; ++d)
            if (value % d == 0) {
                divisors.push_back(d);
                if (value / d != d)
                    divisors.push_back(value / d);
            }
        std::sort(divisors.begin(), divisors.end());
        for (int a : divisors) {
            int b = value / a;
            if (sector != 5 || ((a - b) % 2 == 0 && b != a && b != 3 * a))
                result.push_back(family(sector, a, b));
        }
    }
    return result;
}
int modulo(int a, int n) { return (a % n + n) % n; }
Atlas cap_cylinder(int n, Point z, Point w) {
    Atlas atlas;
    for (int i = 0; i < 2 * n; ++i)
        atlas.patches.push_back({Point{}, z, rotate(z, 1)});
    for (int i = 0; i < n; ++i)
        atlas.patches.push_back({Point{}, z, z + w, w});
    for (int i = 0; i < n; ++i) {
        for (auto seam :
             {std::array{i, 2, (i + 1) % n, 0}, std::array{n + i, 2, n + (i + 1) % n, 0},
              std::array{2 * n + i, 1, 2 * n + (i + 1) % n, 3},
              std::array{i, 1, 2 * n + modulo(-i, n), 0}, std::array{n + i, 1, 2 * n + i, 2}})
            atlas.seams.push_back(seam);
    }
    for (int i = 0; i < n; ++i)
        atlas.rotation.push_back((i + 1) % n);
    for (int i = 0; i < n; ++i)
        atlas.rotation.push_back(n + modulo(i - 1, n));
    for (int i = 0; i < n; ++i)
        atlas.rotation.push_back(2 * n + modulo(i - 1, n));
    return atlas;
}
Atlas double_cap(int n, std::array<Point, 4> wedge) {
    Atlas atlas;
    for (int i = 0; i < 2 * n; ++i) {
        atlas.patches.push_back({wedge[0], wedge[1], wedge[2]});
        atlas.patches.push_back({wedge[0], wedge[2], wedge[3]});
    }
    for (int cap = 0; cap < 2; ++cap)
        for (int j = 0; j < 2 * n; ++j)
            atlas.seams.push_back({2 * n * cap + j, 2, 2 * n * cap + (j + 1) % (2 * n), 0});
    for (int j = 0; j < 2 * n; ++j)
        atlas.seams.push_back({j, 1, 2 * n + modulo(-j, 2 * n), 1});
    for (int j = 0; j < 2 * n; ++j)
        atlas.rotation.push_back((j + 2) % (2 * n));
    for (int j = 0; j < 2 * n; ++j)
        atlas.rotation.push_back(2 * n + modulo(j - 2, 2 * n));
    return atlas;
}
Atlas family_atlas(const Family &f, bool catalogue = false) {
    int n = f.sector == 12 ? 6 : 5;
    if (f.sector != 5 && !catalogue)
        return double_cap(n, {Point{}, Point{0, -f.b}, Point{f.a, -2 * f.a}, Point{f.b, -f.b}});
    const auto &name = f.name;
    int p = f.parameters[0], r = f.parameters.size() == 1 ? p : f.parameters[1];
    if (name == "A1" || name == "F1")
        return cap_cylinder(n, {r, 0}, {0, r});
    if (name == "A2" || name == "F2")
        return cap_cylinder(n, {p, p}, {-p, 2 * p});
    if (name == "A4" || name == "F4")
        return cap_cylinder(n, {r, 0}, {-p, r + 2 * p});
    if (name == "A5" || name == "F5")
        return cap_cylinder(n, {p, p}, {-p - r, 2 * p + r});
    if (name == "H2")
        return cap_cylinder(n, {r, 0}, {-p, 2 * p});
    if (name == "H1")
        return cap_cylinder(n, {p, p}, {-r, r});
    if (name == "I3" || name == "D4")
        return double_cap(n, {Point{}, Point{2 * r, 0}, Point{r, r}, Point{0, 2 * r}});
    Point u, v;
    if (name == "G1") {
        u = {p + r, p};
        v = {p, p + r};
    } else if (name == "I2" || name == "D3") {
        u = {p + r, -2 * (p + r)};
        v = {p + 2 * r, -p - 2 * r};
    } else if (name == "I1" || name == "D2") {
        u = {0, -3 * p - 2 * r};
        v = {p + r, -2 * (p + r)};
    } else
        throw std::runtime_error("unknown cap family");
    return double_cap(n, {Point{}, u, v, rotate(u, 1)});
}
std::vector<int> induced_rotation(const Mesh &mesh) {
    if (mesh.rotation.empty())
        return {};
    std::map<Face, int> labels;
    for (int i = 0; i < static_cast<int>(mesh.faces.size()); ++i)
        labels.emplace(mesh.faces[i], i);
    std::vector<int> result;
    for (auto face : mesh.faces) {
        for (int &v : face)
            v = mesh.rotation.at(v);
        std::sort(face.begin(), face.end());
        result.push_back(labels.at(face));
    }
    auto sorted = result;
    std::sort(sorted.begin(), sorted.end());
    for (int i = 0; i < static_cast<int>(sorted.size()); ++i)
        need(sorted[i] == i, "induced rotation is not a permutation");
    return result;
}
void json_value(std::ostream &out, int x) { out << x; }
void json_value(std::ostream &out, Point p) { out << '[' << p.x << ',' << p.y << ']'; }
template <class Range> void json_value(std::ostream &out, const Range &values) {
    out << '[';
    bool first = true;
    for (const auto &v : values) {
        if (!first)
            out << ',';
        first = false;
        json_value(out, v);
    }
    out << ']';
}
void family_json(std::ostream &out, const Family &f) {
    out << "{\"sector\":" << f.sector << ",\"a\":" << f.a << ",\"b\":" << f.b << ",\"name\":\""
        << f.name << "\",\"point_group\":\"" << f.group << "\",\"vertices\":" << f.vertices()
        << ",\"parameters\":";
    json_value(out, f.parameters);
    out << '}';
}
void atlas_json(std::ostream &out, const Atlas &atlas) {
    out << "{\"patches\":";
    json_value(out, atlas.patches);
    out << ",\"seams\":";
    json_value(out, atlas.seams);
    out << ",\"rotation\":";
    json_value(out, atlas.rotation);
    out << '}';
}
void mesh_json(std::ostream &out, const Mesh &mesh) {
    auto graph = graph_from_faces(mesh.faces);
    std::vector<std::array<int, 2>> edges;
    for (int v = 0; v < static_cast<int>(graph.size()); ++v)
        for (int w : graph[v])
            if (v < w)
                edges.push_back({v, w});
    std::sort(edges.begin(), edges.end());
    out << "{\"faces\":";
    json_value(out, mesh.faces);
    out << ",\"dual_rotation\":";
    json_value(out, mesh.rotation);
    out << ",\"rotation\":";
    json_value(out, induced_rotation(mesh));
    out << ",\"edges\":";
    json_value(out, edges);
    out << '}';
}
Polygon read_polygon(std::istream &in) {
    int n = 0;
    need(bool(in >> n) && n >= 3 && n <= 1000, "invalid polygon size");
    Polygon polygon(n);
    for (auto &p : polygon)
        need(bool(in >> p.x >> p.y) && p.x >= -8000000 && p.x <= 8000000 && p.y >= -8000000 &&
                 p.y <= 8000000,
             "invalid patch coordinate");
    (void)polygon_area(polygon);
    return polygon;
}
Atlas read_atlas(std::istream &in) {
    int nf = 0, ns = 0, nr = 0;
    need(bool(in >> nf >> ns >> nr) && nf > 0 && nf <= 1000 && ns > 0 && ns <= 1000000 && nr == nf,
         "invalid atlas header");
    Atlas atlas;
    for (int f = 0; f < nf; ++f)
        atlas.patches.push_back(read_polygon(in));
    atlas.seams.resize(ns);
    atlas.rotation.resize(nr);
    for (auto &seam : atlas.seams)
        for (int &x : seam)
            need(bool(in >> x), "truncated atlas seams");
    for (int &x : atlas.rotation)
        need(bool(in >> x), "truncated atlas rotation");
    return atlas;
}

// Thurston's shortest-star chart. Eleven copies of a distinguished cone point
// alternate with the other eleven cone points. The ten independent Eisenstein
// coordinates specify the polygon directly; no triangle-domain search is used.
struct Star {
    std::array<Point, 11> tips{};
    Polygon boundary;
    Atlas atlas;
    I area = 0;
};
bool on_segment(Point a, Point b, Point p) {
    return cross(b - a, p - a) == 0 && dot2(p - a, p - b) <= 0;
}
bool segments_meet(Point a, Point b, Point c, Point d) {
    I x = cross(b - a, c - a), y = cross(b - a, d - a);
    I z = cross(d - c, a - c), w = cross(d - c, b - c);
    return (((x < 0 && y > 0) || (x > 0 && y < 0)) &&
            ((z < 0 && w > 0) || (z > 0 && w < 0))) ||
           (x == 0 && on_segment(a, b, c)) || (y == 0 && on_segment(a, b, d)) ||
           (z == 0 && on_segment(c, d, a)) || (w == 0 && on_segment(c, d, b));
}
template<bool nearest_verified>
Star star_boundary(const std::array<Point, 11> &tips, I coordinate_limit, Polygon boundary) {
    Star star;
    star.boundary = std::move(boundary);
    star.boundary.clear();
    star.tips = tips;
    need(tips[0] == Point{}, "the first star tip must be zero");
    for (Point p : tips)
        need(p.x >= -coordinate_limit && p.x <= coordinate_limit &&
                 p.y >= -coordinate_limit && p.y <= coordinate_limit,
             "star coordinate too large");
    star.boundary.reserve(22);
    for (int i = 0; i < 11; ++i) {
        Point p = tips[i], next_tip = tips[(i + 1) % 11];
        Point q = p + rotate(next_tip - p, 1);
        const I radius = norm(next_tip - p);
        need(radius > 0, "repeated star boundary vertex");
        if constexpr (!nearest_verified) {
            bool nearest = true;
            for (Point t : tips) nearest &= norm(q - t) >= radius;
            need(nearest, "star seam is not nearest to its root tips");
        }
        // The positive nearest distance already separates every q from every tip.
        for (int j = 0; j < i; ++j)
            need(p != tips[j] && q != star.boundary[2*j + 1], "repeated star boundary vertex");
        star.boundary.push_back(p);
        star.boundary.push_back(q);
        star.area += cross(p, next_tip) - radius;
    }
    need(star.area > 0, "nonpositive star area");
    // On any open side ending at tip p, |x-t|²-|x-p|² is positive for
    // every other tip t: it is affine, nonnegative at q and positive at p.
    // Thus sides from different tips have disjoint interiors, and a q for
    // another pair cannot lie in a side interior (it has two nearest tips).
    // Only the two adjacent sides from a common tip can overlap. At a q
    // the two sides meet at 60 degrees by construction.
    for (int i = 0; i < 11; ++i) {
        Point p = tips[i], before = star.boundary[(2*i + 21) % 22], after = star.boundary[2*i + 1];
        need(cross(p - before, after - p) != 0 || dot2(p - before, after - p) > 0,
             "overlapping adjacent star sides");
    }
    return star;
}
template<bool nearest_verified = false>
Star star_boundary(const std::array<Point, 11> &tips, I coordinate_limit = 1000000) {
    return star_boundary<nearest_verified>(tips, coordinate_limit, {});
}
Star triangulate_star(Star star) {
    const auto &polygon = star.boundary;
    // Exact ear clipping retains every boundary occurrence, including collinear
    // root tips. A vertex on an ear's boundary blocks that ear as well.
    std::vector<int> remaining(22);
    std::iota(remaining.begin(), remaining.end(), 0);
    std::vector<Face> faces;
    while (remaining.size() > 3) {
        bool clipped = false;
        for (std::size_t i = 0; i < remaining.size(); ++i) {
            int a = remaining[(i + remaining.size() - 1) % remaining.size()];
            int b = remaining[i], c = remaining[(i + 1) % remaining.size()];
            if (cross(polygon[b] - polygon[a], polygon[c] - polygon[a]) <= 0)
                continue;
            bool empty = true;
            for (int v : remaining)
                if (v != a && v != b && v != c &&
                    cross(polygon[b] - polygon[a], polygon[v] - polygon[a]) >= 0 &&
                    cross(polygon[c] - polygon[b], polygon[v] - polygon[b]) >= 0 &&
                    cross(polygon[a] - polygon[c], polygon[v] - polygon[c]) >= 0) {
                    empty = false;
                    break;
                }
            if (!empty)
                continue;
            faces.push_back({a, b, c});
            remaining.erase(remaining.begin() + static_cast<std::ptrdiff_t>(i));
            clipped = true;
            break;
        }
        need(clipped, "star triangulation has no positive ear");
    }
    faces.push_back({remaining[0], remaining[1], remaining[2]});
    std::map<std::pair<int, int>, std::pair<int, int>> darts;
    I area = 0;
    for (int f = 0; f < 20; ++f) {
        Polygon triangle;
        for (int v : faces[f])
            triangle.push_back(polygon[v]);
        area += polygon_area(triangle);
        star.atlas.patches.push_back(std::move(triangle));
        for (int j = 0; j < 3; ++j)
            need(darts.emplace(std::pair(faces[f][j], faces[f][(j + 1) % 3]),
                               std::pair(f, j)).second,
                 "repeated star triangulation dart");
    }
    need(area == star.area, "star triangulation changed area");
    for (const auto &[edge, dart] : darts) {
        auto [a, b] = edge;
        if (auto reverse = darts.find({b, a}); reverse != darts.end() && a < b)
            star.atlas.seams.push_back(
                {dart.first, dart.second, reverse->second.first, reverse->second.second});
    }
    for (int i = 0; i < 11; ++i) {
        auto a = darts.at({2 * i, 2 * i + 1});
        auto b = darts.at({2 * i + 1, (2 * i + 2) % 22});
        star.atlas.seams.push_back({a.first, a.second, b.first, b.second});
    }
    need(star.atlas.seams.size() == 30, "wrong star seam count");
    return star;
}
Star star_chart(const std::array<Point, 11> &tips, I coordinate_limit = 1000000) {
    return triangulate_star(star_boundary(tips, coordinate_limit));
}
std::array<Point, 11> read_star_tips(std::istream &in) {
    std::array<Point, 11> tips{};
    for (int i = 1; i < 11; ++i)
        need(bool(in >> tips[i].x >> tips[i].y), "expected ten Eisenstein star coordinates");
    std::string extra;
    need(!(in >> extra), "extra star coordinates");
    for (Point p : tips)
        need(p.x >= -1000000 && p.x <= 1000000 && p.y >= -1000000 && p.y <= 1000000,
             "star coordinate too large");
    return tips;
}
struct StarReduction {
    std::array<Point, 11> input{}, tips{};
    // kind 0: Hurwitz(index, exponent); kind 1: cyclic start and translation.
    std::vector<std::array<int, 3>> word;
    I initial_energy = 0, final_energy = 0, strict_decreases = 0;
    I primitive_moves = 0, cyclic_changes = 0;
};
I star_energy(const std::array<Point, 11> &p) {
    I result = 0;
    for (int i = 0; i < 11; ++i)
        result += norm(p[(i + 1) % 11] - p[i]);
    return result;
}
I star_area(const std::array<Point, 11> &p) {
    I result = 0;
    for (int i = 0; i < 11; ++i)
        result += cross(p[i], p[(i + 1) % 11]) - norm(p[(i + 1) % 11] - p[i]);
    return result;
}
std::array<I, 3> star_gap(const std::array<Point, 11> &p) {
    std::array<I, 3> best{0, 0, 0};
    for (int i = 0; i < 11; ++i) {
        Point edge = p[(i + 1) % 11] - p[i], q = p[i] + rotate(edge, 1);
        for (int j = 0; j < 11; ++j) {
            I value = norm(q - p[j]) - norm(edge);
            if (value < best[0])
                best = {value, i, j};
        }
    }
    return best;
}
void reduction_bound(const std::array<Point, 11> &p) {
    // This explicit envelope keeps all energy, gap and area arithmetic in int64.
    for (Point x : p)
        need(x.x >= -100000000 && x.x <= 100000000 &&
                 x.y >= -100000000 && x.y <= 100000000,
             "star reduction exceeds the exact integer coordinate envelope");
}
void star_hurwitz(StarReduction &r, int i, int exponent) {
    need(i >= 0 && i < 10 && (exponent == 1 || exponent == 2), "invalid recut");
    for (int k = 0; k < exponent; ++k) {
        r.tips[i + 1] = r.tips[(i + 2) % 11] + rotate(r.tips[i] - r.tips[i + 1], 1);
        reduction_bound(r.tips);
    }
    r.word.push_back({0, i, exponent});
    ++r.primitive_moves;
}
void star_start(StarReduction &r, int start) {
    if (start == 0)
        return;
    auto old = r.tips;
    for (int i = 0; i < 11; ++i)
        r.tips[i] = old[(i + start) % 11] - old[start];
    reduction_bound(r.tips);
    r.word.push_back({1, start, 0});
    ++r.cyclic_changes;
}
int star_bubble(const std::array<Point, 11> &p, int length) {
    for (int i = 0; i < 11; ++i)
        if (p[i] == p[(i + length) % 11])
            return i;
    return -1;
}
void star_escape_two(StarReduction &r, int start) {
    star_start(r, start);
    I energy = star_energy(r.tips);
    for (int attempt = 0; attempt < 2; ++attempt) {
        star_hurwitz(r, 0, 1);
        need(star_energy(r.tips) == energy, "two-bubble recut changed energy");
        if (star_gap(r.tips)[0] < 0)
            return;
    }
    throw std::runtime_error("two-bubble escape failed: mirror-free input is required");
}
StarReduction reduce_star(const std::array<Point, 11> &input) {
    StarReduction r;
    r.input = r.tips = input;
    I area = star_area(input);
    need(area > 0, "star reduction requires a positive period");
    r.initial_energy = star_energy(input);
    while (true) {
        auto [gap, i, j] = star_gap(r.tips);
        if (gap < 0) {
            I energy = star_energy(r.tips);
            if (j < i)
                for (int k = static_cast<int>(i) - 1; k >= j; --k)
                    star_hurwitz(r, k, 1);
            else
                for (int k = static_cast<int>(i); k < j - 1; ++k)
                    star_hurwitz(r, k, 2);
            need(star_energy(r.tips) == energy + gap && star_area(r.tips) == area,
                 "strict recut identity failed");
            ++r.strict_decreases;
            continue;
        }
        need(star_bubble(r.tips, 1) < 0, "collision period has a zero slit");
        if (int start = star_bubble(r.tips, 2); start >= 0) {
            star_escape_two(r, start);
            continue;
        }
        int start = star_bubble(r.tips, 3);
        if (start < 0)
            break;
        star_start(r, start);
        I energy = star_energy(r.tips);
        bool escaped = false;
        for (int attempt = 0; attempt < 5; ++attempt) {
            // Cyclically re-root the closed three-slit bubble. This whole
            // two-move block preserves energy; the intermediate state need not.
            star_hurwitz(r, 1, 1);
            star_hurwitz(r, 0, 1);
            need(star_energy(r.tips) == energy, "three-bubble recut changed energy");
            if (star_gap(r.tips)[0] < 0) {
                escaped = true;
                break;
            }
            if (int two = star_bubble(r.tips, 2); two >= 0) {
                star_escape_two(r, two);
                escaped = true;
                break;
            }
        }
        need(escaped, "three-bubble escape failed: mirror-free input is required");
    }
    r.final_energy = star_energy(r.tips);
    need(star_area(r.tips) == area && r.final_energy <= r.initial_energy,
         "star reduction changed its period norm");
    return r;
}

// Fixed-size specializations keep the short-edge domains in registers. The
// dynamic fallback has the same exact search and imposes no shape-count cutoff.
template <int N> struct Bits {
    std::conditional_t<N == 0, std::vector<U>, std::array<U, N>> a{};
    explicit Bits(int words = 0) {
        if constexpr (N == 0)
            a.resize(words);
    }
    int words() const { return static_cast<int>(a.size()); }
    bool empty() const {
        for (U w : a)
            if (w)
                return false;
        return true;
    }
    bool contains(int k) const { return (a[k / 64] >> (k % 64)) & 1; }
    void set(int k) { a[k / 64] |= U(1) << (k % 64); }
    int count() const {
        int n = 0;
        for (U w : a)
            n += std::popcount(w);
        return n;
    }
    template <class F> void each(F fn) const {
        for (int j = 0; j < words(); ++j)
            for (U w = a[j]; w; w &= w - 1)
                fn(64 * j + std::countr_zero(w));
    }
    int intersection_first(const Bits &b) const {
        for (int j = 0; j < words(); ++j)
            if (U w = a[j] & b.a[j])
                return 64 * j + std::countr_zero(w);
        return -1;
    }
};
struct Bounds {
    std::array<int, 4> lo{}, hi{};
};
struct Counters {
    U nodes = 0, rejected = 0, leaves = 0, arc_revisions = 0, support_searches = 0, duplicates = 0,
      invalid = 0, charts = 0, symmetry_actions = 0, symmetry_rejected = 0;
};
template <int N> class Search {
    using Mask = Bits<N>;
    struct State {
        std::array<Mask, 20> domain;
        std::array<Bounds, 20> bounds;
    };
    struct Range {
        std::vector<int> values;
        std::vector<Mask> below, above;
    };
    const std::vector<Triangle> &ts;
    int D, W, target;
    bool prune, angle_pruning, short_diagonals, symmetry;
    Counters &count;
    std::vector<Mask> supports;
    std::array<Range, 4> ranges;
    std::array<State, 21> states;
    std::vector<int> residue;
    std::array<std::vector<int>, 36> transforms;
    std::array<std::vector<int>, 6> normalized, normal_unit;
    std::vector<Mask> symmetry_support;
    std::vector<ChartAction> actions;
    std::vector<std::pair<int, int>> root_constraints;
    const Chart *chart = nullptr;
    std::array<int, 60> support_offset{};
    std::function<bool(const std::array<Triangle, 20> &)> emit;
    static constexpr U all_arcs = (U(1) << 60) - 1;
    static constexpr unsigned all_constraints = (1U << 13) - 1;
    static constexpr unsigned symmetry_bit = 1U << 13;
    int first(const Mask &mask) const {
        for (int j = 0; j < W; ++j)
            if (mask.a[j])
                return 64 * j + std::countr_zero(mask.a[j]);
        return -1;
    }
    // Compare only a known equal prefix. An unresolved earlier coordinate
    // forbids using a later difference to reject a branch.
    bool lex_leader(State &s, U &arcs, unsigned &constraints) {
        for (const auto &a : actions) {
            int root = a.source[0];
            if (s.domain[root].count() != 1)
                continue;
            int u = mod6(normal_unit[3 * a.reflection + a.shift[0]][first(s.domain[root])] -
                         a.rotation[0]);
            for (int f = 0; f < 20; ++f) {
                int base = 3 * a.reflection + a.shift[f];
                int k = mod6(a.rotation[f] + u), source = a.source[f];
                if (source == f && base == 0 && k == 0)
                    continue;
                const auto &map = transforms[6 * base + k];
                if (s.domain[f].count() != 1 || s.domain[source].count() != 1) {
                    Mask left_keep(W), right_keep(W);
                    int maximum = -1, minimum = first(s.domain[f]);
                    s.domain[source].each([&](int o) { maximum = std::max(maximum, map[o]); });
                    s.domain[f].each([&](int o) {
                        if (source == f ? o <= map[o] : o <= maximum)
                            left_keep.set(o);
                    });
                    if (!intersect(s, f, left_keep, left_keep, arcs, constraints))
                        return false;
                    if (source != f) {
                        s.domain[source].each([&](int o) {
                            if (map[o] >= minimum)
                                right_keep.set(o);
                        });
                        if (!intersect(s, source, right_keep, right_keep, arcs, constraints))
                            return false;
                    }
                    break;
                }
                int left = first(s.domain[f]), right = map[first(s.domain[source])];
                if (left < right)
                    break;
                if (left > right)
                    return false;
            }
        }
        return true;
    }
    bool symmetry_propagate(State &s, U &arcs, unsigned &constraints) {
        for (auto [f, base] : root_constraints) {
            int maximum = -1;
            s.domain[f].each([&](int o) { maximum = std::max(maximum, normalized[base][o]); });
            Mask keep(W);
            s.domain[0].each([&](int o) {
                if (o <= maximum)
                    keep.set(o);
            });
            if (!intersect(s, 0, keep, keep, arcs, constraints))
                return false;
            const auto &support = symmetry_support[base * D + first(s.domain[0])];
            if (!intersect(s, f, support, support, arcs, constraints))
                return false;
        }
        return true;
    }
    int value(int option, int axis, bool upper = false) const {
        return axis == 0 ? static_cast<int>(ts[option].area)
                         : (upper ? ts[option].angle_hi[axis - 1] : ts[option].angle_lo[axis - 1]);
    }
    void bounds(State &s, int f) {
        auto &b = s.bounds[f];
        b.lo.fill(std::numeric_limits<int>::max());
        b.hi.fill(std::numeric_limits<int>::min());
        s.domain[f].each([&](int o) {
            for (int a = 0; a < (angle_pruning ? 4 : 1); ++a) {
                b.lo[a] = std::min(b.lo[a], value(o, a));
                b.hi[a] = std::max(b.hi[a], value(o, a, true));
            }
        });
    }
    void changed(State &s, int f, U &arcs, unsigned &constraints) {
        bounds(s, f);
        constraints |= 1;
        if (!actions.empty())
            constraints |= symmetry_bit;
        for (int d = 3 * f; d < 3 * f + 3; ++d) {
            arcs |= U(1) << chart->paired[d];
            if (angle_pruning)
                constraints |= 1U << (chart->vertex[d] + 1);
        }
    }
    bool intersect(State &s, int f, const Mask &a, const Mask &b, U &arcs, unsigned &constraints) {
        bool different = false, nonempty = false;
        for (int j = 0; j < W; ++j) {
            U v = s.domain[f].a[j] & a.a[j] & b.a[j];
            different |= v != s.domain[f].a[j];
            nonempty |= v != 0;
            s.domain[f].a[j] = v;
        }
        if (!nonempty)
            return false;
        if (different)
            changed(s, f, arcs, constraints);
        return true;
    }
    bool trim(State &s, int f, int axis, int lo, int hi, U &arcs, unsigned &constraints) {
        const auto &r = ranges[axis];
        int l = static_cast<int>(std::lower_bound(r.values.begin(), r.values.end(), lo) -
                                 r.values.begin()),
            h = static_cast<int>(std::upper_bound(r.values.begin(), r.values.end(), hi) -
                                 r.values.begin()) -
                1;
        if (l > h)
            return false;
        return intersect(s, f, r.above[l], r.below[h], arcs, constraints);
    }
    bool propagate(State &s, U arcs, unsigned constraints) {
        while (arcs || constraints) {
            if (arcs) {
                int d = std::countr_zero(arcs);
                arcs &= arcs - 1;
                ++count.arc_revisions;
                int q = chart->paired[d], f = d / 3, g = q / 3;
                bool different = false, nonempty = false;
                for (int j = 0; j < W; ++j) {
                    U old = s.domain[f].a[j], keep = old;
                    for (U bits = old; bits; bits &= bits - 1) {
                        int option = 64 * j + std::countr_zero(bits);
                        const auto &support = supports[support_offset[d] + option];
                        bool valid;
                        if (f == g)
                            valid = support.contains(option);
                        else {
                            int &witness = residue[d * D + option];
                            valid = witness >= 0 && s.domain[g].contains(witness);
                            if (!valid) {
                                ++count.support_searches;
                                witness = support.intersection_first(s.domain[g]);
                                valid = witness >= 0;
                            }
                        }
                        if (!valid)
                            keep &= ~(U(1) << (option % 64));
                    }
                    different |= keep != old;
                    nonempty |= keep != 0;
                    s.domain[f].a[j] = keep;
                }
                if (!nonempty)
                    return false;
                if (different)
                    changed(s, f, arcs, constraints);
            } else {
                int constraint = std::countr_zero(constraints);
                constraints &= constraints - 1;
                if (constraint == 13) {
                    if (!symmetry_propagate(s, arcs, constraints) ||
                        !lex_leader(s, arcs, constraints)) {
                        ++count.symmetry_rejected;
                        return false;
                    }
                } else if (constraint == 0) {
                    int lo = 0, hi = 0;
                    for (auto &b : s.bounds) {
                        lo += b.lo[0];
                        hi += b.hi[0];
                    }
                    if (target < lo || target > hi)
                        return false;
                    if (prune)
                        for (int f = 0; f < 20; ++f)
                            if (!trim(s, f, 0, target - hi + s.bounds[f].hi[0],
                                      target - lo + s.bounds[f].lo[0], arcs, constraints))
                                return false;
                } else if (angle_pruning) {
                    int v = constraint - 1, lo = 0, hi = 0;
                    std::array<std::pair<int, int>, 60> old;
                    int stored = 0;
                    for (int d : chart->corners[v]) {
                        auto &b = s.bounds[d / 3];
                        lo += b.lo[d % 3 + 1];
                        hi += b.hi[d % 3 + 1];
                        old[stored++] = {b.lo[d % 3 + 1], b.hi[d % 3 + 1]};
                    }
                    if (lo > 10 || hi < 10)
                        return false;
                    int index = 0;
                    for (int d : chart->corners[v]) {
                        auto b = old[index++];
                        if (!trim(s, d / 3, d % 3 + 1, 10 - hi + b.second, 10 - lo + b.first, arcs,
                                  constraints))
                            return false;
                    }
                }
            }
        }
        return true;
    }
    bool visit(int depth, U arcs, unsigned constraints) {
        State &s = states[depth];
        ++count.nodes;
        if (!propagate(s, arcs, constraints)) {
            ++count.rejected;
            return false;
        }
        int f = -1, minimum = std::numeric_limits<int>::max(), impact = -1;
        for (int i = 0; i < 20; ++i) {
            int size = s.domain[i].count();
            if (size <= 1)
                continue;
            int degree = 0;
            for (int d = 3 * i; d < 3 * i + 3; ++d)
                degree += s.domain[chart->paired[d] / 3].count() > 1;
            if (size < minimum || (size == minimum && degree > impact)) {
                f = i;
                minimum = size;
                impact = degree;
            }
        }
        if (f < 0) {
            std::array<Triangle, 20> chosen;
            I area = 0;
            for (int i = 0; i < 20; ++i) {
                s.domain[i].each([&](int o) { chosen[i] = ts[o]; });
                area += chosen[i].area;
            }
            need(area == target, "wrong propagated area");
            if (!prune)
                for (int d = 0; d < 60; ++d) {
                    int q = chart->paired[d];
                    if (chosen[d / 3].cot[d % 3] * chosen[q / 3].area +
                            chosen[q / 3].cot[q % 3] * chosen[d / 3].area <
                        0)
                        return false;
                }
            ++count.leaves;
            return emit(chosen);
        }
        bool stop = false;
        s.domain[f].each([&](int option) {
            if (stop)
                return;
            states[depth + 1] = s;
            auto &child = states[depth + 1];
            for (U &w : child.domain[f].a)
                w = 0;
            child.domain[f].set(option);
            U pending = 0;
            unsigned dirty = 0;
            changed(child, f, pending, dirty);
            stop = visit(depth + 1, pending, dirty);
        });
        return stop;
    }

  public:
    Search(const std::vector<Triangle> &shapes, int vertices, bool early, bool shorter,
           bool quotient, Counters &counters)
        : ts(shapes), D(static_cast<int>(ts.size())), W(N == 0 ? (D + 63) / 64 : N),
          target(vertices), prune(early), angle_pruning(early), short_diagonals(early && shorter),
          symmetry(quotient), count(counters), supports(static_cast<std::size_t>(108) * D, Mask(W)),
          residue(static_cast<std::size_t>(60) * D, -1) {
        need(D > 0, "empty triangle domain");
        if (symmetry) {
            std::map<std::array<Point, 3>, int> lookup;
            for (int o = 0; o < D; ++o) {
                if (o)
                    need(ts[o - 1].side < ts[o].side, "triangle domain is not lexicographic");
                lookup.emplace(ts[o].side, o);
            }
            for (int base = 0; base < 6; ++base) {
                normalized[base].resize(D, D);
                normal_unit[base].resize(D);
                for (int k = 0; k < 6; ++k) {
                    auto &map = transforms[base * 6 + k];
                    map.resize(D);
                    for (int o = 0; o < D; ++o) {
                        std::array<Point, 3> sides;
                        for (int j = 0; j < 3; ++j)
                            sides[(base % 3 + (base >= 3 ? -j : j) + 3) % 3] =
                                rotate(base >= 3 ? -conjugate(ts[o].side[j]) : ts[o].side[j], k);
                        map[o] = lookup.at(sides);
                        if (map[o] < normalized[base][o]) {
                            normalized[base][o] = map[o];
                            normal_unit[base][o] = k;
                        }
                    }
                }
                for (int threshold = 0; threshold < D; ++threshold) {
                    Mask mask(W);
                    for (int o = 0; o < D; ++o)
                        if (normalized[base][o] >= threshold)
                            mask.set(o);
                    symmetry_support.push_back(std::move(mask));
                }
            }
        }
        for (auto &s : states)
            for (auto &mask : s.domain)
                mask = Mask(W);
        std::array<std::map<Point, std::vector<int>>, 3> by_edge;
        for (int o = 0; o < D; ++o)
            for (int j = 0; j < 3; ++j)
                by_edge[j][ts[o].side[j]].push_back(o);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 6; ++k)
                    for (int o = 0; o < D; ++o) {
                        auto found = by_edge[j].find(-rotate(ts[o].side[i], k));
                        if (found == by_edge[j].end())
                            continue;
                        auto &mask = supports[((i * 3 + j) * 6 + k) * D + o];
                        auto &short_mask = supports[(54 + (i * 3 + j) * 6 + k) * D + o];
                        for (int b : found->second) {
                            I cot = ts[o].cot[i] * ts[b].area + ts[b].cot[j] * ts[o].area;
                            if (!prune || cot >= 0) {
                                mask.set(b);
                                Point alternative =
                                    ts[o].side[(i + 1) % 3] + rotate(ts[b].side[(j + 2) % 3], -k);
                                if (cot != 0 || norm(alternative) >= norm(ts[o].side[i]))
                                    short_mask.set(b);
                            }
                        }
                    }
        for (int axis = 0; axis < (angle_pruning ? 4 : 1); ++axis) {
            auto &r = ranges[axis];
            for (int o = 0; o < D; ++o) {
                r.values.push_back(value(o, axis));
                r.values.push_back(value(o, axis, true));
            }
            std::sort(r.values.begin(), r.values.end());
            r.values.erase(std::unique(r.values.begin(), r.values.end()), r.values.end());
            for (int v : r.values) {
                Mask lo(W), hi(W);
                for (int o = 0; o < D; ++o) {
                    if (value(o, axis) <= v)
                        lo.set(o);
                    if (value(o, axis, true) >= v)
                        hi.set(o);
                }
                r.below.push_back(std::move(lo));
                r.above.push_back(std::move(hi));
            }
        }
    }
    bool run(const Chart &c, std::function<bool(const std::array<Triangle, 20> &)> callback) {
        chart = &c;
        emit = std::move(callback);
        actions = symmetry ? chart_actions(c) : std::vector<ChartAction>{};
        count.symmetry_actions += actions.size();
        root_constraints.clear();
        for (const auto &a : actions)
            if (a.source[0] != 0)
                root_constraints.emplace_back(a.source[0], 3 * a.reflection + a.shift[0]);
        std::sort(root_constraints.begin(), root_constraints.end());
        root_constraints.erase(std::unique(root_constraints.begin(), root_constraints.end()),
                               root_constraints.end());
        std::fill(residue.begin(), residue.end(), -1);
        bool adjacent[12][12]{};
        for (int d = 0; d < 60; ++d)
            adjacent[c.vertex[d]][c.vertex[c.paired[d]]] = true;
        for (int d = 0; d < 60; ++d) {
            int q = c.paired[d], a = c.vertex[previous(d)], b = c.vertex[previous(q)],
                u = c.vertex[d], v = c.vertex[q];
            // A shorter cocircular flip decreases a finite integer potential.
            // Permit this rejection only when the flipped map remains simple
            // at this diagonal and neither endpoint loses its second corner.
            bool can_flip = short_diagonals && a != b && u != v && a != u && a != v && b != u &&
                            b != v && !adjacent[a][b] && c.corners[u].size() > 2 &&
                            c.corners[v].size() > 2;
            support_offset[d] = ((d % 3 * 3 + q % 3) * 6 + c.phase[d] + 54 * can_flip) * D;
        }
        for (int f = 0; f < 20; ++f) {
            auto &mask = states[0].domain[f];
            std::fill(mask.a.begin(), mask.a.end(), 0);
            for (int o = 0; o < D; ++o) {
                bool allowed = true;
                if (f == 0)
                    for (int k = 1; k < 6; ++k)
                        if (rotate(ts[o].side[0], k) < ts[o].side[0])
                            allowed = false;
                if (f == 0 && allowed)
                    for (const auto &a : actions)
                        if (a.source[0] == 0 && normalized[3 * a.reflection + a.shift[0]][o] < o)
                            allowed = false;
                if (allowed)
                    mask.set(o);
            }
            if (mask.empty())
                return false;
            bounds(states[0], f);
        }
        return visit(0, all_arcs,
                     (angle_pruning ? all_constraints : 1) | (actions.empty() ? 0 : symmetry_bit));
    }
};

enum class InventoryBackend { Auto, Projective, Field, Chain, Dense, Alcove };
const char *inventory_backend_name(InventoryBackend backend) {
    return backend == InventoryBackend::Auto ? "auto" : backend == InventoryBackend::Alcove ? "alcove" : backend == InventoryBackend::Field ? "field" : backend == InventoryBackend::Projective ? "projective" : backend == InventoryBackend::Dense ? "dense" : "chain";
}
struct Options {
    int vertices = 20, bound = 4;
    std::size_t chart_limit = 7595;
    std::string charts, output = "fullerenes.planar", witnesses, replay;
    std::optional<U> expected;
    bool exhaust = false, prune = true, short_diagonals = true, symmetry = true;
    std::optional<Family> cap;
    std::optional<int> family_list;
    std::optional<Point> inventory_divisor;
    std::optional<int> pack_caps_maximum;
    std::string atlas_input, polygon_input, star_input, reduce_star_input;
    std::string cap_ball_input, cap_orbit_input, cap_orbits_input, cusp_normalize_input, cusp_transport_input;
    std::string cusp_competitors_input, cusp_mirror_input, cusp_parabolic_input, cusp_canonical_input;
    std::string decode_identifier_input, cusp_fibers_input;
    bool catalogue_atlas = false, atlas_only = false, symmetric_only = false;
    bool arithmetic = false, identifiers = false, exhaust_heights = false;
    bool cap_cache = true;
    bool cap_chain = true;
    bool full_mirror_checks = false;
    bool fiber_sieve = true;
    bool alcove_tables = false, alcove_table_summary = false, full_alcove_root_checks = false;
    InventoryBackend inventory_backend = InventoryBackend::Auto;
    bool height_pruning = true, central_quotient = true;
    int maximum_height = 0;
};

// Finite cusp coordinates: one hyperbolic plane and two Eisenstein E8 caps.
// All constants below are independently checked in cusp_star_basis.json.
Point multiply(Point a, Point b) {
    return {a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x + a.y * b.y};
}
using Cap = std::array<Point, 4>;
Point cap_pair(const Cap &a, const Cap &b) {
    Point result{};
    for (int j = 0; j < 4; ++j) {
        Point row{3 * b[j].x, 3 * b[j].y};
        if (j)
            row = row - multiply({1, 1}, b[j - 1]);
        if (j < 3)
            row = row - multiply({2, -1}, b[j + 1]);
        result = result + multiply(conjugate(a[j]), row);
    }
    return result;
}
Point divide_theta(Point p) {
    Point q = multiply(p, {2, -1});
    need(q.x % 3 == 0 && q.y % 3 == 0, "nonintegral cap functional");
    return {q.x / 3, q.y / 3};
}
I cap_index(const Cap &a) {
    Point q = cap_pair(a, a);
    need(q.y == 0 && q.x >= 0 && q.x % 3 == 0, "invalid cap norm");
    return q.x / 3;
}
// Checked arithmetic for the identifier primitives. The older graph kernel keeps
// its existing envelope; no unchecked Point operation is used on these inputs.
Wide identifier_add(Wide a, Wide b) {
    Wide c;
    need(!__builtin_add_overflow(a, b, &c), "identifier arithmetic overflow");
    return c;
}
Wide identifier_sub(Wide a, Wide b) {
    Wide c;
    need(!__builtin_sub_overflow(a, b, &c), "identifier arithmetic overflow");
    return c;
}
Wide identifier_mul(Wide a, Wide b) {
    if (a == static_cast<I>(a) && b == static_cast<I>(b))
        return static_cast<Wide>(static_cast<I>(a)) * static_cast<I>(b);
    Wide c;
    need(!__builtin_mul_overflow(a, b, &c), "identifier arithmetic overflow");
    return c;
}
I identifier_narrow(Wide a) {
    need(a >= std::numeric_limits<I>::min() && a <= std::numeric_limits<I>::max(),
         "identifier coordinate overflow");
    return static_cast<I>(a);
}
Wide identifier_floor(Wide a, Wide b) {
    need(b > 0, "identifier denominator must be positive");
    if (a == static_cast<I>(a) && b == static_cast<I>(b)) {
        const I n = static_cast<I>(a), d = static_cast<I>(b);
        return n / d - (n % d < 0);
    }
    return a / b - (a % b < 0);
}
Wide identifier_ceil(Wide a, Wide b) {
    need(b > 0, "identifier denominator must be positive");
    if (a == static_cast<I>(a) && b == static_cast<I>(b)) {
        const I n = static_cast<I>(a), d = static_cast<I>(b);
        return n / d + (n % d > 0);
    }
    return a / b + (a % b > 0);
}
Wide identifier_sqrt(Wide n) {
    need(n >= 0, "negative identifier square root");
    if (n < 2)
        return n;
    if (n <= (Wide{1} << 52)) {
        const auto value = static_cast<U>(n);
        auto x = static_cast<U>(std::sqrt(static_cast<double>(value)));
        // Only the integer inequalities decide the result. Both squares fit in U.
        while (x * x > value) --x;
        while ((x + 1) * (x + 1) <= value) ++x;
        return x;
    }
    // The initial power of two is an upper bound; integer Newton iteration
    // returns the floor without a floating-point estimate, including at 2^126.
    Wide x = 1;
    for (Wide t = n; t >>= 2;)
        x <<= 1;
    x <<= 1;
    while (true) {
        Wide y = (x + n / x) / 2;
        if (y >= x)
            return x;
        x = y;
    }
}
using IdentifierPoint = std::array<Wide, 2>;
IdentifierPoint identifier_point(Point p) { return {p.x, p.y}; }
IdentifierPoint identifier_add(IdentifierPoint a, IdentifierPoint b) {
    return {identifier_add(a[0], b[0]), identifier_add(a[1], b[1])};
}
IdentifierPoint identifier_sub(IdentifierPoint a, IdentifierPoint b) {
    return {identifier_sub(a[0], b[0]), identifier_sub(a[1], b[1])};
}
IdentifierPoint identifier_mul(IdentifierPoint a, IdentifierPoint b) {
    constexpr Wide limit = Wide{1} << 29;
    auto small = [](Wide x) { return x >= -limit && x < limit; };
    if (small(a[0]) && small(a[1]) && small(b[0]) && small(b[1])) {
        // Each resulting coordinate is a sum of at most three products.
        return identifier_point(multiply({static_cast<I>(a[0]), static_cast<I>(a[1])},
                                         {static_cast<I>(b[0]), static_cast<I>(b[1])}));
    }
    return {identifier_sub(identifier_mul(a[0], b[0]), identifier_mul(a[1], b[1])),
            identifier_add(identifier_add(identifier_mul(a[0], b[1]), identifier_mul(a[1], b[0])),
                           identifier_mul(a[1], b[1]))};
}
IdentifierPoint identifier_conjugate(IdentifierPoint a) {
    return {identifier_add(a[0], a[1]), identifier_sub(0, a[1])};
}
Wide identifier_norm(IdentifierPoint a) {
    return identifier_add(identifier_add(identifier_mul(a[0], a[0]), identifier_mul(a[0], a[1])),
                          identifier_mul(a[1], a[1]));
}
Point identifier_narrow(IdentifierPoint a) {
    return {identifier_narrow(a[0]), identifier_narrow(a[1])};
}
IdentifierPoint identifier_cap_pair(const Cap &a, const Cap &b) {
    auto small = [](const Cap &cap) {
        constexpr I limit = I{1} << 24;
        for (Point p : cap)
            if (p.x < -limit || p.x >= limit || p.y < -limit || p.y >= limit)
                return false;
        return true;
    };
    // Each row coefficient is bounded by 9M. Four complex products of
    // a conjugated input (coordinates at most 2M) fit within 216M².
    if (small(a) && small(b)) return identifier_point(cap_pair(a, b));
    IdentifierPoint result{};
    for (int j = 0; j < 4; ++j) {
        IdentifierPoint row{identifier_mul(3, b[j].x), identifier_mul(3, b[j].y)};
        if (j)
            row = identifier_sub(row, identifier_mul({1, 1}, identifier_point(b[j - 1])));
        if (j < 3)
            row = identifier_sub(row, identifier_mul({2, -1}, identifier_point(b[j + 1])));
        result = identifier_add(result,
                                identifier_mul(identifier_conjugate(identifier_point(a[j])), row));
    }
    return result;
}
Wide identifier_cap_index(const Cap &a) {
    auto q = identifier_cap_pair(a, a);
    need(q[1] == 0 && q[0] >= 0 && q[0] % 3 == 0, "invalid identifier cap norm");
    return q[0] / 3;
}
// Fractions occur only while deriving the fixed decomposition, never in the
// enumeration loop. All numerator operations are checked before narrowing.
Rational identifier_fraction_add(Rational a, Rational b) {
    return {identifier_narrow(identifier_add(identifier_mul(a.n, b.d), identifier_mul(b.n, a.d))),
            identifier_narrow(identifier_mul(a.d, b.d))};
}
Rational identifier_fraction_mul(Rational a, Rational b) {
    return {identifier_narrow(identifier_mul(a.n, b.n)), identifier_narrow(identifier_mul(a.d, b.d))};
}
struct CapDecomposition {
    // scale*Q(v) = sum_i weight[i]*(sum_{j>=i} row[i][j]*v[j])^2.
    std::array<std::array<I, 8>, 8> row{};
    std::array<I, 8> weight{};
    I scale = 1;
    CapDecomposition() {
        std::array<Cap, 8> basis{};
        for (int i = 0; i < 8; ++i)
            (i % 2 ? basis[i][i / 2].y : basis[i][i / 2].x) = 1;
        std::array<std::array<Rational, 8>, 8> upper{};
        std::array<Rational, 8> diagonal{}, factors{};
        for (int i = 0; i < 8; ++i) {
            for (int j = i; j < 8; ++j) {
                Cap sum = basis[i];
                for (int k = 0; k < 4; ++k)
                    sum[k] = sum[k] + basis[j][k];
                Rational value(cap_index(sum) - cap_index(basis[i]) - cap_index(basis[j]), 2);
                for (int k = 0; k < i; ++k) {
                    auto term = identifier_fraction_mul(diagonal[k],
                                 identifier_fraction_mul(upper[k][i], upper[k][j]));
                    value = identifier_fraction_add(value, {-term.n, term.d});
                }
                if (j == i) {
                    need(value.n > 0, "cap decomposition is not positive definite");
                    diagonal[i] = value;
                    upper[i][i] = 1;
                } else
                    upper[i][j] = identifier_fraction_mul(value, {diagonal[i].d, diagonal[i].n});
            }
            I denominator = 1;
            for (int j = i; j < 8; ++j)
                denominator = identifier_narrow(identifier_mul(denominator / std::gcd(denominator, upper[i][j].d),
                                                               upper[i][j].d));
            for (int j = i; j < 8; ++j)
                row[i][j] = identifier_narrow(identifier_mul(upper[i][j].n, denominator / upper[i][j].d));
            factors[i] = {diagonal[i].n,
                          identifier_narrow(identifier_mul(diagonal[i].d, identifier_mul(denominator, denominator)))};
            scale = identifier_narrow(identifier_mul(scale / std::gcd(scale, factors[i].d), factors[i].d));
        }
        for (int i = 0; i < 8; ++i)
            weight[i] = identifier_narrow(identifier_mul(factors[i].n, scale / factors[i].d));
        // Independently check every coefficient against the Hermitian form.
        for (int i = 0; i < 8; ++i)
            for (int j = i; j < 8; ++j) {
                Wide coefficient = 0;
                for (int k = 0; k <= i; ++k)
                    coefficient = identifier_add(coefficient,
                                  identifier_mul(weight[k], identifier_mul(row[k][i], row[k][j])));
                auto pairing = identifier_cap_pair(basis[i], basis[j]);
                need(identifier_mul(6, coefficient) == identifier_mul(scale, identifier_add(2 * pairing[0], pairing[1])),
                     "cap decomposition coefficient mismatch");
            }
    }
};
const CapDecomposition &cap_decomposition() {
    static const CapDecomposition decomposition;
    return decomposition;
}
struct CapBall {
    const CapDecomposition &decomposition = cap_decomposition();
    std::array<Wide, 8> offset{}, center{}, residual{}, coordinate{};
    Wide denominator = 1, limit = 0;
    U visited_nodes = 0, count = 0;
    CapBall(Point beta, const Cap &numerator, I center_denominator, I radius_numerator,
            I radius_denominator) {
        need(center_denominator > 0 && radius_denominator > 0,
             "cap-ball denominators must be positive");
        need(radius_numerator >= 0, "cap-ball radius must be nonnegative");
        Wide D = identifier_norm(identifier_point(beta));
        need(D > 0, "cap-ball beta must be nonzero");
        denominator = identifier_mul(center_denominator, D);
        // Q(beta*y-c/q) = Q(q*D*y-conj(beta)*c)/(q^2*D).
        // Cancel rational factors before multiplication to enlarge the exact
        // signed-128 envelope. The integer inequality includes every equality.
        std::array<Wide, 5> factors{decomposition.scale, radius_numerator,
                                  center_denominator, center_denominator, D};
        I divisor = radius_denominator;
        for (auto &factor : factors) {
            I common = std::gcd(static_cast<I>(factor % divisor), divisor);
            factor /= common;
            divisor /= common;
        }
        limit = 1;
        for (Wide factor : factors)
            limit = identifier_mul(limit, factor);
        limit /= divisor;
        for (int j = 0; j < 4; ++j) {
            auto product = identifier_mul(identifier_conjugate(identifier_point(beta)),
                                          identifier_point(numerator[j]));
            for (int k = 0; k < 2; ++k) {
                int i = 2 * j + k;
                offset[i] = identifier_floor(product[k], denominator);
                center[i] = identifier_sub(product[k], identifier_mul(offset[i], denominator));
            }
        }
    }
    Cap current() const {
        Cap value;
        for (int j = 0; j < 4; ++j)
            value[j] = {identifier_narrow(identifier_add(offset[2 * j], coordinate[2 * j])),
                        identifier_narrow(identifier_add(offset[2 * j + 1], coordinate[2 * j + 1]))};
        return value;
    }
    Wide row_tail(int i) const {
        Wide value = identifier_sub(0, identifier_mul(decomposition.row[i][i], center[i]));
        for (int j = i + 1; j < 8; ++j)
            value = identifier_add(value, identifier_mul(decomposition.row[i][j], residual[j]));
        return value;
    }
    Cap babai() {
        for (int i = 7; i >= 0; --i) {
            Wide step = identifier_mul(decomposition.row[i][i], denominator);
            Wide numerator = identifier_sub(0, row_tail(i));
            coordinate[i] = identifier_floor(numerator, step);
            Wide remainder = identifier_sub(numerator, identifier_mul(coordinate[i], step));
            if (remainder > step / 2)
                coordinate[i] = identifier_add(coordinate[i], 1);
            residual[i] = identifier_sub(identifier_mul(denominator, coordinate[i]), center[i]);
        }
        return current();
    }
    template <class Callback> bool enumerate(Callback &&callback) {
        // Orthogonal completed squares are nonnegative. At a node, every
        // omitted coordinate contributes >=0, so the integer interval is a
        // necessary and sufficient bound for this one square. Descending all
        // eight coordinates therefore visits exactly the closed ellipsoid.
        bool continuing = true;
        auto visit = [&](auto &&self, int i, Wide spent) -> void {
            if (!continuing || spent > limit)
                return;
            if (i < 0) {
                need(count < std::numeric_limits<U>::max(), "cap-ball count overflow");
                ++count;
                if constexpr (std::is_same_v<std::invoke_result_t<Callback, const Cap &>, bool>)
                    continuing = callback(current());
                else
                    callback(current());
                return;
            }
            Wide radius = identifier_sqrt((limit - spent) / decomposition.weight[i]);
            Wide tail = row_tail(i), step = identifier_mul(decomposition.row[i][i], denominator);
            Wide low = identifier_ceil(identifier_sub(identifier_sub(0, radius), tail), step);
            Wide high = identifier_floor(identifier_sub(radius, tail), step);
            for (Wide value = low; value <= high;) {
                need(visited_nodes < std::numeric_limits<U>::max(), "cap-ball node count overflow");
                ++visited_nodes;
                coordinate[i] = value;
                residual[i] = identifier_sub(identifier_mul(denominator, value), center[i]);
                Wide linear = identifier_add(identifier_mul(step, value), tail);
                Wide contribution = identifier_mul(decomposition.weight[i], identifier_mul(linear, linear));
                self(self, i - 1, identifier_add(spent, contribution));
                if (!continuing || value == high)
                    break;
                value = identifier_add(value, 1);
            }
        };
        visit(visit, 7, 0);
        return continuing;
    }
};
struct NearestCap {
    Cap representative{}, translation{};
    Wide index = 0;
    U ties = 0, visited_nodes = 0;
};
NearestCap nearest_identifier_cap(const Cap &z, Point beta) {
    CapBall ball(beta, z, 1, 0, 1);
    NearestCap result;
    auto candidate = [&](const Cap &y) {
        Cap representative, translation;
        for (int j = 0; j < 4; ++j) {
            representative[j] = identifier_narrow(identifier_sub(identifier_point(z[j]),
                                  identifier_mul(identifier_point(beta), identifier_point(y[j]))));
            translation[j] = identifier_narrow(identifier_sub(IdentifierPoint{}, identifier_point(y[j])));
        }
        Wide index = identifier_cap_index(representative);
        return std::tuple(index, representative, translation);
    };
    auto [seed_index, seed, seed_translation] = candidate(ball.babai());
    result.index = seed_index;
    result.representative = seed;
    result.translation = seed_translation;
    Wide factor = identifier_mul(ball.decomposition.scale, identifier_norm(identifier_point(beta)));
    ball.limit = identifier_mul(factor, seed_index);
    ball.enumerate([&](const Cap &y) {
        auto [index, representative, translation] = candidate(y);
        if (index < result.index) {
            result.index = index;
            result.ties = 0;
            ball.limit = identifier_mul(factor, index);
        }
        if (index == result.index) {
            ++result.ties;
            if (representative < result.representative || result.ties == 1) {
                result.representative = representative;
                result.translation = translation;
            }
        }
    });
    need(result.ties > 0, "nearest cap enumeration is empty");
    result.visited_nodes = ball.visited_nodes;
    return result;
}
struct NormalizedCusp {
    Point u{}, beta{};
    std::array<NearestCap, 2> caps;
    I central_shift = 0, norm_index = 0;
};
NormalizedCusp normalize_identifier_cusp(Point u, Point beta, const std::array<Cap, 2> &z) {
    Wide D = identifier_norm(identifier_point(beta));
    need(D > 0, "cusp beta must be nonzero");
    NormalizedCusp result;
    result.beta = beta;
    auto original_ratio = identifier_mul(identifier_point(u), identifier_conjugate(identifier_point(beta)));
    Wide original_index = identifier_sub(original_ratio[0],
                                        identifier_add(identifier_cap_index(z[0]), identifier_cap_index(z[1])));
    IdentifierPoint transformed = identifier_point(u);
    for (int i = 0; i < 2; ++i) {
        result.caps[i] = nearest_identifier_cap(z[i], beta);
        const auto &tau = result.caps[i].translation;
        auto pairing = identifier_mul(identifier_cap_pair(tau, z[i]), {2, -1});
        need(pairing[0] % 3 == 0 && pairing[1] % 3 == 0, "nonintegral Heisenberg correction");
        transformed = identifier_add(transformed,
                      identifier_add(identifier_mul(identifier_point(beta), {identifier_cap_index(tau), 0}),
                                     {pairing[0] / 3, pairing[1] / 3}));
    }
    auto ratio = identifier_mul(transformed, identifier_conjugate(identifier_point(beta)));
    Wide center = identifier_add(identifier_add(ratio[0], identifier_mul(2, ratio[1])), D);
    Wide shift = identifier_sub(0, identifier_floor(center, identifier_mul(2, D)));
    result.central_shift = identifier_narrow(shift);
    transformed = identifier_add(transformed, identifier_mul(identifier_point(beta), {0, shift}));
    result.u = identifier_narrow(transformed);
    ratio = identifier_mul(transformed, identifier_conjugate(identifier_point(beta)));
    Wide centered = identifier_add(ratio[0], identifier_mul(2, ratio[1]));
    need(centered >= -D && centered < D, "cusp central interval invariant failed");
    Wide index = identifier_sub(ratio[0], identifier_add(result.caps[0].index, result.caps[1].index));
    need(index == original_index, "Heisenberg norm invariant failed");
    result.norm_index = identifier_narrow(index);
    return result;
}
struct IdentifierCusp {
    Point u{}, beta{};
    std::array<Cap, 2> caps{};
};
Wide identifier_cusp_index(const IdentifierCusp &v) {
    auto ratio = identifier_mul(identifier_point(v.u), identifier_conjugate(identifier_point(v.beta)));
    return identifier_sub(ratio[0], identifier_add(identifier_cap_index(v.caps[0]), identifier_cap_index(v.caps[1])));
}
IdentifierPoint identifier_cusp_pair(const IdentifierCusp &a, const IdentifierCusp &b) {
    auto cross = identifier_add(
        identifier_mul(identifier_mul(identifier_conjugate(identifier_point(a.u)), identifier_point(b.beta)), {2, -1}),
        identifier_mul(identifier_mul(identifier_conjugate(identifier_point(a.beta)), identifier_point(b.u)), {1, 1}));
    return identifier_sub(cross, identifier_add(identifier_cap_pair(a.caps[0], b.caps[0]),
                                                identifier_cap_pair(a.caps[1], b.caps[1])));
}
struct CuspTransportStep {
    std::array<Cap, 2> translations{};
    I central_shift = 0;
    bool swap = false;
};
IdentifierCusp apply_cusp_transport(const IdentifierCusp &v, const CuspTransportStep &step) {
    IdentifierCusp result = v;
    auto u = identifier_point(v.u);
    for (int i = 0; i < 2; ++i) {
        const auto &tau = step.translations[i];
        auto pairing = identifier_mul(identifier_cap_pair(tau, v.caps[i]), {2, -1});
        need(pairing[0] % 3 == 0 && pairing[1] % 3 == 0, "nonintegral cusp transport");
        u = identifier_add(u,
            identifier_add(identifier_mul(identifier_point(v.beta), {identifier_cap_index(tau), 0}),
                           {pairing[0] / 3, pairing[1] / 3}));
        for (int j = 0; j < 4; ++j)
            result.caps[i][j] = identifier_narrow(identifier_add(identifier_point(v.caps[i][j]),
                                                identifier_mul(identifier_point(v.beta), identifier_point(tau[j]))));
    }
    u = identifier_add(u, identifier_mul(identifier_point(v.beta), {0, step.central_shift}));
    result.u = identifier_narrow(u);
    if (step.swap) {
        result.u = identifier_narrow(identifier_mul(identifier_point(v.beta), {1, -1}));
        result.beta = identifier_narrow(u);
    }
    need(identifier_cusp_index(result) == identifier_cusp_index(v), "cusp transport norm invariant failed");
    return result;
}
I identifier_central_shift(Point u, Point beta) {
    Wide D = identifier_norm(identifier_point(beta));
    need(D > 0, "cusp beta must be nonzero");
    auto ratio = identifier_mul(identifier_point(u), identifier_conjugate(identifier_point(beta)));
    Wide center = identifier_add(identifier_add(ratio[0], identifier_mul(2, ratio[1])), D);
    return identifier_narrow(identifier_sub(0, identifier_floor(center, identifier_mul(2, D))));
}
struct CuspTransport {
    IdentifierCusp null, image;
    std::vector<CuspTransportStep> steps;
    std::vector<I> heights;
    Point unit{};
    U tie_escapes = 0;
};
CuspTransport transport_identifier_cusp(IdentifierCusp null, IdentifierCusp image) {
    need(identifier_cusp_index(null) == 0, "cusp transport requires a null vector");
    need(null.u != Point{} || null.beta != Point{} || null.caps != std::array<Cap, 2>{},
         "cusp transport rejects the zero vector");
    CuspTransport result;
    auto initial_pair = identifier_cusp_pair(null, image);
    auto append_step = [&](const CuspTransportStep &step) {
        null = apply_cusp_transport(null, step);
        image = apply_cusp_transport(image, step);
        need(identifier_cusp_pair(null, image) == initial_pair, "cusp transport pairing invariant failed");
        result.steps.push_back(step);
    };
    result.heights.push_back(identifier_narrow(identifier_norm(identifier_point(null.beta))));
    while (null.beta != Point{}) {
        Wide D = identifier_norm(identifier_point(null.beta));
        auto normal = normalize_identifier_cusp(null.u, null.beta, null.caps);
        CuspTransportStep step{{normal.caps[0].translation, normal.caps[1].translation},
                               normal.central_shift, false};
        Wide new_height = identifier_norm(identifier_point(normal.u));
        need(new_height <= D, "centered null cusp increased height");
        if (new_height < D) {
            step.swap = true;
            append_step(step);
        } else {
            // Equality means two E8 deep holes and the lower central endpoint.
            // The nearest differences of either cap generate D8 (index two).
            // The Eisenstein complex structure forces an odd central phase
            // among those differences; its half-period changes the height to
            // 3D/4. Enumerate all ties of cap zero, retaining a strict descent.
            append_step(step);
            CapBall ball(null.beta, null.caps[0], 1, 0, 1);
            ball.limit = identifier_mul(identifier_mul(ball.decomposition.scale, D),
                                        identifier_cap_index(null.caps[0]));
            std::optional<CuspTransportStep> escape;
            ball.enumerate([&](const Cap &nearest) {
                CuspTransportStep alternative;
                for (int j = 0; j < 4; ++j)
                    alternative.translations[0][j] = identifier_narrow(identifier_sub(IdentifierPoint{}, identifier_point(nearest[j])));
                auto moved = apply_cusp_transport(null, alternative);
                alternative.central_shift = identifier_central_shift(moved.u, moved.beta);
                moved = apply_cusp_transport(null, alternative);
                if (identifier_norm(identifier_point(moved.u)) < D && !escape)
                    escape = alternative;
            });
            need(escape.has_value(), "null deep-hole tie did not decrease height");
            escape->swap = true;
            append_step(*escape);
            ++result.tie_escapes;
        }
        Wide height = identifier_norm(identifier_point(null.beta));
        need(height < D, "null cusp descent is not strict");
        result.heights.push_back(identifier_narrow(height));
    }
    need(null.caps == std::array<Cap, 2>{}, "terminal null cusp has nonzero horizontal part");
    need(identifier_norm(identifier_point(null.u)) == 1, "cusp transport requires a primitive null vector");
    result.unit = identifier_narrow(identifier_conjugate(identifier_point(null.u)));
    auto scale = [&](IdentifierCusp &v) {
        v.u = identifier_narrow(identifier_mul(identifier_point(v.u), identifier_point(result.unit)));
        v.beta = identifier_narrow(identifier_mul(identifier_point(v.beta), identifier_point(result.unit)));
        for (auto &cap : v.caps)
            for (auto &p : cap)
                p = identifier_narrow(identifier_mul(identifier_point(p), identifier_point(result.unit)));
    };
    scale(null);
    scale(image);
    need(null.u == Point{1, 0} && null.beta == Point{} && null.caps == std::array<Cap, 2>{},
         "cusp transport did not reach the standard null vector");
    need(identifier_cusp_pair(null, image) == initial_pair, "cusp transport final pairing invariant failed");
    result.null = null;
    result.image = image;
    return result;
}
void identifier_cusp_json(std::ostream &out, const IdentifierCusp &v) {
    out << "{\"u\":";
    json_value(out, v.u);
    out << ",\"beta\":";
    json_value(out, v.beta);
    out << ",\"caps\":";
    json_value(out, v.caps);
    out << '}';
}
IdentifierCusp centered_identifier_cusp(const IdentifierCusp &v) {
    auto normal = normalize_identifier_cusp(v.u, v.beta, v.caps);
    return {normal.u, normal.beta, {normal.caps[0].representative, normal.caps[1].representative}};
}
struct IdentifierFiberStats {
    U nearest_calls = 0, nearest_nodes = 0, interval_tests = 0, rejected_pairs = 0;
};
struct IdentifierNearestCap {
    Wide index = 0;
    IdentifierPoint delta{};
};
IdentifierNearestCap prepare_identifier_nearest(const Cap &cap, Point beta, IdentifierFiberStats &stats) {
    auto nearest = nearest_identifier_cap(cap, beta);
    need(stats.nearest_calls < std::numeric_limits<U>::max(), "nearest call count overflow");
    ++stats.nearest_calls;
    need(stats.nearest_nodes <= std::numeric_limits<U>::max() - nearest.visited_nodes,
         "nearest node count overflow");
    stats.nearest_nodes += nearest.visited_nodes;
    auto pairing = identifier_mul(identifier_cap_pair(nearest.translation, cap), {2, -1});
    need(pairing[0] % 3 == 0 && pairing[1] % 3 == 0, "nonintegral cached nearest correction");
    return {nearest.index,
            identifier_add(identifier_mul(identifier_point(beta), {identifier_cap_index(nearest.translation), 0}),
                           {pairing[0] / 3, pairing[1] / 3})};
}
struct IdentifierFiberIntervals {
    std::vector<std::array<I, 2>> intervals;
    I count = 0;
    Wide position = 0, period = 0;
    void append(Wide low, Wide high) {
        if (low > high)
            return;
        I first = identifier_narrow(low), last = identifier_narrow(high);
        need(first >= 0 && last >= first, "invalid retained fiber interval");
        count = identifier_narrow(identifier_add(count, identifier_add(identifier_sub(last, first), 1)));
        if (!intervals.empty() && identifier_add(intervals.back()[1], 1) == first)
            intervals.back()[1] = last;
        else
            intervals.push_back({first, last});
    }
    I next(I fiber, I total) const {
        for (const auto &interval : intervals)
            if (fiber <= interval[1])
                return std::max(fiber, interval[0]);
        return total;
    }
};
IdentifierFiberIntervals identifier_fiber_intervals(Point u, Point beta, Wide m, I subgroup,
                                                   const IdentifierNearestCap &first,
                                                   const IdentifierNearestCap &second) {
    Wide D = identifier_norm(identifier_point(beta));
    I g = std::gcd(beta.x, beta.y);
    need(D > 0 && g > 0 && subgroup > 0 && g % subgroup == 0,
         "invalid central fiber subgroup");
    IdentifierFiberIntervals result;
    auto scalar = identifier_add(identifier_point(u), identifier_add(first.delta, second.delta));
    auto ratio = identifier_mul(scalar, identifier_conjugate(identifier_point(beta)));
    Wide A = identifier_add(m, identifier_add(first.index, second.index));
    need(ratio[0] == A, "cached nearest corrections lost scalar phase");
    Wide step = identifier_mul(2, D / g);
    result.period = identifier_mul(step, subgroup);
    Wide half = result.period / 2;
    result.position = identifier_add(ratio[0], identifier_mul(2, ratio[1]));
    result.position = identifier_sub(result.position,
        identifier_mul(result.period, identifier_floor(identifier_add(result.position, half), result.period)));
    need(result.position >= -half && result.position < half, "fiber position is not centered");
    Wide gap = identifier_sub(identifier_mul(4, identifier_mul(D, D)), identifier_mul(3, identifier_mul(A, A)));
    if (gap <= 0) {
        result.append(0, subgroup - 1);
        return result;
    }
    // Strict rejection means |t| <= floor(sqrt(gap-1)); equality N(u)=D
    // survives. The original f index is retained, including its actual
    // Heisenberg phase offset; a fresh Bezout solution at A is not equivalent.
    Wide radius = identifier_sqrt(identifier_sub(gap, 1));
    if (radius >= half)
        return result;
    Wide wrap = identifier_ceil(identifier_sub(half, result.position), step);
    need(wrap >= 1 && wrap <= subgroup, "fiber wrap is outside one period");
    auto segment = [&](Wide low, Wide high, Wide base) {
        if (low > high)
            return;
        Wide negative_end = identifier_floor(identifier_sub(identifier_sub(0, identifier_add(radius, 1)), base), step);
        Wide positive_start = identifier_ceil(identifier_sub(identifier_add(radius, 1), base), step);
        result.append(low, std::min(high, negative_end));
        result.append(std::max(low, positive_start), high);
    };
    segment(0, wrap - 1, result.position);
    segment(wrap, subgroup - 1, identifier_sub(result.position, result.period));
    need(result.intervals.size() <= 2 && result.count <= subgroup, "fiber interval decomposition failed");
    return result;
}
bool identifier_old_height_rejection(const IdentifierCusp &raw, I subgroup, IdentifierFiberStats &stats) {
    auto normal = normalize_identifier_cusp(raw.u, raw.beta, raw.caps);
    need(stats.nearest_calls <= std::numeric_limits<U>::max() - 2, "nearest call count overflow");
    stats.nearest_calls += 2;
    for (const auto &cap : normal.caps) {
        need(stats.nearest_nodes <= std::numeric_limits<U>::max() - cap.visited_nodes,
             "nearest node count overflow");
        stats.nearest_nodes += cap.visited_nodes;
    }
    Wide D = identifier_norm(identifier_point(raw.beta));
    I g = std::gcd(raw.beta.x, raw.beta.y);
    need(g > 0 && subgroup > 0 && g % subgroup == 0, "invalid legacy central fiber subgroup");
    Point kappa{-raw.beta.y / g, (raw.beta.x + raw.beta.y) / g};
    auto ratio = identifier_mul(identifier_point(normal.u), identifier_conjugate(identifier_point(raw.beta)));
    Wide position = identifier_add(ratio[0], identifier_mul(2, ratio[1]));
    Wide central_step = identifier_mul(identifier_mul(2, D / g), subgroup);
    Wide shift = identifier_sub(0, identifier_floor(identifier_add(identifier_mul(2, position), central_step),
                                                   identifier_mul(2, central_step)));
    auto scalar = identifier_add(identifier_point(normal.u),
                                identifier_mul(identifier_point(kappa), {identifier_mul(shift, subgroup), 0}));
    return identifier_norm(scalar) < D;
}
std::array<I, 3> bezout(I a, I b);
I positive_mod(I a, I b);
bool identifier_divide(IdentifierPoint numerator, IdentifierPoint denominator, IdentifierPoint &quotient) {
    Wide D = identifier_norm(denominator);
    need(D > 0, "zero identifier divisor");
    auto product = identifier_mul(numerator, identifier_conjugate(denominator));
    if (D == static_cast<I>(D) && product[0] == static_cast<I>(product[0]) && product[1] == static_cast<I>(product[1])) {
        const I x = static_cast<I>(product[0]), y = static_cast<I>(product[1]), d = static_cast<I>(D);
        if (x % d || y % d) return false;
        quotient = {x / d, y / d};
        return true;
    }
    if (product[0] % D || product[1] % D)
        return false;
    quotient = {product[0] / D, product[1] / D};
    return true;
}
IdentifierPoint identifier_gcd(IdentifierPoint a, IdentifierPoint b) {
    while (b != IdentifierPoint{}) {
        Wide D = identifier_norm(b);
        if (D == 1)
            return b;
        auto product = identifier_mul(a, identifier_conjugate(b));
        IdentifierPoint quotient;
        for (int i = 0; i < 2; ++i) {
            quotient[i] = identifier_floor(product[i], D);
            Wide remainder = identifier_sub(product[i], identifier_mul(quotient[i], D));
            if (remainder > D / 2)
                quotient[i] = identifier_add(quotient[i], 1);
        }
        auto remainder = identifier_sub(a, identifier_mul(quotient, b));
        need(identifier_norm(remainder) < D, "Eisenstein gcd did not decrease");
        a = b;
        b = remainder;
    }
    return a;
}
bool primitive_identifier_cusp(const IdentifierCusp &v) {
    auto common = identifier_point(v.beta);
    auto absorb = [&](Point p) {
        if (identifier_norm(common) != 1)
            common = identifier_gcd(identifier_point(p), common);
    };
    absorb(v.u);
    for (const auto &cap : v.caps)
        for (Point p : cap)
            absorb(p);
    return identifier_norm(common) == 1;
}
struct IdentifierEnumerationStats {
    U gamma_count = 0, horizontal_pairs = 0, scalar_candidates = 0, count = 0, visited_nodes = 0;
    bool complete = true;
};
void identifier_count_add(U &a, U b) {
    need(a <= std::numeric_limits<U>::max() - b, "identifier count overflow");
    a += b;
}
Cap identifier_scaled_cap(const Cap &z, Point gamma) {
    Cap result;
    for (int i = 0; i < 4; ++i)
        result[i] = identifier_narrow(identifier_mul(identifier_point(gamma), identifier_point(z[i])));
    return result;
}
Wide identifier_shifted_index(const Cap &y, Point beta, const Cap &center) {
    Cap difference;
    for (int i = 0; i < 4; ++i)
        difference[i] = identifier_narrow(identifier_sub(
            identifier_mul(identifier_point(beta), identifier_point(y[i])), identifier_point(center[i])));
    return identifier_cap_index(difference);
}
template <class Callback>
bool identifier_horizontal_ball(const IdentifierCusp &v, Point gamma, Wide radius,
                                IdentifierEnumerationStats &stats, Callback &&callback) {
    need(radius >= 0, "negative horizontal ball radius");
    auto center_a = identifier_scaled_cap(v.caps[0], gamma);
    auto center_b = identifier_scaled_cap(v.caps[1], gamma);
    Wide D = identifier_norm(identifier_point(v.beta));
    CapBall a(v.beta, center_a, 1, 0, 1);
    Wide factor = identifier_mul(a.decomposition.scale, D);
    a.limit = identifier_mul(factor, radius);
    bool complete = a.enumerate([&](const Cap &first) -> bool {
        Wide remaining = identifier_sub(radius, identifier_shifted_index(first, v.beta, center_a));
        need(remaining >= 0, "first cap exceeds horizontal radius");
        CapBall b(v.beta, center_b, 1, 0, 1);
        b.limit = identifier_mul(factor, remaining);
        bool finished = b.enumerate([&](const Cap &second) -> bool {
            identifier_count_add(stats.horizontal_pairs, 1);
            return callback(std::array<Cap, 2>{first, second});
        });
        identifier_count_add(stats.visited_nodes, b.visited_nodes);
        return finished;
    });
    identifier_count_add(stats.visited_nodes, a.visited_nodes);
    return complete;
}
template <class Callback>
IdentifierEnumerationStats identifier_competitors(const IdentifierCusp &v, Callback &&callback) {
    Wide m = identifier_cusp_index(v), D = identifier_norm(identifier_point(v.beta));
    need(m > 0 && D > 0 && D < identifier_mul(40, m),
         "competing cusps require positive norm and height below 40m");
    Wide four_D_squared = identifier_mul(4, identifier_mul(D, D));
    Wide bound = four_D_squared / identifier_mul(3, identifier_mul(m, m));
    IdentifierEnumerationStats stats;
    IdentifierCusp original;
    original.u = {1, 0};
    stats.count = 1;
    if (!callback(original, D)) {
        stats.complete = false;
        return stats;
    }
    // a>0,b>=0 is one representative per Eisenstein unit class; the
    // half-open wedge keeps both conjugate ideals unless they coincide.
    for (I a = 1; identifier_mul(a, a) <= bound; ++a)
        for (I b = 0;; ++b) {
            Point gamma{a, b};
            Wide n = identifier_norm(identifier_point(gamma));
            if (n > bound)
                break;
            Wide radius = identifier_sub(identifier_sqrt(identifier_mul(four_D_squared, n) / 3),
                                         identifier_mul(m, n));
            if (radius < 0)
                continue;
            identifier_count_add(stats.gamma_count, 1);
            auto [g, s, t] = bezout(a + b, b);
            Point kappa{-b / g, (a + b) / g};
            bool complete = identifier_horizontal_ball(v, gamma, radius, stats,
                [&](const std::array<Cap, 2> &caps) -> bool {
                    Wide Q = identifier_add(identifier_cap_index(caps[0]), identifier_cap_index(caps[1]));
                    if (Q % g)
                        return true;
                    Point A0{identifier_narrow(identifier_mul(s, Q / g)),
                             identifier_narrow(identifier_mul(t, Q / g))};
                    IdentifierCusp w{A0, gamma, caps};
                    auto H0 = identifier_cusp_pair(w, v);
                    auto L = identifier_mul(identifier_mul(identifier_conjugate(identifier_point(kappa)),
                                                           identifier_point(v.beta)), {2, -1});
                    Wide quadratic = identifier_norm(L);
                    Wide linear = identifier_add(identifier_add(identifier_mul(2, identifier_mul(H0[0], L[0])),
                                                               identifier_mul(H0[0], L[1])),
                                                 identifier_add(identifier_mul(H0[1], L[0]),
                                                                identifier_mul(2, identifier_mul(H0[1], L[1]))));
                    Wide constant = identifier_sub(identifier_norm(H0), identifier_mul(3, D));
                    Wide discriminant = identifier_sub(identifier_mul(linear, linear),
                                                        identifier_mul(4, identifier_mul(quadratic, constant)));
                    if (discriminant < 0)
                        return true;
                    Wide root = identifier_sqrt(discriminant), denominator = identifier_mul(2, quadratic);
                    Wide low = identifier_ceil(identifier_sub(identifier_sub(0, root), linear), denominator);
                    Wide high = identifier_floor(identifier_sub(root, linear), denominator);
                    for (Wide j = low; j <= high;) {
                        identifier_count_add(stats.scalar_candidates, 1);
                        w.u = identifier_narrow(identifier_add(identifier_point(A0),
                                               identifier_mul(identifier_point(kappa), {j, 0})));
                        auto H = identifier_cusp_pair(w, v);
                        Wide height = identifier_norm(H);
                        need(identifier_cusp_index(w) == 0 && height % 3 == 0 && height <= identifier_mul(3, D),
                             "competing cusp quadratic invariant failed");
                        if (primitive_identifier_cusp(w)) {
                            identifier_count_add(stats.count, 1);
                            if (!callback(w, height / 3))
                                return false;
                        }
                        if (j == high)
                            break;
                        j = identifier_add(j, 1);
                    }
                    return true;
                });
            if (!complete) {
                stats.complete = false;
                return stats;
            }
        }
    return stats;
}
Cap cap_reflect(Cap a, int j);
const std::vector<Cap> &identifier_cap_roots() {
    static const auto roots = [] {
        Cap first{};
        first[0] = {1, 0};
        std::set<Cap> seen{first};
        std::vector<Cap> result{first};
        for (std::size_t i = 0; i < result.size(); ++i)
            for (int j = 0; j < 4; ++j) {
                auto root = cap_reflect(result[i], j);
                if (seen.insert(root).second)
                    result.push_back(root);
            }
        need(result.size() == 240, "wrong identifier cap root count");
        std::set<Cap> represented;
        std::vector<Cap> lines;
        for (const auto &root : result) {
            if (represented.contains(root))
                continue;
            lines.push_back(root);
            for (int k = 0; k < 6; ++k) {
                Cap other;
                for (int j = 0; j < 4; ++j)
                    other[j] = rotate(root[j], k);
                represented.insert(other);
            }
        }
        need(lines.size() == 40 && represented.size() == 240, "wrong identifier cap root-line count");
        // Keep the first BFS representative of each unit line. This preserves
        // the exact first mirror witness of the original 240-root traversal.
        return lines;
    }();
    return roots;
}
struct IdentifierMirror {
    std::optional<IdentifierCusp> root;
    IdentifierEnumerationStats stats;
    U horizontal_root_checks = 0;
};
IdentifierMirror identifier_mirror(const IdentifierCusp &v, bool cap_mirrors_known_absent = false) {
    Wide m = identifier_cusp_index(v), D = identifier_norm(identifier_point(v.beta));
    need(m > 0 && D > 0 && D < identifier_mul(40, m),
         "mirror test requires positive norm and height below 40m");
    IdentifierMirror result;
    auto test_horizontal = [&](Point gamma, const std::array<Cap, 2> &caps) {
        auto h = identifier_add(identifier_cap_pair(caps[0], v.caps[0]),
                                identifier_cap_pair(caps[1], v.caps[1]));
        auto correction = identifier_mul(identifier_mul(identifier_conjugate(identifier_point(gamma)),
                                                        identifier_point(v.u)), {1, 1});
        IdentifierPoint conjugate_A;
        if (!identifier_divide(identifier_sub(h, correction), identifier_mul(identifier_point(v.beta), {2, -1}), conjugate_A))
            return false;
        IdentifierCusp root{identifier_narrow(identifier_conjugate(conjugate_A)), gamma, caps};
        if (identifier_cusp_index(root) != -1 || identifier_cusp_pair(root, v) != IdentifierPoint{})
            return false;
        result.root = root;
        return true;
    };
    if (!cap_mirrors_known_absent)
        for (int cap = 0; cap < 2; ++cap)
            for (const auto &root : identifier_cap_roots()) {
                std::array<Cap, 2> caps{};
                caps[cap] = root;
                identifier_count_add(result.horizontal_root_checks, 1);
                if (test_horizontal({}, caps))
                    return result;
            }
    Wide bound = D / m;
    for (I a = 1; identifier_mul(a, a) <= bound; ++a)
        for (I b = 0;; ++b) {
            Point gamma{a, b};
            Wide n = identifier_norm(identifier_point(gamma));
            if (n > bound)
                break;
            identifier_count_add(result.stats.gamma_count, 1);
            Wide radius = identifier_sub(D, identifier_mul(m, n));
            identifier_horizontal_ball(v, gamma, radius, result.stats,
                [&](const std::array<Cap, 2> &caps) {
                    return !test_horizontal(gamma, caps);
                });
            if (result.root)
                return result;
        }
    return result;
}
struct IdentifierRing {
    Point beta;
    I D, h11, h12, h22, bezout_s, bezout_t;
    explicit IdentifierRing(Point value) : beta(value) {
        D = identifier_narrow(identifier_norm(identifier_point(beta)));
        need(D > 0, "identifier residue ring requires nonzero beta");
        auto [g, s, t] = bezout(beta.y, identifier_narrow(identifier_add(beta.x, beta.y)));
        h22 = g;
        h11 = D / g;
        Wide residue = identifier_sub(identifier_mul(beta.x, s), identifier_mul(beta.y, t));
        h12 = identifier_narrow(identifier_sub(residue, identifier_mul(identifier_floor(residue, h11), h11)));
        auto scalar_bezout = bezout(identifier_narrow(identifier_add(beta.x, beta.y)), beta.y);
        need(scalar_bezout[0] == g, "inconsistent cusp residue gcd");
        bezout_s = scalar_bezout[1];
        bezout_t = scalar_bezout[2];
    }
    Point reduce(IdentifierPoint p) const {
        Wide quotient = identifier_floor(p[1], h22);
        Wide b = identifier_sub(p[1], identifier_mul(quotient, h22));
        Wide a = identifier_sub(p[0], identifier_mul(quotient, h12));
        a = identifier_sub(a, identifier_mul(identifier_floor(a, h11), h11));
        return {identifier_narrow(a), identifier_narrow(b)};
    }
    Cap reduce(const Cap &cap) const {
        Cap result;
        for (int i = 0; i < 4; ++i)
            result[i] = reduce(identifier_point(cap[i]));
        return result;
    }
    I phase(IdentifierPoint value) const {
        Wide result = identifier_sub(identifier_mul(bezout_s, value[1]), identifier_mul(bezout_t, value[0]));
        return identifier_narrow(identifier_sub(result, identifier_mul(identifier_floor(result, h22), h22)));
    }
};
struct IdentifierCapHash {
    std::size_t operator()(const Cap &cap) const {
        U hash = 0x9e3779b97f4a7c15ULL;
        for (Point p : cap)
            for (I value : {p.x, p.y}) {
                U x = static_cast<U>(value) + 0x9e3779b97f4a7c15ULL;
                x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
                x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
                hash ^= (x ^ (x >> 31)) + (hash << 6) + (hash >> 2);
            }
        return static_cast<std::size_t>(hash);
    }
};
struct IdentifierCapStep {
    Cap cap;
    I phase;
};
IdentifierCapStep identifier_cap_step(const IdentifierRing &ring, const Cap &current, int j) {
    auto reflected = identifier_mul(identifier_point(current[j]), {-1, 1});
    if (j)
        reflected = identifier_add(reflected, identifier_point(current[j - 1]));
    if (j < 3)
        reflected = identifier_add(reflected, identifier_mul(identifier_point(current[j + 1]), {1, -1}));
    Cap next = current;
    next[j] = ring.reduce(reflected);
    I edge_phase = 0;
    if (ring.h22 > 1) {
        IdentifierPoint tau;
        need(identifier_divide(identifier_sub(identifier_point(next[j]), reflected), identifier_point(ring.beta), tau),
             "nonintegral cap orbit translation");
        auto row = identifier_mul({3, 0}, reflected);
        if (j)
            row = identifier_sub(row, identifier_mul({1, 1}, identifier_point(current[j - 1])));
        if (j < 3)
            row = identifier_sub(row, identifier_mul({2, -1}, identifier_point(current[j + 1])));
        auto correction = identifier_mul(identifier_mul(identifier_conjugate(tau), row), {2, -1});
        need(correction[0] % 3 == 0 && correction[1] % 3 == 0,
             "nonintegral cap orbit scalar correction");
        auto delta = identifier_add(identifier_mul(identifier_point(ring.beta), {identifier_norm(tau), 0}),
                                    {correction[0] / 3, correction[1] / 3});
        edge_phase = ring.phase(delta);
    }
    return {next, edge_phase};
}
struct IdentifierCapOrbit {
    struct State { Cap cap; I phase; };
    Point beta;
    std::unordered_map<Cap, std::size_t, IdentifierCapHash> index;
    std::vector<State> states;
    std::size_t minimum = 0;
    I subgroup = 1;
    IdentifierCapOrbit(const IdentifierRing &ring, const Cap &first) : beta(ring.beta), subgroup(ring.h22) {
        index.reserve(155520);
        states.reserve(155520);
        index.emplace(first, 0);
        states.push_back({first, 0});
        for (std::size_t i = 0; i < states.size(); ++i) {
            // A generator changes just one root-chain coefficient. Work in
            // HNF residues here; nearest-point calls are needed only for the
            // eventual minimum, not at every one of the 622080 group edges.
            Cap current = states[i].cap;
            I path_phase = states[i].phase;
            for (int j = 0; j < 4; ++j) {
                auto step = identifier_cap_step(ring, current, j);
                Cap next = step.cap;
                I phase = positive_mod(identifier_narrow(identifier_add(path_phase, step.phase)), ring.h22);
                auto [position, inserted] = index.emplace(next, states.size());
                if (inserted) {
                    need(states.size() < 155520, "cap orbit exceeds its proved group order");
                    states.push_back({next, phase});
                    if (next < states[minimum].cap)
                        minimum = states.size() - 1;
                } else
                    subgroup = std::gcd(subgroup, positive_mod(phase - states[position->second].phase, ring.h22));
            }
        }
        need(155520 % states.size() == 0, "cap orbit size does not divide its group order");
    }
};
struct IdentifierOrbitImage {
    Cap cap;
    I phase = 0, subgroup = 1;
};
struct IdentifierCapTransversal {
    Cap row;
    std::vector<unsigned char> word;
};
using IdentifierCapChain = std::array<std::vector<IdentifierCapTransversal>, 4>;
const IdentifierCapChain &identifier_cap_chain() {
    static const IdentifierCapChain chain = [] {
        IdentifierCapChain result;
        constexpr std::array<std::size_t, 4> orders{240, 27, 8, 3};
        for (int i = 0; i < 4; ++i) {
            Cap first{};
            first[i] = {1, 0};
            auto &level = result[i];
            level.push_back({first, {}});
            std::unordered_map<Cap, std::size_t, IdentifierCapHash> index;
            index.emplace(first, 0);
            for (std::size_t k = 0; k < level.size(); ++k)
                for (int j = i; j < 4; ++j) {
                    Cap row = level[k].row;
                    auto value = identifier_point(row[j]);
                    row[j] = identifier_narrow(identifier_mul(value, {-1, 1}));
                    if (j)
                        row[j - 1] = identifier_narrow(identifier_add(identifier_point(row[j - 1]), value));
                    if (j < 3)
                        row[j + 1] = identifier_narrow(identifier_add(identifier_point(row[j + 1]), identifier_mul(value, {1, -1})));
                    if (index.emplace(row, level.size()).second) {
                        need(level.size() < orders[i], "cap chain row orbit exceeds its certified size");
                        auto word = level[k].word;
                        word.push_back(static_cast<unsigned char>(j));
                        level.push_back({row, std::move(word)});
                    }
                }
            need(level.size() == orders[i], "cap chain row orbit has incorrect size");
        }
        return result;
    }();
    return chain;
}
struct IdentifierChainStats {
    U nodes = 0, row_tests = 0, generator_steps = 0, singleton_shortcuts = 0, duplicate_children = 0;
};
IdentifierOrbitImage identifier_chain_image(const IdentifierRing &ring, const Cap &cap,
                                            IdentifierChainStats &stats, int level = 0) {
    identifier_count_add(stats.nodes, 1);
    if (level == 4)
        return {cap, 0, ring.h22};
    I fixed_subgroup = ring.h22;
    bool fixed = true;
    for (int j = level; j < 4; ++j) {
        auto step = identifier_cap_step(ring, cap, j);
        identifier_count_add(stats.generator_steps, 1);
        if (step.cap != cap) {
            fixed = false;
            break;
        }
        fixed_subgroup = std::gcd(fixed_subgroup, step.phase);
    }
    if (fixed) {
        identifier_count_add(stats.singleton_shortcuts, 1);
        return {cap, 0, fixed_subgroup};
    }
    const auto &transversals = identifier_cap_chain()[level];
    std::vector<std::size_t> survivors;
    std::optional<Point> minimum;
    for (std::size_t k = 0; k < transversals.size(); ++k) {
        IdentifierPoint value{};
        for (int j = 0; j < 4; ++j)
            value = identifier_add(value, identifier_mul(identifier_point(transversals[k].row[j]), identifier_point(cap[j])));
        Point coordinate = ring.reduce(value);
        identifier_count_add(stats.row_tests, 1);
        if (!minimum || coordinate < *minimum) {
            minimum = coordinate;
            survivors.clear();
        }
        if (coordinate == *minimum)
            survivors.push_back(k);
    }
    struct Child { Cap cap; I phase, incoming_subgroup; };
    std::vector<Child> children;
    std::unordered_map<Cap, std::size_t, IdentifierCapHash> child_index;
    for (std::size_t k : survivors) {
        Cap next = cap;
        I phase = 0;
        // Row words act on the right; cap columns therefore use reverse order.
        const auto &word = transversals[k].word;
        for (auto j = word.rbegin(); j != word.rend(); ++j) {
            auto step = identifier_cap_step(ring, next, *j);
            next = step.cap;
            phase = positive_mod(identifier_narrow(identifier_add(phase, step.phase)), ring.h22);
            identifier_count_add(stats.generator_steps, 1);
        }
        need(next[level] == *minimum, "cap chain row and word actions disagree");
        for (int j = 0; j < level; ++j)
            need(next[j] == cap[j], "cap chain changed a fixed prefix coordinate");
        auto [position, inserted] = child_index.emplace(next, children.size());
        if (inserted)
            children.push_back({next, phase, ring.h22});
        else {
            auto &previous = children[position->second];
            previous.incoming_subgroup = std::gcd(previous.incoming_subgroup,
                                                  positive_mod(phase - previous.phase, ring.h22));
            identifier_count_add(stats.duplicate_children, 1);
        }
    }
    std::optional<IdentifierOrbitImage> best;
    for (const auto &child : children) {
        auto image = identifier_chain_image(ring, child.cap, stats, level + 1);
        image.phase = positive_mod(identifier_narrow(identifier_add(child.phase, image.phase)), ring.h22);
        image.subgroup = std::gcd(child.incoming_subgroup, image.subgroup);
        if (!best || image.cap < best->cap)
            best = image;
        else if (image.cap == best->cap)
            best->subgroup = std::gcd(std::gcd(best->subgroup, image.subgroup),
                                       positive_mod(image.phase - best->phase, ring.h22));
    }
    need(best.has_value(), "empty cap chain minimum");
    // Only the complete set of minimum transports is a stabilizer-image coset.
    // Intermediate duplicate-child phase sets are used only through differences
    // of actual paths, never asserted to be closed by themselves.
    best->phase = positive_mod(best->phase, best->subgroup);
    return *best;
}
struct IdentifierParabolicCache {
    // A fixed four-orbit LRU keeps memory independent of the requested size.
    // It stores finite cap group orbits, never fullerene or period graph keys.
    std::vector<IdentifierCapOrbit> orbits;
    U built_states = 0, cache_hits = 0;
    bool use_chain = true;
    IdentifierChainStats chain_stats;
    IdentifierOrbitImage image(const IdentifierRing &ring, const Cap &cap) {
        if (use_chain) {
            auto result = identifier_chain_image(ring, cap, chain_stats);
            built_states = chain_stats.nodes;
            return result;
        }
        std::size_t chosen = orbits.size(), state = 0;
        for (std::size_t i = 0; i < orbits.size(); ++i)
            if (orbits[i].beta == ring.beta) {
                auto found = orbits[i].index.find(cap);
                if (found != orbits[i].index.end()) {
                    chosen = i;
                    state = found->second;
                    identifier_count_add(cache_hits, 1);
                    break;
                }
            }
        if (chosen == orbits.size()) {
            if (orbits.size() == 4)
                orbits.erase(orbits.begin());
            orbits.emplace_back(ring, cap);
            chosen = orbits.size() - 1;
            identifier_count_add(built_states, orbits.back().states.size());
        }
        const auto &orbit = orbits[chosen];
        IdentifierOrbitImage result{orbit.states[orbit.minimum].cap,
                                   positive_mod(orbit.states[orbit.minimum].phase - orbit.states[state].phase, ring.h22),
                                   orbit.subgroup};
        if (chosen + 1 != orbits.size())
            std::rotate(orbits.begin() + chosen, orbits.begin() + chosen + 1, orbits.end());
        return result;
    }
};
void identifier_chain_stats_json(std::ostream &out, const IdentifierParabolicCache &cache) {
    out << ",\"group_backend\":\"" << (cache.use_chain ? "chain" : "bfs")
        << "\",\"chain_row_tests\":" << cache.chain_stats.row_tests
        << ",\"chain_generator_steps\":" << cache.chain_stats.generator_steps
        << ",\"chain_singleton_shortcuts\":" << cache.chain_stats.singleton_shortcuts
        << ",\"chain_duplicate_children\":" << cache.chain_stats.duplicate_children;
}
IdentifierCusp identifier_scale_unit(IdentifierCusp v, Point unit) {
    auto scale = [&](Point p) { return identifier_narrow(identifier_mul(identifier_point(p), identifier_point(unit))); };
    v.u = scale(v.u);
    v.beta = scale(v.beta);
    for (auto &cap : v.caps)
        for (auto &p : cap)
            p = scale(p);
    return v;
}
IdentifierCusp identifier_height_unit(IdentifierCusp v) {
    need(v.beta != Point{}, "parabolic normal form requires nonzero beta");
    for (int i = 0; i < 6; ++i) {
        if (v.beta.x > 0 && v.beta.y >= 0)
            return v;
        v = identifier_scale_unit(v, {0, 1});
    }
    throw std::runtime_error("no unit representative of cusp height");
}
IdentifierCusp identifier_antiunitary(const IdentifierCusp &v) {
    IdentifierCusp result;
    result.u = identifier_narrow(identifier_conjugate(identifier_point(v.u)));
    result.beta = identifier_narrow(identifier_mul({0, 1}, identifier_conjugate(identifier_point(v.beta))));
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 4; ++j)
            result.caps[i][j] = identifier_narrow(identifier_conjugate(identifier_point(v.caps[i][3 - j])));
    need(identifier_cusp_index(result) == identifier_cusp_index(v), "antiunitary norm invariant failed");
    return result;
}
struct IdentifierParabolic {
    IdentifierCusp record;
    std::array<Cap, 2> residue_key;
    I central_subgroup = 1;
    bool antiunitary = false;
};
struct PreparedIdentifierCap {
    Cap residue, nearest;
    IdentifierPoint nearest_delta{};
    Wide residue_index = 0;
    I phase = 0, subgroup = 1;
};
PreparedIdentifierCap prepare_identifier_cap(const Cap &cap, const IdentifierRing &ring,
                                             IdentifierParabolicCache &cache) {
    Cap residue = ring.reduce(cap), tau;
    for (int j = 0; j < 4; ++j) {
        IdentifierPoint quotient;
        need(identifier_divide(identifier_sub(identifier_point(residue[j]), identifier_point(cap[j])),
                               identifier_point(ring.beta), quotient), "nonintegral residue translation");
        tau[j] = identifier_narrow(quotient);
    }
    auto delta_for = [&](const Cap &translation, const Cap &original) {
        auto pairing = identifier_mul(identifier_cap_pair(translation, original), {2, -1});
        need(pairing[0] % 3 == 0 && pairing[1] % 3 == 0, "nonintegral prepared cap correction");
        return identifier_add(identifier_mul(identifier_point(ring.beta), {identifier_cap_index(translation), 0}),
                              {pairing[0] / 3, pairing[1] / 3});
    };
    auto image = cache.image(ring, residue);
    auto nearest = nearest_identifier_cap(image.cap, ring.beta);
    return {image.cap, nearest.representative, delta_for(nearest.translation, image.cap),
            identifier_cap_index(image.cap),
            positive_mod(ring.phase(delta_for(tau, cap)) + image.phase, ring.h22), image.subgroup};
}
IdentifierParabolic combine_identifier_caps(Point u, const IdentifierRing &ring, Wide m,
                                             const PreparedIdentifierCap &first,
                                             const PreparedIdentifierCap &second) {
    I subgroup = std::gcd(first.subgroup, second.subgroup);
    I phase = positive_mod(identifier_narrow(identifier_add(ring.phase(identifier_point(u)),
                                                            identifier_add(first.phase, second.phase))), subgroup);
    Wide R = identifier_add(m, identifier_add(first.residue_index, second.residue_index));
    need(R % ring.h22 == 0, "cap orbit lost scalar solubility");
    Point kappa{-ring.beta.y / ring.h22, (ring.beta.x + ring.beta.y) / ring.h22};
    auto scalar = identifier_add({identifier_mul(ring.bezout_s, R / ring.h22),
                                 identifier_mul(ring.bezout_t, R / ring.h22)},
                                identifier_mul(identifier_point(kappa), {phase, 0}));
    scalar = identifier_add(scalar, identifier_add(first.nearest_delta, second.nearest_delta));
    auto ratio = identifier_mul(scalar, identifier_conjugate(identifier_point(ring.beta)));
    Wide position = identifier_add(ratio[0], identifier_mul(2, ratio[1]));
    Wide step = identifier_mul(identifier_mul(2, ring.D / ring.h22), subgroup);
    Wide low = identifier_ceil(identifier_sub(-Wide(ring.D), position), step);
    Wide high = identifier_sub(identifier_ceil(identifier_sub(ring.D, position), step), 1);
    auto direction = identifier_mul(identifier_point(kappa), {subgroup, 0});
    Wide choice = (direction[0] > 0 || (direction[0] == 0 && direction[1] > 0)) ? low : high;
    need(low <= high, "empty central phase interval");
    scalar = identifier_add(scalar, identifier_mul(direction, {choice, 0}));
    IdentifierCusp v{identifier_narrow(scalar), ring.beta, {first.nearest, second.nearest}};
    std::array<Cap, 2> residue_key{first.residue, second.residue};
    if (residue_key[1] < residue_key[0]) {
        std::swap(residue_key[0], residue_key[1]);
        std::swap(v.caps[0], v.caps[1]);
    }
    need(identifier_cusp_index(v) == m, "parabolic normal form changed the norm");
    return {v, residue_key, subgroup, false};
}
bool identifier_parabolic_less(const IdentifierParabolic &a, const IdentifierParabolic &b) {
    return std::tuple(a.record.beta, a.residue_key, a.record.u) <
           std::tuple(b.record.beta, b.residue_key, b.record.u);
}
bool identifier_parabolic_equal(const IdentifierParabolic &a, const IdentifierParabolic &b) {
    return !identifier_parabolic_less(a, b) && !identifier_parabolic_less(b, a);
}
IdentifierParabolic identifier_unitary_parabolic(IdentifierCusp v, IdentifierParabolicCache &cache) {
    v = identifier_height_unit(v);
    Wide m = identifier_cusp_index(v);
    IdentifierRing ring(v.beta);
    auto first = prepare_identifier_cap(v.caps[0], ring, cache);
    auto second = prepare_identifier_cap(v.caps[1], ring, cache);
    return combine_identifier_caps(v.u, ring, m, first, second);
}
IdentifierParabolic identifier_parabolic(const IdentifierCusp &v, IdentifierParabolicCache &cache) {
    auto unitary = identifier_height_unit(v);
    auto conjugate = identifier_height_unit(identifier_antiunitary(unitary));
    if (conjugate.beta < unitary.beta) {
        auto second = identifier_unitary_parabolic(conjugate, cache);
        second.antiunitary = true;
        return second;
    }
    auto first = identifier_unitary_parabolic(unitary, cache);
    // A larger beta cannot win the declared beta-first comparison order.
    if (conjugate.beta > first.record.beta)
        return first;
    auto second = identifier_unitary_parabolic(conjugate, cache);
    if (identifier_parabolic_less(second, first)) {
        second.antiunitary = true;
        return second;
    }
    return first;
}
struct IdentifierCanonical {
    IdentifierParabolic normal;
    std::optional<IdentifierCusp> root;
    std::vector<I> heights;
    U comparison_cusps = 0, transported_cusps = 0, horizontal_pairs = 0, visited_nodes = 0;
};
Cap cap_reflect(Cap a, int j) {
    a[j] = rotate(a[j], 2) + (j ? a[j - 1] : Point{}) + (j < 3 ? rotate(a[j + 1], -1) : Point{});
    return a;
}
const std::vector<Cap> &cap_root_functionals() {
    static const auto rows = [] {
        Cap first{};
        first[0] = {1, 0};
        std::set<Cap> seen{first};
        std::vector<Cap> queue{first};
        for (std::size_t i = 0; i < queue.size(); ++i)
            for (int j = 0; j < 4; ++j) {
                Cap next_cap = cap_reflect(queue[i], j);
                if (seen.insert(next_cap).second)
                    queue.push_back(next_cap);
            }
        need(seen.size() == 240, "wrong cap root orbit");
        std::set<Cap> lines;
        for (const auto &root : queue) {
            Cap best = root;
            for (int k = 1; k < 6; ++k) {
                Cap other;
                for (int j = 0; j < 4; ++j)
                    other[j] = rotate(root[j], k);
                best = std::min(best, other);
            }
            lines.insert(best);
        }
        need(lines.size() == 40, "wrong cap root-line count");
        std::vector<Cap> result;
        for (const auto &root : lines) {
            Cap row;
            for (int j = 0; j < 4; ++j) {
                Cap unit{};
                unit[j] = {1, 0};
                row[j] = divide_theta(cap_pair(root, unit));
            }
            result.push_back(row);
        }
        return result;
    }();
    return rows;
}
Cap short_cap(Cap z, Point d) {
    static const auto roots = [] {
        Cap first{};
        first[0] = {1, 0};
        std::set<Cap> seen{first};
        std::vector<Cap> result{first};
        for (std::size_t i = 0; i < result.size(); ++i)
            for (int j = 0; j < 4; ++j) {
                Cap r = cap_reflect(result[i], j);
                if (seen.insert(r).second)
                    result.push_back(r);
            }
        need(result.size() == 240, "wrong cap Voronoi root count");
        return result;
    }();
    while (true) {
        I best_norm = cap_index(z);
        Cap best = z;
        for (const auto &r : roots) {
            Cap other;
            for (int j = 0; j < 4; ++j)
                other[j] = z[j] + multiply(d, r[j]);
            I value = cap_index(other);
            if (value < best_norm) {
                best_norm = value;
                best = other;
            }
        }
        if (best == z)
            return z;
        z = best;
    }
}
std::array<I, 3> bezout(I a, I b) {
    I r = a, next = b, s = 1, ns = 0, t = 0, nt = 1;
    while (next) {
        I q = r / next;
        I old = r;
        r = next;
        next = old - q * next;
        old = s;
        s = ns;
        ns = old - q * ns;
        old = t;
        t = nt;
        nt = old - q * nt;
    }
    if (r < 0)
        return {-r, -s, -t};
    return {r, s, t};
}
I positive_mod(I a, I b) {
    I r = a % b;
    return r < 0 ? r + b : r;
}
struct CuspRing {
    Point divisor;
    I size, h11, h12, h22;
    std::vector<Point> elements;
    std::vector<unsigned> addition, omega, inverse_zeta;
    explicit CuspRing(Point d) : divisor(d), size(norm(d)) {
        auto [g, s, t] = bezout(d.y, d.x + d.y);
        h22 = g;
        h11 = size / g;
        h12 = positive_mod(d.x * s - d.y * t, h11);
        need(size > 0 && size <= std::numeric_limits<unsigned>::max(), "cusp ring too large");
        for (I b = 0; b < h22; ++b)
            for (I a = 0; a < h11; ++a)
                elements.push_back({a, b});
        addition.resize(static_cast<std::size_t>(size) * size);
        omega.resize(size);
        inverse_zeta.resize(size);
        for (unsigned i = 0; i < size; ++i) {
            omega[i] = reduce(rotate(elements[i], 2));
            inverse_zeta[i] = reduce(rotate(elements[i], -1));
            for (unsigned j = 0; j < size; ++j)
                addition[i * size + j] = reduce(elements[i] + elements[j]);
        }
        need(reduce(d) == 0 && reduce(rotate(d, 1)) == 0, "wrong cusp residue ring");
    }
    unsigned reduce(Point p) const {
        I b = positive_mod(p.y, h22), q = (p.y - b) / h22;
        return static_cast<unsigned>(positive_mod(p.x - q * h12, h11) + h11 * b);
    }
    unsigned add(unsigned a, unsigned b) const { return addition[a * size + b]; }
    Point centered(Point p) const {
        Point q = multiply(p, conjugate(divisor));
        I a = Rational(2 * q.x + size, 2 * size).floor();
        I b = Rational(2 * q.y + size, 2 * size).floor();
        return p - multiply(divisor, {a, b});
    }
};
struct ProjectiveCapOrbit {
    Cap representative;
    U size, multiplier_order, scalar_cosets;
    bool open;
    Point conductor{};
};
using AlcoveVector = std::array<I, 8>;
using AlcoveMatrix = std::array<std::int16_t, 64>;
struct AlcoveMatrixHash {
    std::size_t operator()(const AlcoveMatrix &matrix) const {
        U value = 1469598103934665603ULL;
        for (auto x : matrix) {
            value ^= static_cast<std::uint16_t>(x);
            value *= 1099511628211ULL;
        }
        return value;
    }
};
std::int16_t alcove_narrow(I x) {
    need(x >= std::numeric_limits<std::int16_t>::min() &&
             x <= std::numeric_limits<std::int16_t>::max(), "alcove table coefficient overflow");
    return static_cast<std::int16_t>(x);
}
AlcoveVector alcove_flat(const Cap &cap) {
    AlcoveVector result{};
    for (int i = 0; i < 4; ++i) {
        result[2 * i] = cap[i].x;
        result[2 * i + 1] = cap[i].y;
    }
    return result;
}
Cap alcove_cap(const AlcoveVector &flat) {
    Cap result{};
    for (int i = 0; i < 4; ++i)
        result[i] = {flat[2 * i], flat[2 * i + 1]};
    return result;
}
I alcove_pair(const Cap &a, const Cap &b) {
    Cap sum{};
    for (int i = 0; i < 4; ++i)
        sum[i] = a[i] + b[i];
    return cap_index(sum) - cap_index(a) - cap_index(b);
}
AlcoveVector alcove_apply(const AlcoveMatrix &matrix, const AlcoveVector &x) {
    AlcoveVector result{};
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
            result[i] = identifier_narrow(identifier_add(result[i], identifier_mul(matrix[8 * i + j], x[j])));
    return result;
}
AlcoveMatrix alcove_right_reflect(const AlcoveMatrix &matrix, const AlcoveVector &root,
                                 const AlcoveVector &row) {
    auto image = alcove_apply(matrix, root);
    AlcoveMatrix result{};
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
            result[8 * i + j] = alcove_narrow(matrix[8 * i + j] - image[i] * row[j]);
    return result;
}
AlcoveMatrix alcove_conjugate(const AlcoveMatrix &matrix, const AlcoveVector &root,
                             const AlcoveVector &row) {
    auto result = alcove_right_reflect(matrix, root, row);
    AlcoveVector image{};
    for (int j = 0; j < 8; ++j)
        for (int i = 0; i < 8; ++i)
            image[j] += row[i] * result[8 * i + j];
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
            result[8 * i + j] = alcove_narrow(result[8 * i + j] - root[i] * image[j]);
    return result;
}
struct AlcoveTables {
    std::array<Cap, 8> simple_roots{}, fundamental_weights{};
    AlcoveMatrix cartan{}, inverse_cartan{}, basis{}, complex_structure{};
    AlcoveVector highest_marks{};
    std::array<AlcoveVector, 9> reflection_roots{}, reflection_rows{};
    std::vector<AlcoveMatrix> conjugates, representatives;
    std::vector<std::array<std::uint16_t, 9>> transitions;
    U edge_checks = 0;
    AlcoveTables() {
        std::set<Cap> positive;
        for (const auto &line : identifier_cap_roots())
            for (int unit = 0; unit < 6; ++unit) {
                Cap root{};
                for (int j = 0; j < 4; ++j)
                    root[j] = rotate(line[j], unit);
                auto flat = alcove_flat(root);
                auto first = std::find_if(flat.begin(), flat.end(), [](I x) { return x != 0; });
                need(first != flat.end(), "zero E8 root");
                if (*first > 0)
                    positive.insert(root);
            }
        need(positive.size() == 120, "incorrect positive E8 root system");
        std::vector<Cap> simple;
        for (const auto &root : positive) {
            bool decomposable = false;
            for (const auto &other : positive) {
                Cap difference{};
                for (int j = 0; j < 4; ++j)
                    difference[j] = root[j] - other[j];
                if (positive.contains(difference)) {
                    decomposable = true;
                    break;
                }
            }
            if (!decomposable)
                simple.push_back(root);
        }
        need(simple.size() == 8, "E8 root system does not have eight simple roots");
        std::copy(simple.begin(), simple.end(), simple_roots.begin());
        std::array<std::array<Rational, 16>, 8> augmented{};
        for (int i = 0; i < 8; ++i) {
            auto flat = alcove_flat(simple[i]);
            for (int j = 0; j < 8; ++j) {
                basis[8 * j + i] = alcove_narrow(flat[j]);
                I value = alcove_pair(simple[i], simple[j]);
                need(value == (i == j ? 2 : 0) || (i != j && value == -1), "incorrect E8 Cartan entry");
                cartan[8 * i + j] = alcove_narrow(value);
                augmented[i][j] = value;
            }
            augmented[i][8 + i] = 1;
        }
        for (int k = 0; k < 8; ++k) {
            int pivot = k;
            while (pivot < 8 && !augmented[pivot][k].n)
                ++pivot;
            need(pivot < 8, "singular E8 Cartan matrix");
            std::swap(augmented[pivot], augmented[k]);
            Rational divisor = augmented[k][k];
            for (auto &entry : augmented[k])
                entry = identifier_fraction_mul(entry, {divisor.d, divisor.n});
            for (int i = 0; i < 8; ++i) {
                if (i == k)
                    continue;
                Rational factor = augmented[i][k];
                for (int j = 0; j < 16; ++j) {
                    auto product = identifier_fraction_mul(factor, augmented[k][j]);
                    augmented[i][j] = identifier_fraction_add(augmented[i][j], {-product.n, product.d});
                }
            }
        }
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j) {
                need(augmented[i][8 + j].d == 1, "E8 Cartan inverse is not integral");
                inverse_cartan[8 * i + j] = alcove_narrow(augmented[i][8 + j].n);
            }
        I highest_height = 0;
        for (const auto &root : positive) {
            AlcoveVector pairings{};
            for (int i = 0; i < 8; ++i)
                pairings[i] = alcove_pair(simple[i], root);
            auto coefficients = alcove_apply(inverse_cartan, pairings);
            need(std::all_of(coefficients.begin(), coefficients.end(), [](I n) { return n >= 0; }),
                 "positive E8 root has negative simple coefficient");
            need(alcove_cap(alcove_apply(basis, coefficients)) == root, "simple-root coordinates do not reconstruct root");
            I height = std::accumulate(coefficients.begin(), coefficients.end(), I(0));
            if (height > highest_height) {
                highest_height = height;
                highest_marks = coefficients;
            }
        }
        need(highest_height == 29, "incorrect E8 highest-root height");
        for (int j = 0; j < 8; ++j) {
            AlcoveVector column{};
            for (int i = 0; i < 8; ++i)
                column[i] = inverse_cartan[8 * i + j];
            fundamental_weights[j] = alcove_cap(alcove_apply(basis, column));
            Cap rotated{};
            for (int i = 0; i < 4; ++i)
                rotated[i] = rotate(simple[j][i], 1);
            for (int i = 0; i < 8; ++i)
                column[i] = alcove_pair(simple[i], rotated);
            column = alcove_apply(inverse_cartan, column);
            need(alcove_cap(alcove_apply(basis, column)) == rotated, "incorrect E8 complex structure");
            for (int i = 0; i < 8; ++i)
                complex_structure[8 * i + j] = alcove_narrow(column[i]);
            reflection_roots[j][j] = 1;
            for (int i = 0; i < 8; ++i)
                reflection_rows[j][i] = cartan[8 * j + i];
        }
        reflection_roots[8] = highest_marks;
        reflection_rows[8] = alcove_apply(cartan, highest_marks);
        conjugates.reserve(4480);
        representatives.reserve(4480);
        transitions.reserve(4480);
        std::unordered_map<AlcoveMatrix, std::uint16_t, AlcoveMatrixHash> index;
        index.reserve(8192);
        index.emplace(complex_structure, 0);
        conjugates.push_back(complex_structure);
        representatives.push_back(basis);
        transitions.push_back({});
        for (std::size_t k = 0; k < conjugates.size(); ++k)
            for (int j = 0; j < 8; ++j) {
                auto next = alcove_conjugate(conjugates[k], reflection_roots[j], reflection_rows[j]);
                auto [position, fresh] = index.emplace(next, static_cast<std::uint16_t>(conjugates.size()));
                if (fresh) {
                    need(conjugates.size() < 4480, "E8 complex-structure orbit exceeds4480");
                    conjugates.push_back(next);
                    representatives.push_back(alcove_right_reflect(representatives[k], reflection_roots[j], reflection_rows[j]));
                    transitions.push_back({});
                }
                transitions[k][j] = position->second;
                ++edge_checks;
            }
        need(conjugates.size() == 4480, "E8 complex-structure orbit is not4480");
        for (std::size_t k = 0; k < conjugates.size(); ++k) {
            auto next = alcove_conjugate(conjugates[k], reflection_roots[8], reflection_rows[8]);
            auto position = index.find(next);
            need(position != index.end(), "highest-root reflection leaves E8 cosets");
            transitions[k][8] = position->second;
            ++edge_checks;
        }
        for (std::size_t k = 0; k < transitions.size(); ++k)
            for (int j = 0; j < 9; ++j)
                need(transitions[transitions[k][j]][j] == k, "alcove coset reflection is not involutive");
    }
};
const AlcoveTables &alcove_tables() {
    static const AlcoveTables result;
    return result;
}
struct AlcoveComponentRecord {
    std::uint16_t coset, component_size;
    U orbit_size;
    bool open;
};
struct AlcovePointRecord {
    AlcoveVector labels;
    std::uint16_t zero_wall_mask;
    U stabilizer_order;
    std::vector<AlcoveComponentRecord> components;
};
// Research extension point of analyze_general_joined_caps.py: a probe may
// replace the exact phase-subgroup computation of the cap inventories. Its
// arguments are the inventory modulus, the local scalar content, the residue
// stabilizer order, the residue ring and the reduced cap. Null (the default)
// keeps the exact chain image.
using InventoryPhaseSubgroup = I (*)(Point full_height, I local_content, I stabilizer,
                                     const IdentifierRing &ring, const Cap &cap,
                                     IdentifierChainStats &stats);
InventoryPhaseSubgroup inventory_phase_subgroup = nullptr;
I inventory_phase(Point full_height, I local_content, I stabilizer, const IdentifierRing &ring,
                  const Cap &cap, IdentifierChainStats &stats) {
    if (inventory_phase_subgroup)
        return inventory_phase_subgroup(full_height, local_content, stabilizer, ring, cap, stats);
    return identifier_chain_image(ring, cap, stats).subgroup;
}
struct CuspCaps {
    std::vector<Cap> representatives;
    U states = 0, open_states = 0, orbits = 0;
    std::vector<I> fibers{};
    U visited_bytes = 0, phase_bytes = 0, phase_edge_checks = 0;
    IdentifierChainStats phase_chain_stats{};
    bool projective = false;
    U projective_states = 0, projective_edges = 0, scalar_lifts = 0, multiplier_entries_peak = 0;
    std::vector<ProjectiveCapOrbit> projective_orbits{};
    bool composite = false;
    U content_strata = 0;
    bool alcove = false;
    U alcove_points = 0, alcove_faces = 0, alcove_face_edges = 0, alcove_component_tests = 0;
    U alcove_table_bytes = 0, alcove_face_bytes = 0;
    U alcove_root_checks = 0, alcove_openness_tests = 0;
    std::vector<AlcovePointRecord> alcove_records{};
};
CuspCaps affine_cusp_caps(Point d, bool dense_phases = false) {
    CuspRing ring(d);
    U D = ring.size;
    need(D <= 65535, "cusp height exceeds native residue indexing envelope");
    std::array<U, 4> power{1, D, D * D, D * D * D};
    U count = power[3] * D;
    need(count / D == power[3], "cusp residue count overflow");
    std::vector<U> visited((count + 63) / 64);
    Point beta = rotate(d, 1);
    auto [g, bs, bt] = bezout(beta.x + beta.y, beta.y);
    std::vector<unsigned char> phases(dense_phases && g > 1 ? count : 0);
    auto mark = [&](U index) {
        U mask = U(1) << (index % 64);
        U &word = visited[index / 64];
        bool fresh = !(word & mask);
        word |= mask;
        return fresh;
    };
    CuspCaps result;
    result.states = count;
    result.visited_bytes = visited.size() * sizeof(U);
    result.phase_bytes = phases.size();
    IdentifierRing phase_ring(beta);
    std::vector<U> queue;
    queue.reserve(155520);
    const auto &rows = cap_root_functionals();
    for (U first = 0; first < count; ++first) {
        if (!mark(first))
            continue;
        ++result.orbits;
        Cap representative;
        U n = first;
        for (int j = 0; j < 4; ++j) {
            representative[j] = ring.elements[n % D];
            n /= D;
        }
        bool open = true;
        for (const auto &row : rows) {
            Point value{};
            for (int j = 0; j < 4; ++j)
                value = value + multiply(row[j], representative[j]);
            if (ring.reduce(value) == 0) {
                open = false;
                break;
            }
        }
        queue.clear();
        queue.push_back(first);
        I fiber_index = g;
        if (dense_phases && g > 1 && open)
            phases[first] = 0;
        for (std::size_t i = 0; i < queue.size(); ++i) {
            std::array<unsigned, 4> x;
            n = queue[i];
            for (int j = 0; j < 4; ++j) {
                x[j] = n % D;
                n /= D;
            }
            for (int j = 0; j < 4; ++j) {
                unsigned y = ring.omega[x[j]];
                if (j)
                    y = ring.add(y, x[j - 1]);
                if (j < 3)
                    y = ring.add(y, ring.inverse_zeta[x[j + 1]]);
                U other = queue[i] - U(x[j]) * power[j] + U(y) * power[j];
                I phase = 0;
                if (dense_phases && g > 1 && open) {
                    identifier_count_add(result.phase_edge_checks, 1);
                    Point rx = rotate(ring.elements[x[j]], 2) +
                               (j ? ring.elements[x[j - 1]] : Point{}) +
                               (j < 3 ? rotate(ring.elements[x[j + 1]], -1) : Point{});
                    Point numerator = multiply(ring.elements[y] - rx, conjugate(beta));
                    need(numerator.x % static_cast<I>(D) == 0 &&
                             numerator.y % static_cast<I>(D) == 0,
                         "nonintegral cap Schreier translation");
                    Point tau{numerator.x / static_cast<I>(D), numerator.y / static_cast<I>(D)};
                    Point row{3 * rx.x, 3 * rx.y};
                    if (j)
                        row = row - multiply({1, 1}, ring.elements[x[j - 1]]);
                    if (j < 3)
                        row = row - multiply({2, -1}, ring.elements[x[j + 1]]);
                    Point delta = multiply(beta, {norm(tau), 0}) +
                                  divide_theta(multiply(conjugate(tau), row));
                    phase = positive_mod(phases[queue[i]] + bs * delta.y - bt * delta.x, g);
                }
                if (mark(other)) {
                    queue.push_back(other);
                    if (dense_phases && g > 1 && open)
                        phases[other] = static_cast<unsigned char>(phase);
                } else if (dense_phases && g > 1 && open)
                    fiber_index = std::gcd(fiber_index, positive_mod(phase - phases[other], g));
            }
        }
        need(queue.size() <= 155520 && 155520 % queue.size() == 0, "wrong cap orbit size");
        if (open) {
            result.open_states += queue.size();
            for (auto &p : representative)
                p = ring.centered(p);
            representative = short_cap(representative, d);
            I stabilizer_order = 155520 / static_cast<I>(queue.size());
            if (!dense_phases && std::gcd(g, stabilizer_order) > 1) {
                // The finite residue traversal proves orbit coverage. Only its
                // stabilizer phase image is needed; the exact fixed chain gets
                // that image once, without a phase byte at every residue state.
                fiber_index = inventory_phase(d, g, stabilizer_order, phase_ring,
                                              phase_ring.reduce(representative), result.phase_chain_stats);
            }
            // The phase-image order divides both g and the residue stabilizer
            // order. A coprime stabilizer therefore leaves the initialized g.
            result.representatives.push_back(representative);
            result.fibers.push_back(fiber_index);
        }
    }
    return result;
}
bool inventory_prime(U n) {
    if (n < 2)
        return false;
    for (U p = 2; p <= n / p; ++p)
        if (n % p == 0)
            return false;
    return true;
}
bool inventory_field(Point d) {
    I D = identifier_narrow(identifier_norm(identifier_point(d)));
    if (inventory_prime(D))
        return true;
    I p = identifier_narrow(identifier_sqrt(D));
    return p * p == D && p % 3 == 2 && inventory_prime(p) &&
           d.x % p == 0 && d.y % p == 0 && norm({d.x / p, d.y / p}) == 1;
}
struct InventoryField {
    CuspRing ring;
    unsigned size, period;
    std::vector<unsigned> powers, logarithms;
    explicit InventoryField(Point d) : ring(d), size(static_cast<unsigned>(ring.size)), period(size - 1) {
        need(inventory_field(d), "projective inventory requires a residue field");
        std::vector<unsigned> factors;
        unsigned remaining = period;
        for (unsigned p = 2; p <= remaining / p; ++p)
            if (remaining % p == 0) {
                factors.push_back(p);
                while (remaining % p == 0)
                    remaining /= p;
            }
        if (remaining > 1)
            factors.push_back(remaining);
        auto raw_multiply = [&](unsigned a, unsigned b) {
            return ring.reduce(multiply(ring.elements[a], ring.elements[b]));
        };
        auto power = [&](unsigned a, unsigned exponent) {
            unsigned result = 1;
            while (exponent) {
                if (exponent & 1)
                    result = raw_multiply(result, a);
                a = raw_multiply(a, a);
                exponent >>= 1;
            }
            return result;
        };
        unsigned generator = 0;
        for (unsigned a = 1; a < size && !generator; ++a) {
            bool primitive = power(a, period) == 1;
            for (unsigned p : factors)
                primitive = primitive && power(a, period / p) != 1;
            if (primitive)
                generator = a;
        }
        need(generator != 0, "residue field has no primitive multiplicative generator");
        logarithms.assign(size, std::numeric_limits<unsigned>::max());
        powers.reserve(period);
        unsigned value = 1;
        for (unsigned exponent = 0; exponent < period; ++exponent) {
            need(value != 0 && logarithms[value] == std::numeric_limits<unsigned>::max(),
                 "invalid residue field multiplicative cycle");
            logarithms[value] = exponent;
            powers.push_back(value);
            value = raw_multiply(value, generator);
        }
        need(value == 1, "residue field multiplicative cycle did not close");
    }
    unsigned mul(unsigned a, unsigned b) const {
        return !a || !b ? 0 : powers[(logarithms[a] + logarithms[b]) % period];
    }
};
CuspCaps projective_cusp_caps(Point d) {
    InventoryField field(d);
    const auto &ring = field.ring;
    U D = field.size;
    need(D <= 65535, "projective cap height exceeds native indexing envelope");
    std::array<U, 5> power{1, D, D * D, D * D * D, D * D * D * D};
    U count = power[3] + power[2] + power[1] + 1;
    using ResidueCap = std::array<unsigned, 4>;
    auto decode = [&](U index) {
        ResidueCap point{};
        for (int first = 0; first < 4; ++first) {
            U block = power[3 - first];
            if (index >= block) {
                index -= block;
                continue;
            }
            point[first] = 1;
            for (int j = first + 1; j < 4; ++j) {
                point[j] = static_cast<unsigned>(index % D);
                index /= D;
            }
            return point;
        }
        throw std::runtime_error("projective cap index outside its domain");
    };
    auto encode = [&](const ResidueCap &point) {
        U offset = 0;
        for (int first = 0; first < 4; ++first) {
            if (!point[first]) {
                offset += power[3 - first];
                continue;
            }
            need(point[first] == 1, "unnormalized projective cap section");
            U tail = 0, place = 1;
            for (int j = first + 1; j < 4; ++j) {
                tail += point[j] * place;
                place *= D;
            }
            return offset + tail;
        }
        throw std::runtime_error("zero projective cap");
    };
    auto lift = [&](const ResidueCap &point, unsigned scalar) {
        Cap cap;
        for (int j = 0; j < 4; ++j)
            cap[j] = ring.elements[field.mul(point[j], scalar)];
        return cap;
    };
    std::vector<U> visited((count + 63) / 64);
    auto mark = [&](U index) {
        need(index < count, "projective cap index overflow");
        U mask = U(1) << (index % 64);
        U &word = visited[index / 64];
        bool fresh = !(word & mask);
        word |= mask;
        return fresh;
    };
    CuspCaps result;
    result.projective = true;
    result.states = power[4];
    result.projective_states = count;
    result.visited_bytes = visited.size() * sizeof(U);
    result.orbits = 1; // The zero affine cap is one closed orbit.
    Point beta = rotate(d, 1);
    IdentifierRing phase_ring(beta);
    I g = phase_ring.h22;
    const auto &root_functionals = cap_root_functionals();
    std::vector<U> queue;
    std::vector<I> path_exponents;
    std::unordered_map<U, std::size_t> orbit_index;
    std::set<unsigned> scalar_units;
    for (int k = 0; k < 6; ++k)
        scalar_units.insert(ring.reduce(rotate({1, 0}, k)));
    U orbit_bound = 155520 / scalar_units.size();
    U reserve = std::min(count, orbit_bound);
    queue.reserve(reserve);
    path_exponents.reserve(reserve);
    orbit_index.reserve(reserve);
    for (U first = 0; first < count; ++first) {
        if (!mark(first))
            continue;
        ResidueCap base = decode(first);
        Cap base_lift = lift(base, 1);
        bool open = true;
        // These are integral coefficients of h(r,z)/theta. Reducing them is
        // valid even in characteristic three; no field division is performed.
        for (const auto &row : root_functionals) {
            Point value{};
            for (int j = 0; j < 4; ++j)
                value = value + multiply(row[j], base_lift[j]);
            if (ring.reduce(value) == 0) {
                open = false;
                break;
            }
        }
        queue.assign(1, first);
        path_exponents.assign(1, 0);
        orbit_index.clear();
        orbit_index.emplace(first, 0);
        I scalar_cosets = field.period;
        for (std::size_t i = 0; i < queue.size(); ++i) {
            ResidueCap current = decode(queue[i]);
            for (int j = 0; j < 4; ++j) {
                ResidueCap next = current;
                next[j] = ring.omega[current[j]];
                if (j)
                    next[j] = ring.add(next[j], current[j - 1]);
                if (j < 3)
                    next[j] = ring.add(next[j], ring.inverse_zeta[current[j + 1]]);
                unsigned pivot = 0;
                for (unsigned value : next)
                    if (value) {
                        pivot = value;
                        break;
                    }
                need(pivot != 0, "reflection killed a projective cap");
                I exponent = positive_mod(path_exponents[i] + field.logarithms[pivot], field.period);
                unsigned inverse = field.powers[(field.period - field.logarithms[pivot]) % field.period];
                for (auto &value : next)
                    value = field.mul(value, inverse);
                U index = encode(next);
                identifier_count_add(result.projective_edges, 1);
                if (mark(index)) {
                    need(queue.size() < orbit_bound, "projective cap orbit exceeds its group order");
                    orbit_index.emplace(index, queue.size());
                    queue.push_back(index);
                    path_exponents.push_back(exponent);
                } else {
                    auto found = orbit_index.find(index);
                    need(found != orbit_index.end(), "projective edge crossed completed orbits");
                    scalar_cosets = std::gcd(scalar_cosets,
                                             positive_mod(exponent - path_exponents[found->second], field.period));
                }
            }
        }
        result.multiplier_entries_peak = std::max<U>(result.multiplier_entries_peak, queue.size());
        U multiplier_order = field.period / scalar_cosets;
        U affine_orbit_size = queue.size() * multiplier_order;
        need(affine_orbit_size > 0 && 155520 % affine_orbit_size == 0,
             "projective scalar lift has incorrect affine orbit size");
        for (int k = 0; k < 6; ++k) {
            unsigned unit = ring.reduce(rotate({1, 0}, k));
            need(unit && field.logarithms[unit] % scalar_cosets == 0,
                 "projective multiplier image omitted an Eisenstein unit");
        }
        identifier_count_add(result.orbits, scalar_cosets);
        result.projective_orbits.push_back({base_lift, queue.size(), multiplier_order,
                                            static_cast<U>(scalar_cosets), open});
        if (!open)
            continue;
        identifier_count_add(result.open_states, queue.size() * field.period);
        for (I exponent = 0; exponent < scalar_cosets; ++exponent) {
            Cap representative = lift(base, field.powers[exponent]);
            for (auto &p : representative)
                p = ring.centered(p);
            representative = short_cap(representative, d);
            I fiber_index = g;
            I stabilizer_order = 155520 / static_cast<I>(affine_orbit_size);
            if (std::gcd(g, stabilizer_order) > 1) {
                fiber_index = inventory_phase(d, g, stabilizer_order, phase_ring,
                                              phase_ring.reduce(representative), result.phase_chain_stats);
            }
            // Scalar lifting only selects a residue vector. It is not an
            // Eisenstein-unit isometry: its norm and central fiber are computed
            // from this integral representative, never from field conjugation.
            identifier_count_add(result.scalar_lifts, 1);
            result.representatives.push_back(representative);
            result.fibers.push_back(fiber_index);
        }
    }
    need(result.projective_edges == 4 * count, "incomplete projective cap inventory");
    return result;
}
std::optional<Point> inventory_quotient(Point numerator, Point denominator) {
    auto product = identifier_mul(identifier_point(numerator),
                                  identifier_conjugate(identifier_point(denominator)));
    Wide D = identifier_norm(identifier_point(denominator));
    need(D > 0, "zero composite inventory divisor");
    if (product[0] % D || product[1] % D)
        return {};
    return identifier_narrow({product[0] / D, product[1] / D});
}
std::vector<Point> inventory_divisors(Point d) {
    I D = norm(d), limit = identifier_narrow(identifier_sqrt(D));
    std::vector<Point> result;
    for (I a = 1; a <= limit; ++a)
        for (I b = 0; b <= limit; ++b) {
            Point candidate{a, b};
            I n = norm(candidate);
            if (n <= D && D % n == 0 && inventory_quotient(d, candidate))
                result.push_back(candidate);
        }
    std::sort(result.begin(), result.end(), [](Point a, Point b) {
        return std::tuple{norm(a), a.x, a.y} < std::tuple{norm(b), b.x, b.y};
    });
    need(!result.empty() && result.front() == Point{1, 0}, "missing unit conductor");
    return result;
}
struct InventoryRing {
    CuspRing ring;
    unsigned size;
    std::vector<unsigned> products, inverses, units;
    explicit InventoryRing(Point d) : ring(d), size(static_cast<unsigned>(ring.size)),
        products(static_cast<std::size_t>(size) * size), inverses(size, size) {
        for (unsigned a = 0; a < size; ++a)
            for (unsigned b = 0; b < size; ++b) {
                unsigned product = ring.reduce(multiply(ring.elements[a], ring.elements[b]));
                products[static_cast<std::size_t>(a) * size + b] = product;
                if (product == 1)
                    inverses[a] = b;
            }
        for (unsigned a = 0; a < size; ++a)
            if (inverses[a] != size)
                units.push_back(a);
        need(size > 1 && !units.empty(), "invalid finite inventory ring");
    }
    unsigned mul(unsigned a, unsigned b) const { return products[static_cast<std::size_t>(a) * size + b]; }
    bool unit(unsigned a) const { return inverses[a] != size; }
};
using InventoryResidueCap = std::array<unsigned, 4>;
struct InventoryLocalChart {
    InventoryRing arithmetic;
    std::vector<unsigned> nonunits, nonunit_index;
    std::array<U, 4> blocks{};
    U count = 0;
    explicit InventoryLocalChart(Point d) : arithmetic(d), nonunit_index(arithmetic.size, arithmetic.size) {
        for (unsigned a = 0; a < arithmetic.size; ++a)
            if (!arithmetic.unit(a)) {
                nonunit_index[a] = static_cast<unsigned>(nonunits.size());
                nonunits.push_back(a);
            }
        U before = 1;
        for (int first = 0; first < 4; ++first) {
            U tail = 1;
            for (int j = first + 1; j < 4; ++j)
                tail *= arithmetic.size;
            blocks[first] = before * tail;
            count += blocks[first];
            before *= nonunits.size();
        }
    }
    InventoryResidueCap decode(U index) const {
        for (int first = 0; first < 4; ++first) {
            if (index >= blocks[first]) {
                index -= blocks[first];
                continue;
            }
            InventoryResidueCap point{};
            for (int j = 0; j < first; ++j) {
                point[j] = nonunits[index % nonunits.size()];
                index /= nonunits.size();
            }
            point[first] = 1;
            for (int j = first + 1; j < 4; ++j) {
                point[j] = static_cast<unsigned>(index % arithmetic.size);
                index /= arithmetic.size;
            }
            need(index == 0, "local projective section decode overflow");
            return point;
        }
        throw std::runtime_error("local projective section index overflow");
    }
    U encode(const InventoryResidueCap &point) const {
        U offset = 0;
        for (int first = 0; first < 4; ++first) {
            if (!arithmetic.unit(point[first])) {
                offset += blocks[first];
                continue;
            }
            need(point[first] == 1, "unnormalized local projective section");
            U index = 0, place = 1;
            for (int j = 0; j < first; ++j) {
                need(nonunit_index[point[j]] < nonunits.size(), "wrong local projective chart");
                index += place * nonunit_index[point[j]];
                place *= nonunits.size();
            }
            for (int j = first + 1; j < 4; ++j) {
                index += place * point[j];
                place *= arithmetic.size;
            }
            need(index < blocks[first], "local projective section encode overflow");
            return offset + index;
        }
        throw std::runtime_error("nonprimitive local projective vector");
    }
};
struct InventoryProjectiveRing {
    InventoryRing arithmetic;
    std::vector<InventoryLocalChart> factors;
    std::vector<std::vector<unsigned>> reductions;
    std::vector<unsigned> crt;
    U count = 1;
    explicit InventoryProjectiveRing(Point d) : arithmetic(d), crt(arithmetic.size, arithmetic.size) {
        U product_norm = 1;
        for (Point prime : inventory_divisors(d)) {
            if (norm(prime) == 1 || !inventory_field(prime))
                continue;
            Point power = prime;
            while (inventory_quotient(d, multiply(power, prime)))
                power = multiply(power, prime);
            factors.emplace_back(power);
            count *= factors.back().count;
            product_norm *= norm(power);
        }
        need(!factors.empty() && product_norm == arithmetic.size,
             "incomplete prime-power factorization for composite inventory");
        for (const auto &factor : factors) {
            std::vector<unsigned> values;
            for (Point p : arithmetic.ring.elements)
                values.push_back(factor.arithmetic.ring.reduce(p));
            reductions.push_back(std::move(values));
        }
        for (unsigned value = 0; value < arithmetic.size; ++value) {
            U key = 0, place = 1;
            for (std::size_t i = 0; i < factors.size(); ++i) {
                key += place * reductions[i][value];
                place *= factors[i].arithmetic.size;
            }
            need(key < crt.size() && crt[key] == arithmetic.size, "nonbijective CRT inventory table");
            crt[key] = value;
        }
    }
    InventoryResidueCap decode(U index) const {
        InventoryResidueCap key{};
        unsigned place = 1;
        for (const auto &factor : factors) {
            auto local = factor.decode(index % factor.count);
            index /= factor.count;
            for (int j = 0; j < 4; ++j)
                key[j] += place * local[j];
            place *= factor.arithmetic.size;
        }
        need(index == 0, "global projective section decode overflow");
        for (auto &value : key)
            value = crt[value];
        return key;
    }
    std::pair<U, unsigned> normalize(const InventoryResidueCap &point) const {
        U index = 0, place = 1;
        unsigned scalar_key = 0, scalar_place = 1;
        for (std::size_t i = 0; i < factors.size(); ++i) {
            const auto &factor = factors[i];
            InventoryResidueCap local{};
            unsigned pivot = factor.arithmetic.size;
            for (int j = 0; j < 4; ++j) {
                local[j] = reductions[i][point[j]];
                if (pivot == factor.arithmetic.size && factor.arithmetic.unit(local[j]))
                    pivot = local[j];
            }
            need(pivot != factor.arithmetic.size, "reflection lost primitive local content");
            unsigned inverse = factor.arithmetic.inverses[pivot];
            for (auto &value : local)
                value = factor.arithmetic.mul(value, inverse);
            index += place * factor.encode(local);
            place *= factor.count;
            scalar_key += scalar_place * pivot;
            scalar_place *= factor.arithmetic.size;
        }
        need(index < count && scalar_key < crt.size(), "global projective normalization overflow");
        unsigned scalar = crt[scalar_key];
        need(arithmetic.unit(scalar), "nonunit projective normalization scalar");
        return {index, scalar};
    }
};
CuspCaps composite_cusp_caps(Point d) {
    I D = norm(d);
    need(D > 0 && D <= 65535, "composite cap height exceeds native indexing envelope");
    CuspCaps result;
    result.composite = true;
    result.states = U(D) * D * D * D;
    result.orbits = 1;
    result.content_strata = 1; // The zero vector has conductor one.
    IdentifierRing target_phase_ring(rotate(d, 1));
    const auto &root_functionals = cap_root_functionals();
    U recovered_states = 1;
    for (Point e : inventory_divisors(d)) {
        if (norm(e) == 1)
            continue;
        ++result.content_strata;
        Point a = *inventory_quotient(d, e);
        InventoryProjectiveRing projective(e);
        const auto &arithmetic = projective.arithmetic;
        const auto &ring = arithmetic.ring;
        IdentifierRing phase_ring(rotate(e, 1));
        need(target_phase_ring.h22 % phase_ring.h22 == 0, "nonintegral content phase scaling");
        I phase_scale = target_phase_ring.h22 / phase_ring.h22;
        std::vector<U> visited((projective.count + 63) / 64);
        result.visited_bytes = std::max<U>(result.visited_bytes, visited.size() * sizeof(U));
        result.projective_states += projective.count;
        auto mark = [&](U index) {
            U mask = U(1) << (index % 64);
            U &word = visited[index / 64];
            bool fresh = !(word & mask);
            word |= mask;
            return fresh;
        };
        std::vector<U> queue;
        std::vector<unsigned> paths;
        std::unordered_map<U, std::size_t> orbit_index;
        for (U first = 0; first < projective.count; ++first) {
            if (!mark(first))
                continue;
            auto base = projective.decode(first);
            Cap base_lift;
            for (int j = 0; j < 4; ++j)
                base_lift[j] = ring.elements[base[j]];
            bool open = true;
            // A root is excluded only when its integral functional vanishes
            // modulo the whole conductor, not separately in some CRT factor.
            for (const auto &row : root_functionals) {
                Point value{};
                for (int j = 0; j < 4; ++j)
                    value = value + multiply(row[j], base_lift[j]);
                if (ring.reduce(value) == 0) {
                    open = false;
                    break;
                }
            }
            queue.assign(1, first);
            paths.assign(1, 1);
            orbit_index.clear();
            orbit_index.emplace(first, 0);
            std::vector<unsigned> subgroup{1};
            std::vector<bool> in_subgroup(arithmetic.size);
            in_subgroup[1] = true;
            auto add_multiplier = [&](unsigned generator) {
                need(arithmetic.unit(generator), "nonunit Schreier multiplier");
                if (in_subgroup[generator])
                    return;
                auto previous = subgroup;
                unsigned power = generator;
                while (!in_subgroup[power]) {
                    for (unsigned h : previous) {
                        unsigned value = arithmetic.mul(h, power);
                        need(!in_subgroup[value], "overlapping multiplier subgroup cosets");
                        in_subgroup[value] = true;
                        subgroup.push_back(value);
                    }
                    power = arithmetic.mul(power, generator);
                }
            };
            for (std::size_t i = 0; i < queue.size(); ++i) {
                auto current = projective.decode(queue[i]);
                for (int j = 0; j < 4; ++j) {
                    auto next = current;
                    next[j] = ring.omega[current[j]];
                    if (j)
                        next[j] = ring.add(next[j], current[j - 1]);
                    if (j < 3)
                        next[j] = ring.add(next[j], ring.inverse_zeta[current[j + 1]]);
                    auto [index, scalar] = projective.normalize(next);
                    unsigned path = arithmetic.mul(paths[i], scalar);
                    ++result.projective_edges;
                    if (mark(index)) {
                        need(queue.size() < 155520, "primitive projective orbit exceeds group order");
                        orbit_index.emplace(index, queue.size());
                        queue.push_back(index);
                        paths.push_back(path);
                    } else {
                        auto found = orbit_index.find(index);
                        need(found != orbit_index.end(), "primitive projective edge crossed orbits");
                        add_multiplier(arithmetic.mul(path, arithmetic.inverses[paths[found->second]]));
                    }
                }
            }
            result.multiplier_entries_peak = std::max<U>(result.multiplier_entries_peak, queue.size());
            U orbit_size = queue.size() * subgroup.size();
            need(orbit_size && 155520 % orbit_size == 0 && arithmetic.units.size() % subgroup.size() == 0,
                 "wrong primitive projective scalar orbit size");
            for (int k = 0; k < 6; ++k)
                need(in_subgroup[ring.reduce(rotate({1, 0}, k))], "missing integral scalar unit in multiplier subgroup");
            U scalar_cosets = arithmetic.units.size() / subgroup.size();
            result.orbits += scalar_cosets;
            recovered_states += queue.size() * arithmetic.units.size();
            result.projective_orbits.push_back({base_lift, queue.size(), subgroup.size(), scalar_cosets, open, e});
            if (!open)
                continue;
            result.open_states += queue.size() * arithmetic.units.size();
            I fiber = phase_ring.h22;
            I stabilizer_order = 155520 / static_cast<I>(orbit_size);
            if (std::gcd(fiber, stabilizer_order) > 1) {
                // Unit scalars modulo e multiply every stabilizer correction
                // by N(c), which is coprime to g_e. Thus the phase-image
                // subgroup is common to all scalar lifts of this orbit.
                fiber = inventory_phase(d, phase_ring.h22, stabilizer_order, phase_ring,
                                        phase_ring.reduce(base_lift), result.phase_chain_stats);
            }
            std::vector<bool> covered(arithmetic.size);
            U emitted = 0;
            for (unsigned scalar : arithmetic.units) {
                if (covered[scalar])
                    continue;
                for (unsigned h : subgroup)
                    covered[arithmetic.mul(scalar, h)] = true;
                Cap representative;
                for (int j = 0; j < 4; ++j)
                    representative[j] = ring.centered(ring.elements[arithmetic.mul(scalar, base[j])]);
                representative = short_cap(representative, e);
                for (auto &p : representative)
                    p = multiply(a, p);
                result.representatives.push_back(representative);
                result.fibers.push_back(phase_scale * fiber);
                ++result.scalar_lifts;
                ++emitted;
            }
            need(emitted == scalar_cosets, "incomplete noncyclic scalar coset lifts");
        }
    }
    need(recovered_states == result.states, "incomplete primitive-content inventory decomposition");
    need(result.projective_edges == 4 * result.projective_states, "incomplete composite projective traversal");
    return result;
}
U alcove_wall_order(const AlcoveTables &tables, unsigned mask) {
    std::array<AlcoveVector, 9> roots = tables.reflection_roots;
    for (auto &x : roots[8])
        x = -x;
    std::array<std::vector<int>, 9> adjacent{};
    for (int i = 0; i < 9; ++i)
        if (mask & (1U << i)) {
            auto row = alcove_apply(tables.cartan, roots[i]);
            for (int j = i + 1; j < 9; ++j)
                if (mask & (1U << j)) {
                    I inner = 0;
                    for (int k = 0; k < 8; ++k)
                        inner += row[k] * roots[j][k];
                    need(inner == 0 || inner == -1, "invalid affine E8 wall diagram");
                    if (inner == -1) {
                        adjacent[i].push_back(j);
                        adjacent[j].push_back(i);
                    }
                }
        }
    auto factorial = [](int n) {
        U result = 1;
        for (int k = 2; k <= n; ++k)
            result *= k;
        return result;
    };
    U order = 1;
    std::array<bool, 9> visited{};
    for (int first = 0; first < 9; ++first) {
        if (!(mask & (1U << first)) || visited[first])
            continue;
        std::vector<int> queue{first};
        visited[first] = true;
        int branch = -1, edge_sum = 0;
        for (std::size_t k = 0; k < queue.size(); ++k) {
            int vertex = queue[k];
            edge_sum += adjacent[vertex].size();
            if (adjacent[vertex].size() == 3) {
                need(branch == -1, "multiple branches in finite E8 wall component");
                branch = vertex;
            }
            need(adjacent[vertex].size() <= 3, "excessive finite E8 wall degree");
            for (int next_vertex : adjacent[vertex])
                if (!visited[next_vertex]) {
                    visited[next_vertex] = true;
                    queue.push_back(next_vertex);
                }
        }
        int n = queue.size();
        need(edge_sum == 2 * (n - 1), "cyclic finite E8 wall component");
        U factor = factorial(n + 1);
        if (branch != -1) {
            std::vector<int> arms;
            for (int first_arm : adjacent[branch]) {
                int length = 1, last = branch, current = first_arm;
                while (adjacent[current].size() == 2) {
                    int next_vertex = adjacent[current][0] == last ? adjacent[current][1] : adjacent[current][0];
                    last = current;
                    current = next_vertex;
                    ++length;
                }
                need(adjacent[current].size() == 1, "finite E8 arm does not end");
                arms.push_back(length);
            }
            std::sort(arms.begin(), arms.end());
            if (arms[0] == 1 && arms[1] == 1)
                factor = (U(1) << (n - 1)) * factorial(n);
            else if (arms == std::vector<int>{1, 2, 2})
                factor = 51840;
            else if (arms == std::vector<int>{1, 2, 3})
                factor = 2903040;
            else {
                need(arms == std::vector<int>{1, 2, 4}, "unknown finite E8 wall type");
                factor = 696729600;
            }
        }
        order *= factor;
    }
    need(order > 0 && 696729600 % order == 0, "finite wall order does not divide E8 Weyl order");
    return order;
}
CuspCaps alcove_cusp_caps(Point d, bool full_root_checks = false) {
    need(d.x >= 1 && d.y == 0 && d.x <= 255,
         "alcove inventory currently requires a positive rational modulus at most255");
    const auto &tables = alcove_tables();
    I p = d.x;
    CuspCaps result;
    result.alcove = true;
    result.states = 1;
    for (int i = 0; i < 8; ++i) {
        need(result.states <= std::numeric_limits<U>::max() / p, "alcove residue count overflow");
        result.states *= p;
    }
    result.alcove_table_bytes = tables.conjugates.size() * sizeof(AlcoveMatrix) +
        tables.representatives.size() * sizeof(AlcoveMatrix) + tables.transitions.size() * sizeof(tables.transitions[0]);
    struct Face {
        struct Component {
            std::uint16_t coset, size;
            std::optional<bool> open;
        };
        U order;
        std::vector<Component> components;
    };
    std::map<unsigned, Face> faces;
    auto face_for = [&](unsigned mask) -> Face & {
        auto existing = faces.find(mask);
        if (existing != faces.end())
            return existing->second;
        Face face{alcove_wall_order(tables, mask), {}};
        std::array<bool, 4480> visited{};
        std::vector<std::uint16_t> queue;
        queue.reserve(4480);
        for (unsigned first = 0; first < 4480; ++first) {
            if (visited[first])
                continue;
            queue.clear();
            queue.push_back(static_cast<std::uint16_t>(first));
            visited[first] = true;
            for (std::size_t k = 0; k < queue.size(); ++k)
                for (int j = 0; j < 9; ++j)
                    if (mask & (1U << j)) {
                        identifier_count_add(result.alcove_face_edges, 1);
                        unsigned next_state = tables.transitions[queue[k]][j];
                        if (!visited[next_state]) {
                            visited[next_state] = true;
                            queue.push_back(next_state);
                        }
                    }
            need(face.order % queue.size() == 0, "alcove component does not divide finite wall order");
            face.components.push_back({static_cast<std::uint16_t>(first), static_cast<std::uint16_t>(queue.size()), std::nullopt});
        }
        result.alcove_face_bytes += face.components.size() * sizeof(face.components[0]);
        ++result.alcove_faces;
        return faces.emplace(mask, std::move(face)).first->second;
    };
    IdentifierRing phase_ring(rotate(d, 1));
    const auto &functionals = cap_root_functionals();
    AlcoveVector labels{};
    U covered_states = 0;
    std::function<void(int, I)> visit = [&](int level, I remaining) {
        if (level < 8) {
            for (I n = 0; n <= remaining / tables.highest_marks[level]; ++n) {
                labels[level] = n;
                visit(level + 1, remaining - n * tables.highest_marks[level]);
            }
            return;
        }
        unsigned mask = remaining == 0 ? (1U << 8) : 0;
        for (int i = 0; i < 8; ++i)
            if (labels[i] == 0)
                mask |= 1U << i;
        auto &face = face_for(mask);
        auto point = alcove_apply(tables.inverse_cartan, labels);
        AlcovePointRecord record{labels, static_cast<std::uint16_t>(mask), face.order, {}};
        record.components.reserve(face.components.size());
        ++result.alcove_points;
        for (auto &component : face.components) {
            U stabilizer = face.order / component.size;
            need(stabilizer > 0 && 155520 % stabilizer == 0, "invalid unitary stabilizer from alcove double coset");
            U orbit_size = 155520 / stabilizer;
            identifier_count_add(covered_states, orbit_size);
            identifier_count_add(result.orbits, 1);
            identifier_count_add(result.alcove_component_tests, 1);
            std::optional<Cap> raw_representative;
            auto cap_image = [&]() -> const Cap & {
                if (!raw_representative)
                    raw_representative = alcove_cap(alcove_apply(tables.representatives[component.coset], point));
                return *raw_representative;
            };
            bool open;
            if (full_root_checks || !component.open) {
                identifier_count_add(result.alcove_openness_tests, 1);
                open = true;
                for (const auto &functional : functionals) {
                    identifier_count_add(result.alcove_root_checks, 1);
                    Point value{};
                    const auto &representative = cap_image();
                    for (int j = 0; j < 4; ++j)
                        value = value + multiply(functional[j], representative[j]);
                    if (value.x % p == 0 && value.y % p == 0) {
                        open = false;
                        break;
                    }
                }
                // The projected affine stabilizer, and hence presence of a
                // root triflection, is constant on each relative-open face.
                // The independent control recomputes every actual root test.
                if (!full_root_checks)
                    component.open = open;
            } else
                open = *component.open;
            record.components.push_back({component.coset, component.size, orbit_size, open});
            if (!open)
                continue;
            identifier_count_add(result.open_states, orbit_size);
            Cap representative = cap_image();
            for (auto &coordinate : representative) {
                coordinate.x -= p * Rational(2 * coordinate.x + p, 2 * p).floor();
                coordinate.y -= p * Rational(2 * coordinate.y + p, 2 * p).floor();
            }
            representative = short_cap(representative, d);
            I fiber = p;
            if (std::gcd(U(p), stabilizer) > 1)
                fiber = inventory_phase(d, p, I(stabilizer), phase_ring, phase_ring.reduce(representative),
                                        result.phase_chain_stats);
            result.representatives.push_back(representative);
            result.fibers.push_back(fiber);
        }
        result.alcove_records.push_back(std::move(record));
    };
    visit(0, p);
    need(covered_states == result.states, "alcove double-coset mass does not equal the complete residue lattice");
    return result;
}
void print_alcove_tables(bool summary = false) {
    const auto &tables = alcove_tables();
    if (summary) {
        std::cout << "{\"all_states\":4480,\"edge_checks\":" << tables.edge_checks
                  << ",\"table_bytes\":" << tables.conjugates.size() * sizeof(AlcoveMatrix) +
                      tables.representatives.size() * sizeof(AlcoveMatrix) + tables.transitions.size() * sizeof(tables.transitions[0])
                  << ",\"complete\":true}\n";
        return;
    }
    auto matrix = [](const AlcoveMatrix &m) {
        std::cout << '[';
        for (int i = 0; i < 8; ++i) {
            if (i)
                std::cout << ',';
            std::cout << '[';
            for (int j = 0; j < 8; ++j) {
                if (j)
                    std::cout << ',';
                std::cout << m[8 * i + j];
            }
            std::cout << ']';
        }
        std::cout << ']';
    };
    std::cout << "{\"all_states\":4480,\"dimension\":8,\"simple_roots_are\":\"rows of four Eisenstein pairs\",\"simple_roots\":";
    json_value(std::cout, tables.simple_roots);
    std::cout << ",\"cartan\":";
    matrix(tables.cartan);
    std::cout << ",\"inverse_cartan\":";
    matrix(tables.inverse_cartan);
    std::cout << ",\"complex_structure\":";
    matrix(tables.complex_structure);
    std::cout << ",\"highest_marks\":[";
    for (int i = 0; i < 8; ++i) {
        if (i)
            std::cout << ',';
        std::cout << tables.highest_marks[i];
    }
    std::cout << "],\"fundamental_weights\":";
    json_value(std::cout, tables.fundamental_weights);
    std::cout << ",\"edge_checks\":" << tables.edge_checks << ",\"transitions\":[";
    for (std::size_t k = 0; k < tables.transitions.size(); ++k) {
        if (k)
            std::cout << ',';
        std::cout << '[';
        for (int j = 0; j < 9; ++j) {
            if (j)
                std::cout << ',';
            std::cout << tables.transitions[k][j];
        }
        std::cout << ']';
    }
    std::cout << "],\"representatives\":[";
    for (std::size_t k = 0; k < tables.representatives.size(); ++k) {
        if (k)
            std::cout << ',';
        matrix(tables.representatives[k]);
    }
    std::cout << "],\"complete\":true}\n";
}
CuspCaps cusp_caps(Point d, InventoryBackend backend = InventoryBackend::Auto, bool full_alcove_root_checks = false) {
    if (backend == InventoryBackend::Auto) {
        // Only literal rational moduli are Weyl-stable here. A nonreal
        // Eisenstein divisor with square norm is not a rational modulus.
        // The measured small-modulus startup crossover is between5and7;
        // retain the projective path below7 and preserve it as a control.
        if (d.y == 0 && d.x >= 7)
            return alcove_cusp_caps(d, full_alcove_root_checks);
        return inventory_field(d) ? projective_cusp_caps(d) : composite_cusp_caps(d);
    }
    if (backend == InventoryBackend::Alcove) {
        if (d.y == 0)
            return alcove_cusp_caps(d, full_alcove_root_checks);
        return inventory_field(d) ? projective_cusp_caps(d) : composite_cusp_caps(d);
    }
    if (backend == InventoryBackend::Projective)
        return inventory_field(d) ? projective_cusp_caps(d) : composite_cusp_caps(d);
    if (backend == InventoryBackend::Field && inventory_field(d))
        return projective_cusp_caps(d);
    return affine_cusp_caps(d, backend == InventoryBackend::Dense);
}
std::array<Point, 11> cusp_star(Point u, Point beta, const Cap &a, const Cap &b) {
    constexpr Cap e{{{0, -1}, {1, -2}, {2, -2}, {2, -1}}};
    constexpr Cap f{{{1, 1}, {0, 2}, {-1, 2}, {-1, 1}}};
    std::array<Point, 11> p{};
    for (int j = 0; j < 4; ++j) {
        p[j + 1] = a[j] + multiply(u, e[j]);
        p[j + 7] = b[j] + multiply(beta, f[j]);
    }
    p[5] = u;
    p[6] = beta;
    reduction_bound(p);
    return p;
}
bool mathematical_star_rejection(std::string_view message) {
    return message == "collision period has a zero slit" ||
           message == "two-bubble escape failed: mirror-free input is required" ||
           message == "three-bubble escape failed: mirror-free input is required" ||
           message == "repeated star boundary vertex" ||
           message == "overlapping adjacent star sides" || message == "star boundary is not simple";
}
} // namespace gf
