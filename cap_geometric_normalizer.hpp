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

// Private mixed cap normal form: strict real-P owners, exact legacy fallback
// on every boundary. This changes the canonical identifier order.
#pragma once
#include "genfullerene.hpp"
#include "cap_polytope_native.hpp"

namespace gf {
struct GeometricNormalizerStats {
    U cap_queries = 0, strict_queries = 0, boundary_fallbacks = 0;
    U facet_steps = 0, maximum_steps = 0, facet_checks = 0, simple_generator_steps = 0;
    U nearest_nodes = 0, nearest_ties = 0;
    void json(std::ostream &out) const {
        out << "{\"geometric_normalizer\":true,\"cap_queries\":" << cap_queries
            << ",\"strict_queries\":" << strict_queries << ",\"boundary_fallbacks\":" << boundary_fallbacks
            << ",\"facet_steps\":" << facet_steps << ",\"maximum_steps\":" << maximum_steps
            << ",\"facet_checks\":" << facet_checks << ",\"simple_generator_steps\":" << simple_generator_steps
            << ",\"nearest_nodes\":" << nearest_nodes << ",\"nearest_ties\":" << nearest_ties << "}\n";
    }
};
struct GeometricFacet {
    std::array<I, 8> normal;
    std::vector<unsigned char> word;
};
// BEGIN CERTIFIED GEOMETRIC FACET WORDS
const std::array<GeometricFacet, 16> &geometric_facet_words() {
    static const std::array<GeometricFacet, 16> result{{
        {{3, 3, -1, -2, 0, 0, 0, 0}, {0}},
        {{3, 0, -2, -1, 0, 0, 0, 0}, {0, 0}},
        {{2, 1, 0, 0, 0, 0, -1, -2}, {0, 1, 2, 3}},
        {{2, 1, 0, 0, -1, -2, -1, 1}, {0, 1, 2, 3, 3}},
        {{2, 1, -1, -2, -1, 1, 0, 0}, {0, 1, 2, 3, 3, 2}},
        {{1, -1, -1, 1, 0, 0, 0, 0}, {0, 1, 2, 3, 3, 2, 1}},
        {{1, 2, 0, 0, 0, 0, 0, 0}, {0, 1, 2, 3, 3, 2, 1, 0}},
        {{1, -1, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 1, 2, 3, 3, 2, 1, 1, 0, 0}},
        {{-2, -1, 3, 3, -1, -2, 0, 0}, {1}},
        {{-1, 1, 3, 0, -2, -1, 0, 0}, {1, 1}},
        {{-2, -1, 2, 1, 0, 0, -2, -1}, {1, 2, 3}},
        {{0, 0, -2, -1, 3, 3, -1, -2}, {2}},
        {{0, 0, -2, -1, 2, 1, -1, 1}, {2, 3}},
        {{0, 0, -1, 1, 3, 0, -2, -1}, {2, 2}},
        {{0, 0, 0, 0, -2, -1, 3, 3}, {3}},
        {{0, 0, 0, 0, -1, 1, 3, 0}, {3, 3}},
    }};
    return result;
}
// END CERTIFIED GEOMETRIC FACET WORDS

struct GeometricPreparedCap {
    PreparedIdentifierCap value;
    bool geometric = false;
};
struct GeometricCapKey {
    bool geometric = false;
    Cap coordinates;
    auto operator<=>(const GeometricCapKey &) const = default;
};
struct GeometricParabolic {
    IdentifierCusp record;
    std::array<GeometricCapKey, 2> residue_key;
    I central_subgroup = 1;
    bool antiunitary = false;
};

IdentifierPoint geometric_translation_delta(const Cap &tau, const Cap &original, Point beta) {
    auto pairing = identifier_mul(identifier_cap_pair(tau, original), {2, -1});
    need(pairing[0] % 3 == 0 && pairing[1] % 3 == 0, "nonintegral geometric Heisenberg correction");
    return identifier_add(identifier_mul(identifier_point(beta), {identifier_cap_index(tau), 0}),
                          {pairing[0] / 3, pairing[1] / 3});
}

GeometricPreparedCap prepare_geometric_cap(const Cap &cap, const IdentifierRing &ring,
                                           IdentifierParabolicCache &cache,
                                           GeometricNormalizerStats *reported = nullptr) {
    GeometricNormalizerStats local;
    auto &stats = reported ? *reported : local;
    identifier_count_add(stats.cap_queries, 1);
    auto nearest = nearest_identifier_cap(cap, ring.beta);
    identifier_count_add(stats.nearest_nodes, nearest.visited_nodes);
    identifier_count_add(stats.nearest_ties, nearest.ties);
    auto fallback = [&]() {
        identifier_count_add(stats.boundary_fallbacks, 1);
        return GeometricPreparedCap{prepare_identifier_cap(cap, ring, cache), false};
    };
    // More than one nearest translate is exactly a Voronoi-boundary event.
    // G preserves nearest multiplicity, so no fold can reach strict P.
    if (nearest.ties != 1)
        return fallback();
    Cap owner = nearest.representative;
    std::array<std::array<Wide, 8>, 16> normals{};
    const auto &facets = geometric_facet_words();
    for (std::size_t i = 0; i < facets.size(); ++i)
        for (int j = 0; j < 4; ++j) {
            I a = facets[i].normal[2 * j], b = facets[i].normal[2 * j + 1];
            normals[i][2 * j] = identifier_sub(identifier_mul(identifier_add(ring.beta.x, ring.beta.y), a),
                                               identifier_mul(ring.beta.y, b));
            normals[i][2 * j + 1] = identifier_add(identifier_mul(ring.beta.y, a), identifier_mul(ring.beta.x, b));
        }
    U steps = 0;
    while (true) {
        bool moved = false;
        for (std::size_t i = 0; i < facets.size(); ++i) {
            identifier_count_add(stats.facet_checks, 1);
            Wide value = 0;
            for (int j = 0; j < 4; ++j)
                value = identifier_add(value, identifier_add(identifier_mul(normals[i][2 * j], owner[j].x),
                                                              identifier_mul(normals[i][2 * j + 1], owner[j].y)));
            if (value >= 0)
                continue;
            for (auto j = facets[i].word.rbegin(); j != facets[i].word.rend(); ++j) {
                auto coordinate = identifier_mul(identifier_point(owner[*j]), {-1, 1});
                if (*j)
                    coordinate = identifier_add(coordinate, identifier_point(owner[*j - 1]));
                if (*j < 3)
                    coordinate = identifier_add(coordinate, identifier_mul(identifier_point(owner[*j + 1]), {1, -1}));
                owner[*j] = identifier_narrow(coordinate);
                identifier_count_add(stats.simple_generator_steps, 1);
            }
            identifier_count_add(steps, 1);
            need(steps < 155520, "geometric cap descent exceeded its finite-group bound");
            moved = true;
            break;
        }
        if (!moved)
            break;
    }
    identifier_count_add(stats.facet_steps, steps);
    stats.maximum_steps = std::max(stats.maximum_steps, steps);
    bool strict = true;
    auto beta_bar = identifier_conjugate(identifier_point(ring.beta));
    std::array<IdentifierPoint, 4> scaled;
    for (int j = 0; j < 4; ++j)
        scaled[j] = identifier_mul(beta_bar, identifier_point(owner[j]));
    for (const auto &row : cap_polytope_projection_rows()[3]) {
        Wide value = identifier_mul(ring.D, row[0]);
        for (int j = 0; j < 4; ++j)
            value = identifier_add(value, identifier_add(identifier_mul(row[2 * j + 1], scaled[j][0]),
                                                          identifier_mul(row[2 * j + 2], scaled[j][1])));
        need(value >= 0, "geometric nearest-cap fold ended outside P");
        strict = strict && value > 0;
    }
    if (!strict)
        return fallback();
    identifier_count_add(stats.strict_queries, 1);
    auto delta = geometric_translation_delta(nearest.translation, cap, ring.beta);
    return {{owner, owner, {}, identifier_cap_index(owner), ring.phase(delta), ring.h22}, true};
}

GeometricParabolic combine_geometric_caps(Point u, const IdentifierRing &ring, Wide m,
                                           const GeometricPreparedCap &a,
                                           const GeometricPreparedCap &b) {
    const auto &first = a.value;
    const auto &second = b.value;
    I subgroup = std::gcd(first.subgroup, second.subgroup);
    I phase = positive_mod(identifier_narrow(identifier_add(ring.phase(identifier_point(u)),
                                                            identifier_add(first.phase, second.phase))), subgroup);
    Wide R = identifier_add(m, identifier_add(first.residue_index, second.residue_index));
    need(R % ring.h22 == 0, "geometric cap normalization lost scalar solubility");
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
    need(low <= high, "empty geometric central phase interval");
    scalar = identifier_add(scalar, identifier_mul(direction, {choice, 0}));
    IdentifierCusp v{identifier_narrow(scalar), ring.beta, {first.nearest, second.nearest}};
    std::array<GeometricCapKey, 2> keys{{{a.geometric, first.residue}, {b.geometric, second.residue}}};
    if (keys[1] < keys[0]) {
        std::swap(keys[0], keys[1]);
        std::swap(v.caps[0], v.caps[1]);
    }
    need(identifier_cusp_index(v) == m, "geometric parabolic normal form changed the norm");
    return {v, keys, subgroup, false};
}

bool geometric_parabolic_less(const GeometricParabolic &a, const GeometricParabolic &b) {
    return std::tuple(a.record.beta, a.residue_key, a.record.u) <
           std::tuple(b.record.beta, b.residue_key, b.record.u);
}
bool geometric_parabolic_equal(const GeometricParabolic &a, const GeometricParabolic &b) {
    return !geometric_parabolic_less(a, b) && !geometric_parabolic_less(b, a);
}
GeometricParabolic geometric_unitary_parabolic(IdentifierCusp v, IdentifierParabolicCache &cache,
                                               GeometricNormalizerStats *stats = nullptr) {
    v = identifier_height_unit(v);
    Wide m = identifier_cusp_index(v);
    IdentifierRing ring(v.beta);
    auto first = prepare_geometric_cap(v.caps[0], ring, cache, stats);
    auto second = prepare_geometric_cap(v.caps[1], ring, cache, stats);
    return combine_geometric_caps(v.u, ring, m, first, second);
}
GeometricParabolic geometric_parabolic(const IdentifierCusp &v, IdentifierParabolicCache &cache,
                                        GeometricNormalizerStats *stats = nullptr) {
    auto unitary = identifier_height_unit(v);
    auto conjugate = identifier_height_unit(identifier_antiunitary(unitary));
    if (conjugate.beta < unitary.beta) {
        auto second = geometric_unitary_parabolic(conjugate, cache, stats);
        second.antiunitary = true;
        return second;
    }
    auto first = geometric_unitary_parabolic(unitary, cache, stats);
    if (conjugate.beta > first.record.beta)
        return first;
    auto second = geometric_unitary_parabolic(conjugate, cache, stats);
    if (geometric_parabolic_less(second, first)) {
        second.antiunitary = true;
        return second;
    }
    return first;
}
} // namespace gf
