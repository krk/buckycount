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

// Private experimental all-modulus cap inventory over the native library.
// Closed-face ownership uses a measured canonical-key dictionary, not a
// dictionary-free boundary rule. Projection constants come from the exact
// cap_polytope.json certificate; certify_cap_polytope_native.py checks them.
#pragma once
#include "genfullerene.hpp"

namespace gf {
struct CapPolytopeStats {
    std::array<U, 4> prefixes{}, integer_rows{};
    U domain_points = 0, boundary_points = 0, boundary_duplicates = 0;
    U canonical_calls = 0, root_checks = 0, dictionary_entries = 0;
    U closed_before_normalization = 0;
    U interior_direct = 0;
    IdentifierChainStats stabilizer_stats{};
};

using CapPolytopeRow = std::array<I, 9>;
// BEGIN CERTIFIED CAP POLYTOPE PROJECTIONS
const std::array<std::vector<CapPolytopeRow>, 4> &cap_polytope_projection_rows() {
    static const std::array<std::vector<CapPolytopeRow>, 4> result{{
        {
            {0, 1, -1, 0, 0, 0, 0, 0, 0},
            {0, 1, 2, 0, 0, 0, 0, 0, 0},
            {1, -1, -1, 0, 0, 0, 0, 0, 0},
            {1, -1, 0, 0, 0, 0, 0, 0, 0},
        },
        {
            {0, -4, -2, 3, 3, 0, 0, 0, 0},
            {0, -2, 2, 3, 0, 0, 0, 0, 0},
            {0, 1, -1, -1, 1, 0, 0, 0, 0},
            {0, 1, -1, 0, 0, 0, 0, 0, 0},
            {0, 1, 2, 0, 0, 0, 0, 0, 0},
            {0, 3, 0, -2, -1, 0, 0, 0, 0},
            {0, 3, 3, -1, -2, 0, 0, 0, 0},
            {1, -1, -1, 0, 0, 0, 0, 0, 0},
            {1, -1, 0, 0, 0, 0, 0, 0, 0},
            {1, 0, -1, -1, 0, 0, 0, 0, 0},
            {1, 0, 1, 0, -1, 0, 0, 0, 0},
        },
        {
            {0, -2, -1, 2, 1, -1, 0, 0, 0},
            {0, -2, -1, 3, 3, -1, -2, 0, 0},
            {0, -1, 1, 3, 0, -2, -1, 0, 0},
            {0, 0, 0, -2, -1, 2, 2, 0, 0},
            {0, 0, 0, -1, 1, 2, 0, 0, 0},
            {0, 1, -1, -1, 1, 0, 0, 0, 0},
            {0, 1, -1, 0, 0, 0, 0, 0, 0},
            {0, 1, 2, 0, 0, 0, 0, 0, 0},
            {0, 2, 1, -1, -2, -1, 1, 0, 0},
            {0, 2, 1, 0, 0, -1, -1, 0, 0},
            {0, 3, 0, -2, -1, 0, 0, 0, 0},
            {0, 3, 3, -1, -2, 0, 0, 0, 0},
            {1, -1, -1, 0, 0, 0, 0, 0, 0},
            {1, -1, 0, 0, 0, 0, 0, 0, 0},
            {1, 0, -1, -1, 0, 0, 0, 0, 0},
            {1, 0, 0, 0, -1, -1, 0, 0, 0},
            {1, 0, 1, 0, -1, 0, 0, 0, 0},
            {3, 0, 0, 0, 0, -1, -2, 0, 0},
        },
        {
            {0, -2, -1, 2, 1, 0, 0, -2, -1},
            {0, -2, -1, 3, 3, -1, -2, 0, 0},
            {0, -1, 1, 3, 0, -2, -1, 0, 0},
            {0, 0, 0, -2, -1, 2, 1, -1, 1},
            {0, 0, 0, -2, -1, 3, 3, -1, -2},
            {0, 0, 0, -1, 1, 3, 0, -2, -1},
            {0, 0, 0, 0, 0, -2, -1, 3, 3},
            {0, 0, 0, 0, 0, -1, 1, 3, 0},
            {0, 1, -1, -1, 1, 0, 0, 0, 0},
            {0, 1, -1, 0, 0, 0, 0, 0, 0},
            {0, 1, 2, 0, 0, 0, 0, 0, 0},
            {0, 2, 1, -1, -2, -1, 1, 0, 0},
            {0, 2, 1, 0, 0, -1, -2, -1, 1},
            {0, 2, 1, 0, 0, 0, 0, -1, -2},
            {0, 3, 0, -2, -1, 0, 0, 0, 0},
            {0, 3, 3, -1, -2, 0, 0, 0, 0},
            {1, -1, -1, 0, 0, 0, 0, 0, 0},
            {1, -1, 0, 0, 0, 0, 0, 0, 0},
            {1, 0, -1, -1, 0, 0, 0, 0, 0},
            {1, 0, 0, 0, -1, -1, 0, 0, 0},
            {1, 0, 0, 0, 0, 0, -1, -1, 0},
            {1, 0, 0, 0, 0, 0, 0, 0, -1},
            {1, 0, 1, 0, -1, 0, 0, 0, 0},
        },
    }};
    return result;
}
// END CERTIFIED CAP POLYTOPE PROJECTIONS

struct CapPolytopePolygonRow { Wide a, b, c; };

template <class Callback>
void cap_polytope_points(Point d, Callback emit, CapPolytopeStats &stats) {
    Wide D = identifier_norm(identifier_point(d));
    need(D > 0 && D <= 65535, "polytope cap norm must be in 1..65535");
    std::array<std::vector<std::array<Wide, 9>>, 4> facets;
    for (int level = 0; level < 4; ++level)
        for (const auto &row : cap_polytope_projection_rows()[level]) {
            std::array<Wide, 9> image{};
            image[0] = identifier_mul(D, row[0]);
            for (int j = 0; j <= level; ++j) {
                image[1 + 2 * j] = identifier_sub(identifier_mul(identifier_add(d.x, d.y), row[1 + 2 * j]),
                                                   identifier_mul(d.y, row[2 + 2 * j]));
                image[2 + 2 * j] = identifier_add(identifier_mul(d.y, row[1 + 2 * j]),
                                                   identifier_mul(d.x, row[2 + 2 * j]));
            }
            facets[level].push_back(image);
        }
    std::array<Wide, 8> point{};
    std::function<void(int)> visit = [&](int level) {
        if (level == 4) {
            Cap cap;
            for (int j = 0; j < 4; ++j)
                cap[j] = {identifier_narrow(point[2 * j]), identifier_narrow(point[2 * j + 1])};
            bool interior = true;
            for (const auto &row : facets[3]) {
                Wide value = row[0];
                for (int j = 0; j < 8; ++j)
                    value = identifier_add(value, identifier_mul(row[j + 1], point[j]));
                need(value >= 0, "polytope iterator emitted outside its exact domain");
                interior = interior && value > 0;
            }
            identifier_count_add(stats.domain_points, 1);
            if (!interior)
                identifier_count_add(stats.boundary_points, 1);
            emit(cap, interior);
            return;
        }
        std::vector<CapPolytopePolygonRow> polygon;
        std::optional<Wide> xmin, xmax;
        bool possible = true;
        auto constrain_x = [&](Wide a, Wide c) {
            if (a == 0) {
                possible = possible && c >= 0;
            } else if (a > 0) {
                Wide value = identifier_floor(c, a);
                xmax = xmax ? std::min(*xmax, value) : value;
            } else {
                Wide value = identifier_ceil(identifier_sub(0, c), identifier_sub(0, a));
                xmin = xmin ? std::max(*xmin, value) : value;
            }
        };
        for (const auto &row : facets[level]) {
            Wide c = row[0];
            for (int j = 0; j < 2 * level; ++j)
                c = identifier_add(c, identifier_mul(row[j + 1], point[j]));
            CapPolytopePolygonRow p{identifier_sub(0, row[2 * level + 1]),
                                      identifier_sub(0, row[2 * level + 2]), c};
            polygon.push_back(p);
            if (p.b == 0)
                constrain_x(p.a, p.c);
        }
        for (const auto &lower : polygon)
            if (lower.b < 0)
                for (const auto &upper : polygon)
                    if (upper.b > 0)
                        constrain_x(identifier_sub(identifier_mul(upper.b, lower.a), identifier_mul(lower.b, upper.a)),
                                    identifier_sub(identifier_mul(upper.b, lower.c), identifier_mul(lower.b, upper.c)));
        need(xmin && xmax, "unbounded polytope coordinate projection");
        if (!possible || *xmin > *xmax)
            return;
        for (Wide x = *xmin; x <= *xmax; ++x) {
            identifier_count_add(stats.integer_rows[level], 1);
            std::optional<Wide> ymin, ymax;
            for (const auto &row : polygon) {
                Wide c = identifier_sub(row.c, identifier_mul(row.a, x));
                if (row.b > 0) {
                    Wide value = identifier_floor(c, row.b);
                    ymax = ymax ? std::min(*ymax, value) : value;
                } else if (row.b < 0) {
                    Wide value = identifier_ceil(identifier_sub(0, c), identifier_sub(0, row.b));
                    ymin = ymin ? std::max(*ymin, value) : value;
                } else
                    need(c >= 0, "polytope row elimination lost a zero-coefficient bound");
            }
            need(ymin && ymax, "unbounded polytope row fiber");
            for (Wide y = *ymin; y <= *ymax; ++y) {
                point[2 * level] = x;
                point[2 * level + 1] = y;
                identifier_count_add(stats.prefixes[level], 1);
                visit(level + 1);
            }
        }
    };
    visit(0);
}

U cap_polytope_stabilizer_order(const IdentifierRing &ring, const Cap &target,
                               const Cap &current, IdentifierChainStats &stats, int level = 0) {
    constexpr std::array<U, 5> suffix_orders{155520, 648, 24, 3, 1};
    identifier_count_add(stats.nodes, 1);
    if (level == 4)
        return current == target ? 1 : 0;
    bool fixed = true;
    for (int j = level; j < 4; ++j) {
        identifier_count_add(stats.generator_steps, 1);
        if (identifier_cap_step(ring, current, j).cap != current) {
            fixed = false;
            break;
        }
    }
    if (fixed) {
        identifier_count_add(stats.singleton_shortcuts, 1);
        return current == target ? suffix_orders[level] : 0;
    }
    std::map<Cap, U> children;
    for (const auto &step : identifier_cap_chain()[level]) {
        IdentifierPoint value{};
        for (int j = 0; j < 4; ++j)
            value = identifier_add(value, identifier_mul(identifier_point(step.row[j]), identifier_point(current[j])));
        identifier_count_add(stats.row_tests, 1);
        if (ring.reduce(value) != target[level])
            continue;
        Cap child = current;
        for (auto j = step.word.rbegin(); j != step.word.rend(); ++j) {
            child = identifier_cap_step(ring, child, *j).cap;
            identifier_count_add(stats.generator_steps, 1);
        }
        need(child[level] == target[level], "stabilizer transversal convention mismatch");
        for (int j = 0; j < level; ++j)
            need(child[j] == target[j], "stabilizer changed an established prefix");
        auto [position, inserted] = children.emplace(child, 0);
        identifier_count_add(position->second, 1);
        if (!inserted)
            identifier_count_add(stats.duplicate_children, 1);
    }
    U result = 0;
    for (const auto &[child, multiplicity] : children) {
        U count = cap_polytope_stabilizer_order(ring, target, child, stats, level + 1);
        need(count <= suffix_orders[level] / multiplicity, "stabilizer transport count overflow");
        identifier_count_add(result, count * multiplicity);
    }
    need(result <= suffix_orders[level], "stabilizer transport count exceeds group order");
    return result;
}

CuspCaps polytope_cusp_caps(Point d, CapPolytopeStats *reported = nullptr,
                           bool audit_total_mass = false) {
    CapPolytopeStats local;
    auto &stats = reported ? *reported : local;
    Wide D = identifier_norm(identifier_point(d));
    need(D > 0 && D <= 65535, "polytope inventory norm must be in 1..65535");
    IdentifierRing ring(identifier_narrow(identifier_mul(identifier_point(d), {0, 1})));
    CuspCaps result;
    result.states = U(ring.D) * ring.D * ring.D * ring.D;
    std::set<Cap> seen;
    U recovered_states = 0;
    const auto &functionals = cap_root_functionals();
    cap_polytope_points(d, [&](const Cap &raw, bool interior) {
        identifier_count_add(stats.canonical_calls, 1);
        auto image = identifier_chain_image(ring, ring.reduce(raw), result.phase_chain_stats);
        if (!seen.insert(image.cap).second) {
            need(!interior, "interior polytope orbit occurred twice");
            identifier_count_add(stats.boundary_duplicates, 1);
            return;
        }
        identifier_count_add(result.orbits, 1);
        bool open = true;
        if (!interior)
            for (const auto &row : functionals) {
                identifier_count_add(stats.root_checks, 1);
                IdentifierPoint value{};
                for (int j = 0; j < 4; ++j)
                    value = identifier_add(value, identifier_mul(identifier_point(row[j]), identifier_point(raw[j])));
                if (ring.reduce(value) == Point{}) {
                    open = false;
                    break;
                }
            }
        // The module has D^4 elements independently of its orbit inventory.
        // Closed-orbit stabilizer counts are needed only by the optional
        // complete-mass audit, not by any returned mathematical field.
        if (!open && !audit_total_mass)
            return;
        U stabilizer = interior ? 1 : cap_polytope_stabilizer_order(ring, image.cap, image.cap,
                                                                  stats.stabilizer_stats);
        need(stabilizer && 155520 % stabilizer == 0, "invalid cap stabilizer order");
        U orbit_size = 155520 / stabilizer;
        if (audit_total_mass)
            identifier_count_add(recovered_states, orbit_size);
        if (!open)
            return;
        if (interior)
            need(image.subgroup == ring.h22, "interior cap has a nontrivial stabilizer phase");
        identifier_count_add(result.open_states, orbit_size);
        // raw is already nearest by construction. Its canonical key and exact
        // phase subgroup are unchanged by replacing it with a legacy lift.
        result.representatives.push_back(raw);
        result.fibers.push_back(image.subgroup);
    }, stats);
    if (audit_total_mass)
        need(recovered_states == result.states, "polytope orbit masses do not cover the residue module");
    stats.dictionary_entries = seen.size();
    // Payload lower bound only: std::set node/allocation overhead is excluded.
    result.visited_bytes = seen.size() * sizeof(Cap);
    return result;
}

// Census-only interface: representatives, fibers, states, phase-chain work and
// dictionary cap payload are exact. The unused all-orbit count (orbits) and
// open-state mass (open_states) are deliberately not computed and remain zero.
// Never use this function to service --cap-inventory or print full metadata.
// KEEP(raw, D) may discard a domain point before any orbit normalization.
template <class Keep>
CuspCaps polytope_census_caps_where(Point d, CapPolytopeStats *reported, Keep keep) {
    CapPolytopeStats local;
    auto &stats = reported ? *reported : local;
    Wide D = identifier_norm(identifier_point(d));
    need(D > 0 && D <= 65535, "polytope census norm must be in 1..65535");
    IdentifierRing ring(identifier_narrow(identifier_mul(identifier_point(d), {0, 1})));
    CuspCaps result;
    result.states = U(ring.D) * ring.D * ring.D * ring.D;
    std::set<Cap> seen;
    const auto &functionals = cap_root_functionals();
    cap_polytope_points(d, [&](const Cap &raw, bool interior) {
        if (!keep(raw, D))
            return;
        // A strict interior has one affine-orbit occurrence and trivial
        // stabilizer. No boundary key or modular normalizer is needed.
        if (interior) {
            identifier_count_add(stats.interior_direct, 1);
            result.representatives.push_back(raw);
            result.fibers.push_back(ring.h22);
            return;
        }
        if (!interior)
            for (const auto &row : functionals) {
                identifier_count_add(stats.root_checks, 1);
                IdentifierPoint value{};
                for (int j = 0; j < 4; ++j)
                    value = identifier_add(value, identifier_mul(identifier_point(row[j]), identifier_point(raw[j])));
                if (ring.reduce(value) == Point{}) {
                    identifier_count_add(stats.closed_before_normalization, 1);
                    return;
                }
            }
        identifier_count_add(stats.canonical_calls, 1);
        auto image = identifier_chain_image(ring, ring.reduce(raw), result.phase_chain_stats);
        if (!seen.insert(image.cap).second) {
            need(!interior, "interior polytope census orbit occurred twice");
            identifier_count_add(stats.boundary_duplicates, 1);
            return;
        }
        result.representatives.push_back(raw);
        result.fibers.push_back(image.subgroup);
    }, stats);
    stats.dictionary_entries = seen.size();
    result.visited_bytes = seen.size() * sizeof(Cap);
    return result;
}
CuspCaps polytope_census_caps(Point d, CapPolytopeStats *reported = nullptr) {
    return polytope_census_caps_where(d, reported, [](const Cap &, Wide) { return true; });
}
} // namespace gf
