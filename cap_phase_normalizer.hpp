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

// Private nearest-lift normal form over the polytope and geometric
// normalizers. The canonical order differs from both older normalizers;
// comparisons must use this convention on both sides.
#pragma once
#include <mutex>
#include "genfullerene.hpp"
#include "cap_geometric_normalizer.hpp"

namespace gf {
struct PhaseNormalizerStats {
    U cap_queries = 0, nearest_nodes = 0, nearest_lifts = 0, tied_queries = 0;
    U integral_queries = 0, reflection_steps = 0, face_tables = 0, face_edges = 0;
    U nontrivial_phases = 0, inventory_boundary = 0, inventory_duplicates = 0, inventory_norm_rejections = 0;
    void json(std::ostream &out) const {
        out << "{\"phase_normalizer\":true,\"cap_queries\":" << cap_queries
            << ",\"nearest_nodes\":" << nearest_nodes << ",\"nearest_lifts\":" << nearest_lifts
            << ",\"tied_queries\":" << tied_queries << ",\"integral_queries\":" << integral_queries
            << ",\"reflection_steps\":" << reflection_steps << ",\"face_tables\":" << face_tables
            << ",\"face_edges\":" << face_edges << ",\"nontrivial_phases\":" << nontrivial_phases
            << ",\"inventory_boundary\":" << inventory_boundary
            << ",\"inventory_duplicates\":" << inventory_duplicates
            << ",\"inventory_norm_rejections\":" << inventory_norm_rejections << "}\n";
    }
};

struct PhaseNearestLift { Cap cap, translation; };
std::vector<PhaseNearestLift> phase_nearest_lifts(const Cap &z, Point beta,
                                                PhaseNormalizerStats &stats) {
    CapBall ball(beta, z, 1, 0, 1);
    auto lift = [&](const Cap &y) {
        PhaseNearestLift result;
        for (int j = 0; j < 4; ++j) {
            result.cap[j] = identifier_narrow(identifier_sub(identifier_point(z[j]),
                identifier_mul(identifier_point(beta), identifier_point(y[j]))));
            result.translation[j] = identifier_narrow(identifier_sub(IdentifierPoint{}, identifier_point(y[j])));
        }
        return result;
    };
    Wide best = identifier_cap_index(lift(ball.babai()).cap);
    Wide factor = identifier_mul(ball.decomposition.scale, identifier_norm(identifier_point(beta)));
    ball.limit = identifier_mul(factor, best);
    std::vector<PhaseNearestLift> result;
    ball.enumerate([&](const Cap &y) {
        auto candidate = lift(y);
        Wide index = identifier_cap_index(candidate.cap);
        if (index < best) {
            best = index;
            ball.limit = identifier_mul(factor, best);
            result.clear();
        }
        if (index == best)
            result.push_back(candidate);
    });
    need(!result.empty(), "empty nearest-lift phase set");
    identifier_count_add(stats.nearest_nodes, ball.visited_nodes);
    identifier_count_add(stats.nearest_lifts, result.size());
    identifier_count_add(stats.tied_queries, result.size() > 1);
    return result;
}

// Integral G-orbits are W(E8) dominant labels plus a zero-wall double coset.
// The tables are derived in the native basis, not imported from the external
// implementation (whose four Eisenstein coordinates have reverse order).
Cap phase_integral_owner(const Cap &cap, PhaseNormalizerStats &stats) {
    identifier_count_add(stats.integral_queries, 1);
    const auto &tables = alcove_tables();
    AlcoveVector labels{};
    for (int i = 0; i < 8; ++i) {
        auto pairing = identifier_cap_pair(tables.simple_roots[i], cap);
        Wide trace = identifier_add(identifier_mul(2, pairing[0]), pairing[1]);
        need(trace % 3 == 0, "nonintegral real E8 pairing");
        labels[i] = identifier_narrow(trace / 3);
    }
    unsigned state = 0;
    U steps = 0;
    while (true) {
        int bad = -1;
        for (int i = 0; i < 8; ++i)
            if (labels[i] < 0) { bad = i; break; }
        if (bad < 0) break;
        I value = labels[bad];
        for (int i = 0; i < 8; ++i)
            labels[i] = identifier_narrow(identifier_sub(labels[i],
                identifier_mul(tables.cartan[8 * i + bad], value)));
        state = tables.transitions[state][bad];
        ++steps;
        // Every move crosses a distinct separating positive-root wall.
        need(steps <= 120, "E8 dominant fold exceeded the positive-root count");
    }
    identifier_count_add(stats.reflection_steps, steps);
    unsigned mask = 0;
    for (int i = 0; i < 8; ++i)
        if (labels[i] == 0) mask |= 1U << i;
    static std::array<std::vector<std::uint16_t>, 256> face_labels;
    static std::array<std::once_flag, 256> initialized;
    auto &components = face_labels[mask];
    std::call_once(initialized[mask], [&] {
        identifier_count_add(stats.face_tables, 1);
        components.assign(4480, 65535);
        std::vector<unsigned> queue;
        for (unsigned first = 0; first < 4480; ++first) {
            if (components[first] != 65535) continue;
            queue.clear(); queue.push_back(first); components[first] = first;
            for (std::size_t k = 0; k < queue.size(); ++k)
                for (int i = 0; i < 8; ++i)
                    if (mask & (1U << i)) {
                        identifier_count_add(stats.face_edges, 1);
                        unsigned next = tables.transitions[queue[k]][i];
                        if (components[next] == 65535) {
                            components[next] = first;
                            queue.push_back(next);
                        }
                    }
        }
    });
    auto coordinates = alcove_apply(tables.inverse_cartan, labels);
    auto owner = alcove_cap(alcove_apply(tables.representatives[components[state]], coordinates));
    need(identifier_cap_index(owner) == identifier_cap_index(cap), "integral cap owner changed norm");
    return owner;
}

GeometricPreparedCap prepare_phase_cap(const Cap &cap, const IdentifierRing &ring,
                                       IdentifierParabolicCache &cache,
                                       PhaseNormalizerStats *reported = nullptr) {
    (void)cache;
    PhaseNormalizerStats local;
    auto &stats = reported ? *reported : local;
    identifier_count_add(stats.cap_queries, 1);
    auto lifts = phase_nearest_lifts(cap, ring.beta, stats);
    std::optional<Cap> owner;
    I phase = 0, subgroup = ring.h22;
    for (const auto &lift : lifts) {
        auto image = phase_integral_owner(lift.cap, stats);
        I correction = ring.phase(geometric_translation_delta(lift.translation, cap, ring.beta));
        if (!owner || image < *owner) {
            owner = image;
            phase = correction;
            subgroup = ring.h22;
        } else if (image == *owner) {
            subgroup = std::gcd(subgroup, positive_mod(correction - phase, ring.h22));
        }
    }
    need(owner.has_value(), "nearest-lift owner is missing");
    identifier_count_add(stats.nontrivial_phases, subgroup < ring.h22);
    // Every lift in the selected integral orbit gives one actual transport.
    // Their phase differences generate exactly the modular stabilizer image;
    // a linear stabilizer fixing a lift has zero scalar correction.
    return {{*owner, *owner, {}, identifier_cap_index(*owner), positive_mod(phase, subgroup), subgroup}, true};
}

GeometricParabolic phase_unitary_parabolic(IdentifierCusp v, IdentifierParabolicCache &cache,
                                          PhaseNormalizerStats *stats = nullptr) {
    v = identifier_height_unit(v);
    IdentifierRing ring(v.beta);
    auto first = prepare_phase_cap(v.caps[0], ring, cache, stats);
    auto second = prepare_phase_cap(v.caps[1], ring, cache, stats);
    return combine_geometric_caps(v.u, ring, identifier_cusp_index(v), first, second);
}
GeometricParabolic phase_parabolic(const IdentifierCusp &v, IdentifierParabolicCache &cache,
                                   PhaseNormalizerStats *stats = nullptr) {
    auto unitary = identifier_height_unit(v);
    auto conjugate = identifier_height_unit(identifier_antiunitary(unitary));
    if (conjugate.beta < unitary.beta) {
        auto result = phase_unitary_parabolic(conjugate, cache, stats);
        result.antiunitary = true;
        return result;
    }
    auto first = phase_unitary_parabolic(unitary, cache, stats);
    if (conjugate.beta > first.record.beta) return first;
    auto second = phase_unitary_parabolic(conjugate, cache, stats);
    if (geometric_parabolic_less(second, first)) {
        second.antiunitary = true;
        return second;
    }
    return first;
}

// Census metadata only, exactly as polytope_census_caps: all-orbit count and
// open-state mass are intentionally omitted. The enumeration domain is the
// unchanged certified 23-facet polytope, not an enclosing lattice ball.
} // namespace gf
