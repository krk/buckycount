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

// Private half-open integral G ownership inside the existing cap polytope,
// over the nearest-lift phase normalizer. The infinitesimal scalar
// 1+epsilon*zeta commutes with G; all16 cone walls have root-line kernels.
#pragma once
#include "genfullerene.hpp"
#include "cap_phase_normalizer.hpp"

namespace gf {
struct BoundaryOwnerStats {
    CapPolytopeStats polytope;
    PhaseNormalizerStats phase;
    U strict_direct = 0, cone_boundary_direct = 0, cone_rejections = 0;
    U voronoi_boundary = 0, horizontal_rejections = 0, duplicate_rejections = 0;
    U norm_rejections = 0, root_checks = 0;
    U preparation_direct = 0, fold_steps = 0, maximum_fold_steps = 0, singular_fallbacks = 0;
    U prefix_rejections = 0, unary_rejections = 0;
    void json(std::ostream &out) const {
        out << "{\"half_open_boundary_owner\":true,\"strict_direct\":" << strict_direct
            << ",\"cone_boundary_direct\":" << cone_boundary_direct
            << ",\"cone_rejections\":" << cone_rejections
            << ",\"voronoi_boundary\":" << voronoi_boundary
            << ",\"horizontal_rejections\":" << horizontal_rejections
            << ",\"duplicate_rejections\":" << duplicate_rejections
            << ",\"norm_rejections\":" << norm_rejections << ",\"root_checks\":" << root_checks
            << ",\"normalizer_queries\":" << phase.cap_queries
            << ",\"preparation_direct\":" << preparation_direct << ",\"fold_steps\":" << fold_steps
            << ",\"maximum_fold_steps\":" << maximum_fold_steps
            << ",\"singular_fallbacks\":" << singular_fallbacks
            << ",\"prefix_rejections\":" << prefix_rejections
            << ",\"unary_rejections\":" << unary_rejections
            << ",\"domain_points\":" << polytope.domain_points << "}\n";
    }
};

struct BoundaryOwnerPosition {
    bool owner = false;
    bool voronoi_boundary = false;
};
BoundaryOwnerPosition boundary_owner_position(const Cap &raw, Point divisor) {
    auto conjugate_divisor = identifier_conjugate(identifier_point(divisor));
    Wide D = identifier_norm(identifier_point(divisor));
    need(D > 0, "half-open cap owner requires a nonzero divisor");
    std::array<IdentifierPoint, 4> scaled;
    for (int j = 0; j < 4; ++j)
        scaled[j] = identifier_mul(conjugate_divisor, identifier_point(raw[j]));
    bool voronoi_boundary = false;
    for (const auto &row : cap_polytope_projection_rows()[3]) {
        Wide value = identifier_mul(D, row[0]);
        Wide derivative = 0;
        for (int j = 0; j < 4; ++j) {
            value = identifier_add(value, identifier_add(identifier_mul(row[2*j+1], scaled[j][0]),
                                                          identifier_mul(row[2*j+2], scaled[j][1])));
            if (row[0] == 0) {
                // zeta*(a+b*zeta) = -b+(a+b)*zeta.
                derivative = identifier_add(derivative,
                    identifier_add(identifier_mul(-Wide(row[2*j+1]), scaled[j][1]),
                                   identifier_mul(row[2*j+2], identifier_add(scaled[j][0], scaled[j][1]))));
            }
        }
        if (value < 0) return {};
        if (row[0] == 0) {
            if (value == 0 && derivative <= 0) return {};
        } else {
            voronoi_boundary = voronoi_boundary || value == 0;
        }
    }
    return {true, voronoi_boundary};
}

std::optional<Cap> boundary_fold_cap(Cap cap, Point beta, BoundaryOwnerStats &stats) {
    const auto beta_bar = identifier_conjugate(identifier_point(beta));
    U steps = 0;
    while (true) {
        std::array<IdentifierPoint, 4> scaled;
        for (int j = 0; j < 4; ++j)
            scaled[j] = identifier_mul(beta_bar, identifier_point(cap[j]));
        bool moved = false, singular = false;
        for (const auto &facet : geometric_facet_words()) {
            Wide value = 0;
            for (int j = 0; j < 4; ++j)
                value = identifier_add(value, identifier_add(identifier_mul(facet.normal[2*j], scaled[j][0]),
                                                              identifier_mul(facet.normal[2*j+1], scaled[j][1])));
            if (value > 0) continue;
            Wide derivative = 0;
            if (value == 0) {
                for (int j = 0; j < 4; ++j)
                    derivative = identifier_add(derivative,
                        identifier_add(identifier_mul(-Wide(facet.normal[2*j]), scaled[j][1]),
                            identifier_mul(facet.normal[2*j+1], identifier_add(scaled[j][0], scaled[j][1]))));
                if (derivative == 0) singular = true;
                if (derivative >= 0) continue;
            }
            for (auto j = facet.word.rbegin(); j != facet.word.rend(); ++j) {
                auto coordinate = identifier_mul(identifier_point(cap[*j]), {-1,1});
                if (*j) coordinate = identifier_add(coordinate, identifier_point(cap[*j-1]));
                if (*j < 3) coordinate = identifier_add(coordinate,
                    identifier_mul(identifier_point(cap[*j+1]), {1,-1}));
                cap[*j] = identifier_narrow(coordinate);
            }
            ++steps;
            need(steps < 155520, "half-open cone fold exceeded the finite G bound");
            moved = true;
            break;
        }
        if (!moved) {
            identifier_count_add(stats.fold_steps, steps);
            stats.maximum_fold_steps = std::max(stats.maximum_fold_steps, steps);
            if (singular) return std::nullopt;
            return cap;
        }
    }
}

GeometricPreparedCap prepare_boundary_inventory_owner(const Cap &owner, const IdentifierRing &ring,
                                                       I subgroup) {
    need(subgroup > 0 && ring.h22 % subgroup == 0, "invalid known cap-owner phase subgroup");
    return {{owner, owner, {}, identifier_cap_index(owner), 0, subgroup}, true};
}
GeometricPreparedCap prepare_boundary_cap(const Cap &cap, const IdentifierRing &ring,
                                          IdentifierParabolicCache &cache,
                                          BoundaryOwnerStats *reported = nullptr) {
    BoundaryOwnerStats local;
    auto &stats = reported ? *reported : local;
    auto position = boundary_owner_position(cap, ring.beta);
    if (position.owner && !position.voronoi_boundary) {
        identifier_count_add(stats.preparation_direct, 1);
        return prepare_boundary_inventory_owner(cap, ring, ring.h22);
    }
    identifier_count_add(stats.phase.cap_queries, 1);
    auto lifts = phase_nearest_lifts(cap, ring.beta, stats.phase);
    std::optional<Cap> owner;
    I phase = 0, subgroup = ring.h22;
    for (const auto &lift : lifts) {
        auto folded = boundary_fold_cap(lift.cap, ring.beta, stats);
        if (!folded) {
            // Presence of a singular nearest lift is invariant under the
            // affine cap action. Tag this disjoint stratum consistently.
            identifier_count_add(stats.singular_fallbacks, 1);
            auto result = prepare_phase_cap(cap, ring, cache, &stats.phase);
            result.geometric = false;
            return result;
        }
        I correction = ring.phase(geometric_translation_delta(lift.translation, cap, ring.beta));
        if (!owner || *folded < *owner) {
            owner = *folded;
            phase = correction;
            subgroup = ring.h22;
        } else if (*folded == *owner) {
            subgroup = std::gcd(subgroup, positive_mod(correction - phase, ring.h22));
        }
    }
    need(owner.has_value(), "half-open nearest owner is missing");
    return {{*owner, *owner, {}, identifier_cap_index(*owner), positive_mod(phase, subgroup), subgroup}, true};
}
GeometricParabolic boundary_unitary_parabolic(IdentifierCusp v, IdentifierParabolicCache &cache,
                                             BoundaryOwnerStats *stats = nullptr) {
    v = identifier_height_unit(v);
    IdentifierRing ring(v.beta);
    auto first = prepare_boundary_cap(v.caps[0], ring, cache, stats);
    auto second = prepare_boundary_cap(v.caps[1], ring, cache, stats);
    return combine_geometric_caps(v.u, ring, identifier_cusp_index(v), first, second);
}
GeometricParabolic boundary_parabolic(const IdentifierCusp &v, IdentifierParabolicCache &cache,
                                      BoundaryOwnerStats *stats = nullptr) {
    auto unitary = identifier_height_unit(v);
    auto conjugate = identifier_height_unit(identifier_antiunitary(unitary));
    if (conjugate.beta < unitary.beta) {
        auto result = boundary_unitary_parabolic(conjugate, cache, stats);
        result.antiunitary = true;
        return result;
    }
    auto first = boundary_unitary_parabolic(unitary, cache, stats);
    if (conjugate.beta > first.record.beta) return first;
    auto second = boundary_unitary_parabolic(conjugate, cache, stats);
    if (geometric_parabolic_less(second, first)) {
        second.antiunitary = true;
        return second;
    }
    return first;
}

// Each homogeneous row of a projection of P is valid on the entire cone C:
// otherwise scale a violating C point into the neighborhood of0 inside V.
// Thus every C+ owner also satisfies the half-open sign rule at each prefix.
// This rejects complete subtrees without modifying the certified projections.
template <class Callback>
void boundary_polytope_points(Point d, Callback emit, BoundaryOwnerStats &reported) {
    auto &stats = reported.polytope;
    Wide D = identifier_norm(identifier_point(d));
    need(D > 0, "half-open polytope requires a nonzero divisor");
    std::array<std::vector<std::array<Wide, 9>>, 4> facets;
    for (int level = 0; level < 4; ++level)
        for (const auto &row : cap_polytope_projection_rows()[level]) {
            std::array<Wide, 9> image{};
            image[0] = identifier_mul(D, row[0]);
            for (int j = 0; j <= level; ++j) {
                image[1+2*j] = identifier_sub(identifier_mul(identifier_add(d.x,d.y), row[1+2*j]),
                                              identifier_mul(d.y, row[2+2*j]));
                image[2+2*j] = identifier_add(identifier_mul(d.y, row[1+2*j]),
                                              identifier_mul(d.x, row[2+2*j]));
            }
            facets[level].push_back(image);
        }
    std::array<Wide, 8> point{};
    std::function<void(int)> visit = [&](int level) {
        if (level == 4) {
            Cap cap;
            for (int j = 0; j < 4; ++j)
                cap[j] = {identifier_narrow(point[2*j]), identifier_narrow(point[2*j+1])};
            bool interior = true;
            for (const auto &row : facets[3]) {
                Wide value = row[0];
                for (int j = 0; j < 8; ++j)
                    value = identifier_add(value, identifier_mul(row[j+1], point[j]));
                need(value >= 0, "half-open iterator emitted outside P");
                interior = interior && value > 0;
            }
            identifier_count_add(stats.domain_points, 1);
            if (!interior) identifier_count_add(stats.boundary_points, 1);
            emit(cap, interior);
            return;
        }
        std::vector<CapPolytopePolygonRow> polygon;
        std::optional<Wide> xmin, xmax;
        bool possible = true;
        auto constrain_x = [&](Wide a, Wide c) {
            if (a == 0) possible = possible && c >= 0;
            else if (a > 0) {
                Wide value = identifier_floor(c,a);
                xmax = xmax ? std::min(*xmax,value) : value;
            } else {
                Wide value = identifier_ceil(identifier_sub(0,c),identifier_sub(0,a));
                xmin = xmin ? std::max(*xmin,value) : value;
            }
        };
        for (const auto &row : facets[level]) {
            Wide c = row[0];
            for (int j = 0; j < 2*level; ++j)
                c = identifier_add(c,identifier_mul(row[j+1],point[j]));
            CapPolytopePolygonRow p{identifier_sub(0,row[2*level+1]),
                                    identifier_sub(0,row[2*level+2]),c};
            polygon.push_back(p);
            if (p.b == 0) constrain_x(p.a,p.c);
        }
        for (const auto &lower : polygon)
            if (lower.b < 0)
                for (const auto &upper : polygon)
                    if (upper.b > 0)
                        constrain_x(identifier_sub(identifier_mul(upper.b,lower.a),identifier_mul(lower.b,upper.a)),
                                    identifier_sub(identifier_mul(upper.b,lower.c),identifier_mul(lower.b,upper.c)));
        need(xmin && xmax, "unbounded half-open coordinate projection");
        if (!possible || *xmin > *xmax) return;
        for (Wide x = *xmin; x <= *xmax; ++x) {
            identifier_count_add(stats.integer_rows[level],1);
            std::optional<Wide> ymin,ymax;
            for (const auto &row : polygon) {
                Wide c = identifier_sub(row.c,identifier_mul(row.a,x));
                if (row.b > 0) {
                    Wide value = identifier_floor(c,row.b);
                    ymax = ymax ? std::min(*ymax,value) : value;
                } else if (row.b < 0) {
                    Wide value = identifier_ceil(identifier_sub(0,c),identifier_sub(0,row.b));
                    ymin = ymin ? std::max(*ymin,value) : value;
                } else need(c >= 0,"half-open row elimination lost a bound");
            }
            need(ymin && ymax,"unbounded half-open row fiber");
            for (Wide y = *ymin; y <= *ymax; ++y) {
                point[2*level]=x; point[2*level+1]=y;
                identifier_count_add(stats.prefixes[level],1);
                bool keep = true;
                for (const auto &row : facets[level]) {
                    if (row[0] != 0) continue;
                    Wide value = 0;
                    for (int j = 0; j <= level; ++j)
                        value = identifier_add(value,identifier_add(identifier_mul(row[2*j+1],point[2*j]),
                                                                    identifier_mul(row[2*j+2],point[2*j+1])));
                    need(value >= 0,"half-open prefix violates a closed projection");
                    if (value != 0) continue;
                    Wide derivative = 0;
                    for (int j = 0; j <= level; ++j)
                        derivative = identifier_add(derivative,
                            identifier_add(identifier_mul(row[2*j+2],point[2*j]),
                                identifier_mul(identifier_sub(row[2*j+2],row[2*j+1]),point[2*j+1])));
                    if (derivative <= 0) { keep = false; break; }
                }
                if (keep) visit(level+1);
                else identifier_count_add(reported.prefix_rejections,1);
            }
        }
    };
    visit(0);
}

// Walk geometric owners directly. The residue space has D^4 elements, but
// the iterator only needs D, the eight coordinates, and the stabilizer step.
template <class Callback>
void visit_boundary_caps(const IdentifierRing &ring, Callback emit, BoundaryOwnerStats &stats,
                         Wide minimum_norm = 0,
                         const std::function<bool(const Cap &, Point)> &eligible = {}) {
    const Point d = ring.beta;
    IdentifierParabolicCache unused;
    const auto &functionals = cap_root_functionals();
    boundary_polytope_points(d, [&](const Cap &raw, bool interior) {
        if (identifier_cap_index(raw) < minimum_norm) {
            identifier_count_add(stats.norm_rejections, 1);
            return;
        }
        if (eligible && !eligible(raw,d)) {
            identifier_count_add(stats.unary_rejections,1);
            return;
        }
        if (interior) {
            emit(raw, ring.h22);
            identifier_count_add(stats.strict_direct, 1);
            return;
        }
        auto position = boundary_owner_position(raw, d);
        if (!position.owner) {
            identifier_count_add(stats.cone_rejections, 1);
            return;
        }
        if (!position.voronoi_boundary) {
            // The perturbed point lies in int(C), hence has no linear root
            // isotropy. Unique nearest translation then excludes every
            // affine root stabilizer and every nonzero stabilizer phase.
            emit(raw, ring.h22);
            identifier_count_add(stats.cone_boundary_direct, 1);
            return;
        }
        identifier_count_add(stats.voronoi_boundary, 1);
        for (const auto &row : functionals) {
            identifier_count_add(stats.root_checks, 1);
            IdentifierPoint value{};
            for (int j = 0; j < 4; ++j)
                value = identifier_add(value, identifier_mul(identifier_point(row[j]), identifier_point(raw[j])));
            if (ring.reduce(value) == Point{}) {
                identifier_count_add(stats.horizontal_rejections, 1);
                return;
            }
        }
        auto image = prepare_boundary_cap(raw, ring, unused, &stats);
        need(image.geometric, "horizontal-open cap reached the singular boundary fallback");
        if (image.value.residue != raw) {
            identifier_count_add(stats.duplicate_rejections, 1);
            return;
        }
        need(image.value.phase == 0, "selected cap owner has a nonzero canonical transport phase");
        emit(raw, image.value.subgroup);
    }, stats);
}

// This diagnostic reports the complete residue-space cardinality in a U.
CuspCaps boundary_census_caps(Point d, BoundaryOwnerStats *reported = nullptr, Wide minimum_norm = 0,
                             const std::function<bool(const Cap &, Point)> &eligible = {}) {
    BoundaryOwnerStats local;
    auto &stats = reported ? *reported : local;
    IdentifierRing ring(d);
    const Wide square = identifier_mul(ring.D, ring.D);
    need(square <= std::numeric_limits<U>::max() / square, "cap residue state count exceeds uint64");
    CuspCaps result;
    result.states = static_cast<U>(square * square);
    visit_boundary_caps(ring, [&](const Cap &raw, I step) {
        result.representatives.push_back(identifier_scaled_cap(raw, {0,1}));
        result.fibers.push_back(step);
    }, stats, minimum_norm, eligible);
    return result;
}
} // namespace gf
