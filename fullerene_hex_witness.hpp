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

namespace gf::fullerene_rank {

inline Wide hex_height(Point beta) {
    Wide a = beta.x, b = beta.y, c = a + b;
    return std::max({a, -a, b, -b, c, -c});
}

// Write T=H(beta)^2. An H competitor has norm height <=T, giving
// 3m^2 N(gamma) <= 4DT by projection to the negative complement.
inline Wide hex_gamma_bound(Point beta, I m) {
    Wide h = hex_height(beta), d = identifier_norm(identifier_point(beta));
    need(d > 0 && m > 0, "hexagonal cusp comparison needs positive height and index");
    return identifier_mul(identifier_mul(4, d), identifier_mul(h, h)) /
           identifier_mul(3, identifier_mul(m, m));
}

// The horizontal identity bounds the sum over the two cap differences by
// (2/3)Re(gamma*conj(beta)*bar(theta)*q)-m*N(gamma), with H(q)<=H(beta).
// The hexagon's vertices are H(beta) times the six units. For z=x+y*zeta,
// (2/3)Re(bar(theta)*z)=x+y, whose six rotations have maximum H(z).
inline Wide hex_cap_radius(Point beta, I m, Point gamma) {
    Wide n = identifier_norm(identifier_point(gamma));
    need(beta != Point{} && m > 0 && n > 0, "invalid hexagonal cusp radius parameters");
    auto z = identifier_narrow(
        identifier_mul(identifier_point(gamma), identifier_conjugate(identifier_point(beta))));
    return identifier_sub(identifier_mul(hex_height(beta), hex_height(z)), identifier_mul(m, n));
}

template <class Visit>
bool visit_hex_cap_neighbors(Point beta, Point gamma, const Cap &cap, Wide radius, Visit &&visit) {
    if (radius < 0)
        return true;
    auto center = identifier_scaled_cap(cap, gamma);
    CapBall ball(beta, center, 1, 0, 1);
    ball.limit = identifier_mul(
        identifier_mul(ball.decomposition.scale, identifier_norm(identifier_point(beta))), radius);
    return ball.enumerate([&](const Cap &y) {
        Wide energy = identifier_shifted_index(y, beta, center);
        need(energy <= radius, "cap exceeds hexagonal competitor radius");
        return visit(y, energy);
    });
}

// Intersect the three closed facet strips on the scalar line b + j*s.
inline std::optional<std::pair<Wide, Wide>> hex_scalar_interval(IdentifierPoint b, IdentifierPoint s,
                                                                Wide height) {
    need(s != IdentifierPoint{} && height >= 0, "invalid hexagonal scalar line");
    std::array<Wide, 3> base{b[0], b[1], identifier_add(b[0], b[1])};
    std::array<Wide, 3> slope{s[0], s[1], identifier_add(s[0], s[1])};
    Wide low = std::numeric_limits<Wide>::min(), high = std::numeric_limits<Wide>::max();
    for (std::size_t i = 0; i < base.size(); ++i) {
        Wide c = base[i], d = slope[i];
        if (d < 0) {
            c = identifier_sub(0, c);
            d = identifier_sub(0, d);
        }
        if (d == 0) {
            if (c < -height || c > height)
                return std::nullopt;
        } else {
            low = std::max(low, identifier_ceil(identifier_sub(-height, c), d));
            high = std::min(high, identifier_floor(identifier_sub(height, c), d));
        }
    }
    if (low > high)
        return std::nullopt;
    return std::pair{low, high};
}

// One unit representative per non-original primitive null line w with
// H(h(w,v)/bar(theta)) <= H(beta). A false callback stops the search.
template <class Visit> bool visit_hex_competitors(const IdentifierCusp &v, Visit &&visit) {
    I m = identifier_narrow(identifier_cusp_index(v));
    Wide height = hex_height(v.beta), bound = hex_gamma_bound(v.beta, m);
    for (I a = 1; identifier_mul(a, a) <= bound; ++a) {
        for (I b = 0; identifier_norm({a, b}) <= bound; ++b) {
            Point gamma{a, b};
            Wide radius = hex_cap_radius(v.beta, m, gamma);
            if (radius < 0)
                continue;
            auto [g, p, q] = bezout(a + b, b);
            Point kappa{-b / g, (a + b) / g};
            auto slope =
                identifier_mul(identifier_conjugate(identifier_point(kappa)), identifier_point(v.beta));
            IdentifierEnumerationStats stats;
            bool complete =
                identifier_horizontal_ball(v, gamma, radius, stats, [&](const std::array<Cap, 2> &caps) {
                    Wide energy =
                        identifier_add(identifier_cap_index(caps[0]), identifier_cap_index(caps[1]));
                    if (energy % g)
                        return true;
                    IdentifierPoint initial{identifier_mul(p, energy / g), identifier_mul(q, energy / g)};
                    IdentifierCusp w{identifier_narrow(initial), gamma, caps};
                    IdentifierPoint base;
                    need(identifier_divide(identifier_cusp_pair(w, v), {2, -1}, base),
                         "nonintegral hexagonal cusp coefficient");
                    auto range = hex_scalar_interval(base, slope, height);
                    if (!range)
                        return true;
                    for (Wide j = range->first; j <= range->second; ++j) {
                        w.u = identifier_narrow(
                            identifier_add(initial, identifier_mul(identifier_point(kappa), {j, 0})));
                        if (primitive_identifier_cusp(w)) {
                            auto coefficient = identifier_add(base, identifier_mul(slope, {j, 0}));
                            if (!visit(w, hex_height(identifier_narrow(coefficient))))
                                return false;
                        }
                    }
                    return true;
                });
            if (!complete)
                return false;
        }
    }
    return true;
}

} // namespace gf::fullerene_rank
