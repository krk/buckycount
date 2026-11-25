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

// For two raw caps in a fixed chart, the nearest representatives and both
// Heisenberg corrections are independent of the scalar phase. The remaining
// central translation uses exactly the half-open strip of normalize_identifier_cusp.
class CensusCentralStrip {
    Point beta_;
    Wide d_, twice_d_, height_, x_weight_, y_weight_;

public:
    explicit CensusCentralStrip(Point beta)
        : beta_(beta), d_(identifier_norm(identifier_point(beta))), twice_d_(identifier_mul(2, d_)),
          height_(Wide(beta.x) + beta.y), x_weight_(Wide(beta.x) - beta.y),
          y_weight_(2 * Wide(beta.x) + beta.y) {
        need(d_ > 0, "cannot center a zero height");
    }

    IdentifierCusp center(Point u, const Cap& first, const Cap& second,
                          IdentifierPoint first_delta, IdentifierPoint second_delta) const {
        const auto shifted = identifier_add(identifier_point(u), identifier_add(first_delta, second_delta));
        // If w = shifted * conjugate(beta), w.x + 2*w.y is one dot
        // product. The central translation direction is beta*zeta.
        const auto projection = identifier_add(identifier_mul(x_weight_, shifted[0]),
                                               identifier_mul(y_weight_, shifted[1]));
        const auto wraps = identifier_floor(identifier_add(d_, projection), twice_d_);
        const IdentifierPoint centered{
            identifier_add(shifted[0], identifier_mul(beta_.y, wraps)),
            identifier_sub(shifted[1], identifier_mul(height_, wraps))};
        return {identifier_narrow(centered), beta_, {first, second}};
    }
};

inline IdentifierCusp census_centered_period(Point u, Point beta, const Cap& first, const Cap& second,
                                             IdentifierPoint first_delta, IdentifierPoint second_delta) {
    return CensusCentralStrip(beta).center(u, first, second, first_delta, second_delta);
}

using CensusTips = std::array<Point, 11>;

class CensusTipReduction {
    // For |coordinates|<=512, q=p_i+zeta*edge has coordinates <=1536
    // and q-p_j <=2048. Thus 128*(Q(q-p_j)-Q(edge)) + 11*i+j and the
    // individual terms in its expanded dot product fit signed 32-bit words.
    static constexpr int limit = 512, key_base = 128;
    alignas(64) std::array<int, 16> xs_{}, ys_{}, keys_{};

    static int norm32(int x, int y) { return x * x + x * y + y * y; }
    void put(std::size_t i, int x, int y) {
        xs_[i] = x;
        ys_[i] = y;
        keys_[i] = key_base * norm32(x, y) + static_cast<int>(i);
    }
    int gap() const {
        int best = 0;
        for (int i = 0; i < 11; ++i) {
            auto next = i == 10 ? 0 : i + 1;
            int dx = xs_[next] - xs_[i], dy = ys_[next] - ys_[i];
            int qx = xs_[i] - dy, qy = ys_[i] + dx + dy;
            int off = key_base * (norm32(qx, qy) - norm32(dx, dy)) + 11 * i;
            int ax = key_base * (2 * qx + qy), ay = key_base * (qx + 2 * qy);
            int row = std::numeric_limits<int>::max();
            // Five padding lanes copy tip zero with larger tie indices.
            // The fixed trip count lets compilers vectorize the exact minimum.
            for (int j = 0; j < 16; ++j)
                row = std::min(row, keys_[j] - ax * xs_[j] - ay * ys_[j]);
            best = std::min(best, off + row);
        }
        return best;
    }
    bool recut(int k, bool twice) {
        int n = k + 2 == 11 ? 0 : k + 2;
        int a, b;
        if (twice) {
            // Two positive Hurwitz moves at k equal
            // p[k] + zeta^(-1) * (p[k+2]-p[k+1]).
            int dx = xs_[n] - xs_[k + 1], dy = ys_[n] - ys_[k + 1];
            a = xs_[k] + dx + dy;
            b = ys_[k] - dx;
        } else {
            int dx = xs_[k] - xs_[k + 1], dy = ys_[k] - ys_[k + 1];
            a = xs_[n] - dy;
            b = ys_[n] + dx + dy;
        }
        if (std::max(std::abs(a), std::abs(b)) > limit) return false;
        put(k + 1, a, b);
        return true;
    }

  public:
    CensusTips reduce(const CensusTips& input) {
        for (auto v : input)
            if (v.x > limit || v.x < -limit || v.y > limit || v.y < -limit) return reduce_star(input).tips;
        auto area = star_area(input), energy = star_energy(input);
        need(area > 0, "star reduction requires a positive period");
        for (int j = 0; j < 16; ++j) {
            auto v = input[j < 11 ? j : 0];
            put(j, static_cast<int>(v.x), static_cast<int>(v.y));
        }
        for (;;) {
            int score = gap();
            if (!score) break;
            int ij = score % key_base;
            if (ij < 0) ij += key_base;
            energy += (score - ij) / key_base;
            const int i = ij / 11, j = ij % 11;
            if (j < i) {
                for (int k = i - 1; k >= j; --k)
                    if (!recut(k, false)) return reduce_star(input).tips;
            } else {
                for (int k = i; k < j - 1; ++k)
                    if (!recut(k, true)) return reduce_star(input).tips;
            }
        }
        CensusTips tips;
        for (int j = 0; j < 11; ++j)
            tips[j] = {xs_[j], ys_[j]};
        need(star_energy(tips) == energy && star_area(tips) == area, "stream recuts changed period");
        // The zero-gap collision and bubble rules agree with the generic
        // reducer. On these uncommon inputs reuse its exact escape logic.
        if (star_bubble(tips, 1) >= 0 || star_bubble(tips, 2) >= 0 || star_bubble(tips, 3) >= 0)
            return reduce_star(tips).tips;
        return tips;
    }
};

inline CensusTips census_reduced_tips(const IdentifierCusp& centered) {
    return CensusTipReduction{}.reduce(cusp_star(centered.u, centered.beta, centered.caps[0], centered.caps[1]));
}

} // namespace gf::fullerene_rank
