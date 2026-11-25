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

#include "fullerene_hex_witness.hpp"

namespace gf::fullerene_rank {

class CensusContactPhases {
    struct Target {
        Wide energy;
        Point unit;
    };
    struct Direction {
        Point gamma;
        IdentifierPoint scalar_factor;
        Wide radius = 0;
        std::vector<Target> targets;
        std::map<Point, std::uint64_t> phase_masks;
    };
    struct Neighbor {
        Wide energy;
        Point projection;
    };

    IdentifierRing ring_;
    std::vector<Direction> directions_;
    std::uint64_t valid_offsets_ = 0;

    static IdentifierPoint cap_projection(const Cap& y, const Cap& z) {
        IdentifierPoint result;
        need(identifier_divide(identifier_cap_pair(y, z), {2, -1}, result),
             "nonintegral contact cap projection");
        return result;
    }

public:
    struct Prepared {
        unsigned horizontal = 0;
        std::vector<std::vector<Neighbor>> neighbors;
    };

    CensusContactPhases(I m, Point beta) : ring_(beta) {
        const I h = beta.x + beta.y;
        need(m > 0 && beta.x > 0 && (beta.y == 0 || beta.y >= beta.x)
             && Wide(h) * h < Wide(6) * m, "contact family outside ordinary height bound");
        const Point step{-beta.y / ring_.h22, h / ring_.h22};
        for (I offset = 0; offset < 64; offset += step.y)
            valid_offsets_ |= std::uint64_t{1} << offset;

        // A contact satisfies D/n >= m/4 and has one gamma in the half-open
        // unit wedge. Its normalized pairing equals (2-zeta)*unit.
        const Wide bound = 4 * Wide(ring_.D) / m;
        for (I a = 1; Wide(a) * a <= bound; ++a)
            for (I b = 0; identifier_norm({a, b}) <= bound; ++b) {
                const Point gamma{a, b};
                const Wide n = identifier_norm({a, b});
                const auto coeff = identifier_mul(identifier_conjugate(identifier_point(beta)),
                                                  identifier_point(gamma));
                Direction direction{gamma, identifier_mul({0, 1}, identifier_conjugate({a, b})), 0, {}, {}};
                for (int k = 0; k < 6; ++k) {
                    const Point unit = rotate({1, 0}, k);
                    const auto trace = identifier_mul(coeff, identifier_point(multiply({2, -1}, unit)));
                    const Wide numerator = 2 * trace[0] + trace[1];
                    if (numerator % 3) continue;
                    const Wide energy = ring_.D - Wide(m) * n + numerator / 3;
                    if (energy < 0) continue;
                    direction.targets.push_back({energy, unit});
                    direction.radius = std::max(direction.radius, energy);
                }
                if (direction.targets.empty()) continue;
                const auto slope = identifier_mul(direction.scalar_factor, identifier_point(step));
                for (I t = 0; t * step.y < 64; ++t) {
                    const auto residue = ring_.reduce(identifier_mul(slope, {t, 0}));
                    direction.phase_masks[residue] |= std::uint64_t{1} << (t * step.y);
                }
                directions_.push_back(std::move(direction));
            }
    }

    Prepared prepare(const Cap& cap) const {
        Prepared result;
        // Fixed root per cap-unit line. Each soluble pairing specifies a
        // distinct scalar A, hence a distinct horizontal contact root line.
        for (const auto& root : identifier_cap_roots()) {
            const auto projection = cap_projection(root, cap);
            for (int k = 0; k < 6; ++k)
                if (ring_.reduce(identifier_add(projection, identifier_point(rotate({1, 0}, k)))) == Point{})
                    ++result.horizontal;
        }
        for (const auto& dir : directions_) {
            auto& neighbors = result.neighbors.emplace_back();
            visit_hex_cap_neighbors(ring_.beta, dir.gamma, cap, dir.radius, [&](const Cap& y, Wide energy) {
                neighbors.push_back({energy, ring_.reduce(cap_projection(y, cap))});
                return true;
            });
            std::sort(neighbors.begin(), neighbors.end(), [](auto& a, auto& b) { return a.energy < b.energy; });
        }
        return result;
    }

    // A set bit j represents raw scalar q=base.y+j, with p determined by
    // the size equation. Both caps must be prepared in this coefficient.
    std::array<std::uint64_t, 31> histogram(const Prepared& first, const Prepared& second,
                                            Point base, std::uint64_t mask) const {
        std::array<std::uint64_t, 31> result{};
        if (!mask) return result;
        need(!(mask & ~valid_offsets_), "contact mask contains a nonsoluble scalar offset");
        const auto horizontal = first.horizontal + second.horizontal;
        need(horizontal < result.size() && first.neighbors.size() == directions_.size()
             && second.neighbors.size() == directions_.size(), "invalid prepared contact caps");
        result[horizontal] = std::popcount(mask);
        std::array<std::uint64_t, 5> bits{};
        std::uint64_t touched = 0;
        for (std::size_t i = 0; i < directions_.size(); ++i) {
            const auto& dir = directions_[i];
            const auto& left = first.neighbors[i];
            const auto& right = second.neighbors[i];
            if (left.empty() || right.empty() || left.front().energy + right.front().energy > dir.radius)
                continue;
            const Point scalar = ring_.reduce(identifier_mul(dir.scalar_factor, identifier_point(base)));
            for (const auto& target : dir.targets)
                for (const auto& x : left) {
                    const auto remaining = target.energy - x.energy;
                    if (remaining < right.front().energy) break;
                    auto y = std::lower_bound(right.begin(), right.end(), remaining,
                                              [](auto& row, Wide e) { return row.energy < e; });
                    for (; y != right.end() && y->energy == remaining; ++y) {
                        // Solving h(r,v)=(2-zeta)*unit fixes conjugate(A).
                        // The trace energy equality already enforces Q(r)=-1.
                        // Divisibility by beta is the last condition, linear
                        // in the phase and independent of the root lift.
                        auto residue = ring_.reduce(identifier_point(x.projection + y->projection + target.unit - scalar));
                        const auto it = dir.phase_masks.find(residue);
                        if (it == dir.phase_masks.end()) continue;
                        auto carry = mask & it->second;
                        touched |= carry;
                        for (auto& digit : bits) {
                            if (!carry) break;
                            const auto old = digit;
                            digit ^= carry;
                            carry &= old;
                        }
                        need(!carry, "more than thirty contact lines");
                    }
                }
        }
        // Recover each distinct exceptional contact count in one mask
        // operation per binary digit, independent of phase multiplicity.
        result[horizontal] -= std::popcount(touched);
        while (touched) {
            const auto one = touched & (0 - touched);
            auto same = touched;
            unsigned count = horizontal;
            for (unsigned i = 0; i < bits.size(); ++i) {
                if (bits[i] & one) {
                    count += 1U << i;
                    same &= bits[i];
                } else {
                    same &= ~bits[i];
                }
            }
            need(count < result.size(), "more than thirty pentagon contacts");
            result[count] += std::popcount(same);
            touched ^= same;
        }
        return result;
    }
};

} // namespace gf::fullerene_rank
