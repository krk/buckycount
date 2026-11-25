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
#include "census_small_height.hpp"

namespace gf::fullerene_rank {

// Compile strict nonhorizontal witnesses to products of unary posting lists.
// A bit is cleared only when the same witnessed product rules out every
// scalar residue. Equality is deliberately left to the ordinary tie test.
template<class HeightLine> class CensusWitnessProducts {
    struct Feature {
        I norm, cost, energy;
        IdentifierPoint projection;
        auto operator<=>(const Feature&) const = default;
    };
    struct Posting {
        Feature key;
        std::vector<std::size_t> caps;
    };
    const IdentifierRing& ring_;
    I m_, radius_;
    std::size_t size_, words_, first_ = 0, end_ = 0;
    I divisor_;
    IdentifierPoint scalar_factor_, energy_factor_, increment_;
    std::optional<HeightLine> line_;
    std::vector<Posting> postings_;
    std::vector<U> rejected_;
    static constexpr std::size_t width = std::numeric_limits<U>::digits;
    static constexpr std::size_t tile_words = 1U << 22;

    bool excludes(const Feature& a, const Feature& b) const {
        const Wide total = Wide(m_) + a.norm + b.norm;
        const Wide energy = Wide(a.energy) + b.energy;
        if (energy % divisor_ || total % ring_.h22) return false;
        const Wide n = total / ring_.h22, e = energy / divisor_;
        auto value = identifier_sub(identifier_mul(scalar_factor_, {ring_.bezout_s*n, ring_.bezout_t*n}),
                                    identifier_add(a.projection, b.projection));
        value = identifier_add(value, {energy_factor_[0]*e, energy_factor_[1]*e});
        for (I phase = 0; phase < ring_.h22; ++phase) {
            if (line_->visit(value, [](Point, Wide) { return true; })) return false;
            value = identifier_add(value, increment_);
        }
        return true;
    }

    void build_tile(std::size_t first) {
        first_ = first;
        end_ = first + std::min(size_ - first, std::max<std::size_t>(1, tile_words / words_));
        rejected_.assign((end_ - first) * words_, 0);
        std::vector<U> bad(words_);
        for (const auto& a : postings_) {
            auto lo = std::lower_bound(a.caps.begin(), a.caps.end(), first);
            if (lo == a.caps.end() || *lo >= end_) continue;
            std::fill(bad.begin(), bad.end(), 0);
            bool any = false;
            for (const auto& b : postings_) {
                if (Wide(a.key.cost) + b.key.cost > radius_) break;
                if (!excludes(a.key, b.key)) continue;
                any = true;
                for (auto j : b.caps) bad[j / width] |= U{1} << (j % width);
            }
            if (!any) continue;
            for (; lo != a.caps.end() && *lo < end_; ++lo) {
                auto row = rejected_.data() + (*lo - first) * words_;
                for (std::size_t k = *lo / width; k < words_; ++k) row[k] |= bad[k];
            }
        }
    }

  public:
    CensusWitnessProducts(I m, const IdentifierRing& ring, const std::vector<CensusHeightCap>& caps,
                          Point gamma)
        : ring_(ring), m_(m), radius_(identifier_narrow(hex_cap_radius(ring.beta, m, gamma))),
          size_(caps.size()), words_((size_ + width - 1) / width) {
        if (radius_ < 0 || caps.empty()) return;
        auto [d, p, q] = bezout(gamma.x + gamma.y, gamma.y);
        divisor_ = d;
        const IdentifierPoint kappa{-gamma.y/d, (gamma.x+gamma.y)/d};
        const auto beta = identifier_point(ring.beta);
        scalar_factor_ = identifier_mul({0, 1}, identifier_conjugate(identifier_point(gamma)));
        energy_factor_ = identifier_mul(identifier_conjugate({p,q}), beta);
        increment_ = identifier_mul(scalar_factor_, {-ring.beta.y/ring.h22, (ring.beta.x+ring.beta.y)/ring.h22});
        line_.emplace(identifier_mul(identifier_conjugate(kappa), beta), ring.beta);
        std::map<Feature, std::vector<std::size_t>> members;
        for (std::size_t i = 0; i < size_; ++i) {
            const auto& z = caps[i].geometric.value.nearest;
            visit_hex_cap_neighbors(ring.beta, gamma, z, radius_, [&](const Cap& y, Wide cost) {
                IdentifierPoint projection;
                need(identifier_divide(identifier_cap_pair(y, z), {2,-1}, projection), "invalid witness projection");
                Feature key{identifier_narrow(caps[i].nearest_norm), identifier_narrow(cost),
                            identifier_narrow(identifier_cap_index(y)), projection};
                auto& ids = members[key];
                if (ids.empty() || ids.back() != i) ids.push_back(i);
                return true;
            });
        }
        for (auto& [key, ids] : members) postings_.push_back({key, std::move(ids)});
        std::sort(postings_.begin(), postings_.end(),
                  [](const Posting& a, const Posting& b) { return a.key.cost < b.key.cost; });
        if (!postings_.empty() && 2*Wide(postings_.front().key.cost) > radius_) postings_.clear();
    }

    bool active() const { return !postings_.empty(); }
    std::span<const U> rejected(std::size_t i) {
        if (!active()) return {};
        if (i < first_ || i >= end_) build_tile(i);
        return {rejected_.data() + (i - first_) * words_, words_};
    }
};

} // namespace gf::fullerene_rank
