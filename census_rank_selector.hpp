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

#include "census_family_table.hpp"
#include "census_point_owner.hpp"
#include "census_phase_group.hpp"
#include "census_orbit_family.hpp"
#include "census_centered_stream.hpp"
#include <memory_resource>

namespace gf::fullerene_rank {
unsigned long ordinary_rank_max_index();

// Keep the low, frequently decremented part native. A refill is necessary
// only when a skipped block crosses a machine-word boundary.
class CensusRankCursor {
    mpz_class high_;
    unsigned long low_;

    void refill() {
        low_ = high_.fits_ulong_p() ? high_.get_ui() : ULONG_MAX;
        high_ -= low_;
    }
public:
    explicit CensusRankCursor(mpz_class rank) : high_(std::move(rank)) { refill(); }
    bool skip(unsigned long count) {
        if (low_ < count && high_ != 0) {
            high_ += low_;
            refill();
        }
        if (low_ < count) return false;
        low_ -= count;
        return true;
    }
    unsigned long offset() const { return low_; }
    mpz_class residual() const { return high_ + low_; }
};

class CensusRankSelector {
    I m;
    mpz_class expected;
    // Consume one family's mass or return the graph period at its local rank.
    std::optional<IdentifierCusp> take_family(Point beta, mpz_class &rank) const {
        const auto mass = CENSUS_HEAP_EVAL(family_mass, census_orbit_family_mass(m, beta));
        if (mass) {
            if (rank >= *mass) {
                rank -= *mass;
                return std::nullopt;
            }
            auto family = CENSUS_HEAP_EVAL(family_build, CensusOrbitFamilyIndex(m, beta));
            need(family.total() == *mass, "orbit directory disagrees with family mass");
            return centered_identifier_cusp(CENSUS_HEAP_EVAL(family_select, family.select(rank)));
        }
        CENSUS_HEAP_SCOPE(owner);
        const auto caps = census_small_height_inventory(beta);
        if (caps.empty()) return std::nullopt;
        CensusPointOwner owner(m, beta, caps);
        CensusRankCursor remaining(std::move(rank));
        struct Selected { IdentifierCusp centered; };
        try {
            owner.emit<true>([&](const auto &window, const auto &x, const auto &y) {
                if (remaining.skip(static_cast<unsigned long>(window.size()))) return;
                const auto scalar = window[static_cast<I>(remaining.offset())];
                throw Selected{census_centered_period(scalar, beta, x.nearest, y.nearest, x.to_nearest, y.to_nearest)};
            });
        } catch (const Selected &answer) { return answer.centered; }
        rank = remaining.residual();
        return std::nullopt;
    }
public:
    CensusRankSelector(I size, mpz_class total) : m(size), expected(std::move(total)) {}
    IdentifierCusp select_centered(mpz_class rank) const {
        if (rank < 0 || rank >= expected)
            throw std::out_of_range("rank outside progressive census");
        std::optional<IdentifierCusp> result;
        visit_census_families(m, [&](Point beta) {
            result = take_family(beta, rank);
            return !result;
        });
        if (result) return *result;
        throw std::logic_error("family masses did not cover the requested rank");
    }
    IdentifierCusp select_family_centered(Point beta, mpz_class local_rank) const {
        if (local_rank >= 0)
            if (auto period = take_family(beta, local_rank)) return *period;
        throw std::out_of_range("rank outside census family");
    }
};

} // namespace gf::fullerene_rank
