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

struct CensusLocalHeight {
    bool operator==(const CensusLocalHeight &) const = default;
};
template<class Height> struct CensusPhaseWindow {
    [[no_unique_address]] Height beta;
    std::size_t first, second;
    I low_q;
    std::uint32_t phases;
    unsigned long end;
};

using CensusPhaseGroup = CensusPhaseWindow<Point>;
using CensusLocalPhaseGroup = CensusPhaseWindow<CensusLocalHeight>;

// For one cap pair, q arrives in strictly increasing order. Opening a new
// window retains the cumulative rank even across a gap wider than the mask.
template<class Bank, class Height>
inline void append_census_phase(Bank& groups, Height beta,
                                std::size_t first, std::size_t second, I q) {
    const unsigned long count = groups.empty() ? 0 : groups.back().end;
    need(count < std::numeric_limits<unsigned long>::max(), "phase group rank overflow");
    if (!groups.empty()) {
        auto& last = groups.back();
        if (last.first == first && last.second == second && last.beta == beta) {
            const Wide offset = Wide(q) - last.low_q;
            need(offset >= 0, "phase coordinates must increase");
            if (offset < std::numeric_limits<std::uint32_t>::digits) {
                const auto flag = std::uint32_t{1} << static_cast<unsigned>(offset);
                need(last.phases < flag, "duplicate or decreasing phase coordinate");
                last.phases |= flag;
                ++last.end;
                return;
            }
        }
    }
    groups.push_back({beta, first, second, q, 1, count + 1});
}

} // namespace gf::fullerene_rank
