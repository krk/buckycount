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

// Exact ordinary and isolated-pentagon counts.
#pragma once

#include "prefix_engine_compressed.hpp"

namespace {

// Twelve vertex-disjoint pentagons require at least sixty vertices.
constexpr unsigned long first_ipr_index = 30;

#include "ipr_hecke_arithmetic.hpp"
#include "ipr_real_affine_arithmetic.hpp"

Row evaluate_index(const Integer& m, bool include_ipr = false) {
    Row row = evaluate_ordinary_index(m);
    if (include_ipr) {
        Rational ipr_ordinary, ipr_oriented;
        if (m >= first_ipr_index) {
            const Rational positive = ipr_arithmetic::oriented_negative(m);
            Rational negative = ipr_real::negative_split(m);
            negative += ipr_real::negative_antipodal(m) + ipr_arithmetic::rotary_negative(m);
            ipr_ordinary = Rational(row.ordinary) - (positive + negative) / 2 -
                ipr_arithmetic::ordinary_cusp(m);
            ipr_oriented = Rational(row.oriented) - positive -
                ipr_arithmetic::oriented_cusp(m);
        }
        if (ipr_ordinary.get_den() != 1 || ipr_oriented.get_den() != 1 ||
            ipr_ordinary < 0 || ipr_oriented < 0 ||
            ipr_ordinary > row.ordinary || ipr_oriented > row.oriented)
            fail("the IPR formula violates count normalization at m=" + m.get_str());
        IprCounts counts{ipr_ordinary.get_num(), ipr_oriented.get_num(), 0};
        counts.achiral = 2 * counts.ordinary - counts.oriented;
        if (counts.achiral < 0 || counts.achiral > counts.oriented ||
            counts.achiral > row.achiral)
            fail("the IPR formula violates mirror normalization at m=" + m.get_str());
        row.ipr = std::move(counts);
    }
    return row;
}

inline void evaluate_batch(std::vector<Row>& rows, const Integer& first,
                           bool ipr, unsigned workers) {
    parallel_jobs(rows.size(), workers, [&](std::size_t i) {
        rows[i] = evaluate_index(first+i, ipr);
    });
}
} // namespace

#include "prefix_engine_cli.hpp"
