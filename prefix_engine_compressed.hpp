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
#include "prefix_engine_formula.hpp"
namespace {
// The seven sector rows, their independent achiral expansion, integrality,
// the oriented bound and Burnside parity, without the IPR counts.
inline Row evaluate_ordinary_index(const Integer& m) {
    if (m == 0) fail("the index must be a positive integer");
    prepare_shared_weights(m);
    const Expansions& rows = expansions();
    Row row;
    row.index = m;
    row.oriented = oriented_count(m);
    row.split = evaluate_rational_row(rows.split, m);
    row.antipodal = evaluate_rational_row(rows.antipodal, m);
    row.q2 = q2_open(m);
    row.q3 = q3_open(m);
    row.q4 = q4_open(m);
    row.order10 = order10_open(m);
    row.order12 = order12_open(m);
    const Rational sectors =
        row.split + row.antipodal + row.q2 + row.q3 + row.q4 + row.order10 + row.order12;
    const Rational expansion = evaluate_rational_row(rows.achiral, m);
    if (mpq_equal(sectors.get_mpq_t(), expansion.get_mpq_t()) == 0) {
        fail(
            "the seven sector rows give a_m = " + rational_text(sectors) +
            " but the single expansion gives " + rational_text(expansion) + " at m=" + m.get_str()
        );
    }
    if (sectors.get_den() != 1 || sgn(sectors) < 0) {
        fail("the achiral correction is not a nonnegative integer at m=" + m.get_str());
    }
    row.achiral = sectors.get_num();
    if (row.achiral > row.oriented) fail("the achiral correction exceeds the oriented count at m=" + m.get_str());
    const Integer sum = row.oriented + row.achiral;
    if (mpz_odd_p(sum.get_mpz_t()) != 0) fail("Burnside parity fails at m=" + m.get_str());
    row.ordinary = sum / 2;
    return row;
}

} // namespace
