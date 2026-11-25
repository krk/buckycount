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

// Exact runtime for generated IPR oriented/rotary tables; GMP arithmetic only.
#pragma once
namespace ipr_arithmetic {
inline Integer affine_floor(const Integer& m, unsigned p, unsigned b, unsigned d) {
    if (m <= 0 || !p || !d) fail("invalid IPR affine argument");
    return (m*p+b)/d;
}

inline std::optional<Integer> affine_argument(const Integer& m, unsigned p, unsigned b, unsigned d) {
    const Integer value = m*p+b;
    if (value%d != 0) return std::nullopt;
    return value/d;
}
inline void require_arithmetic_index(const Integer& m) {
    if (m <= 0) fail("IPR arithmetic requires a positive index");
}

inline Rational oriented_negative(const Integer& m) {
    require_arithmetic_index(m);
    static const std::vector<Integer> numerators = [] {
        std::vector<Integer> values;
        values.reserve(std::size(oriented_terms));
        for (const auto& term : oriented_terms) values.emplace_back(term.numerator);
        return values;
    }();
    static const Integer denominator(oriented_denominator);
    Integer total;
    std::size_t index = 0;
    for (const auto& term : oriented_terms) {
        const Integer& coefficient = numerators[index++];
        if (m % term.modulus != term.residue) continue;
        const auto n = affine_argument(m, term.p, term.b, term.d);
        if (n) total += coefficient * two_character_sum(term.weight, term.left, term.right, *n);
    }
    return make_rational(total, denominator);
}

inline Integer affine_tau(const Integer& m, unsigned long p, unsigned long b, unsigned long d) {
    const auto n = affine_argument(m, p, b, d);
    return n ? Integer(divisor_count(*n)) : Integer(0);
}

inline Integer affine_same(const Integer& m, unsigned long b) {
    const auto n = affine_argument(m, 1, b, 2);
    return n ? Integer(same_parity_pairs(*n)) : Integer(0);
}

inline Rational q2_negative(const Integer& m) {
    require_arithmetic_index(m);
    const Integer n = affine_floor(m, 1, 2, 1);
    return scaled(2 * affine_tau(m, 1, 2, 6), 3)
        + scaled(affine_tau(m, 1, 2, 2) - Integer(restricted_divisor_count(6, 4, n)), 6);
}

inline Rational q3_binary_negative(const Integer& m) {
    require_arithmetic_index(m);
    Rational total = -scaled(affine_tau(m, 3, 2, 4) + affine_tau(m, 3, 4, 4), 12)
        - scaled(affine_tau(m, 3, 6, 4) - affine_tau(m, 1, 2, 4), 2);
    total -= scaled(3 * affine_tau(m, 1, 2, 4) + affine_tau(m, 1, 4, 4)
                    + 2 * affine_tau(m, 1, 6, 4) + 2 * affine_tau(m, 1, 10, 4), 2);
    total += scaled(-2 * affine_tau(m, 1, 2, 2) + affine_same(m, 2)
                    - 2 * affine_same(m, 4) - 2 * affine_tau(m, 1, 6, 2)
                    + affine_same(m, 6), 2);
    return total;
}

inline Rational q3_negative(const Integer& m) {
    require_arithmetic_index(m);
    const Integer n1 = affine_floor(m, 9, 12, 1);
    const Integer n2 = affine_floor(m, 9, 6, 1);
    const Rational first = scaled(-s_at(n1, 2) + 3 * t_at(n1, 2), 32)
        - scaled(s_at(n1, 4) + 3 * t_at(n1, 4), 48);
    const Rational second = scaled(-s_at(n2, 4) + 3 * t_at(n2, 4), 48);
    return Rational(4, 5) * (first + second) + q3_binary_negative(m);
}

inline Rational pell_mass(unsigned field, const Integer& n, bool negative = false) {
    if (field == 3) {
        const Integer second = two_character_sum(0, -3, -4, n);
        return scaled(two_character_sum(0, 1, 12, n) + (negative ? -second : second), 4);
    }
    if (field == 6 && !negative)
        return scaled(two_character_sum(0, 1, 24, n) + two_character_sum(0, -8, -3, n), 4);
    fail("unknown IPR Pell genus");
}

inline Rational b27_zero(const Integer& n) {
    return n % 3 == 0 ? pell_mass(3, n / 3) : Rational();
}

inline Rational q4_negative(const Integer& m) {
    require_arithmetic_index(m);
    const Integer n1 = affine_floor(m, 1, 1, 1);
    const Integer n2 = affine_floor(m, 1, 2, 1);
    Rational total = t9_mass(affine_floor(m, 3, 3, 1)) - t27_mass(n1);
    total += Rational(2, 3) * (pell_mass(6, affine_floor(m, 4, 18, 1))
                              - pell_mass(6, affine_floor(m, 4, 6, 1)));
    total -= Rational(2, 3) * (pell_mass(3, n1, true) - b27_zero(n1));
    total -= Rational(1, 3) * (pell_mass(3, n2, true) + b27_zero(n2));
    total += scaled(2 * Integer(restricted_divisor_count(6, 3, n2)), 3);
    total += scaled(Integer(same_parity_pairs(n2)) - Integer(restricted_divisor_count(6, 1, n2)), 6);
    const auto square = [&m](unsigned long divisor, unsigned long shift) -> Integer {
        return Integer(square_indicator(divisor, affine_floor(m, 1, shift, 1)));
    };
    total -= square(1, 6);
    total += scaled(square(3, 2) - 4 * square(3, 5) + square(3, 6)
                    - 2 * square(3, 9) - 2 * square(3, 14), 3);
    total -= scaled(square(1, 1) - square(9, 1), 18);
    total -= scaled(square(9, 2) + square(1, 3) - square(9, 3), 3);
    return total;
}

inline Rational order10_negative(const Integer& m) {
    require_arithmetic_index(m);
    if (m % 5 != 0) return Rational();
    return scaled(2 * Integer(square_indicator(1, affine_floor(m / 5, 3, 1, 1))), 5);
}

inline Rational rotary_negative(const Integer& m) {
    return q2_negative(m) + q3_negative(m) + q4_negative(m) + order10_negative(m);
}

inline Integer ordinary_cusp(const Integer& m) {
    require_arithmetic_index(m);
    return Integer(ordinary_cusp_values[residue(m, std::size(ordinary_cusp_values))]);
}

inline Integer oriented_cusp(const Integer& m) {
    require_arithmetic_index(m);
    return Integer(oriented_cusp_values[residue(m, std::size(oriented_cusp_values))]);
}

} // namespace ipr_arithmetic
