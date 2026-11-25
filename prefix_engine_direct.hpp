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

// Included inside the arithmetic engine namespace. Compute the coefficients
// needed by a single index without allocating tables through that index.
#pragma once

namespace direct_arithmetic {

using Factors = std::vector<std::pair<Integer, unsigned long>>;

inline Factors factor(const Integer& n) {
    if (n <= 0) fail("factorization needs a positive argument");
    return ExactFactorization(n).primes();
}

inline const Factors& cached_factors(const Integer& n) {
    thread_local std::map<Integer, Factors> recent;
    if (recent.size() >= 2048) recent.clear();
    auto it = recent.find(n);
    if (it == recent.end()) it = recent.emplace(n, factor(n)).first;
    return it->second;
}
// At p^e sum chi_left(p)^(e-j) (chi_right(p) p^w)^j for j=0..e.
inline Integer divisor_sum(unsigned weight, long left, long right, const Integer& n) {
    Integer total = 1;
    for (const auto& [p, exponent] : cached_factors(n)) {
        Integer base;
        mpz_pow_ui(base.get_mpz_t(), p.get_mpz_t(), weight);
        base *= kronecker(right, p);
        const int chi = kronecker(left, p);
        Integer sum = 1, power = 1;
        for (unsigned long j = 0; j < exponent; ++j) {
            power *= base;
            sum = chi * sum + power;
        }
        total *= sum;
    }
    return total;
}

enum class QuadraticSeries { hurwitz, even_cohen, odd_cohen };
template<class Visit>
void quadratic_factors(const Integer& D, QuadraticSeries series, Visit visit) {
    const bool hurwitz = series == QuadraticSeries::hurwitz;
    const bool odd = series == QuadraticSeries::odd_cohen;
    const Integer last = integer_square_root(hurwitz ? Integer(D / 3) : D);
    for (Integer r = odd ? 0 : residue(D, 2); r <= last; r += odd ? 1 : 2) {
        Integer n;
        if (hurwitz) n = (D + r*r)/4;
        else n = (D-r*r)/(odd ? 1 : 4);
        visit(r, n, n == 0 ? Factors{} : factor(n));
    }
}
inline int chi_minus_four(const Integer& n) {
    const auto r = residue(n, 4);
    return r == 1 ? 1 : r == 3 ? -1 : 0;
}
// Enumerate reduced forms (a,b,c) of discriminant -D. For fixed |b|,
// a*c=(D+b^2)/4, so only its divisors with |b| <= a <= c can contribute.
inline Integer fundamental_hurwitz(const Integer& D) {
    Integer total = 0;
    std::vector<Integer> divisors;
    quadratic_factors(D, QuadraticSeries::hurwitz, [&](const Integer& b, const Integer& product, const Factors& factors) {
        const auto last_a = integer_square_root(product);
        divisors.assign(1, 1);
        for (const auto& [p, exponent] : factors) {
            const auto count = divisors.size();
            Integer power = 1;
            for (unsigned long j = 0; j < exponent; ++j) {
                power *= p;
                if (power > last_a) break;
                const Integer bound = last_a / power;
                for (std::size_t i = 0; i < count; ++i)
                    if (divisors[i] <= bound) divisors.push_back(divisors[i] * power);
            }
        }
        Integer subtotal = 0;
        for (auto a : divisors) {
            if (a < b) continue;
            const bool equal = a * a == product;
            // Interior forms occur at both signs of b. The elliptic forms
            // at D=3,4 have stabilizers of order three and two respectively.
            if (equal && a == b) subtotal += 4;
            else if (equal && !b) subtotal += 6;
            else if (!b || a == b || equal) subtotal += 12;
            else subtotal += 24;
        }
        total += subtotal;
    });
    return total;
}

using PrimitiveValues = std::map<std::pair<unsigned, Integer>, Integer>;
inline PrimitiveValues& primitive_values() {
    thread_local PrimitiveValues cache;
    return cache;
}

// Cohen, "Computing L-Functions of Quadratic Characters at Negative
// Integers", Prop. 6.1 and section 6.4: U(4) preserves M_{7/2}(Gamma_0(4)),
// which has dimension two. Write B(n)=sum chi_-4(d)*d^2 and
// A(n)=sum chi_-4(n/d)*d^2, with B(0)=-1/4, A(0)=0.
// For n>0 the coefficient of q^n in F=theta*(4B-2A) is
// 252H_3(n). U(4)F=theta*(4B-128A): the two basis
// series theta*B and theta*A are independent at q^0,q^1, and both
// sides have constant -1 and coefficient -126 at q^1.
// With squarefree radicand s, each signed fundamental discriminant has
// magnitude s or 4s. U(4) puts the odd-weight sums at s. The even-weight
// sums also use s-r^2, with a division by four when +s is fundamental.
// Hence all four weights can share one factorization per shift.
using WeightSet = std::array<bool, 5>;

inline Integer prime_power_sum(const Integer& base, unsigned long exponent) {
    Integer sum = 1 + base;
    for (unsigned long j = 1; j < exponent; ++j) sum = sum * base + 1;
    return sum;
}

inline PrimitiveValues fundamental_values(const Integer& s, WeightSet need) {
    PrimitiveValues values;
    const Integer positive = s % 4 == 1 ? s : 4 * s;
    const Integer negative = s % 4 == 3 ? s : 4 * s;
    // Gauss: r3(s)=12H(4s) for s=1,2 mod 4; r3(s)=24H(s) for s=3 mod 8.
    // For s=7 mod 8 the three-square count vanishes.
    if (need[1] && negative % 8 == 7) {
        values.emplace(std::pair{1U, negative}, fundamental_hurwitz(negative));
        need[1] = false;
    }
    if (!need[1] && !need[2] && !need[3] && !need[4]) return values;

    const bool even_only = !need[1] && !need[3];
    const auto series = even_only ? QuadraticSeries::even_cohen : QuadraticSeries::odd_cohen;
    std::array<Integer, 5> sums{};
    quadratic_factors(even_only ? positive : s, series, [&](const Integer& r, const Integer& n, const Factors& factors) {
        // even_only already supplies (positive-r^2)/4. In the shared
        // traversal n=s-r^2, so remove two factors of 2 if positive=s.
        const bool divide_four = !even_only && positive == s;
        const bool even_term = n && (!divide_four || n % 4 == 0);
        const bool want2 = need[2] && even_term, want4 = need[4] && even_term;
        Integer two_squares = n ? 4 : 1, sigma1 = 1, B = 1, sigma3 = 1;
        for (const auto& [p, exponent] : factors) {
            if (p != 2) {
                const int chi = chi_minus_four(p);
                if (need[1]) {
                    if (chi == 1) two_squares *= exponent + 1;
                    else if (exponent % 2) two_squares = 0;
                }
                if (need[3]) {
                    Integer base = Integer(p) * p;
                    if (chi < 0) base = -base;
                    B *= prime_power_sum(base, exponent);
                }
            }
            if (want2 || want4) {
                const auto e = exponent - (divide_four && p == 2 ? 2 : 0);
                if (e && want2) sigma1 *= prime_power_sum(Integer(p), e);
                if (e && want4) sigma3 *= prime_power_sum(Integer(p) * p * p, e);
            }
        }
        auto add = [&](unsigned k, const Integer& term) {
            sums[k] += term;
            if (r) sums[k] += term;
        };
        if (need[1]) add(1, two_squares);
        if (want2) add(2, sigma1);
        if (need[3]) {
            if (n == 0) B = -1;
            else {
                const auto a = mpz_scan1(n.get_mpz_t(), 0);
                Integer scale = negative == s ? 2 : 128;
                scale <<= 2 * a;
                if (chi_minus_four(n >> a) < 0) scale = -scale;
                B *= 4 - scale;
            }
            add(3, B);
        }
        if (want4) add(4, sigma3);
    });
    if (need[1] && negative == s) {
        if (mpz_odd_p(sums[1].get_mpz_t())) fail("odd r3 in the Hurwitz reconstruction");
        sums[1] /= 2;
    }
    sums[2] *= -24;
    sums[4] *= 240;
    // Include the theta constant and the weight-two square correction.
    if (s == 1) { sums[2] -= 10; sums[4] += 2; }
    for (unsigned k = 1; k <= 4; ++k)
        if (need[k]) values.emplace(std::pair{k, k % 2 ? negative : positive}, std::move(sums[k]));
    return values;
}

inline Integer squarefree_part(const Factors& factors) {
    Integer part = 1;
    for (const auto& [p, exponent] : factors)
        if (exponent % 2) part *= p;
    return part;
}

// Write (-1)^k n = D f^2. The multiplier at p^e || f is
// sigma_{2k-1}(p^e) - (D/p) p^(k-1) sigma_{2k-1}(p^(e-1)).
inline Integer half_integer_numerator(unsigned k, const Integer& n) {
    if (n == 0 || (n % 4 != 0 && n % 4 != (k % 2 ? 3 : 1))) return Integer(0);
    const auto& factors = cached_factors(n);
    const auto s = squarefree_part(factors);
    Integer D = s;
    const bool scale_by_four = D % 4 != (k % 2 ? 3UL : 1UL);
    if (scale_by_four) D *= 4;
    const Integer signed_D = k % 2 ? -D : D;
    auto& cache = primitive_values();
    const auto key = std::pair{k, D};
    if (cache.size() >= 2048) cache.clear();
    if (!cache.contains(key)) {
        Integer primitive;
        if (!prefix_cache().get(k, D, primitive)) {
            WeightSet need{};
            need[k] = true;
            primitive = fundamental_values(s, need).at(key);
            prefix_cache().put(k, D, primitive);
        }
        cache.emplace(key, std::move(primitive));
    }
    Integer value = cache.at(key);
    for (const auto& [p, exponent] : factors) {
        const auto e = exponent / 2 - (p == 2 && scale_by_four);
        if (!e) continue;
        Integer power, twist;
        mpz_pow_ui(power.get_mpz_t(), p.get_mpz_t(), 2 * k - 1);
        mpz_pow_ui(twist.get_mpz_t(), p.get_mpz_t(), k - 1);
        Integer factor = 1;
        Integer increment = power - kronecker(signed_D, p) * twist;
        for (unsigned long j = 0; j < e; ++j) {
            factor += increment;
            if (j + 1 < e) increment *= power;
        }
        value *= factor;
    }
    return value;
}

} // namespace direct_arithmetic
