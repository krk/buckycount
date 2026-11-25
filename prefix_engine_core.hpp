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

// Arithmetic shared by the standalone counting and ranking commands.
#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cctype>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_factor.h>
#include <cerrno>
#include <fcntl.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <memory>
#include <system_error>
#include <unistd.h>

namespace {

using Integer = mpz_class;
using Rational = mpq_class;

[[noreturn]] inline void fail(const std::string& message) {
    throw std::runtime_error(message);
}

inline std::string rational_text(const Rational& value) { return value.get_str(); }

inline Rational make_rational(const Integer& numerator, const Integer& denominator) {
    Rational value(numerator, denominator);
    value.canonicalize();
    return value;
}

inline Rational parse_rational(const std::string& text) {
    Rational value;
    if (mpq_set_str(value.get_mpq_t(), text.c_str(), 10) != 0) {
        fail("cannot parse the rational constant " + text);
    }
    value.canonicalize();
    return value;
}

inline unsigned worker_count(unsigned requested) {
    return requested ? requested : std::max(1U, std::thread::hardware_concurrency());
}

inline unsigned long ceil_div(unsigned long n, unsigned long d) {
    return n / d + static_cast<unsigned long>(n % d != 0);
}

template <typename Function>
void parallel_jobs(std::size_t count, unsigned workers, Function job) {
    if (count == 0) return;
    workers = static_cast<unsigned>(std::min<std::size_t>(worker_count(workers), count));
    if (workers == 1) {
        for (std::size_t i = 0; i < count; ++i) job(i);
        return;
    }
    std::atomic<std::size_t> next{0};
    std::exception_ptr error;
    std::mutex error_mutex;
    std::vector<std::thread> threads;
    try {
        for (unsigned worker = 0; worker < workers; ++worker) {
            threads.emplace_back([&] {
                try {
                    for (;;) {
                        const std::size_t i = next.fetch_add(1, std::memory_order_relaxed);
                        if (i >= count) return;
                        job(i);
                    }
                } catch (...) {
                    const std::lock_guard<std::mutex> lock(error_mutex);
                    if (!error) error = std::current_exception();
                    next.store(count, std::memory_order_relaxed);
                }
            });
        }
    } catch (...) {
        next.store(count, std::memory_order_relaxed);
        for (auto& thread : threads) thread.join();
        throw;
    }
    for (auto& thread : threads) thread.join();
    if (error) std::rethrow_exception(error);
}


inline unsigned long residue(const Integer& n, unsigned long modulus) {
    return mpz_fdiv_ui(n.get_mpz_t(), modulus);
}

inline int kronecker(const Integer& d, const Integer& n) {
    return mpz_kronecker(d.get_mpz_t(), n.get_mpz_t());
}

inline unsigned long two_adic_valuation(const Integer& n, Integer& odd) {
    if (n <= 0) fail("expected a positive dyadic argument");
    const auto v = mpz_scan1(n.get_mpz_t(), 0);
    odd = n >> v;
    return v;
}

inline Integer integer_square_root(const Integer& n) {
    Integer root;
    mpz_sqrt(root.get_mpz_t(), n.get_mpz_t());
    return root;
}
inline Integer power(const Integer& n, unsigned long exponent) {
    Integer result;
    mpz_pow_ui(result.get_mpz_t(), n.get_mpz_t(), exponent);
    return result;
}

constexpr std::array<unsigned, 5> cohen_denominators{0, 12, 120, 252, 240};
#include "prefix_engine_flint.hpp"
#include "prefix_cache.hpp"
#include "prefix_engine_direct.hpp"

inline Integer two_character_sum(unsigned weight, long left, long right, const Integer& n) {
    return direct_arithmetic::divisor_sum(weight, left, right, n);
}
inline Rational hurwitz_class_number(const Integer& n) {
    if (n == 0) return Rational(-1, 12);
    return make_rational(direct_arithmetic::half_integer_numerator(1, n), 12);
}
inline Rational cohen_number(unsigned k, const Integer& n) {
    return make_rational(direct_arithmetic::half_integer_numerator(k, n), cohen_denominators.at(k));
}
inline Rational t9_mass(const Integer& m) {
    return hurwitz_class_number(36*m)/6 - hurwitz_class_number(9*m)/3
        - hurwitz_class_number(4*m)/2 + hurwitz_class_number(m);
}
inline Integer character_convolution(long left, long right, const Integer& n) {
    const Integer value = two_character_sum(0, left, right, n);
    return value;
}

inline Integer divisor_count(const Integer& n) {
    const Integer value = two_character_sum(0, 1, 1, n);
    return value;
}

inline Integer divisor_power_sum(const Integer& n, unsigned weight) {
    return two_character_sum(weight, 1, 1, n);
}

inline Integer same_parity_pairs(const Integer& n) {
    if (n % 2 != 0) return divisor_count(n);
    return n % 4 == 0 ? divisor_count(n / 4) : 0;
}

inline unsigned long square_indicator(unsigned long d, const Integer& n) {
    if (n % d != 0) return 0;
    const Integer root = integer_square_root(n / d);
    return root * root == n / d;
}

inline Rational pell_count(unsigned long field, const Integer& n) {
    if (field == 3) return make_rational(
        Integer(character_convolution(1, 12, n)) + Integer(character_convolution(-3, -4, n)), Integer(2));
    if (field == 6) return make_rational(
        Integer(character_convolution(1, 24, n)) + Integer(character_convolution(-8, -3, n)), Integer(2));
    fail("unknown Pell character family");
}

inline Integer twisted_s(const Integer& n) {
    return two_character_sum(1, 1, -3, n);
}
inline Integer twisted_t(const Integer& n) {
    return two_character_sum(1, -3, 1, n);
}

inline Rational t27_mass(const Integer& m) {
    if (m % 3 == 0) return t9_mass(m / 3);
    if (m % 3 == 1) {
        return t9_mass(3 * m) / 2;
    }
    return Rational();
}

inline Integer restricted_divisor_count(unsigned long modulus, unsigned long residue,
                                       const Integer& n) {
    if (n == 0) return 0;
    if (modulus == 0 || residue >= modulus)
        fail("unprepared restricted-divisor argument");
    const auto tau = [&n](unsigned long d) -> Integer {
        return n % d == 0 ? Integer(divisor_count(n / d)) : Integer(0);
    };
    const auto same = [](const Integer& m) -> Integer {
        if (m % 2 != 0) return Integer(divisor_count(m));
        return m % 4 == 0 ? Integer(divisor_count(m / 4)) : Integer(0);
    };
    Integer value;
    if (residue == 0) value = tau(modulus);
    else if (modulus == 4 && residue == 2) {
        // (chi_4*chi_4)(v) is tau(v) for odd v, and zero for even v.
        value = -tau(2) + 2 * tau(4);
        if (n % 2 == 0 && (n / 2) % 2 != 0) value += 2 * tau(2);
    } else if (modulus == 6 && residue == 3) value = tau(3) - 2 * tau(6) + 2 * tau(12);
    else if (modulus == 6 && residue == 1) {
        if (n % 3 == 1) value += same(n);
        if (n % 9 == 0) value += same(n / 9);
    } else if (modulus == 6 && residue == 4) {
        if (n % 3 == 1) value += tau(2);
        if (n % 9 == 0) value += tau(18);
    } else fail("unsupported restricted-divisor family");
    if (sgn(value) < 0)
        fail("restricted divisor count exceeds its exact word range");
    return value;
}


} // namespace
