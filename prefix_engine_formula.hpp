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

// The explicit mixed-weight Grand Formula shared by the theta/Hecke prefix
// engines: named arithmetic functions, the three literal expansion rows, the
// five rotary rows, the oriented count, and the complete sector assembly with
// its achiral, integrality and Burnside checks.
//
// The literal rows and normalization implement Paper I, Appendix A. The public
// coefficient verifier compares them with the printed TeX and structured export;
// the chamber-return verifier independently checks the split Hurwitz block.
#pragma once

#include "prefix_engine_core.hpp"

namespace {

enum class Kind {
    divisor_sum,      // sigma_w(m/j) or E_w[psi,phi](m/j)
    cohen,            // H_k(p m/q), optionally twisted by chi_c on the argument
    hurwitz,          // H(p m/q), optionally restricted to a dyadic class of m
    convolution,      // tau(m/j) or (chi_a*chi_b)(m/j)
    square_indicator, // I_d(m)
    restricted_divisors,  // d_{k,r}(m)
    pell,             // N_3^+(m/j) or N_6^+(m/j)
    same_parity,      // d_same(m/j)
    t9,               // T_9(m)
    t27,              // T_27(m)
};

struct NamedFunction {
    std::string name;
    std::size_t slot = 0;
    Kind kind = Kind::divisor_sum;
    unsigned weight = 0;
    long left = 1;
    long right = 1;
    unsigned long numerator = 1;    // the argument is numerator * m / denominator
    unsigned long denominator = 1;
    int parity = -1;                // v_2(m) mod 2 restriction, or -1
    int unit = -1;                  // odd part of m mod 8 restriction, or -1
    long twist = 1;                 // Kronecker twist on a Cohen argument
    unsigned long parameter = 0;    // d for I_d, k for d_{k,r}, field for N^+
    unsigned long residue = 0;      // r for d_{k,r}
};

struct Term {
    const char* name;
    const char* coefficient;
};

struct ParsedTerm {
    NamedFunction function;
    Rational coefficient;
};

struct ParsedRow {
    std::vector<ParsedTerm> terms;
};

// Parse "m", "m/q", "pm" or "pm/q".
inline void parse_argument(const std::string& text, NamedFunction& function, const std::string& whole) {
    std::size_t position = 0;
    unsigned long p = 1;
    if (position < text.size() && std::isdigit(static_cast<unsigned char>(text[position]))) {
        p = 0;
        while (position < text.size() && std::isdigit(static_cast<unsigned char>(text[position]))) {
            p = 10 * p + static_cast<unsigned long>(text[position] - '0');
            ++position;
        }
    }
    if (position >= text.size() || text[position] != 'm') fail("cannot parse the argument of " + whole);
    ++position;
    unsigned long q = 1;
    if (position < text.size()) {
        if (text[position] != '/') fail("cannot parse the argument of " + whole);
        ++position;
        q = 0;
        while (position < text.size() && std::isdigit(static_cast<unsigned char>(text[position]))) {
            q = 10 * q + static_cast<unsigned long>(text[position] - '0');
            ++position;
        }
        if (position != text.size() || q == 0) fail("cannot parse the argument of " + whole);
    }
    if (p == 0) fail("cannot parse the argument of " + whole);
    function.numerator = p;
    function.denominator = q;
}

inline long parse_signed(const std::string& text, const std::string& whole) {
    if (text.empty()) fail("an empty integer in " + whole);
    std::size_t position = 0;
    bool negative = false;
    if (text[position] == '-') {
        negative = true;
        ++position;
    }
    if (position >= text.size()) fail("cannot parse an integer in " + whole);
    long value = 0;
    for (; position < text.size(); ++position) {
        if (!std::isdigit(static_cast<unsigned char>(text[position]))) {
            fail("cannot parse an integer in " + whole);
        }
        if (value > (std::numeric_limits<long>::max() - 9) / 10) fail("an integer overflowed in " + whole);
        value = 10 * value + (text[position] - '0');
    }
    return negative ? -value : value;
}

// The argument text between the parentheses that follow `head`.
inline std::string parenthesized(const std::string& name, std::size_t start, std::size_t& end) {
    if (start >= name.size() || name[start] != '(') fail("cannot parse the named function " + name);
    end = name.find(')', start);
    if (end == std::string::npos) fail("cannot parse the named function " + name);
    return name.substr(start + 1, end - start - 1);
}

inline NamedFunction parse_named_function(const std::string& name) {
    NamedFunction function;
    function.name = name;
    std::size_t end = 0;
    if (name.rfind("sigma_", 0) == 0) {
        const std::size_t open = name.find('(');
        function.kind = Kind::divisor_sum;
        function.weight = static_cast<unsigned>(parse_signed(name.substr(6, open - 6), name));
        parse_argument(parenthesized(name, open, end), function, name);
    } else if (name.rfind("E_", 0) == 0) {
        const std::size_t bracket = name.find('[');
        const std::size_t comma = name.find(',', bracket);
        const std::size_t close = name.find(']', comma);
        function.kind = Kind::divisor_sum;
        function.weight = static_cast<unsigned>(parse_signed(name.substr(2, bracket - 2), name));
        function.left = parse_signed(name.substr(bracket + 1, comma - bracket - 1), name);
        function.right = parse_signed(name.substr(comma + 1, close - comma - 1), name);
        parse_argument(parenthesized(name, close + 1, end), function, name);
    } else if (name.rfind("H_", 0) == 0) {
        const std::size_t open = name.find('(');
        function.kind = Kind::cohen;
        function.weight = static_cast<unsigned>(parse_signed(name.substr(2, open - 2), name));
        parse_argument(parenthesized(name, open, end), function, name);
    } else if (name.rfind("chi_", 0) == 0) {
        // chi_c(N) H_2(N): the twist acts on the full Cohen argument N.
        const std::size_t open = name.find('(');
        function.kind = Kind::cohen;
        function.twist = parse_signed(name.substr(4, open - 4), name);
        const std::string twist_argument = parenthesized(name, open, end);
        const std::string rest = name.substr(end + 1);
        if (rest.rfind(" H_", 0) != 0) fail("cannot parse the twisted Cohen term " + name);
        const std::size_t open_again = rest.find('(');
        function.weight = static_cast<unsigned>(parse_signed(rest.substr(3, open_again - 3), name));
        std::size_t end_again = 0;
        const std::string cohen_argument = parenthesized(rest, open_again, end_again);
        if (cohen_argument != twist_argument) fail("a twisted Cohen term has two arguments: " + name);
        parse_argument(cohen_argument, function, name);
        end = name.size() - 1;
    } else if (name.rfind("H(", 0) == 0) {
        function.kind = Kind::hurwitz;
        parse_argument(parenthesized(name, 1, end), function, name);
        const std::string rest = name.substr(end + 1);
        if (!rest.empty()) {
            // " [2:v%2=P,u=U|3:all]": a dyadic square class of the index m.
            if (rest.rfind(" [2:v%2=", 0) != 0 || rest.back() != ']') {
                fail("cannot parse the class restriction of " + name);
            }
            const std::size_t comma = rest.find(",u=");
            const std::size_t bar = rest.find("|3:all]");
            if (comma == std::string::npos || bar == std::string::npos) {
                fail("cannot parse the class restriction of " + name);
            }
            function.parity = static_cast<int>(parse_signed(rest.substr(8, comma - 8), name));
            function.unit = static_cast<int>(parse_signed(rest.substr(comma + 3, bar - comma - 3), name));
            if (function.parity < 0 || function.parity > 1 || function.unit % 2 == 0 ||
                function.unit < 1 || function.unit > 7) {
                fail("an impossible dyadic class in " + name);
            }
        }
    } else if (name.rfind("tau(", 0) == 0) {
        function.kind = Kind::convolution;
        parse_argument(parenthesized(name, 3, end), function, name);
    } else if (name.rfind("(chi_", 0) == 0) {
        const std::size_t star = name.find("*chi_");
        const std::size_t close = name.find(')', star);
        function.kind = Kind::convolution;
        function.left = parse_signed(name.substr(5, star - 5), name);
        function.right = parse_signed(name.substr(star + 5, close - star - 5), name);
        parse_argument(parenthesized(name, close + 1, end), function, name);
    } else if (name.rfind("I_", 0) == 0) {
        const std::size_t open = name.find('(');
        function.kind = Kind::square_indicator;
        function.parameter = static_cast<unsigned long>(parse_signed(name.substr(2, open - 2), name));
        parse_argument(parenthesized(name, open, end), function, name);
    } else if (name.rfind("d_{", 0) == 0) {
        const std::size_t comma = name.find(',');
        const std::size_t close = name.find('}');
        function.kind = Kind::restricted_divisors;
        function.parameter = static_cast<unsigned long>(parse_signed(name.substr(3, comma - 3), name));
        function.residue = static_cast<unsigned long>(parse_signed(name.substr(comma + 1, close - comma - 1), name));
        parse_argument(parenthesized(name, close + 1, end), function, name);
    } else if (name.rfind("N_", 0) == 0) {
        const std::size_t caret = name.find("^+");
        function.kind = Kind::pell;
        function.parameter = static_cast<unsigned long>(parse_signed(name.substr(2, caret - 2), name));
        parse_argument(parenthesized(name, caret + 2, end), function, name);
    } else if (name.rfind("d_same(", 0) == 0) {
        function.kind = Kind::same_parity;
        parse_argument(parenthesized(name, 6, end), function, name);
    } else if (name == "T_9(m)") {
        function.kind = Kind::t9;
    } else if (name == "T_27(m)") {
        function.kind = Kind::t27;
    } else {
        fail("an unknown named arithmetic function " + name);
    }
    return function;
}

// The argument p m / q of a named function, or nothing when q does not divide p m.
inline std::optional<Integer> scaled_argument(const NamedFunction& function, const Integer& m) {
    const Integer product = function.numerator * m;
    if (product % function.denominator != 0) return std::nullopt;
    return product / function.denominator;
}

inline bool hurwitz_class_applies(const NamedFunction& function, const Integer& m) {
    if (function.parity < 0) return true;
    Integer odd_part;
    const auto valuation = two_adic_valuation(m, odd_part);
    return static_cast<int>(valuation % 2) == function.parity &&
           static_cast<int>(residue(odd_part, 8)) == function.unit;
}

inline Rational evaluate(const NamedFunction& function, const Integer& m) {
    switch (function.kind) {
        case Kind::divisor_sum: {
            const auto n = scaled_argument(function, m);
            if (!n) return Rational();
            return Rational(two_character_sum(function.weight, function.left, function.right, *n));
        }
        case Kind::cohen: {
            const auto n = scaled_argument(function, m);
            if (!n) return Rational();
            Rational value = cohen_number(function.weight, *n);
            if (function.twist != 1 && sgn(value) != 0) value *= Rational(kronecker(function.twist, *n));
            return value;
        }
        case Kind::hurwitz: {
            if (!hurwitz_class_applies(function, m)) return Rational();
            const auto n = scaled_argument(function, m);
            if (!n) return Rational();
            return hurwitz_class_number(*n);
        }
        case Kind::convolution: {
            const auto n = scaled_argument(function, m);
            if (!n) return Rational();
            return Rational(character_convolution(function.left, function.right, *n));
        }
        case Kind::square_indicator:
            return Rational(square_indicator(function.parameter, m));
        case Kind::restricted_divisors:
            return Rational(restricted_divisor_count(function.parameter, function.residue, m));
        case Kind::pell: {
            const auto n = scaled_argument(function, m);
            if (!n) return Rational();
            return pell_count(function.parameter, *n);
        }
        case Kind::same_parity: {
            const auto n = scaled_argument(function, m);
            if (!n) return Rational();
            return Rational(same_parity_pairs(*n));
        }
        case Kind::t9:
            return t9_mass(m);
        case Kind::t27:
            return t27_mass(m);
    }
    fail("an unknown named function kind");
}

// ---------------------------------------------------------------------------
// The three recorded expansions of the explicit mixed-weight Grand Formula.
// Each line is one named function and its rational coefficient, grouped by
// the weight of the signature-(1,q) Siegel series it belongs to.

// M_split^open(m): the split involutive row, 171 terms.
constexpr Term split_row[] = {
    // weight 5: 3 terms
    {"sigma_4(m)", "697/622080"},
    {"sigma_4(m/2)", "-697/311040"},
    {"sigma_4(m/4)", "697/19440"},
    // weight 9/2: 6 terms
    {"H_4(4m)", "-41/5184"},
    {"H_4(m)", "-41/324"},
    {"H_4(12m)", "-1/10368"},
    {"H_4(3m)", "-1/648"},
    {"H_4(4m/3)", "-1/128"},
    {"H_4(m/3)", "-1/8"},
    // weight 4: 12 terms
    {"E_3[-4,1](m)", "1/24"},
    {"E_3[1,-4](m)", "7/2592"},
    {"E_3[-8,1](m)", "1/81"},
    {"E_3[1,-8](m)", "1/5184"},
    {"E_3[-3,1](m)", "5/288"},
    {"E_3[-3,1](m/2)", "-7/576"},
    {"E_3[-3,1](m/4)", "17/36"},
    {"E_3[1,-3](m)", "17/15552"},
    {"E_3[1,-3](m/2)", "-7/15552"},
    {"E_3[1,-3](m/4)", "5/486"},
    {"E_3[-4,1](m/3)", "1/24"},
    {"E_3[1,-4](m/3)", "1/192"},
    // weight 7/2: 14 terms
    {"H_3(4m)", "-13/144"},
    {"H_3(m)", "-131/72"},
    {"H_3(8m)", "7/288"},
    {"H_3(m/2)", "-14/9"},
    {"H_3(12m)", "1/216"},
    {"H_3(3m)", "13/216"},
    {"H_3(4m/3)", "-3/32"},
    {"H_3(m/3)", "-15/8"},
    {"H_3(24m)", "1/1152"},
    {"H_3(3m/2)", "-1/18"},
    {"H_3(8m/3)", "3/128"},
    {"H_3(m/6)", "-3/2"},
    {"H_3(36m)", "1/864"},
    {"H_3(9m)", "1/216"},
    // weight 3: 18 terms
    {"sigma_2(m)", "5393/20736"},
    {"sigma_2(m/4)", "-1933/2592"},
    {"sigma_2(m/12)", "-73/36"},
    {"sigma_2(m/2)", "-1415/10368"},
    {"sigma_2(m/3)", "-37/288"},
    {"sigma_2(m/6)", "23/72"},
    {"E_2[1,8](m)", "-5/576"},
    {"E_2[1,8](m/3)", "-1/64"},
    {"E_2[8,1](m)", "1/4"},
    {"E_2[8,1](m/3)", "1/4"},
    {"E_2[-3,-4](m)", "1/48"},
    {"E_2[1,12](m)", "-1/864"},
    {"E_2[12,1](m)", "1/24"},
    {"sigma_2(m/8)", "-5/6"},
    {"E_2[-3,-8](m)", "-1/96"},
    {"E_2[-8,-3](m)", "1/54"},
    {"E_2[1,24](m)", "-1/1728"},
    {"E_2[24,1](m)", "1/12"},
    // weight 5/2: 22 terms
    {"H_2(4m)", "913/576"},
    {"H_2(m)", "-367/72"},
    {"H_2(8m)", "-1/12"},
    {"H_2(m/2)", "-17/3"},
    {"H_2(12m)", "359/3456"},
    {"H_2(3m)", "-161/432"},
    {"H_2(4m/3)", "189/128"},
    {"H_2(m/3)", "-75/16"},
    {"H_2(m/4)", "-10/3"},
    {"chi_8(m) H_2(m)", "-5/6"},
    {"H_2(24m)", "7/288"},
    {"H_2(8m/3)", "-3/32"},
    {"H_2(3m/2)", "4/9"},
    {"H_2(m/6)", "-6"},
    {"H_2(36m)", "-1/216"},
    {"H_2(9m)", "1/27"},
    {"H_2(3m/4)", "-1/3"},
    {"H_2(m/12)", "-3"},
    {"chi_8(3m) H_2(3m)", "-1/12"},
    {"chi_8(m/3) H_2(m/3)", "-3/4"},
    {"H_2(72m)", "1/72"},
    {"H_2(9m/2)", "1/9"},
    // weight 2: 23 terms
    {"E_1[-4,1](m)", "19/24"},
    {"E_1[1,-4](m)", "-215/864"},
    {"E_1[-8,1](m)", "35/54"},
    {"E_1[1,-8](m)", "-263/1728"},
    {"E_1[-3,1](m)", "13/288"},
    {"E_1[-3,1](m/4)", "-295/288"},
    {"E_1[1,-3](m)", "2039/15552"},
    {"E_1[1,-3](m/4)", "431/7776"},
    {"E_1[-3,1](m/2)", "439/576"},
    {"E_1[1,-3](m/2)", "-4591/15552"},
    {"E_1[1,-4](m/2)", "-1/6"},
    {"E_1[-4,1](m/2)", "1/4"},
    {"E_1[-24,1](m)", "1/6"},
    {"E_1[-3,8](m)", "1/48"},
    {"E_1[1,-24](m)", "-1/72"},
    {"E_1[-4,1](m/3)", "19/24"},
    {"E_1[1,-4](m/3)", "-199/576"},
    {"E_1[-3,1](m/8)", "3/2"},
    {"E_1[1,-3](m/8)", "-1/3"},
    {"E_1[-8,1](m/3)", "1/18"},
    {"E_1[1,-8](m/3)", "-13/144"},
    {"E_1[-4,1](m/6)", "1/4"},
    {"E_1[1,-4](m/6)", "-1/4"},
    // weight 3/2: 55 terms
    {"H(4m)", "-31/216"},
    {"H(m)", "-217/216"},
    {"H(2m)", "-3/16"},
    {"H(8m)", "-43/96"},
    {"H(m/2)", "-7/3"},
    {"H(12m)", "-335/2592"},
    {"H(3m)", "137/1296"},
    {"H(4m/3)", "-31/216"},
    {"H(m/3)", "-217/216"},
    {"H(m/4) [2:v%2=0,u=1|3:all]", "-3/4"},
    {"H(16m) [2:v%2=0,u=3|3:all]", "9/16"},
    {"H(4m) [2:v%2=0,u=3|3:all]", "-21/16"},
    {"H(m/4) [2:v%2=0,u=5|3:all]", "-3/4"},
    {"H(16m) [2:v%2=0,u=7|3:all]", "-3/32"},
    {"H(m/4) [2:v%2=1,u=1|3:all]", "-3/4"},
    {"H(m/4) [2:v%2=1,u=3|3:all]", "-3/4"},
    {"H(m/4) [2:v%2=1,u=5|3:all]", "-3/4"},
    {"H(m/4) [2:v%2=1,u=7|3:all]", "-3/4"},
    {"H(24m)", "-101/1152"},
    {"H(m/6)", "-11/6"},
    {"H(8m/3)", "-133/384"},
    {"H(3m/2)", "-7/18"},
    {"H(6m)", "-3/64"},
    {"H(2m/3)", "-9/64"},
    {"H(36m)", "-335/2592"},
    {"H(9m)", "137/1296"},
    {"H(48m) [2:v%2=0,u=1|3:all]", "-1/16"},
    {"H(12m) [2:v%2=0,u=1|3:all]", "7/48"},
    {"H(16m/3) [2:v%2=0,u=1|3:all]", "9/16"},
    {"H(4m/3) [2:v%2=0,u=1|3:all]", "-21/16"},
    {"H(m/12) [2:v%2=0,u=3|3:all]", "-3/4"},
    {"H(3m/4) [2:v%2=0,u=3|3:all]", "1/12"},
    {"H(48m) [2:v%2=0,u=5|3:all]", "1/96"},
    {"H(16m/3) [2:v%2=0,u=5|3:all]", "-3/32"},
    {"H(m/12) [2:v%2=0,u=7|3:all]", "-3/4"},
    {"H(3m/4) [2:v%2=0,u=7|3:all]", "1/12"},
    {"H(m/12) [2:v%2=1,u=1|3:all]", "-3/4"},
    {"H(3m/4) [2:v%2=1,u=1|3:all]", "1/12"},
    {"H(m/12) [2:v%2=1,u=3|3:all]", "-3/4"},
    {"H(3m/4) [2:v%2=1,u=3|3:all]", "1/12"},
    {"H(m/12) [2:v%2=1,u=5|3:all]", "-3/4"},
    {"H(3m/4) [2:v%2=1,u=5|3:all]", "1/12"},
    {"H(m/12) [2:v%2=1,u=7|3:all]", "-3/4"},
    {"H(3m/4) [2:v%2=1,u=7|3:all]", "1/12"},
    {"H(72m)", "1/72"},
    {"H(9m/2)", "1/9"},
    {"H(9m/4) [2:v%2=0,u=1|3:all]", "1/12"},
    {"H(144m) [2:v%2=0,u=3|3:all]", "-1/16"},
    {"H(36m) [2:v%2=0,u=3|3:all]", "7/48"},
    {"H(9m/4) [2:v%2=0,u=5|3:all]", "1/12"},
    {"H(144m) [2:v%2=0,u=7|3:all]", "1/96"},
    {"H(9m/4) [2:v%2=1,u=1|3:all]", "1/12"},
    {"H(9m/4) [2:v%2=1,u=3|3:all]", "1/12"},
    {"H(9m/4) [2:v%2=1,u=5|3:all]", "1/12"},
    {"H(9m/4) [2:v%2=1,u=7|3:all]", "1/12"},
    // weight 1: 14 terms
    {"tau(m/2)", "-259931/622080"},
    {"tau(m/6)", "-1217/7776"},
    {"tau(m/4)", "286633/311040"},
    {"(chi_1*chi_12)(m/2)", "1/24"},
    {"(chi_1*chi_12)(m)", "155/2592"},
    {"(chi_1*chi_24)(m)", "251/5184"},
    {"(chi_-8*chi_-3)(m)", "157/2592"},
    {"tau(m/12)", "635/1944"},
    {"(chi_1*chi_8)(m)", "85/576"},
    {"tau(m)", "144073/622080"},
    {"(chi_4*chi_4)(m/2)", "17/48"},
    {"(chi_1*chi_8)(m/3)", "85/576"},
    {"tau(m/3)", "865/7776"},
    {"(chi_4*chi_4)(m/6)", "1/12"},
    // weight 1/2: 4 terms
    {"I_1(m)", "-74233/622080"},
    {"I_2(m)", "-23/288"},
    {"I_3(m)", "-74233/622080"},
    {"I_6(m)", "-23/288"},
};

// M_antipodal^open(m): the antipodal involutive row, 13 terms.
constexpr Term antipodal_row[] = {
    // weight 5: 1 terms
    {"sigma_4(m/2)", "41/311040"},
    // weight 4: 2 terms
    {"E_3[-3,1](m/2)", "-1/576"},
    {"E_3[1,-3](m/2)", "-1/15552"},
    // weight 3: 3 terms
    {"sigma_2(m/4)", "-1/96"},
    {"sigma_2(m/2)", "85/10368"},
    {"sigma_2(m/6)", "-1/72"},
    // weight 2: 4 terms
    {"E_1[-3,1](m/4)", "-1/32"},
    {"E_1[1,-3](m/4)", "1/96"},
    {"E_1[-3,1](m/2)", "-11/576"},
    {"E_1[1,-3](m/2)", "53/15552"},
    // weight 1: 3 terms
    {"tau(m/2)", "2849/311040"},
    {"tau(m/6)", "7/1944"},
    {"tau(m/4)", "1/256"},
};

// a_m: the whole achiral correction as one expansion, 181 terms.  The lower
// sectors enter through T_27(m), N_3^+(m/3), d_{6,1}, d_{6,4}, tau(m/10),
// d_same(m/5) and the square indicators I_5, I_9, I_15.
constexpr Term achiral_row[] = {
    // weight 5: 3 terms
    {"sigma_4(m)", "697/622080"},
    {"sigma_4(m/2)", "-41/19440"},
    {"sigma_4(m/4)", "697/19440"},
    // weight 9/2: 6 terms
    {"H_4(4m)", "-41/5184"},
    {"H_4(m)", "-41/324"},
    {"H_4(12m)", "-1/10368"},
    {"H_4(3m)", "-1/648"},
    {"H_4(4m/3)", "-1/128"},
    {"H_4(m/3)", "-1/8"},
    // weight 4: 12 terms
    {"E_3[-4,1](m)", "1/24"},
    {"E_3[1,-4](m)", "7/2592"},
    {"E_3[-8,1](m)", "1/81"},
    {"E_3[1,-8](m)", "1/5184"},
    {"E_3[-3,1](m)", "5/288"},
    {"E_3[-3,1](m/2)", "-1/72"},
    {"E_3[-3,1](m/4)", "17/36"},
    {"E_3[1,-3](m)", "17/15552"},
    {"E_3[1,-3](m/2)", "-1/1944"},
    {"E_3[1,-3](m/4)", "5/486"},
    {"E_3[-4,1](m/3)", "1/24"},
    {"E_3[1,-4](m/3)", "1/192"},
    // weight 7/2: 14 terms
    {"H_3(4m)", "-13/144"},
    {"H_3(m)", "-131/72"},
    {"H_3(8m)", "7/288"},
    {"H_3(m/2)", "-14/9"},
    {"H_3(12m)", "1/216"},
    {"H_3(3m)", "13/216"},
    {"H_3(4m/3)", "-3/32"},
    {"H_3(m/3)", "-15/8"},
    {"H_3(24m)", "1/1152"},
    {"H_3(3m/2)", "-1/18"},
    {"H_3(8m/3)", "3/128"},
    {"H_3(m/6)", "-3/2"},
    {"H_3(36m)", "1/864"},
    {"H_3(9m)", "1/216"},
    // weight 3: 18 terms
    {"sigma_2(m)", "5393/20736"},
    {"sigma_2(m/4)", "-535/648"},
    {"sigma_2(m/12)", "-73/36"},
    {"sigma_2(m/2)", "-557/5184"},
    {"sigma_2(m/3)", "-37/288"},
    {"sigma_2(m/6)", "11/36"},
    {"E_2[1,8](m)", "-5/576"},
    {"E_2[1,8](m/3)", "-1/64"},
    {"E_2[8,1](m)", "1/4"},
    {"E_2[8,1](m/3)", "1/4"},
    {"E_2[-3,-4](m)", "1/48"},
    {"E_2[1,12](m)", "-1/864"},
    {"E_2[12,1](m)", "1/24"},
    {"sigma_2(m/8)", "-5/6"},
    {"E_2[-3,-8](m)", "-1/96"},
    {"E_2[-8,-3](m)", "1/54"},
    {"E_2[1,24](m)", "-1/1728"},
    {"E_2[24,1](m)", "1/12"},
    // weight 5/2: 22 terms
    {"H_2(4m)", "913/576"},
    {"H_2(m)", "-367/72"},
    {"H_2(8m)", "-1/12"},
    {"H_2(m/2)", "-17/3"},
    {"H_2(12m)", "359/3456"},
    {"H_2(3m)", "-161/432"},
    {"H_2(4m/3)", "189/128"},
    {"H_2(m/3)", "-75/16"},
    {"H_2(m/4)", "-10/3"},
    {"chi_8(m) H_2(m)", "-5/6"},
    {"H_2(24m)", "7/288"},
    {"H_2(8m/3)", "-3/32"},
    {"H_2(3m/2)", "4/9"},
    {"H_2(m/6)", "-6"},
    {"H_2(36m)", "-1/216"},
    {"H_2(9m)", "1/27"},
    {"H_2(3m/4)", "-1/3"},
    {"H_2(m/12)", "-3"},
    {"chi_8(3m) H_2(3m)", "-1/12"},
    {"chi_8(m/3) H_2(m/3)", "-3/4"},
    {"H_2(72m)", "1/72"},
    {"H_2(9m/2)", "1/9"},
    // weight 2: 23 terms
    {"E_1[-4,1](m)", "19/24"},
    {"E_1[1,-4](m)", "-215/864"},
    {"E_1[-8,1](m)", "35/54"},
    {"E_1[1,-8](m)", "-263/1728"},
    {"E_1[-3,1](m)", "7/96"},
    {"E_1[-3,1](m/4)", "-7/6"},
    {"E_1[1,-3](m)", "2471/15552"},
    {"E_1[1,-3](m/4)", "5/486"},
    {"E_1[-3,1](m/2)", "35/48"},
    {"E_1[1,-3](m/2)", "-2809/7776"},
    {"E_1[1,-4](m/2)", "-1/6"},
    {"E_1[-4,1](m/2)", "1/4"},
    {"E_1[-24,1](m)", "1/6"},
    {"E_1[-3,8](m)", "1/48"},
    {"E_1[1,-24](m)", "-1/72"},
    {"E_1[-4,1](m/3)", "19/24"},
    {"E_1[1,-4](m/3)", "-199/576"},
    {"E_1[-3,1](m/8)", "3/2"},
    {"E_1[1,-3](m/8)", "-1/3"},
    {"E_1[-8,1](m/3)", "1/18"},
    {"E_1[1,-8](m/3)", "-13/144"},
    {"E_1[-4,1](m/6)", "1/4"},
    {"E_1[1,-4](m/6)", "-1/4"},
    // weight 3/2: 56 terms
    {"H(4m)", "77/216"},
    {"H(m)", "-433/216"},
    {"H(2m)", "-3/16"},
    {"H(8m)", "-43/96"},
    {"H(m/2)", "-7/3"},
    {"H(12m)", "-335/2592"},
    {"H(3m)", "137/1296"},
    {"H(4m/3)", "-31/216"},
    {"H(m/3)", "-217/216"},
    {"H(m/4) [2:v%2=0,u=1|3:all]", "-3/4"},
    {"H(16m) [2:v%2=0,u=3|3:all]", "9/16"},
    {"H(4m) [2:v%2=0,u=3|3:all]", "-21/16"},
    {"H(m/4) [2:v%2=0,u=5|3:all]", "-3/4"},
    {"H(16m) [2:v%2=0,u=7|3:all]", "-3/32"},
    {"H(m/4) [2:v%2=1,u=1|3:all]", "-3/4"},
    {"H(m/4) [2:v%2=1,u=3|3:all]", "-3/4"},
    {"H(m/4) [2:v%2=1,u=5|3:all]", "-3/4"},
    {"H(m/4) [2:v%2=1,u=7|3:all]", "-3/4"},
    {"H(24m)", "-101/1152"},
    {"H(m/6)", "-11/6"},
    {"H(8m/3)", "-133/384"},
    {"H(3m/2)", "-7/18"},
    {"H(6m)", "-3/64"},
    {"H(2m/3)", "-9/64"},
    {"H(36m)", "-767/2592"},
    {"H(9m)", "569/1296"},
    {"H(48m) [2:v%2=0,u=1|3:all]", "-1/16"},
    {"H(12m) [2:v%2=0,u=1|3:all]", "7/48"},
    {"H(16m/3) [2:v%2=0,u=1|3:all]", "9/16"},
    {"H(4m/3) [2:v%2=0,u=1|3:all]", "-21/16"},
    {"H(m/12) [2:v%2=0,u=3|3:all]", "-3/4"},
    {"H(3m/4) [2:v%2=0,u=3|3:all]", "1/12"},
    {"H(48m) [2:v%2=0,u=5|3:all]", "1/96"},
    {"H(16m/3) [2:v%2=0,u=5|3:all]", "-3/32"},
    {"H(m/12) [2:v%2=0,u=7|3:all]", "-3/4"},
    {"H(3m/4) [2:v%2=0,u=7|3:all]", "1/12"},
    {"H(m/12) [2:v%2=1,u=1|3:all]", "-3/4"},
    {"H(3m/4) [2:v%2=1,u=1|3:all]", "1/12"},
    {"H(m/12) [2:v%2=1,u=3|3:all]", "-3/4"},
    {"H(3m/4) [2:v%2=1,u=3|3:all]", "1/12"},
    {"H(m/12) [2:v%2=1,u=5|3:all]", "-3/4"},
    {"H(3m/4) [2:v%2=1,u=5|3:all]", "1/12"},
    {"H(m/12) [2:v%2=1,u=7|3:all]", "-3/4"},
    {"H(3m/4) [2:v%2=1,u=7|3:all]", "1/12"},
    {"H(72m)", "1/72"},
    {"H(9m/2)", "1/9"},
    {"H(9m/4) [2:v%2=0,u=1|3:all]", "1/12"},
    {"H(144m) [2:v%2=0,u=3|3:all]", "-1/16"},
    {"H(36m) [2:v%2=0,u=3|3:all]", "7/48"},
    {"H(9m/4) [2:v%2=0,u=5|3:all]", "1/12"},
    {"H(144m) [2:v%2=0,u=7|3:all]", "1/96"},
    {"H(9m/4) [2:v%2=1,u=1|3:all]", "1/12"},
    {"H(9m/4) [2:v%2=1,u=3|3:all]", "1/12"},
    {"H(9m/4) [2:v%2=1,u=5|3:all]", "1/12"},
    {"H(9m/4) [2:v%2=1,u=7|3:all]", "1/12"},
    {"T_27(m)", "-1"},
    // weight 1: 20 terms
    {"tau(m/2)", "-293113/622080"},
    {"tau(m/6)", "-325/7776"},
    {"tau(m/4)", "43811/38880"},
    {"(chi_1*chi_12)(m/2)", "1/24"},
    {"(chi_1*chi_12)(m)", "371/2592"},
    {"(chi_-4*chi_-3)(m)", "1/12"},
    {"(chi_1*chi_24)(m)", "683/5184"},
    {"(chi_-8*chi_-3)(m)", "373/2592"},
    {"tau(m/12)", "905/1944"},
    {"(chi_1*chi_8)(m)", "85/576"},
    {"tau(m)", "144073/622080"},
    {"(chi_4*chi_4)(m/2)", "23/48"},
    {"(chi_1*chi_8)(m/3)", "85/576"},
    {"tau(m/3)", "1/7776"},
    {"(chi_4*chi_4)(m/6)", "1/12"},
    {"d_{6,4}(m)", "-1/18"},
    {"N_3^+(m/3)", "1/6"},
    {"d_{6,1}(m)", "-1/18"},
    {"tau(m/10)", "2/5"},
    {"d_same(m/5)", "2/5"},
    // weight 1/2: 7 terms
    {"I_1(m)", "-143353/622080"},
    {"I_2(m)", "-23/288"},
    {"I_3(m)", "-108793/622080"},
    {"I_6(m)", "-23/288"},
    {"I_9(m)", "1/18"},
    {"I_5(m)", "-2/5"},
    {"I_15(m)", "-2/5"},
};

template <std::size_t Count>
ParsedRow parse_row(const Term (&row)[Count], std::map<std::string, std::size_t>& function_slots) {
    ParsedRow result;
    result.terms.reserve(Count);
    for (const Term& term : row) {
        result.terms.push_back({parse_named_function(term.name), parse_rational(term.coefficient)});
        const auto [position, inserted] = function_slots.try_emplace(term.name, function_slots.size());
        static_cast<void>(inserted);
        auto& parsed = result.terms.back();
        parsed.function.slot = position->second;
        if (sgn(parsed.coefficient) == 0) fail(std::string("a zero coefficient on ") + term.name);
    }
    return result;
}

struct Expansions {
    ParsedRow split;
    ParsedRow antipodal;
    ParsedRow achiral;
    std::size_t slot_count;
};

inline const Expansions& expansions() {
    static const Expansions parsed = [] {
        std::map<std::string, std::size_t> slots;
        return Expansions{parse_row(split_row, slots), parse_row(antipodal_row, slots),
                          parse_row(achiral_row, slots), slots.size()};
    }();
    return parsed;
}

// Group the half-integral terms used by this index before evaluating the rows.
inline void prepare_shared_weights(const Integer& m) {
    using namespace direct_arithmetic;
    static const auto functions = [] {
        std::vector<NamedFunction> list;
        std::set<std::size_t> seen;
        const auto& rows = expansions();
        for (const auto* row : {&rows.split, &rows.antipodal, &rows.achiral})
            for (const auto& term : row->terms) {
                const auto& f = term.function;
                if ((f.kind == Kind::hurwitz || f.kind == Kind::cohen) && seen.insert(f.slot).second)
                    list.push_back(f);
            }
        return list;
    }();
    std::map<Integer, WeightSet> requested;
    for (const auto& f : functions) {
        if (!hurwitz_class_applies(f, m)) continue;
        const auto n = scaled_argument(f, m);
        const auto k = f.kind == Kind::hurwitz ? 1U : f.weight;
        if (!n || (*n % 4 != 0 && *n % 4 != (k % 2 ? 3 : 1))) continue;
        requested[squarefree_part(cached_factors(*n))][k] = true;
    }
    auto& cache = primitive_values();
    if (cache.size() >= 2048) cache.clear();
    for (auto& [s, weights] : requested) {
        if (std::count(weights.begin(), weights.end(), true) < 2) continue;
        for (unsigned k = 1; k <= 4; ++k) {
            if (!weights[k]) continue;
            const Integer d = s % 4 == (k % 2 ? 3 : 1) ? s : 4 * s;
            const auto key = std::pair{k, d};
            Integer value;
            if (cache.contains(key)) weights[k] = false;
            else if (prefix_cache().get(k, d, value)) {
                cache.emplace(key, std::move(value));
                weights[k] = false;
            }
        }
        if (std::count(weights.begin(), weights.end(), true) < 2) continue;
        for (auto& [key, value] : fundamental_values(s, weights)) {
            prefix_cache().put(key.first, key.second, value);
            cache.emplace(key, std::move(value));
        }
    }
}

// Rational row assembly. Each named function is evaluated at most once per
// index and worker thread, then combined with the literal coefficients.
inline Rational evaluate_rational_row(const ParsedRow& row, const Integer& m) {
    thread_local std::vector<Rational> values(expansions().slot_count);
    thread_local std::vector<Integer> stamps(expansions().slot_count, 0);
    Rational total;
    for (const ParsedTerm& term : row.terms) {
        const std::size_t slot = term.function.slot;
        if (stamps[slot] != m) {
            values[slot] = evaluate(term.function, m);
            stamps[slot] = m;
        }
        const Rational& value = values[slot];
        if (sgn(value) == 0) continue;
        total += term.coefficient * value;
    }
    return total;
}

// ---------------------------------------------------------------------------
// The five rotary rows, as printed in the Exact all-index Grand Formula with
// the closed Q_2, Q_3, Q_4 ambient masses of the fixed half-lattices
// U(3)+A_2(-1), U(2)+D_4(-1) and <1,-3,-3,-3>.

inline Rational scaled(const Integer& value, unsigned long divisor) {
    return make_rational(value, Integer(divisor));
}

inline Integer s_at(const Integer& m, unsigned long j) {
    return m % j == 0 ? twisted_s(m / j) : Integer(0);
}
inline Integer t_at(const Integer& m, unsigned long j) {
    return m % j == 0 ? twisted_t(m / j) : Integer(0);
}
inline Integer sigma2_at(const Integer& m, unsigned long j) {
    return m % j == 0 ? divisor_power_sum(m / j, 2) : Integer(0);
}
inline Integer tau_at(const Integer& m, unsigned long j) {
    return m % j == 0 ? Integer(divisor_count(m / j)) : Integer(0);
}

// M_2^amb(m) = (T(m/2) - S(m/2)) / 18.
inline Rational q2_ambient(const Integer& m) {
    return scaled(t_at(m, 2) - s_at(m, 2), 18);
}

// M_2^open(m) = M_2^amb(m) - tau(m/6)/9 - d_{6,4}(m)/18.
inline Rational q2_open(const Integer& m) {
    return q2_ambient(m) - scaled(tau_at(m, 6), 9)
        - scaled(Integer(restricted_divisor_count(6, 4, m)), 18);
}

// M_3^amb(m) = (3 sigma_2(m/2) - 10 sigma_2(m/4)) / 144.
inline Rational q3_ambient(const Integer& m) {
    return scaled(3 * sigma2_at(m, 2) - 10 * sigma2_at(m, 4), 144);
}

// M_3^open(m) = M_3^amb(m) - R_3^(4)(m) + I_3^(2)(m) with
// R_3^(4)(m) = (3T(m/2) - S(m/2))/24 - S(m/4)/18 and
// I_3^(2)(m) = (11/144) d_{4,0}(m) + (1/16) d_{4,2}(m) + (1/36) d_{12,0}(m).
inline Rational q3_open(const Integer& m) {
    const Rational rank_four = scaled(3 * t_at(m, 2) - s_at(m, 2), 24) - scaled(s_at(m, 4), 18);
    const Rational rank_two = scaled(Integer(11 * restricted_divisor_count(4, 0, m)), 144) +
                              scaled(Integer(restricted_divisor_count(4, 2, m)), 16) +
                              scaled(Integer(restricted_divisor_count(12, 0, m)), 36);
    return q3_ambient(m) - rank_four + rank_two;
}

// M_4^amb(m) = (S(m) + T(m) - 2S(m/2) + 2T(m/2) - 4S(m/4) - 4T(m/4)) / 36.
inline Rational q4_ambient(const Integer& m) {
    return scaled(
        s_at(m, 1) + t_at(m, 1) - 2 * s_at(m, 2) + 2 * t_at(m, 2)
        - 4 * s_at(m, 4) - 4 * t_at(m, 4), 36
    );
}

// M_4^open(m) = M_4^amb(m) + (N_3^+(m) + N_6^+(m) + [3 | m] N_3^+(m/3))/6
//               - d_{6,1}(m)/18 - d_{6,3}(m)/9 - T_9(m) - T_27(m)
//               - I_1(m)/9 - I_3(m)/18 + I_9(m)/18.
inline Rational q4_open(const Integer& m) {
    Rational pell = pell_count(3, m) + pell_count(6, m);
    if (m % 3 == 0) pell += pell_count(3, m / 3);
    Rational total = q4_ambient(m) + pell / Rational(6);
    total -= scaled(Integer(restricted_divisor_count(6, 1, m)), 18);
    total -= scaled(Integer(restricted_divisor_count(6, 3, m)), 9);
    total -= t9_mass(m);
    total -= t27_mass(m);
    total -= scaled(Integer(square_indicator(1, m)), 9);
    total -= scaled(Integer(square_indicator(3, m)), 18);
    total += scaled(Integer(square_indicator(9, m)), 18);
    return total;
}

// M_10^open(m) = (2/5) (tau(m/10) + d_same(m/5) - I_5(m) - I_15(m)).
inline Rational order10_open(const Integer& m) {
    Integer total = tau_at(m, 10);
    if (m % 5 == 0) total += Integer(same_parity_pairs(m / 5));
    total -= Integer(square_indicator(5, m) + square_indicator(15, m));
    return scaled(2 * total, 5);
}

// M_12^open(m) = tau(m/12) / 3.
inline Rational order12_open(const Integer& m) {
    return scaled(tau_at(m, 12), 3);
}

// ---------------------------------------------------------------------------
// The oriented count BB_m, from Theorem 3.26 of Engel, Smillie and Goedgebeur.
// E_{k,chi,psi}(q^t) contributes sum_{d | m/t} chi((m/t)/d) psi(d) d^(k-1) at
// q^m; chi_3 is the nontrivial character modulo three.  The printed E_2(q^9)
// and E_2(q^12) coefficients 1/64 and 2/27 are transcription defects of the
// displayed formula; the authors' archived Sage basis forces 1/18 and 1/48,
// which the values BB_9 = 0 and BB_12 = 1 confirm.

enum class Character { trivial, chi_three };

struct OrientedTerm {
    long numerator;
    unsigned two;    // the denominator is 2^two 3^three 5^five * extra
    unsigned three;
    unsigned five;
    unsigned long extra;
    unsigned weight;
    unsigned long dilation;
    Character chi;   // on (m/t)/d
    Character psi;   // on d
};

constexpr Character triv = Character::trivial;
constexpr Character chi3 = Character::chi_three;

constexpr OrientedTerm oriented_terms[] = {
    {809, 15, 13, 2, 1, 10, 1, triv, triv},
    {-1, 16, 11, 2, 1, 9, 1, triv, chi3},
    {-1, 16, 3, 2, 1, 9, 1, chi3, triv},
    {515, 13, 11, 0, 1, 8, 1, triv, triv},
    {-1, 6, 9, 1, 1, 8, 2, triv, triv},
    {-1, 11, 4, 1, 1, 8, 3, triv, triv},
    {47, 14, 12, 0, 1, 7, 1, triv, chi3},
    {1, 8, 9, 1, 7, 7, 2, triv, chi3},
    {-7, 14, 3, 1, 1, 7, 1, chi3, triv},
    {-1, 8, 3, 1, 7, 7, 2, chi3, triv},
    {2138939, 14, 11, 2, 1, 6, 1, triv, triv},
    {3571, 9, 8, 1, 1, 6, 2, triv, triv},
    {-29, 11, 6, 1, 1, 6, 3, triv, triv},
    {1, 6, 3, 1, 1, 6, 6, triv, triv},
    {274421, 15, 11, 2, 1, 5, 1, triv, chi3},
    {-8231, 11, 9, 1, 1, 5, 2, triv, chi3},
    {-2393, 15, 2, 2, 1, 5, 1, chi3, triv},
    {-203, 11, 2, 1, 1, 5, 2, chi3, triv},
    {12074417, 13, 12, 0, 1, 4, 1, triv, triv},
    {287281, 10, 9, 1, 1, 4, 2, triv, triv},
    {-2672137, 11, 9, 1, 1, 4, 3, triv, triv},
    {-1, 2, 3, 0, 1, 4, 4, triv, triv},
    {-4289, 8, 6, 0, 1, 4, 6, triv, triv},
    {3343037, 14, 11, 0, 1, 3, 1, triv, chi3},
    {227779, 12, 8, 1, 1, 3, 2, triv, chi3},
    {1, 3, 3, 0, 1, 3, 3, triv, chi3},
    {1, 5, 2, 0, 1, 3, 4, triv, chi3},
    {-33103, 14, 3, 1, 1, 3, 1, chi3, triv},
    {-713, 12, 1, 1, 1, 3, 2, chi3, triv},
    {1, 3, 1, 0, 1, 3, 3, chi3, triv},
    {-1, 5, 1, 0, 1, 3, 4, chi3, triv},
    {4314057521, 15, 11, 2, 1, 2, 1, triv, triv},
    {-59136661, 14, 9, 1, 1, 2, 2, triv, triv},
    {-32541151, 11, 10, 1, 1, 2, 3, triv, triv},
    {1, 5, 1, 0, 1, 2, 4, triv, triv},
    {2, 0, 1, 1, 1, 2, 5, triv, triv},
    {5023687, 10, 8, 1, 1, 2, 6, triv, triv},
    {1, 1, 2, 0, 1, 2, 9, triv, triv},    // resolved: 1/18, printed as 1/64
    {1, 4, 1, 0, 1, 2, 12, triv, triv},   // resolved: 1/48, printed as 2/27
    {-24019585289, 15, 13, 2, 1, 1, 1, triv, chi3},
    {-3212143, 11, 9, 1, 1, 1, 2, triv, chi3},
    {-1, 2, 4, 0, 1, 1, 3, triv, chi3},
    {-7, 7, 3, 0, 1, 1, 4, triv, chi3},
    {-2, 0, 1, 1, 1, 1, 5, triv, chi3},
    {-1, 1, 2, 0, 1, 1, 6, triv, chi3},
};

inline long character_discriminant(Character character) { return character == triv ? 1 : -3; }

inline Integer oriented_count(const Integer& m) {
    struct Weights {
        Integer denominator = 1;
        std::array<Integer, std::size(oriented_terms)> numerators;
    };
    static const Weights weights = [] {
        Weights result;
        std::array<Integer, std::size(oriented_terms)> denominators;
        for (std::size_t i = 0; i < denominators.size(); ++i) {
            const OrientedTerm& term = oriented_terms[i];
            Integer& denominator = denominators[i];
            Integer power;
            mpz_ui_pow_ui(denominator.get_mpz_t(), 2, term.two);
            mpz_ui_pow_ui(power.get_mpz_t(), 3, term.three);
            denominator *= power;
            mpz_ui_pow_ui(power.get_mpz_t(), 5, term.five);
            denominator *= power;
            denominator *= term.extra;
            mpz_lcm(result.denominator.get_mpz_t(), result.denominator.get_mpz_t(), denominator.get_mpz_t());
        }
        for (std::size_t i = 0; i < denominators.size(); ++i) {
            Integer& numerator = result.numerators[i];
            mpz_divexact(numerator.get_mpz_t(), result.denominator.get_mpz_t(), denominators[i].get_mpz_t());
            numerator *= oriented_terms[i].numerator;
        }
        return result;
    }();
    Integer total;
    for (std::size_t i = 0; i < weights.numerators.size(); ++i) {
        const OrientedTerm& term = oriented_terms[i];
        if (m % term.dilation != 0) continue;
        const Integer value = two_character_sum(term.weight - 1,
            character_discriminant(term.chi), character_discriminant(term.psi), m / term.dilation);
        mpz_addmul(total.get_mpz_t(), weights.numerators[i].get_mpz_t(), value.get_mpz_t());
    }
    if (sgn(total) < 0 || !mpz_divisible_p(total.get_mpz_t(), weights.denominator.get_mpz_t())) {
        fail("the oriented formula did not produce a nonnegative integer at m=" + m.get_str());
    }
    mpz_divexact(total.get_mpz_t(), total.get_mpz_t(), weights.denominator.get_mpz_t());
    return total;
}

// ---------------------------------------------------------------------------
// One coefficient of the theorem.

struct IprCounts {
    Integer ordinary;
    Integer oriented;
    Integer achiral;
};


struct Row {
    Integer index = 0;
    Integer oriented;
    Integer achiral;
    Integer ordinary;
    Rational split;
    Rational antipodal;
    Rational q2;
    Rational q3;
    Rational q4;
    Rational order10;
    Rational order12;
    std::optional<IprCounts> ipr;
};

} // namespace
