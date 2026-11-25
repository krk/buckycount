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
#include "census_tied_cap_prefix.hpp"
#include "census_small_height.hpp"
#include <gmpxx.h>

namespace gf::fullerene_rank {

// Every vector avoiding all forty linear cap mirrors has Q >= 10.
// The finite cap-ball certificate, including its first witness, is replayed
// by check_short_shells in verify_cap_ball.hpp. For an open cap z and (gamma,beta)=1,
// beta*y-gamma*z is open modulo beta, so the bound applies to every y.
template<bool include_unit = true>
inline bool census_orbit_no_rival(I m, Point beta) {
    const auto d = identifier_norm(identifier_point(beta));
    const Wide bound = std::max(hex_gamma_bound(beta, m), d / m);
    for (I a = 1; Wide(a) * a <= bound; ++a)
        for (I b = 0; identifier_norm({a, b}) <= bound; ++b) {
            const Point gamma{a, b};
            if constexpr (!include_unit)
                if (gamma == Point{1, 0}) continue;
            const auto root = d - Wide(m) * identifier_norm({a, b});
            const auto radius = std::max(root, hex_cap_radius(beta, m, gamma));
            if (radius < 0) continue;
            if (radius >= 20 ||
                identifier_norm(identifier_gcd(identifier_point(gamma), identifier_point(beta))) != 1)
                return false;
        }
    return true;
}

namespace orbit_mass_detail {
// A row has at most size_t::max partners, each with at most I::max phases.
inline Wide row_offset(const mpz_class &rank, const mpz_class &begin) {
    static_assert(std::numeric_limits<std::size_t>::digits + std::numeric_limits<I>::digits <= 127);
    // The nonnegative difference fits 127 bits even when both operands are
    // larger. Subtraction modulo 2^128 yields that exact difference.
    auto low = [](const mpz_class &value) {
        __uint128_t bits = 0;
        for (unsigned shift = 0; shift < 128; shift += GMP_NUMB_BITS)
            bits |= __uint128_t(mpz_getlimbn(value.get_mpz_t(), shift / GMP_NUMB_BITS)) << shift;
        return bits;
    };
    const auto offset = low(rank) - low(begin);
    need(offset <= __uint128_t(std::numeric_limits<Wide>::max()), "orbit row offset overflow");
    return static_cast<Wide>(offset);
}

inline mpz_class unordered_pairs(const mpz_class& count) { return count * (count + 1) / 2; }

// Primitive beta has neither dyadic nor 5-primary content since 2 and 5
// are inert. Its theta-valuation is at most one. At the other primes the
// forty-mirror complement has characteristic polynomial chi_F; Steinberg
// and averaging imply the finite unitary action on it is free.
inline mpz_class primitive_cap_orbits(I d) {
    const mpz_class q(static_cast<long>(d % 3 == 0 ? d / 3 : d));
    if (d % 3 == 0)
        return (q - 1) * (q - 3) * (q - 7) * (3 * q - 7) / 5760;
    return (q - 1) * (q - 7) * (q - 13) * (q - 19) / 155520;
}

struct Type {
    I norm, step, conjugate_phase;
    auto operator<=>(const Type&) const = default;
};
using Histogram = std::map<Type, unsigned long>;

inline I phase_mass(const IdentifierRing& ring, I m, const Type& a, const Type& b, bool fixed) {
    const Wide rhs = Wide(m) + a.norm + b.norm;
    if (rhs % ring.h22) return 0;
    const I step = std::gcd(a.step, b.step);
    if (!fixed) return step;
    // t parametrizes the integral solutions modulo step with phase(t)=t.
    // Conjugation fixes the size and reverses the homogeneous solution.
    // The affine offset is its phase on the phase-zero Bezout lift together
    // with the two single-cap transports.
    const I divisor = std::gcd(I{2}, step);
    if (divisor == 1) return 1;
    // Here 2 divides step and hence h22. Reduction modulo h22,
    // the integral mirror and the phase functional all preserve parity.
    const Wide quotient = (rhs / ring.h22) % 2;
    auto u = identifier_narrow({(ring.bezout_s % 2) * quotient, (ring.bezout_t % 2) * quotient});
    auto mirror = rotate(conjugate(u), ring.beta.y ? 0 : -1);
    const I offset = ring.phase(identifier_point(mirror)) + a.conjugate_phase + b.conjugate_phase;
    return offset % divisor == 0 ? divisor : 0;
}

inline mpz_class pairs(const IdentifierRing& ring, I m, const Histogram& histogram, bool fixed) {
    mpz_class result = 0;
    for (auto i = histogram.begin(); i != histogram.end(); ++i)
        for (auto j = i; j != histogram.end(); ++j) {
            const I weight = phase_mass(ring, m, i->first, j->first, fixed);
            if (weight == 0) continue;
            const mpz_class a(i->second);
            result += static_cast<long>(weight) *
                      (i == j ? unordered_pairs(a) : mpz_class(a * j->second));
        }
    return result;
}

// Keep the conjugation action on caps and its central scalar translation.
// These are unary data. An unordered pair has step gcd(step_i,step_j).
struct CapTypes {
    IdentifierRing ring;
    bool reversible;
    std::vector<Type> type;
    std::vector<std::size_t> mate;
    std::vector<I> delta;

    CapTypes(const IdentifierRing& r, const std::vector<Cap>& caps, const std::vector<I>& steps)
        : ring(r), reversible(r.beta.y == 0 || r.beta.x == r.beta.y) {
        std::map<GeometricCapKey, std::size_t> ids;
        IdentifierParabolicCache transports;
        std::vector<I> phases;
        for (std::size_t i = 0; i < caps.size(); ++i) {
            type.push_back({identifier_narrow(identifier_cap_index(caps[i]) % (2 * ring.h22)), steps[i], 0});
            if (!reversible) continue;
            auto prepared = prepare_boundary_cap(caps[i], ring, transports);
            need(ids.emplace(GeometricCapKey{prepared.geometric, prepared.value.residue}, i).second,
                 "duplicate horizontal cap orbit");
            phases.push_back(prepared.value.phase);
        }
        if (!reversible) return;
        for (std::size_t i = 0; i < caps.size(); ++i) {
            Cap reflected;
            for (unsigned j = 0; j < 4; ++j)
                reflected[j] = rotate(conjugate(caps[i][3 - j]), ring.beta.y ? 0 : -1);
            const auto image = prepare_boundary_cap(reflected, ring, transports);
            auto k = ids.at({image.geometric, image.value.residue});
            mate.push_back(k);
            delta.push_back(image.value.phase - phases[k]);
            type[i].conjugate_phase = positive_mod(delta.back(), 2);
        }
        for (std::size_t i = 0; i < caps.size(); ++i)
            need(mate[mate[i]] == i, "cap conjugation is not an involution");
    }

    mpz_class total(I m) const {
        Histogram all, fixed;
        mpz_class swapped = 0;
        for (std::size_t i = 0; i < type.size(); ++i) {
            ++all[{type[i].norm, type[i].step, 0}];
            if (!reversible) continue;
            if (i == mate[i]) ++fixed[type[i]];
            if (i < mate[i])
                swapped += static_cast<long>(phase_mass(ring, m, type[i], type[mate[i]], true));
        }
        auto unitary = pairs(ring, m, all, false);
        if (!reversible) return unitary;
        mpz_class doubled = unitary + swapped + pairs(ring, m, fixed, true);
        need(mpz_even_p(doubled.get_mpz_t()), "nonintegral ordinary family mass");
        return doubled / 2;
    }
};

inline mpz_class parabolic_mass(I m, const IdentifierRing& ring) {
    std::vector<Cap> raw;
    std::vector<I> steps;
    BoundaryOwnerStats stats;
    visit_boundary_caps(ring, [&](const Cap& cap, I step) {
        raw.push_back(cap);
        steps.push_back(step);
    }, stats);
    return CapTypes(ring, raw, steps).total(m);
}
} // namespace orbit_mass_detail

// Return an exact count for the whole family, or no value if a nonhorizontal
// competitor still needs checking. This entry point takes no reference
// prefixes or graph records. Primitive nonsymmetric heights need no inventory.
inline std::optional<mpz_class> census_orbit_family_mass(I m, Point beta) {
    if (!small_height_detail::admissible(m, beta)) return mpz_class(0);
    if (identifier_norm(identifier_point(beta)) > std::numeric_limits<I>::max()) return std::nullopt;
    const IdentifierRing ring(beta);
    if (ring.h22 == 1) {
        const auto count = orbit_mass_detail::primitive_cap_orbits(ring.D);
        if (count == 0) return mpz_class(0);
        if (beta.x != beta.y && beta.y && census_orbit_no_rival(m, beta))
            return orbit_mass_detail::unordered_pairs(count);
    }
    if (!census_orbit_no_rival(m, beta)) return std::nullopt;
    return orbit_mass_detail::parabolic_mass(m, ring);
}

// A selected family retains one entry per cap and one suffix count per cap.
// No pair is built to locate the first cap. A row scan then selects the
// second cap and phase in the native pooled lexicographic order.
class CensusOrbitFamilyIndex {
    using Type = orbit_mass_detail::Type;
    using Histogram = orbit_mass_detail::Histogram;
    I m_;
    Point beta_;
    std::vector<CensusHeightCap> caps_;
    orbit_mass_detail::CapTypes info_;
    std::vector<mpz_class> row_begin_;

    std::map<I, std::map<Type, CensusTiedCapPrefix>> partners_;

    using PartnerBank = std::map<Type, CensusTiedCapPrefix>;
    Wide prefix_mass(std::size_t i, std::size_t last, const PartnerBank &bank) const {
        CENSUS_HEAP_SCOPE(family_prefix);
        // At most K partners, each of weight at most the 63-bit phase step.
        Wide mass = 0;
        auto x = info_.mate[i];
        for (const auto &[b, index] : bank) {
            auto step = orbit_mass_detail::phase_mass(info_.ring, m_, info_.type[i], b, false);
            if (!step) continue;
            auto kept_count = index.count(i, x, last);
            mass += Wide(kept_count) * step;
            // Equality means either two fixed caps or the single transposed pair.
            auto fixed_count = x == i ? index.fixed_count(i, last)
                                      : std::size_t(i < x && x <= last && info_.type[x] == b);
            if (fixed_count) {
                auto fixed_phases = orbit_mass_detail::phase_mass(info_.ring, m_, info_.type[i], b, true);
                mass -= Wide(fixed_count) * ((step - fixed_phases) / 2);
            }
        }
        return mass;
    }

    static std::vector<CensusHeightCap> checked_inventory(I m, Point beta) {
        need(small_height_detail::admissible(m, beta) &&
             identifier_norm(identifier_point(beta)) <= std::numeric_limits<I>::max() &&
             census_orbit_no_rival(m, beta), "orbit selector requires a certified representable height");
        return census_small_height_inventory(beta);
    }

    static orbit_mass_detail::CapTypes cap_types(Point beta, const std::vector<CensusHeightCap>& caps) {
        std::vector<Cap> raw;
        std::vector<I> steps;
        for (const auto& c : caps) {
            raw.push_back(c.pooled);
            steps.push_back(c.phase_step);
        }
        return {IdentifierRing(beta), raw, steps};
    }

    void insert_type(const Type& a, Histogram& bins, bool fixed, mpz_class& sum) {
        using orbit_mass_detail::phase_mass;
        sum += static_cast<long>(phase_mass(info_.ring, m_, a, a, fixed));
        for (const auto& [b, count] : bins)
            sum += mpz_class(count) * static_cast<long>(phase_mass(info_.ring, m_, a, b, fixed));
        ++bins[a];
    }

    bool kept(std::size_t i, std::size_t j) const {
        if (!info_.reversible) return true;
        auto [a, b] = std::minmax(info_.mate[i], info_.mate[j]);
        return std::tie(i, j) <= std::tie(a, b);
    }

    bool fixed(std::size_t i, std::size_t j) const {
        if (!info_.reversible) return false;
        auto [a, b] = std::minmax(info_.mate[i], info_.mate[j]);
        return i == a && j == b;
    }

  public:
    explicit CensusOrbitFamilyIndex(I m, Point beta)
        : m_(m), beta_(beta), caps_(checked_inventory(m, beta)), info_(cap_types(beta, caps_)),
          row_begin_(caps_.size() + 1) {
        if (info_.reversible) {
            std::map<Type, std::vector<std::pair<long, long>>> banks;
            for (std::size_t j = 0; j < caps_.size(); ++j)
                banks[info_.type[j]].emplace_back(j, info_.mate[j]);
            for (const auto &[type, bank] : banks) partners_[type.norm % info_.ring.h22].try_emplace(type, bank);
        }
        Histogram all, fixed_bins;
        mpz_class unitary = 0, invariant = 0;
        for (std::size_t k = caps_.size(); k-- > 0;) {
            // Retain exactly the stable set {c: c>=k and mate(c)>=k}.
            // An ordinary class lies in its unordered square precisely when
            // its native minimal pair has both entries at least k.
            const auto j = info_.reversible ? info_.mate[k] : k;
            if (k <= j) {
                auto insert = [&](std::size_t i) {
                    const auto& a = info_.type[i];
                    insert_type({a.norm, a.step, 0}, all, false, unitary);
                };
                insert(k);
                if (k < j) {
                    insert(j);
                    invariant += static_cast<long>(orbit_mass_detail::phase_mass(
                        info_.ring, m, info_.type[k], info_.type[j], true));
                } else if (info_.reversible) {
                    insert_type(info_.type[k], fixed_bins, true, invariant);
                }
            }
            row_begin_[k] = unitary;
            if (info_.reversible) {
                row_begin_[k] += invariant;
                need(mpz_even_p(row_begin_[k].get_mpz_t()), "nonintegral ordinary suffix mass");
                row_begin_[k] /= 2;
            }
        }
        const mpz_class total_mass = row_begin_.front();
        for (auto &suffix : row_begin_) suffix = total_mass - suffix;
    }

    const mpz_class& total() const { return row_begin_.back(); }

    IdentifierCusp select(const mpz_class& rank) const {
        need(rank >= 0 && rank < total(), "orbit family rank out of range");
        auto next = std::upper_bound(row_begin_.begin(), row_begin_.end(), rank);
        const auto i = static_cast<std::size_t>(next - row_begin_.begin() - 1);
        Wide local = orbit_mass_detail::row_offset(rank, row_begin_[i]);
        const auto& ring = info_.ring;
        // Primitive skew heights give one phase for every unordered pair.
        const bool triangular = ring.h22 == 1 && !info_.reversible;
        auto first_cap = triangular ? i + static_cast<std::size_t>(local) : i;
        if (triangular) local = 0;
        auto limit = caps_.size();
        if (info_.reversible) {
            const I residue = positive_mod(identifier_narrow(-(Wide(m_) + info_.type[i].norm) % ring.h22), ring.h22);
            const auto &bank = partners_.at(residue);
            // A prefix query already costs O(types * log K). Resolve a final
            // interval of at most that width directly, with the same bound.
            auto width = bank.size() * std::bit_width(caps_.size());
            auto low = i, high = caps_.size() - 1;
            while (high - low > width) {
                auto mid = low + (high - low) / 2;
                if (prefix_mass(i, mid, bank) <= local) low = mid + 1;
                else high = mid;
            }
            first_cap = low;
            limit = high + 1;
            if (low > i) local -= prefix_mass(i, low - 1, bank);
        }
        for (std::size_t j = first_cap; j < limit; ++j) {
            if (!kept(i, j)) continue;
            const auto& a = info_.type[i];
            const auto& b = info_.type[j];
            const I step = orbit_mass_detail::phase_mass(ring, m_, a, b, false);
            if (!step) continue;
            const bool self = fixed(i, j);
            const I choices = self ? (step + orbit_mass_detail::phase_mass(ring, m_, a, b, true)) / 2 : step;
            if (local >= static_cast<long>(choices)) {
                local -= static_cast<long>(choices);
                continue;
            }
            const Wide rhs = m_ + caps_[i].norm + caps_[j].norm;
            const Wide h = hex_height(beta_), g = ring.h22, dq = h / g;
            const Wide q0 = ring.bezout_t * (rhs / g);
            const Wide shift = beta_.y ? 0 : rhs / h;
            const Wide first = beta_.y ? identifier_floor(identifier_floor(h - 1 - shift, 2) - q0, dq) - step + 1 : identifier_ceil(identifier_ceil(-h - shift, 2) - q0, dq);
            for (I k = 0; k < step; ++k) {
                const Wide q = q0 + (first + k) * dq;
                const Point u = identifier_narrow({(rhs - beta_.y * q) / h, q});
                if (self && step > 1) {
                    auto anti = rotate(conjugate(u), beta_.y ? 0 : -1);
                    auto delta = ring.phase(identifier_point(anti)) + info_.delta[i] + info_.delta[j] -
                                 ring.phase(identifier_point(u));
                    const auto mate_k = positive_mod(k + delta, step);
                    if ((!beta_.y && mate_k < k) || (beta_.y && mate_k > k)) continue;
                }
                if (local == 0) return {u, beta_, {caps_[i].pooled, caps_[j].pooled}};
                --local;
            }
            throw std::logic_error("ordinary phase count and selection disagree");
        }
        throw std::logic_error("ordinary row count and selection disagree");
    }
};

} // namespace gf::fullerene_rank
