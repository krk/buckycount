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
#include "census_heap_scopes.hpp"

#include <bit>
#include <span>

#include "cap_boundary_owner.hpp"
#include "fullerene_hex_witness.hpp"

namespace gf::fullerene_rank {

// At bounded height the root-open cap residues have a small finite inventory.
// Choose their lexicographic representatives in the *pooled* [-H,H] box;
// the nearest representatives of the original inventory use a different order.
struct CensusHeightCap {
    Cap pooled;
    I phase_step;
    Wide norm;
    Wide nearest_norm;
    I nearest_phase;
    Cap nearest;
    IdentifierPoint to_nearest;
    GeometricPreparedCap geometric;
};

namespace small_height_detail {
class PooledCapReducer {
    struct Row {
        I x_shift, offset;
    };
    I radius_, h11_;
    std::vector<Row> rows_;
    std::vector<I> columns_;
    std::vector<Point> pooled_;

    // Each operand is already in [-H,H]^2. Bound both coordinates of
    // every tested linear form, including the one-coordinate cap reflection.
    static I coefficient_bound() {
        static const I bound = [] {
            I maximum = 5;
            for (const auto& level : identifier_cap_chain())
                for (const auto& rep : level) {
                    I real = 0, imag = 0;
                    for (auto p : rep.row) {
                        real += std::abs(p.x) + std::abs(p.y);
                        imag += std::abs(p.y) + std::abs(p.x + p.y);
                    }
                    maximum = std::max({maximum, real, imag});
                }
            return maximum;
        }();
        return bound;
    }

    std::size_t index(Point p) const {
        const auto& row = rows_[radius_ + p.y];
        return row.offset + columns_[radius_ + p.x + row.x_shift];
    }

    using Path = std::array<const IdentifierCapTransversal*, 4>;
    Cap minimum(const Cap& cap, std::size_t level, Path& path) const {
        if (level == 4) return cap;
        std::array<const IdentifierCapTransversal*, 240> tied;
        std::size_t count = 0;
        Point least{std::numeric_limits<I>::max(), 0};
        for (const auto& rep : identifier_cap_chain()[level]) {
            Point sum{};
            for (std::size_t j = 0; j < cap.size(); ++j) {
                const auto a = rep.row[j], b = cap[j];
                sum.x += a.x * b.x - a.y * b.y;
                sum.y += a.y * b.x + (a.x + a.y) * b.y;
            }
            const auto q = pooled_[index(sum)];
            if (q < least) {
                least = q;
                count = 0;
            }
            if (q == least) tied[count++] = &rep;
        }
        std::optional<Cap> best;
        std::vector<Cap> seen;
        for (std::size_t k = 0; k < count; ++k) {
            const auto& word = tied[k]->word;
            auto child = cap;
            for (auto i = word.rbegin(); i != word.rend(); ++i) {
                child = cap_reflect(child, *i);
                child[*i] = pooled_[index(child[*i])];
            }
            if (count > 1) {
                if (std::find(seen.begin(), seen.end(), child) != seen.end()) continue;
                seen.push_back(child);
            }
            Path candidate{};
            auto result = minimum(child, level + 1, candidate);
            if (!best || result < *best) {
                best = result;
                candidate[level] = tied[k];
                path = candidate;
            }
        }
        return best.value();
    }

  public:
    explicit PooledCapReducer(const IdentifierRing& ring)
        : radius_(coefficient_bound() * (ring.beta.x + ring.beta.y)), h11_(ring.h11), rows_(2 * radius_ + 1),
          columns_(2 * radius_ + ring.h11), pooled_(ring.D, {std::numeric_limits<I>::max(), 0}) {
        // In HNF, reducing (x,y) first takes (0,y) to (s,t), then takes
        // (x+s) modulo h11. Both lookups have bounded arguments.
        for (I y = -radius_; y <= radius_; ++y) {
            auto p = ring.reduce(IdentifierPoint{0, y});
            rows_[radius_ + y] = {p.x, p.y * ring.h11};
        }
        for (I x = -radius_; x < radius_ + ring.h11; ++x)
            columns_[radius_ + x] = positive_mod(x, ring.h11);
        const I h = ring.beta.x + ring.beta.y;
        for (I x = -h; x <= h; ++x)
            for (I y = -h; y <= h; ++y) {
                auto& p = pooled_[index({x, y})];
                p = std::min(p, Point{x, y});
            }
        // A centered half-open fundamental parallelogram spanned by beta and
        // zeta*beta fits in the pooled box of the ordinary height wedge.
        need(std::none_of(pooled_.begin(), pooled_.end(), [](Point p) {
            return p.x == std::numeric_limits<I>::max();
        }), "pooled box misses a height residue");
    }

    Cap operator()(Cap residue) const {
        for (auto& q : residue) q = pooled_[q.y * h11_ + q.x];
        Path path{};
        return minimum(residue, 0, path);
    }

    // Lift the winning finite group word to the original nearest owner.
    // Its image is a nearest translate of the pooled cap and carries the
    // exact scalar phase taking that cap back to the geometric owner.
    std::pair<Cap, Cap> from_owner(const IdentifierRing& ring, const Cap& owner) const {
        auto residue = ring.reduce(owner);
        for (auto& q : residue) q = pooled_[q.y * h11_ + q.x];
        Path path{};
        auto pooled = minimum(residue, 0, path);
        auto nearest = owner;
        for (auto rep : path)
            for (auto i = rep->word.rbegin(); i != rep->word.rend(); ++i)
                nearest = cap_reflect(nearest, *i);
        return {pooled, nearest};
    }
};
} // namespace small_height_detail

inline std::vector<CensusHeightCap> census_small_height_inventory(Point beta) {
    CENSUS_HEAP_SCOPE(inventory);
    need(beta.x > 0 && (beta.y == 0 || beta.y >= beta.x), "height outside ordinary wedge");
    IdentifierRing ring(beta);
    std::optional<small_height_detail::PooledCapReducer> minimum_cap;
    std::vector<CensusHeightCap> caps;
    BoundaryOwnerStats stats;
    visit_boundary_caps(ring, [&](const Cap &owner, I step) {
        if (!minimum_cap) minimum_cap.emplace(ring);
        auto [pooled, nearest_image] = minimum_cap->from_owner(ring, owner);
        Cap translation;
        for (unsigned j = 0; j < 4; ++j) {
            IdentifierPoint shift;
            need(identifier_divide(identifier_sub(identifier_point(nearest_image[j]), identifier_point(pooled[j])),
                                   identifier_point(beta), shift),
                 "incorrect nearest-orbit transport");
            translation[j] = identifier_narrow(shift);
        }
        auto image_delta = geometric_translation_delta(translation, pooled, beta);
        need(step > 0 && ring.h22 % step == 0, "invalid cap phase subgroup");
        auto nearest = nearest_image;
        auto correction = image_delta;
        if (boundary_owner_position(owner, beta).voronoi_boundary) {
            auto lex = nearest_identifier_cap(pooled, beta);
            nearest = lex.representative;
            correction = geometric_translation_delta(lex.translation, pooled, beta);
        }
        const auto norm = identifier_cap_index(owner);
        GeometricPreparedCap geometric{{owner, owner, {}, norm, positive_mod(ring.phase(image_delta), step), step}, true};
        caps.push_back(
            {pooled, step, identifier_cap_index(pooled), norm, ring.phase(correction), nearest, correction, geometric});
    }, stats);
    std::sort(caps.begin(), caps.end(), [](const auto& a, const auto& b) { return a.pooled < b.pooled; });
    need(std::adjacent_find(caps.begin(), caps.end(), [](const auto& a, const auto& b) { return a.pooled == b.pooled; }) == caps.end(), "duplicate cap orbit in height inventory");
    return caps;
}

inline std::vector<Cap> census_small_height_caps(Point beta) {
    std::vector<Cap> caps;
    for (auto& c : census_small_height_inventory(beta))
        caps.push_back(c.pooled);
    return caps;
}

// Every ordinary owner at this height has two individually minimum caps
// in nondecreasing order. The size equation and the half-open central strip
// then leave only the listed scalars. Other cusps and phase ties still need
// the pooled affine owner predicate.
namespace small_height_detail {
struct SwapWindow {
    Wide first = 0, size = 0;
};

// Project the open swap hexagon onto its scalar size and phase. For fixed
// R=H*p+b*q, all integral scalars are base+t*kappa; the three strips give
// one interval of t. An entire cap pair is excluded when this interval
// covers every residue of its joint phase subgroup.
inline std::vector<SwapWindow> swap_windows(const IdentifierRing& ring, I m) {
    const Wide h = identifier_add(ring.beta.x, ring.beta.y), g = ring.h22;
    const Wide maximum = identifier_mul(h, h - 1) / g;
    if (identifier_ceil(m, g) > maximum) return {};
    need(maximum >= 0 && maximum < static_cast<Wide>(std::vector<SwapWindow>().max_size()), "projected hexagon is too large");
    std::vector<SwapWindow> windows(static_cast<std::size_t>(maximum + 1));
    const IdentifierPoint kappa{-ring.beta.y / g, h / g};
    for (Wide r = identifier_ceil(m, g); r <= maximum; ++r) {
        const IdentifierPoint base{identifier_mul(ring.bezout_s, r), identifier_mul(ring.bezout_t, r)};
        if (auto hit = hex_scalar_interval(base, kappa, h - 1))
            windows[static_cast<std::size_t>(r)] = {hit->first, hit->second - hit->first + 1};
    }
    return windows;
}

inline bool admissible(I m, Point beta) {
    need(m > 0, "small-height census requires a positive size");
    need(beta.x > 0 && (beta.y == 0 || beta.y >= beta.x), "height outside ordinary wedge");
    const Wide h = identifier_add(beta.x, beta.y);
    return identifier_mul(h, h) < identifier_mul(6, m);
}
// A fully covered swap window rejects every pair of two cap types. Store the
// surviving partner indices as words, so iterating their bits retains the
// pooled inventory order, including diagonal pairs. Partition first-cap rows
// into batches with at most one machine word's worth of distinct types:
// even if every cap has a different norm, temporary storage is linear.
class SwapPartners {
    using Word = std::uint64_t;
    static constexpr std::size_t width = std::numeric_limits<Word>::digits;
    struct Block { std::size_t word; Word bits; };
    struct Type {
        Wide norm;
        I step;
        std::vector<Block> members;
    };
    std::vector<Type> types_;
    std::vector<std::size_t> row_type_;
    std::vector<std::optional<std::vector<Block>>> survivors_;
    std::vector<std::size_t> batch_;
    std::vector<Word> words_;
    I m_, content_;
    const std::vector<SwapWindow>& windows_;

    const std::vector<Block>& partners(std::size_t id) {
        if (survivors_[id]) return *survivors_[id];
        if (batch_.size() == width) {
            for (auto old : batch_) survivors_[old].reset();
            batch_.clear();
        }
        batch_.push_back(id);
        auto& result = survivors_[id].emplace();
        const auto& first = types_[id];
        for (const auto& second : types_) {
            const Wide total = m_ + first.norm + second.norm;
            if (total % content_) continue;
            const Wide index = total / content_;
            if (index < static_cast<Wide>(windows_.size()) &&
                windows_[static_cast<std::size_t>(index)].size >= std::gcd(first.step, second.step))
                continue;
            for (auto b : second.members) words_[b.word] |= b.bits;
        }
        for (std::size_t w = 0; w < words_.size(); ++w) {
            if (words_[w]) result.push_back({w, words_[w]});
            words_[w] = 0;
        }
        return result;
    }
  public:
    SwapPartners(I m, I content, const std::vector<CensusHeightCap>& caps,
                 const std::vector<SwapWindow>& windows)
        : words_((caps.size() + width - 1) / width), m_(m), content_(content), windows_(windows) {
        std::map<std::pair<Wide, I>, std::size_t> ids;
        row_type_.reserve(caps.size());
        for (std::size_t i = 0; i < caps.size(); ++i) {
            const std::pair key{caps[i].nearest_norm, caps[i].phase_step};
            const auto [it, added] = ids.try_emplace(key, types_.size());
            if (added) types_.push_back({key.first, key.second, {}});
            row_type_.push_back(it->second);
            auto& members = types_[it->second].members;
            if (members.empty() || members.back().word != i / width) members.push_back({i / width, 0});
            members.back().bits |= Word{1} << (i % width);
        }
        survivors_.resize(types_.size());
    }

    template<class Visit> void visit_words(std::size_t i, Visit&& emit, std::span<const U> rejected = {}) {
        const auto& blocks = partners(row_type_[i]);
        auto it = std::lower_bound(blocks.begin(), blocks.end(), i / width,
                                   [](Block b, std::size_t word) { return b.word < word; });
        for (; it != blocks.end(); ++it) {
            Word bits = rejected.empty() ? it->bits : it->bits & ~rejected[it->word];
            if (it->word == i / width) bits &= ~Word{0} << (i % width);
            if (bits) emit(it->word, bits);
        }
    }

    template<class Visit> void visit(std::size_t i, Visit&& emit, std::span<const U> rejected = {}) {
        visit_words(i, [&](std::size_t word, U bits) {
            while (bits) {
                emit(word * width + std::countr_zero(bits));
                bits &= bits - 1;
            }
        }, rejected);
    }
};

// After removing the cyclic swap window, the surviving ordinal sequence
// consists of at most two increasing runs. Index it without expanding points.
struct ScalarSelection {
    Wide p, q, dp, dq;
    I step, residue, excluded;
    I size() const { return step - excluded; }
    Point operator[](I index) const {
        I k = index;
        if (residue < excluded) k += excluded - residue;
        else if (k >= step - residue) k += excluded;
        return identifier_narrow({p - k * dp, q + k * dq});
    }
    // Intersect an ordinal interval with the complement of the cyclic swap
    // exclusion. Each emitted range is contiguous in the original progression.
    template<class Visit> void intersect(I lo, I hi, Visit visit) const {
        auto emit = [&](I a, I b) {
            a = std::max(a, lo);
            b = std::min(b, hi);
            if (a < b) visit(a, b);
        };
        if (!excluded) emit(0, step);
        else if (residue < excluded) emit(excluded - residue, step - residue);
        else {
            emit(0, step - residue);
            emit(step - residue + excluded, step);
        }
    }
    // Reflection pairs k with (c-k) mod step. Its affine formula changes
    // at k=c+1. On either side the canonical-order comparison is one bound.
    template<class Visit> void reflected(I c, bool prefer_greater, Visit visit) const {
        if (prefer_greater) {
            intersect((c + 1) / 2, c + 1, visit);
            intersect(std::max(c + 1, (step + c + 1) / 2), step, visit);
        } else {
            intersect(0, c / 2 + 1, visit);
            intersect(c + 1, (step + c) / 2 + 1, visit);
        }
    }
    ScalarSelection run(I lo, I hi) const {
        return {p - lo * dp, q + lo * dq, dp, dq, hi - lo, 0, 0};
    }
    struct Iterator {
        const ScalarSelection &window;
        I index;
        Point operator*() const { return window[index]; }
        void operator++() { ++index; }
        bool operator!=(const Iterator &other) const { return index != other.index; }
    };
    Iterator begin() const { return {*this, 0}; }
    Iterator end() const { return {*this, size()}; }
};

template <bool stabilize, bool selecting = false, class Visit, class Witness = std::nullptr_t>
void groups(I m, Point beta, const std::vector<CensusHeightCap>& caps, Visit visit, Witness witness = nullptr) {
    const Wide h = identifier_add(beta.x, beta.y);
    IdentifierRing ring(beta);
    const auto content = ring.h22;
    const auto q_step = h / content;
    std::vector<SwapWindow> windows;
    if constexpr (stabilize) windows = swap_windows(ring, m);
    // Split Q(z)=H*k+r. For fixed r_i+r_j the scalar strip is unchanged on
    // skew heights when k_i+k_j is added to p. The real strip then moves
    // its lower q endpoint by half this translation.
    struct NormParts { Wide quotient, nearest; I residue; };
    std::vector<NormParts> parts;
    parts.reserve(caps.size());
    for (const auto& cap : caps) {
        I r = identifier_narrow(cap.norm % h);
        // Q(nearest) and Q(pooled) agree modulo content.
        parts.push_back({cap.norm / h, (cap.nearest_norm - r) / content, r});
    }
    struct ScalarBase { Wide p, q, t, size; };
    std::vector<ScalarBase> bases(static_cast<std::size_t>(2 * h - 1));
    const Wide upper = (h - 1) / 2;
    const Wide dp = beta.y / content;
    for (std::size_t r = 0; r < bases.size(); ++r) {
        const Wide n = identifier_add(m, r);
        if (n % content) continue;
        const Wide size = n / content;
        const Wide q0 = identifier_mul(ring.bezout_t, size);
        const Wide t = identifier_floor(upper - q0, q_step);
        const Wide q = q0 + t * q_step;
        bases[r] = {(n - beta.y * q) / h, q, t, size};
    }

    std::optional<SwapPartners> partners;
    if (!windows.empty() || witness != nullptr) partners.emplace(m, content, caps, windows);
    std::vector<Point> scalars;
    // At row i each list contains exactly the still-available j >= i,
    // in inventory order. The norm equation selects one whole list.
    std::vector<std::size_t> heads(content, caps.size()), next(caps.size());
    for (std::size_t i = caps.size(); i-- > 0;) {
        auto& head = heads[parts[i].residue % content];
        next[i] = head;
        head = i;
    }
    for (std::size_t i = 0; i < caps.size(); ++i) {
        const I residue_i = parts[i].residue % content;
        const I required = positive_mod(-(m % content) - residue_i, content);
        auto first_match = heads[required];
        heads[residue_i] = next[i];
        auto emit_pair = [&](std::size_t j) {
            const I step = stabilize ? std::gcd(caps[i].phase_step, caps[j].phase_step) : content;
            auto window = [&]() -> ScalarSelection {
                const auto& base = bases[parts[i].residue + parts[j].residue];
                SwapWindow excluded;
                if constexpr (stabilize) {
                    const Wide id = base.size + parts[i].nearest + parts[j].nearest;
                    if (step > 1 && id < static_cast<Wide>(windows.size())) excluded = windows[static_cast<std::size_t>(id)];
                }
                const Wide quotient = parts[i].quotient + parts[j].quotient;
                Wide p = base.p + quotient;
                Wide q, t;
                if (beta.y == 0) {
                    q = identifier_ceil(-h - p, 2);
                    t = q; // bezout(H,0)=(H,1,0), and the scalar step is (0,1).
                } else {
                    // The skew lexicographic minimum in each stabilizer class
                    // occurs in the last step positions, taken in ascending q.
                    q = base.q - (step - 1) * q_step;
                    p += (step - 1) * dp;
                    t = base.t - ring.bezout_t * quotient - (step - 1);
                }
                I residue = 0;
                if constexpr (stabilize)
                    if (excluded.size)
                        residue = positive_mod(identifier_narrow((t - excluded.first) % step) + caps[i].nearest_phase + caps[j].nearest_phase, step);
                return {p, q, dp, q_step, step, residue, static_cast<I>(std::min<Wide>(excluded.size, step))};
            };
            if constexpr (selecting) {
                // The owner can reject the cap pair before constructing its scalar window.
                visit(i, j, step, window);
            } else {
                scalars.clear();
                for (Point u : window()) scalars.push_back(u);
                if (!scalars.empty()) visit(i, j, scalars);
            }
        };
        std::span<const U> rejected;
        if constexpr (!std::is_same_v<Witness, std::nullptr_t>)
            if (witness) rejected = witness->rejected(i);
        if (partners) partners->visit(i, emit_pair, rejected);
        else for (std::size_t j = first_match; j < caps.size(); j = next[j]) emit_pair(j);
    }
}
} // namespace small_height_detail

template <class Visit>
void census_small_height_candidates(I m, Point beta, Visit visit) {
    if (!small_height_detail::admissible(m, beta)) return;
    auto inventory = census_small_height_inventory(beta);
    small_height_detail::groups<false>(m, beta, inventory, [&](auto i, auto j, const auto& scalars) {
        for (auto u : scalars)
            visit(IdentifierCusp{u, beta, {inventory[i].pooled, inventory[j].pooled}});
    });
}

// Each retained cap pair appears once. Scalars are in increasing q order.
// Every removed phase has an explicit strict competitor: either a smaller
// scalar at identical caps or a bounded-height swap after nearest translation.
template <class Visit>
void census_small_height_groups(I m, Point beta, Visit visit) {
    if (!small_height_detail::admissible(m, beta)) return;
    auto caps = census_small_height_inventory(beta);
    small_height_detail::groups<true>(m, beta, caps, [&](auto i, auto j, const auto& scalars) {
        visit(std::array<Cap, 2>{caps[i].pooled, caps[j].pooled}, scalars);
    });
}
template <class Visit>
void census_small_height_groups(I m, Point beta, const std::vector<CensusHeightCap>& caps, Visit visit) {
    if (!small_height_detail::admissible(m, beta)) return;
    small_height_detail::groups<true>(m, beta, caps, [&](auto i, auto j, const auto& scalars) {
        visit(std::array<Cap, 2>{caps[i].pooled, caps[j].pooled}, scalars);
    });
}
} // namespace gf::fullerene_rank
