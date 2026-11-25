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
#include "census_witness_products.hpp"

namespace gf::fullerene_rank {

struct CensusHexLine {
    IdentifierPoint slope;
    unsigned facet;
    Wide direction;
    explicit CensusHexLine(IdentifierPoint s) : slope(s), facet(0), direction(s[0]) {
        const std::array<Wide, 3> projections{s[0], s[1], identifier_add(s[0], s[1])};
        for (unsigned i = 1; i < 3; ++i)
            if (std::max(projections[i], -projections[i]) > std::max(direction, -direction)) {
                facet = i;
                direction = projections[i];
            }
        need(direction != 0, "zero null-line direction");
    }
    template <class Visit>
    bool visit(IdentifierPoint base, Wide h, Visit accept) const {
        Wide v = facet == 2 ? identifier_add(base[0], base[1]) : base[facet];
        if (direction < 0) v = -v;
        Wide step = std::max(direction, -direction);
        // The strongest facet selects at most three scalar lifts since
        // H(slope)>=|beta|>=sqrt(3)H(beta)/2. Check the two remaining
        // facets directly, avoiding their divisions.
        Wide j = identifier_ceil(identifier_sub(-h, v), step);
        Wide projected = identifier_add(v, identifier_mul(j, step));
        if (projected > h) return true;
        auto c = identifier_add(base, {identifier_mul(j, slope[0]), identifier_mul(j, slope[1])});
        for (; projected <= h; projected += step, ++j, c = identifier_add(c, slope)) {
            auto coefficient = identifier_narrow(c);
            auto height = hex_height(coefficient);
            if (height <= h && !accept(coefficient, j, height)) return false;
        }
        return true;
    }
};

// Exact complex-linear action of a fixed null transport. The descent
// supplies each integral generator. Its Heisenberg scalar correction is
// linear in the transported cap and beta, so compose the generators once.
class CensusTransport {
    struct Term {
        unsigned column;
        Point coefficient;
    };
    using Column = std::array<Point, 10>;
    std::array<std::vector<Term>, 10> rows_;
    static Column pack(const IdentifierCusp& v) {
        Column a{v.u, v.beta};
        for (unsigned i = 0; i < 8; ++i)
            a[i + 2] = v.caps[i / 4][i % 4];
        return a;
    }
    static IdentifierCusp unpack(const Column& a) {
        IdentifierCusp v{a[0], a[1], {}};
        for (unsigned i = 0; i < 8; ++i)
            v.caps[i / 4][i % 4] = a[i + 2];
        return v;
    }

  public:
    explicit CensusTransport(const CuspTransport& transform) {
        for (unsigned j = 0; j < 10; ++j) {
            Column column{};
            column[j] = {1, 0};
            auto v = unpack(column);
            for (const auto& step : transform.steps)
                v = apply_cusp_transport(v, step);
            column = pack(identifier_scale_unit(v, transform.unit));
            for (unsigned i = 0; i < 10; ++i)
                if (column[i] != Point{}) rows_[i].push_back({j, column[i]});
        }
    }
    IdentifierCusp operator()(const IdentifierCusp& v) const {
        auto in = pack(v);
        Column out{};
        for (unsigned i = 0; i < 10; ++i) {
            IdentifierPoint value{};
            for (auto [j, c] : rows_[i]) {
                auto p = in[j];
                // Signed 64-bit pair products fit in Wide.
                value = identifier_add(value, {identifier_sub(Wide(c.x) * p.x, Wide(c.y) * p.y),
                                               identifier_add(Wide(c.y) * p.x, identifier_mul(identifier_add(c.x, c.y), p.y))});
            }
            out[i] = identifier_narrow(value);
        }
        return unpack(out);
    }
};

// For gamma=1, translate the witness caps by -y and its scalar by
// -(A.y+Q(y))*zeta. Nullity gives A.x=Q(y), so the result is (0,1,0).
// Swap and multiply by zeta to send it to e. Its transported beta is
// exactly h(w,v)/(2-zeta), already computed by the null-line test.
inline IdentifierCusp census_unit_transport(const IdentifierCusp& v, const std::array<Cap, 2>& y, Point coefficient) {
    IdentifierCusp image{v.beta, coefficient, {}};
    for (unsigned side = 0; side < 2; ++side)
        for (unsigned j = 0; j < 4; ++j) {
            auto shifted = identifier_sub(identifier_point(v.caps[side][j]),
                                          identifier_mul(identifier_point(v.beta), identifier_point(y[side][j])));
            image.caps[side][j] = identifier_narrow({-shifted[1], identifier_add(shifted[0], shifted[1])});
        }
    return image;
}

// For a fixed cap pair both the norm-solving Bezout base and its nearest
// translation are fixed. Only the phase and the final central shift vary.
struct CensusNearestFiber {
    const IdentifierRing& ring;
    I subgroup, phase_offset;
    IdentifierPoint base, kappa;
    Wide position, phase_increment;
    std::array<Cap, 2> nearest;

    CensusNearestFiber(const IdentifierRing& r, I m, const GeometricPreparedCap& a,
                       const GeometricPreparedCap& b) : ring(r), subgroup(std::gcd(a.value.subgroup, b.value.subgroup)),
                                                        phase_offset(a.value.phase + b.value.phase), kappa{-r.beta.y / r.h22, (Wide(r.beta.x) + r.beta.y) / r.h22},
                                                        phase_increment(2 * Wide(r.D) / r.h22), nearest{a.value.nearest, b.value.nearest} {
        Wide rhs = identifier_add(m, identifier_add(a.value.residue_index, b.value.residue_index));
        need(rhs % r.h22 == 0, "nonintegral nearest fiber");
        auto n = rhs / r.h22;
        base = identifier_add({identifier_mul(r.bezout_s, n), identifier_mul(r.bezout_t, n)},
                              identifier_add(a.value.nearest_delta, b.value.nearest_delta));
        auto ratio = identifier_mul(base, identifier_conjugate(identifier_point(r.beta)));
        position = identifier_add(ratio[0], identifier_mul(2, ratio[1]));
        if (GeometricCapKey{b.geometric, b.value.residue} < GeometricCapKey{a.geometric, a.value.residue})
            std::swap(nearest[0], nearest[1]);
    }
    Point scalar_at_phase(I input_phase) const {
        const auto phase = positive_mod(input_phase + phase_offset, subgroup);
        const Wide p = position + phase_increment * phase, step = phase_increment * subgroup;
        Wide turn = kappa[0] == 0 ? identifier_ceil(-Wide(ring.D) - p, step) : identifier_floor(Wide(ring.D) - 1 - p, step);
        auto t = identifier_add(phase, identifier_mul(turn, subgroup));
        return identifier_narrow(identifier_add(base, {identifier_mul(t, kappa[0]), identifier_mul(t, kappa[1])}));
    }
};

inline Point census_ordinary_beta(Point q) {
    for (int i = 0; i < 6; ++i) {
        if (q.x > 0 && q.y >= 0) return q.y == 0 || q.x <= q.y ? q : Point{q.y, q.x};
        q = rotate(q, 1);
    }
    throw std::logic_error("vanishing height of a positive period");
}

// Quotient the closed height hexagon by the scalar progression of a null
// witness. A coset contains either a strictly smaller height/key, a short
// list of equal keys, or no competitor. The Euclidean division also gives
// the lift index required to recover the original null witness at a tie.
class CensusHeightQuotient {
    struct Hit {
        bool smaller = false;
        std::vector<std::pair<Point, I>> ties;
    };
    unsigned facet_;
    I step_, transverse_, lo_, hi_, orientation_;
    std::vector<Hit> hits_;

    std::tuple<Wide, Wide, Wide> project(IdentifierPoint value) const {
        Wide along = facet_ == 2 ? identifier_add(value[0], value[1]) : value[facet_];
        Wide shift = identifier_floor(along, step_);
        return {along - shift * step_, identifier_sub(value[facet_ == 1 ? 0 : 1], identifier_mul(shift, transverse_)),
                shift};
    }
    std::size_t index(I residue, I transverse) const { return static_cast<std::size_t>(transverse - lo_) * step_ + residue; }

  public:
    CensusHeightQuotient(IdentifierPoint slope, Point beta) {
        CensusHexLine line(slope);
        facet_ = line.facet;
        orientation_ = line.direction < 0 ? -1 : 1;
        step_ = identifier_narrow(orientation_ * line.direction);
        transverse_ = identifier_narrow(orientation_ * slope[facet_ == 1 ? 0 : 1]);
        const I h = beta.x + beta.y;
        // The strongest facet has |transverse| <= step. Hence reducing a
        // point of the hexagon moves its transverse coordinate by <=h+step.
        lo_ = -2 * h - step_;
        hi_ = -lo_;
        hits_.resize(static_cast<std::size_t>(hi_ - lo_ + 1) * step_);
        for (I x = -h; x <= h; ++x)
            for (I y = std::max(-h, -h - x); y <= std::min(h, h - x); ++y) {
                Point q{x, y};
                I height = hex_height(q);
                const auto [r, t, j] = project({x, y});
                auto& hit = hits_[index(identifier_narrow(r), identifier_narrow(t))];
                if (height < h || census_ordinary_beta(q) < beta)
                    hit.smaller = true;
                else if (census_ordinary_beta(q) == beta)
                    hit.ties.emplace_back(q, identifier_narrow(j));
            }
    }
    template <class Visit> bool visit(IdentifierPoint value, Visit visit_tie) const {
        const auto [r, t, j] = project(value);
        if (t < lo_ || t > hi_) return true;
        const auto& hit = hits_[index(static_cast<I>(r), static_cast<I>(t))];
        if (hit.smaller) return false;
        for (const auto& [q, k] : hit.ties)
            if (!visit_tie(q, orientation_ * (k - j))) return false;
        return true;
    }
};

// The gamma=1 nulls are exactly swaps following independent translations
// of the two caps, up to central translation. If the lifted cap norms are
// n1,n2, the scalar before the swap has size key H*p+b*q=m+n1+n2.
// Its phase is the sum of the input and both translation phases, modulo
// the joint cap stabilizer. Project every lift at once onto these keys.
class CensusUnitProjection {
  public:
    struct UnitImage { I index = std::numeric_limits<I>::max(); U phases = 0; };
    struct Mask {
        U reject = 0, tie = 0;
    };

  private:
    struct Lift {
        I norm, phase;
        std::vector<Cap> caps;
        std::unique_ptr<std::array<UnitImage, 6>> images;
    };
    I m_, h_;
    bool projectable_;
    std::vector<std::vector<Mask>> rows_;
    std::vector<std::vector<Lift>> lifts_;
    std::vector<std::vector<std::pair<Point, I>>> ties_;
  public:
    CensusUnitProjection(I m, const IdentifierRing& ring, const std::vector<CensusHeightCap>& inventory,
                         bool primitive_batch = false)
        : m_(m), h_(identifier_narrow(hex_height(ring.beta))),
          projectable_(h_ <= 63 || (primitive_batch && ring.h22 == 1)) {
        need(m > 0 && ring.beta.x > 0 && (ring.beta.y == 0 || ring.beta.y >= ring.beta.x),
             "invalid unit projection size or height");
        // Primitive batches encode sizes in words and have only one phase.
        // Other users retain the scalar-line traversal beyond a word of phases.
        if (!projectable_ || inventory.empty() || m > identifier_mul(h_, h_))
            return;
        auto minimum = std::min_element(inventory.begin(), inventory.end(), [](const auto& a, const auto& b) {
                           return a.nearest_norm < b.nearest_norm;
                       })->nearest_norm;
        const Wide square = identifier_mul(h_, h_);
        if (identifier_add(m, identifier_mul(2, minimum)) > square)
            return;
        need(square < std::numeric_limits<I>::max(), "unit projection rows exceed native indices");
        rows_.resize(ring.h22 + 1);
        ties_.resize(h_ * h_ + 1);
        for (I g = 1; g <= ring.h22; ++g)
            if (ring.h22 % g == 0)
                rows_[g].resize(h_ * h_ + 1);
        const Point root = rotate(ring.beta, -1);
        for (I x = -h_; x <= h_; ++x)
            for (I y = -h_; y <= h_; ++y) {
                Point u{x, y};
                auto H = hex_height(u);
                const Wide R = Wide(h_) * x + Wide(ring.beta.y) * y;
                if (R < m || H > h_)
                    continue;
                auto simple = census_ordinary_beta(u);
                bool smaller = H < h_ || simple < ring.beta;
                bool tie = !smaller && simple == ring.beta;
                // The Coxeter word (r0 r1 r2 r3)^5 rotates each cap by zeta.
                // For u in the unit orbit of beta, the swap therefore preserves
                // both horizontal orbits. Its scalars before and after the swap
                // have zero phase modulo the content of beta, so this tie is
                // the same parabolic class. Reflected ties still need comparison.
                if (tie)
                    for (int k = 0; k < 6; ++k)
                        if (rotate(u, k) == ring.beta) {
                            tie = false;
                            break;
                        }
                if (!smaller && !tie && u != root)
                    continue;

                I p = ring.phase(identifier_point(u));
                if (tie && u != root)
                    ties_[R].emplace_back(u, p);
                for (I g = 1; g <= ring.h22; ++g)
                    if (!rows_[g].empty()) {
                        auto& cell = rows_[g][R];
                        U bit = U{1} << (p % g);
                        if (smaller || u == root)
                            cell.reject |= bit;
                        if (tie)
                            cell.tie |= bit;
                    }
            }
        I limit = identifier_narrow(h_ * h_ - m - minimum);
        for (const auto& cap : inventory) {
            auto& signatures = lifts_.emplace_back();
            if (cap.nearest_norm > limit)
                continue;
            std::map<std::pair<I, I>, std::vector<Cap>> by_signature;
            CapBall ball(ring.beta, cap.pooled, 1, limit, 1);
            ball.enumerate([&](const Cap& y) {
                Cap tau;
                for (unsigned j = 0; j < 4; ++j)
                    tau[j] = identifier_narrow(identifier_sub(IdentifierPoint{}, identifier_point(y[j])));
                auto delta = geometric_translation_delta(tau, cap.pooled, ring.beta);
                Wide norm = cap.norm + Wide(h_) * delta[0] + Wide(ring.beta.y) * delta[1];
                need(norm >= 0 && norm <= limit, "lifted cap outside unit projection");
                Cap lifted;
                for (unsigned j = 0; j < 4; ++j)
                    lifted[j] = identifier_narrow(identifier_sub(identifier_point(cap.pooled[j]),
                                                                 identifier_mul(identifier_point(ring.beta), identifier_point(y[j]))));
                by_signature[{identifier_narrow(norm), ring.phase(delta)}].push_back(lifted);
            });
            for (auto& [key, caps] : by_signature)
                signatures.push_back({key.first, key.second, std::move(caps), nullptr});
        }
    }
    // These unary keys determine mask() for every partner and subgroup.
    // Tied image comparisons still require the actual lifted caps.
    std::vector<std::pair<I,I>> mask_key(std::size_t i) const {
        std::vector<std::pair<I,I>> key;
        if (!rows_.empty())
            for (const auto& lift : lifts_[i]) key.emplace_back(lift.norm, lift.phase);
        return key;
    }
    Mask mask(std::size_t i, std::size_t j, I g) const {
        if (rows_.empty())
            return {};
        Mask result;
        U all = (U{1} << g) - 1;
        auto rotate_back = [&](U bits, I turns) { return ((bits >> turns) | (bits << (g - turns))) & all; };
        for (const auto& a : lifts_[i])
            for (const auto& b : lifts_[j]) {
                const Wide r = Wide(m_) + a.norm + b.norm;
                if (r > h_ * h_)
                    break;
                auto w = rows_[g][r];
                auto shift = (g == 1 ? 0 : (a.phase + b.phase) % g);
                result.reject |= rotate_back(w.reject, shift);
                if (result.reject == all)
                    return result;
                result.tie |= rotate_back(w.tie, shift);
            }
        return result;
    }
    bool projected() const { return projectable_; }
    std::span<const Mask> primitive_rows() const {
        return rows_.empty() ? std::span<const Mask>{} : std::span<const Mask>(rows_[1]);
    }
    // The tables are indexed by m + lift norms. Increasing m can only
    // discard lifts and low rows built at the original size.
    void advance_size(I m) {
        need(m >= m_, "unit projection sizes must increase");
        m_ = m;
    }
    struct TieRule { Point beta; std::size_t left, right; };
    template <class Visit> bool visit_tie_rules_at(I m, std::size_t i, std::size_t j, I g, I phase, Visit visit) const {
        if (ties_.empty())
            return true;
        for (std::size_t x = 0; x < lifts_[i].size(); ++x)
            for (std::size_t y = 0; y < lifts_[j].size(); ++y) {
                const auto& a = lifts_[i][x];
                const auto& b = lifts_[j][y];
                const Wide r = Wide(m) + a.norm + b.norm;
                if (r > h_ * h_)
                    break;
                if (ties_[r].empty())
                    continue;
                auto lifted = (g == 1 ? 0 : (phase + a.phase + b.phase) % g);
                for (const auto& [u, p] : ties_[r])
                    if (p % g == lifted && !visit(TieRule{u,x,y})) return false;
            }
        return true;
    }
    template <class Visit> bool visit_tie_rules(std::size_t i, std::size_t j, I g, I phase, Visit visit) const {
        return visit_tie_rules_at(m_, i, j, g, phase, visit);
    }
    std::vector<TieRule> tie_rules(std::size_t i, std::size_t j, I g, I phase) const {
        std::vector<TieRule> rules;
        visit_tie_rules(i,j,g,phase,[&](TieRule rule){rules.push_back(rule);return true;});
        return rules;
    }
    static int reflected_turn(Point u, Point beta) {
        Point image = rotate(conjugate(u), 1);
        for (int t = 0; t < 6; ++t, image = rotate(image, 1))
            if (image == beta) return t;
        throw std::logic_error("unit tie has no reflected orientation");
    }

    // For independent lift lists the least sorted pair consists of their
    // two least cap indices. Only equal minima can contribute a phase tie.
    template<class Canonical> UnitImage reflected_minimum(std::size_t i, std::size_t lift, int turn,
                                                          I content, I input_phase, Canonical canonical) {
        auto &group = lifts_[i][lift];
        if (!group.images) group.images = std::make_unique<std::array<UnitImage, 6>>();
        auto &best = (*group.images)[turn];
        if (best.index != std::numeric_limits<I>::max()) return best;
        UnitImage result;
        for (const auto &cap : group.caps) {
            Cap image;
            for (unsigned k = 0; k < 4; ++k) image[k] = rotate(conjugate(cap[3 - k]), turn);
            const auto [id, phase] = canonical(image);
            if (id < result.index) result = {id, 0};
            if (id == result.index)
                result.phases |= U{1} << positive_mod(phase - input_phase, content);
        }
        best = result;
        return best;
    }

};

// The finite cap group permutes the 40 root lines and rotates each
// functional by a unit. Thus the multiset of their values modulo beta
// and the six units is a cap-orbit invariant. A zero value is exactly
// a horizontal mirror. A unique inventory signature certifies the
// minimum pooled cap; ambiguous signatures use the exact normalizer.
class CensusRootSignature {
    const IdentifierRing& ring_;
    I height_, radius_ = 0, stride_ = 0;
    std::vector<U> labels_;
    std::vector<Point> lifts_;
    std::unordered_map<U, const CensusHeightCap*> classes_;
    static constexpr std::size_t root_count = 40;
    // Flatten the two coordinates of every root value directly into the
    // label grid. Columns are contiguous so the 40 independent linear
    // forms can be evaluated together, with no intermediate cap vectors.
    std::array<std::array<std::uint32_t, root_count>, 8> grid_columns_{};

    static U label(I id) {
        U value = static_cast<U>(id) + 0x9e3779b97f4a7c15ULL;
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }
    U signature(Cap cap) const {
        for (auto& p : cap)
            if (p.x < -height_ || p.x > height_ || p.y < -height_ || p.y > height_) {
                auto q = ring_.reduce(identifier_point(p));
                p = lifts_[q.y * ring_.h11 + q.x];
            }
        U hash = 0;
        if (labels_.size() <= std::numeric_limits<std::uint32_t>::max()) {
            std::array<std::uint32_t, root_count> indices;
            indices.fill(static_cast<std::uint32_t>(radius_ * (stride_ + 1)));
            for (unsigned j = 0; j < 4; ++j) {
                const auto x = static_cast<std::uint32_t>(cap[j].x);
                const auto y = static_cast<std::uint32_t>(cap[j].y);
                for (std::size_t k = 0; k < root_count; ++k)
                    indices[k] += x * grid_columns_[2*j][k] + y * grid_columns_[2*j+1][k];
            }
            // Unsigned intermediates can wrap. Each final grid index is
            // nonnegative and fits in 32 bits by the radius construction.
            for (auto index : indices) {
                const U value = labels_[index];
                if (!value) return 0;
                hash += value;
            }
            return hash ? hash : 1;
        }
        for (const auto& row : cap_root_functionals()) {
            Point q{};
            for (unsigned j = 0; j < 4; ++j)
                q = q + multiply(row[j], cap[j]);
            U u = labels_[(q.y + radius_) * stride_ + q.x + radius_];
            if (!u) return 0;
            hash += u;
        }
        return hash ? hash : 1;
    }

  public:
    CensusRootSignature(const IdentifierRing& ring, const std::vector<CensusHeightCap>& caps)
        : ring_(ring), height_(ring.beta.x + ring.beta.y), lifts_(ring.D) {
        if (caps.empty()) return;
        std::vector<U> residues(ring.D);
        for (I x = -height_; x <= height_; ++x)
            for (I y = -height_; y <= height_; ++y) {
                auto r = ring.reduce(IdentifierPoint{x, y});
                lifts_[r.y * ring.h11 + r.x] = {x, y};
            }
        for (I i = 0; i < ring.D; ++i) {
            Point p = lifts_[i];
            I id = i;
            for (int j = 1; j < 6; ++j) {
                p = rotate(p, 1);
                auto r = ring.reduce(identifier_point(p));
                id = std::min(id, r.y * ring.h11 + r.x);
            }
            residues[i] = id ? std::max(U{1}, label(id)) : 0;
        }
        for (const auto& row : cap_root_functionals()) {
            I rx = 0, ry = 0;
            for (auto p : row) {
                rx += std::abs(p.x) + std::abs(p.y);
                ry += std::abs(p.y) + std::abs(p.x + p.y);
            }
            radius_ = std::max({radius_, rx * height_, ry * height_});
        }
        stride_ = 2 * radius_ + 1;
        labels_.resize(static_cast<std::size_t>(stride_) * stride_);
        for (I y = -radius_; y <= radius_; ++y) {
            auto q = ring.reduce(IdentifierPoint{-radius_, y});
            for (I x = -radius_; x <= radius_; ++x) {
                labels_[(y + radius_) * stride_ + x + radius_] = residues[q.y * ring.h11 + q.x];
                if (++q.x == ring.h11) q.x = 0;
            }
        }
        const auto& forms = cap_root_functionals();
        need(forms.size() == root_count, "wrong root signature dimension");
        for (std::size_t k = 0; k < root_count; ++k)
            for (unsigned j = 0; j < 4; ++j) {
                const auto [a, b] = forms[k][j];
                grid_columns_[2*j][k] = static_cast<std::uint32_t>(a + stride_ * b);
                grid_columns_[2*j+1][k] = static_cast<std::uint32_t>(stride_ * (a+b) - b);
            }
        classes_.reserve(caps.size());
        for (const auto& cap : caps) {
            auto id = signature(cap.pooled);
            need(id != 0, "root mirror in open cap inventory");
            auto [it, inserted] = classes_.emplace(id, &cap);
            if (!inserted) it->second = nullptr;
        }
    }
    std::pair<bool, const CensusHeightCap*> classify(const Cap& cap) const {
        auto id = signature(cap);
        if (!id) return {true, nullptr};
        auto pos = classes_.find(id);
        need(pos != classes_.end(), "open cap signature absent from inventory");
        return {pos->second != nullptr, pos->second};
    }
};

// Numeric owner on the bounded-height base case. The cap inventory gives one
// minimum raw cap per open affine-cap orbit. Its phase stabilizer, together
// with the norm equation, recovers the minimum RAW scalar for a transported
// parabolic class. Thus ties compare in precisely the same order as the
// pooled affine predicate, even though the witness search uses nearest caps.
class CensusPointOwner {
    using Inventory = std::vector<CensusHeightCap>;
    struct Neighbor {
        Cap y;
        Wide cost, q;
        IdentifierPoint projection;
    };
    struct Gamma {
        Point gamma, kappa;
        I divisor, p, q;
        Wide root_radius, null_radius;
        CensusHeightQuotient line;
        IdentifierPoint scalar_factor, energy_factor;
        std::vector<std::optional<std::vector<Neighbor>>> neighbors;
    };
    I m_;
    IdentifierRing ring_;
    IdentifierParabolicCache parabolic_;
    std::map<Cap, GeometricPreparedCap> preparations_;
    std::map<GeometricCapKey, const CensusHeightCap*> minima_;
    const Inventory& inventory_;
    CensusUnitProjection unit_;
    CensusRootSignature signature_;
    struct ConjugateCap {
        const CensusHeightCap* minimum;
        I phase;
    };
    std::vector<ConjugateCap> conjugates_;
    std::vector<Gamma> gammas_;
    using TransportKey = std::tuple<Point, Point, std::array<Cap, 2>>;
    std::map<TransportKey, CensusTransport> transports_;

    const GeometricPreparedCap& prepare(const Cap& cap) {
        const auto known = preparations_.find(cap);
        if (known != preparations_.end()) return known->second;
        return preparations_.emplace(cap, prepare_boundary_cap(cap, ring_, parabolic_)).first->second;
    }
    static GeometricCapKey cap_key(const GeometricPreparedCap& cap) {
        return {cap.geometric, cap.value.residue};
    }
    const CensusHeightCap *minimum(const Cap &cap) {
        auto [known, p] = signature_.classify(cap);
        if (known) return p;
        auto it = minima_.find(cap_key(prepare(cap)));
        return it == minima_.end() ? nullptr : it->second;
    }
    bool unitary_less(const IdentifierCusp& image, const std::array<Cap, 2> &caps, I input_phase, I ordinal) {
        auto first = minimum(image.caps[0]);
        if (!first) return true;
        auto second = minimum(image.caps[1]);
        if (!second) return true;
        if (second->pooled < first->pooled) std::swap(first, second);
        if (first->pooled != caps[0]) return first->pooled < caps[0];
        if (second->pooled != caps[1]) return second->pooled < caps[1];
        const auto g = std::gcd(first->phase_step, second->phase_step);
        if (g == 1) return false;
        const auto& a = prepare(image.caps[0]);
        const auto& b = prepare(image.caps[1]);
        const auto& x = first->geometric;
        const auto& y = second->geometric;
        need(g == std::gcd(x.value.subgroup, y.value.subgroup), "transport changed scalar stabilizer");
        I phase = positive_mod(ring_.phase(identifier_point(image.u)) + a.value.phase + b.value.phase - x.value.phase - y.value.phase, g);
        const I shift = phase - input_phase;
        return positive_mod(ordinal + (ring_.beta.y ? -shift : shift), g) < ordinal;
    }

    bool unit_phase_owned(I m, std::size_t i, std::size_t j, I step, I phase, I ordinal) {
        return unit_.visit_tie_rules_at(m, i, j, step, phase, [&](CensusUnitProjection::TieRule rule) {
            const int turn = CensusUnitProjection::reflected_turn(rule.beta, ring_.beta);
            auto image = [&](std::size_t id, std::size_t lift) {
                return unit_.reflected_minimum(id, lift, turn, ring_.h22, inventory_[id].geometric.value.phase,
                    [&](const Cap &cap) -> std::pair<I,I> {
                        const auto *p = minimum(cap);
                        if (!p) return {-1, 0};
                        const I phase = inventory_[id].phase_step == 1 ? 0 : prepare(cap).value.phase;
                        return {p - inventory_.data(), phase};
                    });
            };
            const auto x = image(i, rule.left), y = image(j, rule.right);
            const auto [lo, hi] = std::minmax(x.index, y.index);
            if (lo != static_cast<I>(i)) return lo > static_cast<I>(i);
            if (hi != static_cast<I>(j)) return hi > static_cast<I>(j);
            if (step == 1) return true;
            const auto scalar = rotate(conjugate(rotate(ring_.beta, -1)), turn);
            const I offset = ring_.phase(identifier_point(scalar)) - phase;
            for (U a = x.phases; a; a &= a - 1)
                for (U b = y.phases; b; b &= b - 1) {
                    const I shift = offset + std::countr_zero(a) + std::countr_zero(b);
                    if (positive_mod(ordinal + (ring_.beta.y ? -shift : shift), step) < ordinal) return false;
                }
            return true;
        });
    }

    bool rival(const IdentifierCusp& moved, const std::array<Cap, 2> &caps, I phase, I ordinal) {
        for (bool anti : {false, true}) {
            Point beta = anti ? rotate(conjugate(moved.beta), 1) : moved.beta;
            for (int turn = 0; turn < 6; ++turn) {
                if (beta == ring_.beta) {
                    auto orient = [&](Point p) { return rotate(anti ? conjugate(p) : p, turn); };
                    IdentifierCusp v{orient(moved.u), ring_.beta, {}};
                    for (unsigned side = 0; side < 2; ++side)
                        for (unsigned i = 0; i < 4; ++i)
                            v.caps[side][i] = orient(moved.caps[side][anti ? 3 - i : i]);
                    if (unitary_less(v, caps, phase, ordinal)) return true;
                    break;
                }
                beta = rotate(beta, 1);
            }
        }
        return false;
    }

    IdentifierCusp transport(const IdentifierCusp& w, const IdentifierCusp& v) {
        TransportKey key{w.u, w.beta, w.caps};
        auto it = transports_.find(key);
        if (it == transports_.end()) {
            auto move = transport_identifier_cusp(w, {});
            it = transports_.try_emplace(std::move(key), move).first;
        }
        return it->second(v);
    }

    const std::vector<Neighbor>& neighbors(Gamma& g, std::size_t index) {
        auto& slot = g.neighbors[index];
        if (slot) return *slot;
        const auto& cap = inventory_[index].geometric.value.nearest;
        std::vector<Neighbor> out;
        visit_hex_cap_neighbors(ring_.beta, g.gamma, cap, std::max(g.root_radius, g.null_radius),
                                [&](const Cap& y, Wide cost) {
                                    IdentifierPoint projection;
                                    need(identifier_divide(identifier_cap_pair(y, cap), {2, -1}, projection),
                                         "nonintegral horizontal cusp projection");
                                    out.push_back({y, cost, identifier_cap_index(y), projection});
                                    return true;
                                });
        std::stable_sort(out.begin(), out.end(), [](const auto& a, const auto& b) { return a.cost < b.cost; });
        return slot.emplace(std::move(out));
    }

    bool owns_phase(const std::array<Cap, 2> &caps, I phase, I ordinal, const CensusNearestFiber &fiber, std::size_t i, std::size_t j) {
        struct Tie {
            IdentifierCusp w;
            Point coefficient;
        };
        const Point scalar_point = fiber.scalar_at_phase(phase);
        std::optional<Tie> first_tie;
        std::vector<Tie> extra_ties;
        for (auto& g : gammas_) {
            const auto& first = neighbors(g, i);
            if (first.empty()) continue;
            const auto& second = neighbors(g, j);
            if (second.empty()) continue;
            const auto scalar = identifier_mul(g.scalar_factor, identifier_point(scalar_point));
            const Wide radius = std::max(g.root_radius, g.null_radius);
            for (const auto& x : first) {
                if (x.cost + second.front().cost > radius) break;
                for (const auto& y : second) {
                    const Wide cost = x.cost + y.cost;
                    if (cost > radius) break;
                    // Divide the Hermitian pairing by (2-zeta) before
                    // pairing the lists. Both horizontal quotients are
                    // integral, and (1+zeta)/(2-zeta)=zeta. The same offset
                    // serves both the root equation and the null line.
                    const auto offset = identifier_sub(scalar, identifier_add(x.projection, y.projection));
                    const Wide energy = x.q + y.q;
                    // The horizontal identity at zero pairing and root index -1
                    // gives cost = N(beta) - m*N(gamma).
                    if (cost == g.root_radius) {
                        IdentifierPoint abar;
                        if (identifier_divide(identifier_sub(IdentifierPoint{}, offset), identifier_point(ring_.beta), abar)) {
                            auto value = identifier_mul(identifier_conjugate(abar), identifier_conjugate(identifier_point(g.gamma)));
                            if (value[0] - energy == -1) return false;
                        }
                    }
                    if (cost > g.null_radius || energy % g.divisor != 0) continue;
                    auto e = energy / g.divisor;
                    auto base = identifier_add(offset, {identifier_mul(g.energy_factor[0], e), identifier_mul(g.energy_factor[1], e)});
                    if (!g.line.visit(base, [&](Point coefficient, Wide j) {
                            IdentifierPoint u{identifier_mul(g.p, e), identifier_mul(g.q, e)};
                            IdentifierCusp w{identifier_narrow(identifier_add(u, identifier_mul(identifier_point(g.kappa), {j, 0}))),
                                             g.gamma,
                                             {x.y, y.y}};
                            // If w=d*w0 for a nonunit d, the primitive coefficient has
                            // H at most 2*H(w)/(sqrt(3)*|d|) <= 2*H(w)/3.
                            if (g.gamma != Point{1, 0} && !primitive_identifier_cusp(w)) return false;
                            if (!first_tie) first_tie.emplace(w, coefficient);
                            else extra_ties.push_back({w, coefficient});
                            return true;
                        })) return false;
                }
            }
        }
        if (!first_tie) return true;
        const IdentifierCusp v{scalar_point, ring_.beta, fiber.nearest};
        auto better = [&](const Tie &tie) {
            const auto &w = tie.w;
            auto moved = w.beta == Point{1, 0} ? census_unit_transport(v, w.caps, tie.coefficient) : transport(w, v);
            return rival(moved, caps, phase, ordinal);
        };
        return !better(*first_tie) && std::none_of(extra_ties.begin(), extra_ties.end(), better);
    }

  public:
    CensusPointOwner(I m, Point beta, const Inventory& inventory, bool primitive_batch = false)
        : m_(m), ring_(beta), inventory_(inventory), unit_(m, ring_, inventory, primitive_batch),
          signature_(ring_, inventory) {
        for (const auto& c : inventory_) {
            const auto& info = c.geometric;
            need(info.value.subgroup == c.phase_step, "cap inventory and phase subgroup disagree");
            need(minima_.emplace(cap_key(info), &c).second, "cap minima have duplicate parabolic keys");
        }
        // Conjugation induces a permutation of the horizontal orbit
        // inventory on the two fixed rays, with an additive scalar phase.
        // Thus an entire unordered cap pair can be rejected before any
        // scalar or competing cusp is examined.
        if (beta.y == 0 || beta.x == beta.y) {
            for (const auto& c : inventory_) {
                Cap mate;
                for (unsigned i = 0; i < 4; ++i)
                    mate[i] = rotate(conjugate(c.pooled[3 - i]), beta.y ? 0 : -1);
                const auto& image = prepare(mate);
                auto k = minima_.find(cap_key(image));
                need(k != minima_.end() && k->second->phase_step == c.phase_step, "conjugation left open cap inventory");
                conjugates_.push_back({k->second, image.value.phase - k->second->geometric.value.phase});
            }
        }
        const auto bound = std::max(hex_gamma_bound(beta, m), Wide(ring_.D / m));
        for (I a = 1; identifier_mul(a, a) <= bound; ++a)
            for (I b = 0; identifier_norm({a, b}) <= bound; ++b) {
                Point gamma{a, b};
                if (gamma == Point{1, 0} && unit_.projected()) continue;
                auto [d, p, q] = bezout(a + b, b);
                Point kappa{-b / d, (a + b) / d};
                auto slope = identifier_mul(identifier_conjugate(identifier_point(kappa)), identifier_point(beta));
                Wide root = ring_.D - identifier_mul(m, identifier_norm({a, b}));
                Wide null = hex_cap_radius(beta, m, gamma);
                auto scalar_factor = identifier_mul({0, 1}, identifier_conjugate(identifier_point(gamma)));
                auto energy_factor = identifier_mul(identifier_conjugate({p, q}), identifier_point(beta));
                if (std::max(root, null) >= 0)
                    gammas_.push_back({gamma, kappa, d, p, q, root, null, CensusHeightQuotient(slope, beta), scalar_factor,
                                       energy_factor, std::vector<std::optional<std::vector<Neighbor>>>(inventory.size())});
            }
    }

    // At fixed beta both search radii decrease by (m'-m) N(gamma).
    // Previously generated neighbour lists remain sorted supersets; emit()
    // already intersects them with the current radius. The cap orbits,
    // height quotients, tie transports and unit mask rows are size independent.
    void set_size(I m) {
        need(m >= m_, "point owner sizes must increase");
        if (m == m_) return;
        const auto bound = std::max(hex_gamma_bound(ring_.beta, m), Wide(ring_.D / m));
        for (auto &g : gammas_) {
            const Wide n = identifier_norm(identifier_point(g.gamma));
            const Wide decrease = (Wide(m) - m_) * n;
            g.root_radius -= decrease;
            g.null_radius -= decrease;
        }
        std::erase_if(gammas_, [&](const Gamma &g) {
            return std::max(g.root_radius, g.null_radius) < 0 || identifier_norm(identifier_point(g.gamma)) > bound;
        });
        unit_.advance_size(m);
        m_ = m;
    }

    bool has_unit_projection() const { return unit_.projected(); }
    std::vector<Wide> primitive_unit_masses(std::span<const I> sizes);

    template <bool selecting = false, class Visit> void emit(Visit visit) {
        if (!small_height_detail::admissible(m_, ring_.beta))
            return;
        // The first nonunit norm in the Eisenstein lattice is three.
        CensusWitnessProducts<CensusHeightQuotient> products(m_, ring_, inventory_, {1,1});
        small_height_detail::groups<true, true>(
            m_, ring_.beta, inventory_, [&](std::size_t i, std::size_t j, I phase_step, const auto &make_window) {
                std::array<Cap, 2> caps{inventory_[i].pooled, inventory_[j].pooled};
                bool fixed = false;
                I delta = 0, g = 1;
                if (!conjugates_.empty()) {
                    const auto& a = conjugates_[i];
                    const auto& b = conjugates_[j];
                    auto [first, second] = std::minmax(a.minimum, b.minimum);
                    if (std::tie(first->pooled, second->pooled) < std::tie(caps[0], caps[1]))
                        return;
                    fixed = first->pooled == caps[0] && second->pooled == caps[1];
                    if (fixed) {
                        delta = a.phase + b.phase;
                        g = std::gcd(first->phase_step, second->phase_step);
                    }
                }
                auto projection = unit_.mask(i, j, phase_step);
                if (phase_step < 64 && projection.reject == ((U{1} << phase_step) - 1))
                    return;
                const auto scalars = make_window();
                if (!scalars.size()) return;
                std::size_t first = i, second = j;
                std::optional<CensusNearestFiber> fiber;
                auto publish = [&](const auto &window) {
                    if constexpr (selecting) visit(window, inventory_[i], inventory_[j]);
                    else for (Point u : window) {
                        IdentifierCusp v{u, ring_.beta, caps};
                        if constexpr (std::is_invocable_v<Visit&, const IdentifierCusp&, const CensusHeightCap&, const CensusHeightCap&>)
                            visit(v, inventory_[i], inventory_[j]);
                        else visit(v);
                    }
                };
                bool clear = true;
                for (auto &gamma : gammas_) {
                    const auto &left = neighbors(gamma, i);
                    if (left.empty()) continue;
                    const auto &right = neighbors(gamma, j);
                    if (!right.empty() && left.front().cost + right.front().cost <=
                            std::max(gamma.root_radius, gamma.null_radius)) {
                        clear = false;
                        break;
                    }
                }
                if (clear && (!fixed || g == 1) && !(projection.reject || projection.tie)) {
                    publish(scalars);
                    return;
                }
                auto surviving = [&](I lo, I hi, bool tied, I phase) {
                    if (tied && !unit_phase_owned(m_, i, j, phase_step, phase, ring_.beta.y ? phase_step - 1 - lo : lo))
                        return;
                    if (clear) {
                        publish(scalars.run(lo, hi));
                        return;
                    }
                    if (!fiber) {
                        if (inventory_[second].geometric.value.residue < inventory_[first].geometric.value.residue)
                            std::swap(first, second);
                        fiber.emplace(ring_, m_, inventory_[first].geometric, inventory_[second].geometric);
                    }
                    // Decide rival inequalities in phase coordinates. Only
                    // accepted ordinal intervals reach the period constructor.
                    I start = lo;
                    for (I k = lo; k < hi; ++k) {
                        const I ordinal = ring_.beta.y ? phase_step - 1 - k : k;
                        if (!owns_phase(caps, phase, ordinal, *fiber, first, second)) {
                            if (start < k) publish(scalars.run(start, k));
                            start = k + 1;
                        }
                        if (++phase == phase_step) phase = 0;
                    }
                    if (start < hi) publish(scalars.run(start, hi));
                };
                // Phase increases by one along the primitive scalar step.
                // Rotate the unit masks into ordinal order once for the pair.
                const auto phase0 = phase_step == 1 || (clear && !fixed && !(projection.reject || projection.tie))
                    ? 0 : ring_.phase({scalars.p, scalars.q}) % phase_step;
                auto projected = [&](I lo, I hi) {
                    if (!(projection.reject || projection.tie)) {
                        surviving(lo, hi, false, phase_step == 1 ? 0 : (phase0 + lo) % phase_step);
                        return;
                    }
                    need(phase_step < 64, "projected unit mask exceeds word");
                    const U all = (U{1} << phase_step) - 1;
                    auto ordinal_mask = [&](U bits) { return ((bits >> phase0) | (bits << (phase_step - phase0))) & all; };
                    U allowed = ((U{1} << hi) - (U{1} << lo)) & ~ordinal_mask(projection.reject);
                    const U tied = ordinal_mask(projection.tie);
                    while (allowed) {
                        const I start = std::countr_zero(allowed);
                        const bool tie = (tied >> start) & 1;
                        const I length = tie ? 1 : std::countr_one((allowed & ~tied) >> start);
                        surviving(start, start + length, tie, (phase0 + start) % phase_step);
                        allowed &= ~(((U{1} << length) - 1) << start);
                    }
                };
                if (!fixed || g == 1) {
                    scalars.intersect(0, phase_step, projected);
                    return;
                }
                // Reflection sends ordinal k to (c-k) mod g. Retain the
                // smaller ordinal on the real ray and the larger on a skew ray.
                const auto mirrored = rotate(conjugate(identifier_narrow({scalars.p, scalars.q})), ring_.beta.y ? 0 : -1);
                const I c = positive_mod(ring_.phase(identifier_point(mirrored)) + delta - phase0, g);
                scalars.reflected(c, ring_.beta.y != 0, projected);
            }, products.active() ? &products : nullptr);
    }
};

template <class Visit> void census_point_owners(I m, Point beta, Visit visit) {
    auto caps = census_small_height_inventory(beta);
    if (caps.empty()) return;
    CensusPointOwner owner(m, beta, caps);
    owner.emit(std::move(visit));
}

} // namespace gf::fullerene_rank
