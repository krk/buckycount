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
#include "census_period_graph.hpp"
#include "census_ordered_stream.hpp"
#include "census_family_table.hpp"
#include "census_pentagon_contacts.hpp"
#include "census_contact_phases.hpp"
#include "census_rank_selector.hpp"
#include "star_dense_decoder.hpp"
#include <deque>
#include <charconv>
#include <chrono>
#include <ctime>
#include <thread>

#include "census_heap_scopes.hpp"


namespace {
using namespace gf;
using namespace gf::fullerene_rank;

// The owner iterator emits a cap pair's accepted scalars consecutively.
template<class Emit> CensusContactHistogram visit_family_contacts(I m, Point beta, Emit emit) {
    CENSUS_HEAP_SCOPE(contact);
    CensusContactHistogram histogram{};
    if (!small_height_detail::admissible(m, beta)) return histogram;
    const auto caps = census_small_height_inventory(beta);
    if (caps.empty()) return histogram;
    CensusPointOwner owner(m, beta, caps);
    CensusContactPhases phases(m, beta);
    std::vector<std::optional<CensusContactPhases::Prepared>> prepared(caps.size());
    auto cap_data = [&](std::size_t i) -> const auto & {
        if (!prepared[i]) prepared[i].emplace(phases.prepare(caps[i].pooled));
        return *prepared[i];
    };
    std::size_t first = 0, second = 0;
    Point base{};
    std::uint64_t mask = 0;
    auto flush = [&] {
        if (!mask) return;
        const auto bins = phases.histogram(cap_data(first), cap_data(second), base, mask);
        for (std::size_t i = 0; i < bins.size(); ++i)
            if (bins[i]) histogram[i] += static_cast<unsigned long>(bins[i]);
        mask = 0;
    };
    owner.emit([&](const auto &v, const auto &x, const auto &y) {
        const auto i = static_cast<std::size_t>(&x - caps.data());
        const auto j = static_cast<std::size_t>(&y - caps.data());
        if (!mask || i != first || j != second || Wide(v.u.y) - base.y >= 64) {
            flush();
            first = i;
            second = j;
            base = v.u;
        }
        const auto offset = Wide(v.u.y) - base.y;
        need(offset >= 0 && offset < 64, "out-of-order contact phase");
        const auto bit = std::uint64_t{1} << static_cast<unsigned>(offset);
        need(!(mask & bit), "duplicate contact phase");
        mask |= bit;
        emit(v, x, y);
    });
    flush();
    return histogram;
}

mpz_class family_count(I m, Point beta) {
    if (const auto mass = CENSUS_HEAP_EVAL(family_mass, census_orbit_family_mass(m, beta))) return *mass;
    CENSUS_HEAP_SCOPE(owner);
    const auto caps = census_small_height_inventory(beta);
    if (caps.empty()) return 0;
    CensusPointOwner owner(m, beta, caps);
    mpz_class count = 0;
    owner.emit<true>([&](const auto &window, const auto &, const auto &) {
        count += static_cast<unsigned long>(window.size());
    });
    return count;
}

using FamilyEnds = std::vector<std::pair<Point, mpz_class>>;

void fill_family_table(I m, const mpz_class &expected, CensusFamilyTable *report, FamilyEnds *ends) {
    mpz_class total = 0;
    if (expected != 0)
        visit_census_families(m, [&](Point beta) {
            const mpz_class begin = total;
            if (report) {
                const auto histogram = visit_family_contacts(m, beta, [](const auto &, const auto &, const auto &) {});
                for (const auto &count : histogram) total += count;
                report->record_cumulative(beta, total, histogram);
            } else total += family_count(m, beta);
            if (ends && total > begin) ends->emplace_back(beta, total);
            return total != expected;
        });
    need(total == expected, "family counts disagree with ordinary isomer count");
}

mpz_class parse_rank(std::string_view s) {
    CENSUS_HEAP_SCOPE(parse_rank);
    if (s.empty() || s.find_first_not_of("0123456789") != std::string_view::npos)
        throw std::invalid_argument("expected nonnegative decimal rank");
    return mpz_class(std::string(s), 10);
}

I parse_vertices(std::string_view s) {
    const auto n = parse_rank(s);
    need(n >= 20 && n <= 2UL * ordinary_rank_max_index() && mpz_even_p(n.get_mpz_t()),
         "expected an even vertex count of at least 20 within the arithmetic range");
    return static_cast<I>(n.get_ui());
}

Point parse_family(std::string_view a_text, std::string_view b_text, I vertices) {
    const auto a = parse_rank(a_text), b = parse_rank(b_text);
    const auto maximum = identifier_sqrt(3 * Wide(vertices) - 1);
    need(a > 0 && (b == 0 || a <= b) && a + b <= static_cast<unsigned long>(maximum),
         "height outside ordinary wedge");
    return {static_cast<I>(a.get_ui()), static_cast<I>(b.get_ui())};
}

void write_period(std::ostream &output, I n, const IdentifierCusp &period, star_dense::Workspace &work) {
    CENSUS_HEAP_SCOPE(encode);
    write_graph(output, decode_period(n, period, work));
}

void write_planar_header(I vertices, std::ostream &out = std::cout) {
    out << (vertices <= 255 ? ">>planar_code<<" : ">>planar_code le<<");
}

std::size_t graph_batch_size(I n) {
    const I bytes = n <= 255 ? 1 + 4 * n : 3 + 8 * n;
    return static_cast<std::size_t>(std::max<I>(1, 131072 / bytes));
}


} // namespace
