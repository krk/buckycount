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

#include "genfullerene.hpp"
#include <cerrno>
#include <cstdlib>
#include <gmpxx.h>
#include <system_error>
#include <sys/stat.h>
#include <unistd.h>

namespace gf::fullerene_rank {

// Each of twelve pentagons has five sides; a contact uses two of those sides.
using CensusContactHistogram = std::array<mpz_class, 31>;

// Ordinary coefficient order is by hexagonal height, with the real ray last.
class CensusFamilyCursor {
    I maximum_, h_ = 1, a_ = 1;
public:
    explicit CensusFamilyCursor(I m) : maximum_(identifier_narrow(identifier_sqrt(6 * Wide(m) - 1))) {}
    std::optional<Point> next() {
        if (h_ > maximum_) return std::nullopt;
        if (a_ <= h_ / 2) {
            const I a = a_++;
            return Point{a, h_ - a};
        }
        a_ = 1;
        return Point{h_++, 0};
    }
};

template<class Visit> void visit_census_families(I m, Visit visit) {
    CensusFamilyCursor families(m);
    while (const auto beta = families.next())
        if (!visit(*beta)) break;
}

class CensusFamilyTable {
    struct Observation {
        mpz_class end;
        CensusContactHistogram contacts;
    };
    I vertices_;
    mpz_class expected_;
    std::vector<Point> labels_;
    std::vector<Observation> observed_;

public:
    CensusFamilyTable(I m, mpz_class expected)
        : vertices_(identifier_narrow(2 * Wide(m))), expected_(std::move(expected)) {
        need(expected_ >= 0, "negative expected census total");
        visit_census_families(m, [&](Point beta) {
            labels_.push_back(beta);
            return true;
        });
    }

    void record_cumulative(Point beta, const mpz_class& end, const CensusContactHistogram& contacts) {
        need(observed_.size() < labels_.size() && labels_[observed_.size()] == beta,
             "out-of-order or duplicate census family");
        mpz_class cumulative = observed_.empty() ? mpz_class(0) : observed_.back().end;
        for (const auto& count : contacts) {
            need(count >= 0, "negative contact bin");
            cumulative += count;
        }
        need(end == cumulative && end <= expected_, "family contact bins disagree with the cumulative count");
        observed_.push_back({end, contacts});
    }

    void write_csv(std::ostream& out) const {
        need(observed_.empty() ? expected_ == 0 : observed_.back().end == expected_, "incomplete family table");
        out << "vertices,height,a,b,graphs,rank_begin,rank_end,source";
        const CensusContactHistogram zero{};
        for (std::size_t k = 0; k < zero.size(); ++k) out << ",contacts_" << std::to_string(k);
        out << '\n';
        mpz_class begin = 0;
        for (std::size_t i = 0; i < labels_.size(); ++i) {
            const bool observed = i < observed_.size();
            const auto& end = observed ? observed_[i].end : expected_;
            const mpz_class size = end - begin;
            const auto beta = labels_[i];
            out << std::to_string(vertices_) << ',' << std::to_string(beta.x + beta.y) << ','
                << std::to_string(beta.x) << ',' << std::to_string(beta.y) << ',' << size.get_str()
                << ',' << begin.get_str() << ',' << end.get_str() << ','
                << (observed ? "observed" : "zero_from_total");
            for (const auto& count : observed ? observed_[i].contacts : zero)
                out << ',' << count.get_str();
            out << '\n';
            begin = end;
        }
        out.flush();
        need(bool(out), "family table output failed");
    }

};

} // namespace gf::fullerene_rank
