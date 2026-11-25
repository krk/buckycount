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
#include <algorithm>
#include <cstddef>
#include <bit>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gf::fullerene_rank {
// Exact fixed-witness prefix on a product with separated transported cap keys.
// Indices are the original pooled cap ordinals. A negative image is singular.
class CensusTiedCapPrefix {
    using Key = long;
    std::vector<Key> ids_;
    std::vector<std::size_t> self_;
    std::vector<Key> alphabet_, fixed_ids_;
    struct Level {
        std::vector<std::uint64_t> bits;
        std::vector<std::size_t> prefix;
        std::size_t zeros;
        std::size_t ones(std::size_t pos) const {
            auto word = pos / 64, bit = pos % 64;
            return prefix[word] + (bit ? std::popcount(bits[word] & ((std::uint64_t{1} << bit) - 1)) : 0);
        }
    };
    std::vector<Level> levels_;

    std::size_t less(std::size_t lo, std::size_t hi, std::size_t bound) const {
        if (bound >= alphabet_.size()) return hi - lo;
        std::size_t result = 0, depth = levels_.size();
        for (const auto &level : levels_) {
            auto l = level.ones(lo), r = level.ones(hi);
            if ((bound >> --depth) & 1) {
                result += (hi - lo) - (r - l);
                lo = level.zeros + l;
                hi = level.zeros + r;
            } else {
                lo -= l;
                hi -= r;
            }
        }
        return result;
    }

  public:
    explicit CensusTiedCapPrefix(const std::vector<std::pair<Key,Key>>& right) : self_{0} {
        for (auto [j,y] : right) {
            if (j < 0 || (!ids_.empty() && ids_.back() >= j))
                throw std::invalid_argument("expected distinct increasing cap indices");
            alphabet_.push_back(y);
            ids_.push_back(j);
            self_.push_back(self_.back() + (y >= j));
            if (j == y) fixed_ids_.push_back(j);
        }
        std::sort(alphabet_.begin(), alphabet_.end());
        alphabet_.erase(std::unique(alphabet_.begin(), alphabet_.end()), alphabet_.end());
        if (alphabet_.empty()) return;
        std::vector<std::size_t> values, next(right.size());
        for (auto [j,y] : right)
            values.push_back(std::lower_bound(alphabet_.begin(), alphabet_.end(), y) - alphabet_.begin());
        for (auto depth = std::bit_width(alphabet_.size() - 1); depth > 0; --depth) {
            auto &level = levels_.emplace_back();
            level.bits.resize((values.size() + 63) / 64);
            for (std::size_t j = 0; j < values.size(); ++j)
                if ((values[j] >> (depth - 1)) & 1)
                    level.bits[j / 64] |= std::uint64_t{1} << (j % 64);
            level.prefix.push_back(0);
            for (auto word : level.bits) level.prefix.push_back(level.prefix.back() + std::popcount(word));
            level.zeros = values.size() - level.prefix.back();
            std::size_t zero = 0, one = level.zeros;
            for (auto value : values) {
                if ((value >> (depth - 1)) & 1) next[one++] = value;
                else next[zero++] = value;
            }
            values.swap(next);
        }
    }
    std::size_t fixed_count(Key i, Key last) const {
        if (last < i) return 0;
        return std::upper_bound(fixed_ids_.begin(), fixed_ids_.end(), last) -
               std::lower_bound(fixed_ids_.begin(), fixed_ids_.end(), i);
    }
    std::size_t count(Key i, Key x, Key last) const {
        if (i < 0) throw std::invalid_argument("negative first cap index");
        if (x < i || last < i) return 0;
        auto begin = std::lower_bound(ids_.begin(), ids_.end(), i);
        auto end = std::upper_bound(begin, ids_.end(), last);
        const auto lo = static_cast<std::size_t>(begin - ids_.begin());
        const auto hi = static_cast<std::size_t>(end - ids_.begin());
        if (x == i) return self_[hi] - self_[lo];
        const auto bounded = std::upper_bound(begin, end, x) - ids_.begin();
        const auto lower = std::lower_bound(alphabet_.begin(), alphabet_.end(), i);
        const std::size_t k = lower - alphabet_.begin();
        const auto after = k + (lower != alphabet_.end() && *lower == i);
        return hi - lo - less(lo, bounded, k) - less(bounded, hi, after);
    }
};
} // namespace gf::fullerene_rank
