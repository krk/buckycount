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
#include <filesystem>
#include <functional>
#include <ostream>
#include <iostream>
#include <gmpxx.h>
#include <string>
#include <string_view>
#include <vector>

struct Request {
    std::string command, size, selection, family, format, output = "-";
    std::filesystem::path cache = ".cache";
    bool help = false, estimates = true, use_cache = true, ipr = false, details = false, families = false;
    std::vector<std::string> inputs;
    unsigned workers = 0;
};

int buckycount_main(int argc, char** argv);
mpz_class graph_count(const Request&, const mpz_class& index);
void graph_request(const Request&, std::ostream&);
void decode_request(const Request&);
void atomic_bytes(const std::filesystem::path&, const std::string&);
void write_atomic(const std::filesystem::path&, const std::function<void(std::ostream&)>&);

inline std::string disk_size(const mpz_class& bytes, bool upper_bound = false) {
    if (bytes < 1000) return bytes.get_str() + (bytes == 1 ? " byte" : " bytes");
    constexpr std::string_view units[] = {"kB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB", "RB", "QB"};
    std::size_t unit = 0;
    mpz_class scale = 1000;
    while (unit + 1 < std::size(units) && bytes >= 1000 * scale) {
        ++unit;
        scale *= 1000;
    }
    // Round upper bounds up so the displayed disk allowance stays sufficient.
    const mpz_class bias = upper_bound ? mpz_class(scale - 1) : mpz_class(scale / 2);
    const mpz_class tenths = (10 * bytes + bias) / scale;
    auto digits = tenths.get_str();
    if (digits.back() == '0') digits.pop_back();
    else digits.insert(digits.size() - 1, ".");
    return digits + " " + std::string(units[unit]);
}
