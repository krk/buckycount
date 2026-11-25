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

#include "commands.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <fstream>
#include <set>
#include <thread>
#include <unistd.h>
#include <iostream>
#include <stdexcept>

namespace fs = std::filesystem;

static void decode(const std::string& input, std::ostream& out) {
    std::ifstream file;
    if (input != "-") {
        file.open(input, std::ios::binary);
        if (!file) throw std::runtime_error("cannot open " + input);
    }
    std::istream& in = input == "-" ? std::cin : file;
    std::string header;
    while (header.size() < 19 && !header.ends_with("<<")) {
        int c = in.get();
        if (c == EOF) break;
        header += static_cast<char>(c);
    }
    const bool big_endian = header == ">>planar_code be<<";
    if (!big_endian && header != ">>planar_code<<" && header != ">>planar_code le<<")
        throw std::runtime_error("invalid planar-code header");
    for (std::size_t graph = 1;; ++graph) {
        int n = in.get();
        if (n == EOF && in.eof() && !in.bad()) break;
        auto error = [&](const char* reason) { return std::runtime_error("graph " + std::to_string(graph) + ": " + reason); };
        auto byte = [&] {
            int value = in.get();
            if (value == EOF) throw error("truncated frame");
            return value;
        };
        const bool wide = n == 0;
        auto number = [&] {
            int first = byte();
            if (!wide) return first;
            int second = byte();
            return big_endian ? 256 * first + second : first + 256 * second;
        };
        if (wide) n = number();
        if (n <= 0) throw error("invalid vertex count");
        std::vector<std::array<int, 3>> rows(n);
        for (int v = 1; v <= n; ++v) {
            auto& row = rows[v - 1];
            for (auto it = row.begin(); it != row.end(); ++it) {
                int w = number();
                if (!w || w > n || w == v || std::find(row.begin(), it, w) != it)
                    throw error("invalid cubic adjacency");
                *it = w;
            }
            if (number()) throw error("expected three neighbors");
        }
        for (int v = 1; v <= n; ++v)
            for (int w : rows[v - 1])
                if (std::find(rows[w - 1].begin(), rows[w - 1].end(), v) == rows[w - 1].end())
                    throw error("asymmetric adjacency");
        out << "graph " << graph << ", " << n << " vertices\n";
        for (int v = 1; v <= n; ++v) {
            out << v << ':';
            for (int w : rows[v - 1])
                out << ' ' << w;
            out << '\n';
        }
        if (!out) throw error("output write failed");
    }
}

void decode_request(const Request& r) {
    const bool batch = r.inputs.size() > 1;
    std::set<fs::path> inputs, destinations;
    std::vector<fs::path> outputs;
    for (const auto& name : r.inputs) {
        if (batch && name == "-") throw std::runtime_error("stdin requires a single input");
        if (name != "-") inputs.insert(fs::weakly_canonical(fs::absolute(name)));
    }

    for (const auto& name : r.inputs) {
        fs::path path = r.output;
        if (batch) {
            path = fs::path(name).replace_extension(".adj");
            if (r.output != "-") path = fs::path(r.output) / path.filename();
        }
        if (path != "-") {
            const auto key = fs::weakly_canonical(fs::absolute(path));
            if (inputs.contains(key) || !destinations.insert(key).second)
                throw std::runtime_error("conflicting output: " + path.string());
        }
        outputs.push_back(path);
        std::cerr << "output: " << (path == "-" ? "stdout" : fs::absolute(path).string()) << '\n';
    }
    if (r.estimates) {
        std::cerr << "cache estimate: 0 bytes\n";
    }
    if (!batch && outputs.front() == "-") {
        decode(r.inputs.front(), std::cout);
        if (!std::cout.flush()) throw std::runtime_error("output write failed");
        return;
    }

    std::atomic_size_t next{0};
    std::vector<std::string> errors(outputs.size());
    auto run = [&] {
        for (;;) {
            const auto i = next.fetch_add(1);
            if (i >= outputs.size()) return;
            try {
                write_atomic(outputs[i], [&](std::ostream& out) { decode(r.inputs[i], out); });
            } catch (const std::exception& error) {
                errors[i] = r.inputs[i] + ": " + error.what();
            }
        }
    };
    {
        std::vector<std::jthread> workers;
        for (std::size_t i = 1; i < std::min<std::size_t>(r.workers, outputs.size()); ++i)
            workers.emplace_back(run);
        run();
    }
    std::string failures;
    for (const auto& error : errors)
        if (!error.empty()) failures += error + '\n';
    if (!failures.empty()) throw std::runtime_error(failures);
}
