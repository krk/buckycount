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

// Shared options, row evaluation and serialization for exact prefix calculators.
#pragma once

#include "commands.hpp"
#include "prefix_engine_formula.hpp"
#include <sstream>

namespace {

struct Options {
    Integer first = 1, last = 0;
    unsigned workers = 0;
    std::string output, format = "csv";
    bool ipr = false;
};

inline void print_count_cache_estimate(const Options& options, bool disk) {
    if (!disk) {
        std::cerr << "cache estimate: 0 bytes (disabled)\n";
        return;
    }
    const Integer rows = options.last - options.first + 1;
    // Estimate distinct discriminants across a range. Through m=500000 the ordinary
    // prefix uses 911898 discriminants; a single row typically uses six.
    Integer primitives = std::min(Integer(6 * rows), Integer((1824 * options.last + 999) / 1000));
    const Integer first_ipr = std::max(options.first, Integer(30));
    if (options.ipr && first_ipr <= options.last) {
        const Integer ipr_rows = options.last - first_ipr + 1;
        // Shifted IPR arguments overlap across adjacent rows. Short cold
        // runs use about 250 keys at one row, and long prefixes about 11 per row.
        primitives = std::max(primitives, std::min(Integer(250 + 50 * (ipr_rows - 1)),
                                                   Integer(11 * options.last + 120)));
    }
    const Integer width = (mpz_sizeinbase(options.last.get_mpz_t(), 2) + 7) / 8;
    const Integer bytes = primitives * (24 + 2 * width) + ((primitives + 127) / 128) * 20 + 24;
    std::cerr << "cache estimate: " << disk_size(bytes) << " shared coefficients\n";
}

inline unsigned long parse_number(const std::string& text, const std::string& flag, bool zero = false) {
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
        fail(flag + " needs an unsigned decimal integer");
    Integer value;
    if (mpz_set_str(value.get_mpz_t(), text.c_str(), 10) != 0 ||
        !mpz_fits_ulong_p(value.get_mpz_t()) || (!zero && value == 0))
        fail(flag + " is outside the supported integer range");
    return value.get_ui();
}

inline void set_count_interval(Options& options, const std::string& selector) {
    const auto endpoint = [&](const std::string& text, const std::string& label) {
        if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
            fail(label + " needs an unsigned decimal integer");
        const Integer value(text, 10);
        if (value < 1) fail(label + " must be positive");
        return value;
    };
    const std::size_t dots = selector.find("..");
    if (dots == std::string::npos) {
        options.last = endpoint(selector, "endpoint");
        options.first = options.last;
    } else {
        options.first = endpoint(selector.substr(0, dots), "range start");
        options.last = endpoint(selector.substr(dots + 2), "range end");
    }
    if (options.first > options.last) fail("range is empty");
    if (options.first % 2 != 0 || options.last % 2 != 0) fail("vertex endpoints must be even");
    options.first /= 2;
    options.last /= 2;
}

inline void write_row(std::ostream& out, const Row& row, bool json, const Integer& last, char sep) {
    if (!json) {
        out << row.index.get_str() << sep << Integer(2 * row.index).get_str() << sep << row.ordinary.get_str() << sep
            << row.oriented.get_str() << sep << row.achiral.get_str() << sep
            << rational_text(row.split) << sep << rational_text(row.antipodal) << sep
            << rational_text(row.q2) << sep << rational_text(row.q3) << sep
            << rational_text(row.q4) << sep << rational_text(row.order10) << sep
            << rational_text(row.order12);
        if (row.ipr)
            out << sep << row.ipr->ordinary.get_str() << sep
                << row.ipr->oriented.get_str() << sep << row.ipr->achiral.get_str();
        out << '\n';
        return;
    }
    out << "  {\"m\": \"" << row.index.get_str() << "\", \"vertices\": \"" << Integer(2 * row.index).get_str()
        << "\", \"iso\": \"" << row.ordinary.get_str() << '"'
        << ", \"BB_m\": \"" << row.oriented.get_str() << '"'
        << ", \"a_m\": \"" << row.achiral.get_str() << '"'
        << ", \"split\": \"" << rational_text(row.split) << "\""
        << ", \"antipodal\": \"" << rational_text(row.antipodal) << "\""
        << ", \"Q2\": \"" << rational_text(row.q2) << "\""
        << ", \"Q3\": \"" << rational_text(row.q3) << "\""
        << ", \"Q4\": \"" << rational_text(row.q4) << "\""
        << ", \"order10\": \"" << rational_text(row.order10) << "\""
        << ", \"order12\": \"" << rational_text(row.order12) << "\"";
    if (row.ipr)
        out << ", \"ipr\": \"" << row.ipr->ordinary.get_str() << '"'
            << ", \"ipr_oriented\": \"" << row.ipr->oriented.get_str() << '"'
            << ", \"ipr_achiral\": \"" << row.ipr->achiral.get_str() << '"';
    out << (row.index == last ? "}\n" : "},\n");
}

inline void write_prefix(std::ostream& out, const Options& options,
                         bool plain = false, bool range = false) {
    const bool json = options.format == "json";
    const char sep = ',';
    if (plain) {
        if (range) out << "vertices," << (options.ipr ? "ipr" : "iso") << '\n';
    } else if (json) out << "[\n";
    else {
        out << "m";
        for (const auto* column : {"vertices", "iso", "BB_m", "a_m", "split", "antipodal",
                                   "Q2", "Q3", "Q4", "order10", "order12"}) out << sep << column;
        if (options.ipr) out << sep << "ipr" << sep << "ipr_oriented" << sep << "ipr_achiral";
        out << '\n';
    }

    if (!out.flush()) fail("cannot write output " + options.output);
    const Integer length = options.last - options.first + 1;
    const unsigned window = length < 512 ? length.get_ui() : 512;
    const unsigned workers = std::min(worker_count(options.workers), window);
    const unsigned batch = window / workers;
    Integer next = options.first;
    std::mutex mutex;
    std::condition_variable ready;
    std::exception_ptr error;
    auto stop = [&](std::exception_ptr cause) {
        const std::lock_guard lock(mutex);
        if (!error) error = cause;
        ready.notify_all();
    };
    // Each lane retains its arithmetic caches across batches. At most one
    // window is in flight, and only the lane owning the next rows can write.
    auto run = [&](unsigned lane) {
        try {
            std::vector<Row> rows;
            rows.reserve(batch);
            for (Integer first = options.first + lane * batch; first <= options.last;
                 first += workers * batch) {
                rows.clear();
                for (unsigned i = 0; i < batch && first + i <= options.last; ++i) {
                    rows.push_back(evaluate_index(first + i, options.ipr));
                }
                std::unique_lock lock(mutex);
                ready.wait(lock, [&] { return error || next == first; });
                if (error) return;
                for (const auto& row : rows) {
                    if (plain) {
                        if (range) out << Integer(2 * row.index).get_str() << ',';
                        out << (options.ipr ? row.ipr->ordinary : row.ordinary).get_str() << '\n';
                    } else write_row(out, row, json, options.last, sep);
                }
                if (!out.flush()) fail("cannot write output " + options.output);
                next += rows.size();
                ready.notify_all();
            }
        } catch (...) {
            stop(std::current_exception());
        }
    };
    {
        std::vector<std::jthread> threads;
        try {
            for (unsigned lane = 1; lane < workers; ++lane) threads.emplace_back(run, lane);
            run(0);
        } catch (...) {
            stop(std::current_exception());
        }
    }
    if (error) std::rethrow_exception(error);
    if (json) out << "]\n";
}

} // namespace
