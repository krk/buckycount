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

#include "prefix_engine_buckycount.hpp"
#include "calculator.hpp"
#include "commands.hpp"
#include <iomanip>
#include <set>

namespace fs = std::filesystem;

namespace {
const char* help = R"(Usage:
  buckycount [count] N|A..B [options]
  buckycount unrank N R|A..B|all|- -o FILE.pc [options]
  buckycount unrank N cache [options]
  buckycount decode FILE.pc [FILE.pc ...] [options]

N is an even vertex count. Ranks start at 0. Ranges include both endpoints.

Count options:
  --ipr              Count only isolated-pentagon cages
  --details          Show the count breakdown
  --format csv|json  Write counts as CSV or JSON
  --families         List families and their rank intervals for one vertex count

Unrank options:
  --family A,B       Use ranks within family A,B when emitting graphs

Count, unrank, and decode options:
  -o, --output PATH  Count or decode to PATH (directory for multiple inputs). Required for graph output.
  --workers W        Limit workers for count, unrank all, or multi-file decode
  --no-estimates     Skip size estimates
  -h, --help         Show help for the program or the selected command

Count and unrank options:
  --cache-dir DIR    Set the count or unrank cache directory (default: .cache)
  --no-cache         Disable caching for count or unrank (except unrank N cache)

Counting formulas and rank order: https://doi.org/10.5281/zenodo.24289727
)";

void command_help(const std::string& command) {
    if (command == "count")
        std::cout << R"(Usage: buckycount [count] N|A..B [options]

Options:
  --ipr              Count isolated-pentagon cages
  --details          Show the count breakdown
  --format csv|json  Select the count output format
  --families         List families and rank intervals for one vertex count
  -o, --output FILE  Write to FILE (default: stdout)
  --workers W        Set the worker limit
)";
    else if (command == "unrank")
        std::cout << R"(Usage: buckycount unrank N R|A..B|all|- -o FILE.pc [options]
       buckycount unrank N cache [--cache-dir DIR] [--no-estimates]

Ranks start at 0. Use - to read ranks from stdin, or cache to build the rank index.

Options:
  -o, --output FILE  Write binary planar code (required when emitting graphs)
  --family A,B       Select ranks within one family
  --workers W        Set the worker limit for unrank all
)";
    else if (command == "decode")
        std::cout << R"(Usage: buckycount decode FILE.pc [FILE.pc ...] [options]

Use - as the sole input to read planar code from stdin.

Options:
  -o, --output PATH  Write adjacency text to a file (directory for multiple inputs)
  --workers W        Set the worker limit for multiple input files
)";
    if (command == "count" || command == "unrank") {
        std::cout << "  --cache-dir DIR    Choose the cache directory (default: .cache)\n";
        std::cout << (command == "unrank"
                          ? "  --no-cache         Disable disk caching when emitting graphs\n"
                          : "  --no-cache         Disable disk caching\n");
    }
    std::cout << "  --no-estimates     Skip size estimates\n";
    std::cout << "  -h, --help         Show this help\n";
}

Request parse_request(std::string command, int start, int argc, char** argv) {
    Request r;
    r.command = std::move(command);
    std::vector<std::string> positional;
    std::set<std::string> seen;
    bool options = true;
    const bool count = r.command == "count";
    const bool caches = count || r.command == "unrank";
    for (int i = start; i < argc; ++i) {
        std::string arg = argv[i];
        if (options && arg == "--") {
            options = false;
            continue;
        }
        if (options && r.command != "decode" && arg.size() > 1 && arg[0] == '-' && arg[1] >= '0' && arg[1] <= '9')
            fail(positional.empty() ? "vertex count must be positive" : "rank must be nonnegative");
        if (!options || arg == "-" || !arg.starts_with("-")) {
            positional.push_back(arg);
            continue;
        }
        const auto equal = arg.starts_with("--") ? arg.find('=') : std::string::npos;
        std::string value;
        if (equal != std::string::npos) {
            value = arg.substr(equal + 1);
            arg.resize(equal);
        }
        if (arg == "-o") arg = "--output";
        if (arg == "-h") arg = "--help";
        const bool flag = arg == "--help" || (caches && arg == "--no-cache") ||
                          arg == "--no-estimates" ||
                          (count && (arg == "--ipr" || arg == "--details" || arg == "--families"));
        const bool takes_value = arg == "--output" || (caches && arg == "--cache-dir") ||
                                 arg == "--workers" ||
                                 (count && arg == "--format") || (r.command == "unrank" && arg == "--family");
        if (!flag && !takes_value) fail("unknown option for " + r.command + ": " + arg + " (see --help)");
        if (!seen.insert(arg).second) fail("duplicate option " + arg);
        if (flag && equal != std::string::npos) fail(arg + " does not take a value");
        if (arg == "--help") r.help = true;
        else if (arg == "--no-cache") r.use_cache = false;
        else if (arg == "--no-estimates") r.estimates = false;
        else if (arg == "--ipr") r.ipr = true;
        else if (arg == "--details") r.details = true;
        else if (arg == "--families") r.families = true;
        else {
            if (equal == std::string::npos) {
                if (i + 1 == argc || (argv[i + 1][0] == '-' && argv[i + 1][1] != '\0'))
                    fail(arg + " requires a value");
                value = argv[++i];
            }
            if (value.empty()) fail(arg + " requires a nonempty value");
            if (arg == "--output") r.output = value;
            else if (arg == "--cache-dir")
                r.cache = value;
            else if (arg == "--format") {
                if (value != "csv" && value != "json") fail("--format expects csv or json");
                r.format = value;
            } else if (arg == "--family")
                r.family = value;
            else {
                const auto w = parse_number(value, "workers", true);
                if (w > std::numeric_limits<unsigned>::max()) fail("worker count too large");
                r.workers = static_cast<unsigned>(w);
            }
        }
    }
    if (r.help) return r;
    const unsigned online = std::max(1u, std::thread::hardware_concurrency());
    r.workers = r.workers ? std::min(r.workers, online) : online;
    if (r.command == "decode") {
        r.inputs = std::move(positional);
        if (r.inputs.empty()) fail("decode requires at least one input file (use - for stdin)");
        if (r.inputs.size() > 1 && seen.contains("--output") && !fs::is_directory(r.output))
            fail("multiple inputs require an output directory");
    } else {
        if (positional.empty()) fail(r.command + " requires a vertex count (for example: " + r.command + (count ? " 60)" : " 60 9 -o cage.pc)"));
        if (!count && positional.size() == 1) fail("unrank requires a rank, an inclusive range, all, cache, or -");
        if (positional.size() > (count ? 1u : 2u)) fail("unexpected argument for " + r.command + ": " + positional[count ? 1 : 2]);
        r.size = positional[0];
        if (r.families && (r.ipr || r.details || r.format == "json" || r.size.find("..") != std::string::npos))
            fail("--families requires one size and cannot be combined with --ipr, --details, or --format json");
        if (r.command == "unrank") {
            r.selection = positional[1];
            if (r.selection == "cache") {
                if (!r.use_cache || seen.contains("--output") || seen.contains("--family"))
                    fail("unrank cache builds the full family index; omit --no-cache, --output and --family");
            } else if (r.output == "-")
                fail("unrank requires an output filename for binary planar code (use -o FILE.pc)");
            if (seen.contains("--workers") && r.selection != "all") fail("--workers requires unrank all");
        }
    }
    return r;
}

Options count_options(const Request& r) {
    Options opts;
    opts.workers = r.workers;
    opts.ipr = r.ipr;
    opts.output = r.output;
    opts.format = r.format.empty() ? "csv" : r.format;
    set_count_interval(opts, r.size);
    return opts;
}

void count_request(std::ostream& out, const Request& r, const Options& opts) {
    const bool range = r.size.find("..") != std::string::npos;
    const Integer rows = opts.last - opts.first + 1;
    if (r.estimates) {
        print_count_cache_estimate(opts, r.use_cache);
        std::cerr << "output estimate: " << rows.get_str() << (rows == 1 ? " row" : " rows") << "\n";
    }

    write_prefix(out, opts, !r.details && r.format.empty(), range);
}
} // namespace

mpz_class graph_count(const Request& r, const mpz_class& index) {
    if (r.use_cache) prefix_cache().open(r.cache / "coeff", r.output);
    if (r.estimates) {
        Options opts;
        opts.first = opts.last = index;
        print_count_cache_estimate(opts, r.use_cache);
    }
    const auto total = evaluate_ordinary_index(index).ordinary;
    if (r.selection == "cache") prefix_cache().finish(false);
    return total;
}

void write_atomic(const fs::path& path, const std::function<void(std::ostream&)>& write) {
    OutputFile file(path.string());
    file.publish(path.string(), write);
}

void atomic_bytes(const fs::path& path, const std::string& bytes) {
    write_atomic(path, [&](std::ostream& out) { out.write(bytes.data(), bytes.size()); });
}

int buckycount_main(int argc, char** argv) {
    try {
        const std::string command = argc > 1 ? argv[1] : "";
        if (argc == 1 || (argc == 2 && (command == "--help" || command == "-h"))) {
            std::cout << help;
            return 0;
        }
        const bool named = command == "count" || command == "unrank" || command == "decode";
        if (!named && !command.starts_with("-") && command.find_first_of("0123456789.") != 0)
            fail("unknown command '" + command + "' (use count, unrank, or decode; see --help)");
        const auto r = parse_request(named ? command : "count", named ? 2 : 1, argc, argv);
        if (r.help) {
            command_help(r.command);
            return 0;
        }
        std::optional<Options> count;
        if (r.command == "count" && !r.families) count = count_options(r);
        std::cerr << "request:";
        for (int i = 1; i < argc; ++i)
            std::cerr << ' ' << std::quoted(argv[i]);
        std::cerr << '\n';
        std::cerr << "cache: " << fs::absolute(r.cache).lexically_normal().string()
                  << (r.command == "decode" ? " (unused)" : r.use_cache ? "" : " (disabled)") << '\n';
        OutputCleanup cleanup;
        if (r.command == "decode") {
            decode_request(r);
            return 0;
        }
        const auto destination = count ? count->output : r.output;
        if (r.selection != "cache") std::cerr << "output: " << (destination == "-" ? "stdout" : fs::absolute(destination).lexically_normal().string()) << '\n';
        if (count && r.use_cache) prefix_cache().open(r.cache / "coeff", destination);
        write_atomic(destination, [&](std::ostream& out) {
            if (count) count_request(out, r, *count);
            else graph_request(r, out);
            if (!out.flush()) fail("output write failed");
            if (r.selection != "cache") prefix_cache().finish(destination == "-");
        });
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "buckycount: " << e.what() << '\n';
        return 1;
    }
}
