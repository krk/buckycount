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

#include "census_graph_output.hpp"
#include "commands.hpp"
#include <fstream>

namespace {
// The cursor visits floor(h/2) interior labels and one real-ray label at height h.
unsigned long families_through(I h) {
    return static_cast<unsigned long>(Wide(h) * h / 4 + h);
}

unsigned long family_id(Point beta) {
    const I h = beta.x + beta.y;
    return families_through(h - 1) + (beta.y ? beta.x : h / 2 + 1);
}

Point family_label(unsigned long id, I max_height) {
    I lo = 1, hi = max_height;
    while (lo < hi) {
        const I mid = lo + (hi - lo) / 2;
        if (families_through(mid) < id) lo = mid + 1;
        else hi = mid;
    }
    const I a = static_cast<I>(id - families_through(lo - 1));
    return a <= lo / 2 ? Point{a, lo - a} : Point{lo, 0};
}

// Little-endian base-128 integers. GMP handles ranks beyond machine width.
void put_uint(std::string& bytes, const mpz_class& value) {
    const auto offset = bytes.size();
    const auto size = (mpz_sizeinbase(value.get_mpz_t(), 2) + 6) / 7;
    bytes.resize(offset + size);
    mpz_export(bytes.data() + offset, nullptr, -1, 1, 0, 1, value.get_mpz_t());
    for (auto i = offset; i + 1 < bytes.size(); ++i) bytes[i] |= 0x80;
}

mpz_class get_positive(std::string_view& bytes) {
    for (std::size_t i = 0; i < bytes.size(); ++i) {
        const auto byte = static_cast<unsigned char>(bytes[i]);
        if (byte >= 128) continue;
        need(byte != 0, "noncanonical family index");
        mpz_class value;
        mpz_import(value.get_mpz_t(), i + 1, -1, 1, 0, 1, bytes.data());
        bytes.remove_prefix(i + 1);
        return value;
    }
    throw std::runtime_error("truncated family index");
}

std::uint64_t index_checksum(I n, std::string_view bytes) {
    std::uint64_t hash = 14695981039346656037ULL ^ static_cast<std::uint64_t>(n);
    for (unsigned char c : bytes) hash = (hash ^ c) * 1099511628211ULL;
    return hash;
}

void save_index(const std::filesystem::path& path, I n, const FamilyEnds& ends) {
    std::string bytes;
    unsigned long previous = 0;
    mpz_class begin = 0;
    for (const auto& [beta, end] : ends) {
        const auto id = family_id(beta);
        put_uint(bytes, id - previous);
        put_uint(bytes, end - begin);
        previous = id;
        begin = end;
    }
    auto hash = index_checksum(n, bytes);
    for (int i = 0; i < 8; ++i, hash >>= 8) bytes.push_back(static_cast<char>(hash & 255));
    std::filesystem::create_directories(path.parent_path());
    atomic_bytes(path, bytes);
}

FamilyEnds load_index(const std::filesystem::path& path, I n, I max_height, const mpz_class& total) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    const auto size = in.tellg();
    need(size >= 8, "missing or truncated family index");
    std::string storage(static_cast<std::size_t>(size), '\0');
    in.seekg(0);
    need(bool(in.read(storage.data(), size)), "cannot read family index");
    std::string_view bytes(storage.data(), storage.size() - 8);
    auto hash = index_checksum(n, bytes);
    for (std::size_t i = bytes.size(); i < storage.size(); ++i, hash >>= 8)
        if (static_cast<unsigned char>(storage[i]) != (hash & 255))
            throw std::runtime_error("damaged family index: " + path.string());
    FamilyEnds ends;
    unsigned long id = 0;
    mpz_class end = 0;
    while (!bytes.empty()) {
        const auto gap = get_positive(bytes);
        need(gap <= families_through(max_height) - id, "family index exceeds height bound");
        id += gap.get_ui();
        end += get_positive(bytes);
        need(end <= total, "family index exceeds census");
        ends.emplace_back(family_label(id, max_height), end);
    }
    need(end == total, "incomplete family index");
    return ends;
}

mpz_class interval_begin(const FamilyEnds& ends, FamilyEnds::const_iterator it) {
    return it == ends.begin() ? mpz_class(0) : std::prev(it)->second;
}
} // namespace

void graph_request(const Request& r, std::ostream& out) {
    const I n = parse_vertices(r.size), m = n / 2;
    const bool index_only = r.selection == "cache";
    if (!r.families && !index_only) need(n <= 65534, "planar code supports at most 65534 vertices");
    const auto path = r.cache / ("family-ranks-v1-C" + std::to_string(n) + ".bin");
    if (r.use_cache && r.output != "-") need(std::filesystem::weakly_canonical(path) != std::filesystem::weakly_canonical(r.output), "output collides with cache entry");
    const bool hit = r.use_cache && std::filesystem::exists(path);
    if (index_only) std::cerr << "output: " << std::filesystem::absolute(path).lexically_normal().string() << '\n';
    const bool build_table = index_only || r.families || (r.selection == "all" && r.family.empty());
    std::optional<Point> beta;
    if (!r.family.empty()) {
        const auto comma = r.family.find(',');
        need(comma != 0 && comma != std::string::npos && comma + 1 < r.family.size() &&
             r.family.find_first_not_of("0123456789,") == std::string::npos && r.family.find(',', comma + 1) == std::string::npos,
             "--family expects two nonnegative integers A,B");
        beta = parse_family(r.family.substr(0, comma), r.family.substr(comma + 1), n);
    }
    const bool batch = r.selection == "-", all = r.selection == "all";
    mpz_class first = 0, last = 0;
    if (!all && !batch && !r.families && !index_only) {
        const auto dots = r.selection.find("..");
        if (dots != std::string::npos)
            need(dots != 0 && dots + 2 < r.selection.size() &&
                 r.selection.find_first_not_of("0123456789", dots + 2) == std::string::npos,
                 "expected an inclusive rank range A..B");
        first = parse_rank(r.selection.substr(0, dots));
        last = dots == std::string::npos ? first : parse_rank(r.selection.substr(dots + 2));
        need(first <= last, "reversed rank range");
    }
    const auto total = graph_count(r, static_cast<unsigned long>(m));
    const I max_height = static_cast<I>(identifier_sqrt(3 * Wide(n) - 1));
    if (r.estimates) {
        // Each coefficient row contains at most 38 integers bounded by N or the total.
        const mpz_class h = static_cast<unsigned long>(max_height);
        const mpz_class bound = 1024 + h * (h + 1) / 2 * (100 + 38 * std::max(r.size.size(), total.get_str().size()));
        const mpz_class labels = families_through(max_height);
        const auto row_bytes = (mpz_sizeinbase(labels.get_mpz_t(), 2) + 6) / 7 +
                               (mpz_sizeinbase(total.get_mpz_t(), 2) + 6) / 7;
        const mpz_class index_bound = 8 + std::min(labels, total) * row_bytes;
        if (!r.use_cache || !build_table || hit)
            std::cerr << "family cache estimate: 0 additional bytes (" << (hit ? "hit" : "unused") << ")\n";
        else
            std::cerr << "family cache estimate: at most " << disk_size(index_bound, true) << " additional\n";
        if (r.families)
            std::cerr << "output estimate: at most " << disk_size(bound, true) << '\n';
        else if (!index_only) {
            const unsigned long bytes = n <= 255 ? 1 + 4 * n : 3 + 8 * n;
            if (batch)
                std::cerr << "output estimate: " << disk_size(bytes) << " per graph plus stream header; count depends on stdin\n";
            else {
                const mpz_class count = all ? total : mpz_class(last - first + 1);
                std::cerr << "output estimate: " << (all && beta ? "at most " : "") << disk_size(count * bytes + (n <= 255 ? 15 : 18), all && beta)
                          << " (" << (all && beta ? "up to " : "") << count.get_str() << (count == 1 ? " graph)\n" : " graphs)\n");
            }
        }
    }
    FamilyEnds ends;
    if (hit) ends = load_index(path, n, max_height, total);
    if (index_only) {
        if (!hit) {
            fill_family_table(m, total, nullptr, &ends);
            save_index(path, n, ends);
        }
        return;
    }
    if (r.families) {
        CensusFamilyTable report(m, total);
        FamilyEnds computed;
        fill_family_table(m, total, &report, &computed);
        if (hit) need(ends == computed, "family index disagrees with report");
        else if (r.use_cache) save_index(path, n, computed);
        report.write_csv(out);
        return;
    }
    CensusRankSelector selector(m, total);
    if (beta && !hit && !all && !batch && first == last) {
        need(first < total && small_height_detail::admissible(m, *beta), "rank outside selected family or census");
        const auto period = selector.select_family_centered(*beta, first);
        star_dense::Workspace work;
        write_planar_header(n, out);
        write_period(out, n, period, work);
        return;
    }
    mpz_class available = total;
    if (beta) {
        if (hit) {
            const auto it = std::lower_bound(ends.begin(), ends.end(), family_id(*beta),
                [](const auto& entry, unsigned long id) { return family_id(entry.first) < id; });
            available = it != ends.end() && it->first == *beta ? mpz_class(it->second - interval_begin(ends, it)) : mpz_class(0);
        } else if (!all) {
            available = small_height_detail::admissible(m, *beta) ? family_count(m, *beta) : mpz_class(0);
        }
    }
    if (!all && !batch) need(last < available, "rank outside selected family or census");
    if (all) {
        unsigned workers = r.workers;
        if (available < workers) workers = std::max(1UL, available.get_ui());
        write_planar_header(n, out);
        if (available != 0) {
            std::vector<star_dense::Workspace> work(workers);
            auto produce = [&](auto emit) {
                if (beta) {
                    visit_family_periods(m, *beta, emit);
                    return;
                }
                mpz_class count = 0;
                visit_census_families(m, [&](Point family) {
                    const auto begin = count;
                    visit_family_periods(m, family, [&](const IdentifierCusp& period) {
                        emit(period);
                        ++count;
                    });
                    need(count <= total, "family stream exceeds the ordinary isomer count");
                    if (r.use_cache && !hit && count != begin) ends.emplace_back(family, count);
                    return count != total;
                });
                need(count == total, "family stream disagrees with the ordinary isomer count");
            };
            auto batch_size = graph_batch_size(n);
            const mpz_class target = (available + 4UL * workers - 1) / (4UL * workers);
            if (target < batch_size) batch_size = std::max(1UL, target.get_ui());
            census_write_ordered(out, workers, batch_size, produce,
                [&](auto& stream, const auto& period, unsigned lane) { write_period(stream, n, period, work[lane]); });
        }
        if (!beta && r.use_cache && !hit) save_index(path, n, ends);
        return;
    }
    star_dense::Workspace work;
    if (!batch && first != last) {
        mpz_class skip = first, remaining = last - first + 1;
        write_planar_header(n, out);
        auto stream_family = [&](Point family) {
            if (!small_height_detail::admissible(m, family)) return true;
            if (skip > 0 && !hit && !beta) {
                const auto mass = census_orbit_family_mass(m, family);
                if (mass && skip >= *mass) {
                    skip -= *mass;
                    return true;
                }
            }
            const auto caps = census_small_height_inventory(family);
            if (caps.empty()) return true;
            CensusPointOwner owner(m, family, caps);
            struct Finished {};
            try {
                owner.emit<true>([&](const auto& window, const auto& x, const auto& y) {
                    const auto length = static_cast<unsigned long>(window.size());
                    if (skip >= length) {
                        skip -= length;
                        return;
                    }
                    const auto start = skip.get_ui();
                    skip = 0;
                    for (auto i = start; i < length; ++i) {
                        const auto period = census_centered_period(window[i], family, x.nearest,
                            y.nearest, x.to_nearest, y.to_nearest);
                        write_period(out, n, period, work);
                        if (--remaining == 0) throw Finished{};
                    }
                });
            } catch (const Finished&) { return false; }
            return true;
        };
        if (beta) stream_family(*beta);
        else if (hit) {
            auto it = std::upper_bound(ends.begin(), ends.end(), first,
                [](const mpz_class& value, const auto& entry) { return value < entry.second; });
            skip -= interval_begin(ends, it);
            while (it != ends.end() && stream_family(it->first)) ++it;
        } else visit_census_families(m, stream_family);
        need(remaining == 0, "family stream ended before last requested rank");
        return;
    }
    bool started = false;
    auto select = [&](const mpz_class& rank) {
        need(rank < available, "rank outside selected family or census");
        if (!started) {
            write_planar_header(n, out);
            started = true;
        }
        if (beta)
            write_period(out, n, selector.select_family_centered(*beta, rank), work);
        else if (hit) {
            auto it = std::upper_bound(ends.begin(), ends.end(), rank,
                                       [](const mpz_class& value, const auto& entry) { return value < entry.second; });
            write_period(out, n, selector.select_family_centered(it->first, rank - interval_begin(ends, it)), work);
        } else
            write_period(out, n, selector.select_centered(rank), work);
    };
    if (batch) {
        std::string token;
        while (std::cin >> token)
            select(parse_rank(token));
        need(std::cin.eof() && !std::cin.bad(), "rank input failed");
        if (!started) write_planar_header(n, out);
    } else select(first);
}
