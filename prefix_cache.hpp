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

// One immutable pack per million discriminants. Each footer index entry
// addresses the first of up to 128 variable-length records.
#pragma once

class PrefixCache {
    using Values = std::array<std::optional<Integer>, 4>;
    using Path = std::filesystem::path;
    static constexpr unsigned long span = 1000000;
    static constexpr unsigned stride = 128;
    static constexpr unsigned entry_size = 12; // uint32 key and uint64 byte offset
    Path root_;
    std::mutex mutex_;
    std::map<Integer, Values> pending_;

    static std::uint64_t hash(std::string_view bytes, std::uint64_t h = 14695981039346656037ULL) {
        for (unsigned char byte : bytes)
            h = (h ^ byte) * 1099511628211ULL;
        return h;
    }
    Path path(const Integer &base) const {
        auto dir = root_ / "primitive";
        auto digits = base.get_str();
        while (digits.size() > 100) {
            dir /= digits.substr(0, 90);
            digits.erase(0, 90);
        }
        // The last six digits are zero, including when earlier digits form directories.
        const Integer end = Integer(digits, 10) + span - 1;
        return dir / (digits + "-" + end.get_str() + ".bin");
    }
    static void append_word(std::string &data, std::uint64_t value, unsigned width) {
        for (unsigned i = 0; i < width; ++i)
            data.push_back(static_cast<char>(value >> (8 * i)));
    }
    static std::uint64_t word(std::string_view data, unsigned width) {
        if (data.size() < width) fail("truncated coefficient pack");
        std::uint64_t value = 0;
        for (unsigned i = 0; i < width; ++i)
            value |= std::uint64_t(static_cast<unsigned char>(data[i])) << (8 * i);
        return value;
    }
    static std::string checksum(const Integer &discriminant, std::string_view data) {
        auto h = hash(data, hash(discriminant.get_str(), hash("primitive-gmp")));
        std::string bytes;
        for (unsigned i = 0; i < 8; ++i)
            bytes.push_back(static_cast<char>(h >> (8 * i)));
        return bytes;
    }
    // The first byte carries the sign and magnitude length. 127 introduces
    // an eight-byte little-endian length. Magnitudes are least-significant byte first.
    static void append_integer(std::string &data, const Integer &value) {
        const std::size_t size = value == 0 ? 0 : (mpz_sizeinbase(value.get_mpz_t(), 2) + 7) / 8;
        data.push_back(static_cast<char>((value < 0 ? 128 : 0) | std::min<std::size_t>(size, 127)));
        if (size >= 127) {
            const std::uint64_t length = size;
            for (unsigned i = 0; i < 8; ++i)
                data.push_back(static_cast<char>(length >> (8 * i)));
        }
        if (size) {
            const auto start = data.size();
            data.resize(start + size);
            mpz_export(data.data() + start, nullptr, -1, 1, 0, 0, value.get_mpz_t());
        }
    }
    static Integer read_integer(std::string_view &data) {
        if (data.empty())
            fail("missing coefficient cache integer");
        const auto header = static_cast<unsigned char>(data.front());
        data.remove_prefix(1);
        std::uint64_t size = header & 127;
        if (size == 127) {
            if (data.size() < 8)
                fail("truncated coefficient cache length");
            size = 0;
            for (unsigned i = 0; i < 8; ++i)
                size |= std::uint64_t(static_cast<unsigned char>(data[i])) << (8 * i);
            data.remove_prefix(8);
            if (size < 127)
                fail("invalid coefficient cache length");
        }
        if (size > data.size() || (size && data[size - 1] == 0) || (!size && (header & 128)))
            fail("invalid coefficient cache integer");
        Integer value;
        mpz_import(value.get_mpz_t(), size, -1, 1, 0, 0, data.data());
        if (header & 128)
            mpz_neg(value.get_mpz_t(), value.get_mpz_t());
        data.remove_prefix(size);
        return value;
    }
    static std::string encode(const Integer &discriminant, const Values &values) {
        std::string data(1, 0);
        for (unsigned i = 0; i < values.size(); ++i) {
            if (!values[i])
                continue;
            data[0] |= static_cast<char>(1U << i);
            append_integer(data, *values[i]);
        }
        return data + checksum(discriminant, data);
    }
    static Values decode(const Integer &discriminant, std::string_view bytes) {
        if (bytes.size() < 10)
            fail("truncated coefficient cache record");
        const std::string_view data(bytes.data(), bytes.size() - 8);
        if (checksum(discriminant, data) != std::string_view(bytes).substr(data.size()))
            fail("coefficient cache checksum mismatch at D" + discriminant.get_str());
        const auto mask = static_cast<unsigned char>(data[0]);
        if (!mask || mask > 15)
            fail("invalid coefficient cache mask");
        Values values;
        auto payload = data.substr(1);
        for (unsigned i = 0; i < values.size(); ++i) {
            if (mask & (1U << i))
                values[i] = read_integer(payload);
        }
        if (!payload.empty())
            fail("unexpected trailing coefficient data");
        return values;
    }
    static bool merge(Values &values, const Values &updates) {
        bool changed = false;
        for (unsigned i = 0; i < updates.size(); ++i) {
            if (!updates[i])
                continue;
            if (values[i] && values[i] != updates[i])
                fail("inconsistent cached coefficient");
            if (!values[i]) {
                values[i] = updates[i];
                changed = true;
            }
        }
        return changed;
    }
    // A record holds its local key, encoded payload length and checksummed values.
    static void record(std::string &out, unsigned key, std::string_view payload) {
        append_word(out, key, 4);
        if (payload.size() < 255) out.push_back(static_cast<char>(payload.size()));
        else {
            out.push_back(static_cast<char>(255));
            append_word(out, payload.size(), 8);
        }
        out.append(payload);
    }
    struct Record { unsigned key; std::string_view values; };
    static Record next(std::string_view &data) {
        const auto key = word(data, 4);
        data.remove_prefix(4);
        auto size = word(data, 1);
        data.remove_prefix(1);
        if (size == 255) {
            size = word(data, 8);
            data.remove_prefix(8);
        }
        if (key >= span || size < 10 || size > data.size()) fail("invalid coefficient pack record");
        const auto values = data.substr(0, size);
        data.remove_prefix(size);
        return {static_cast<unsigned>(key), values};
    }

    struct Descriptor {
        int fd;
        explicit Descriptor(int value) : fd(value) {}
        Descriptor(const Descriptor &) = delete;
        ~Descriptor() { if (fd >= 0) ::close(fd); }
    };
    static void io_error(const char *operation) {
        throw std::system_error(errno, std::generic_category(), operation);
    }
    struct Mapping {
        void *address = MAP_FAILED;
        std::size_t size = 0;
        ~Mapping() { if (address != MAP_FAILED) ::munmap(address, size); }
    };
    struct Pack {
        Mapping mapping;
        ino_t inode = 0;
        std::string_view records, index;
        std::unique_ptr<std::once_flag[]> checked;

        Pack(const Path &file, const Integer &base) {
            Descriptor input(::open(file.c_str(), O_RDONLY | O_CLOEXEC));
            if (input.fd < 0) {
                if (errno == ENOENT) return;
                io_error("opening coefficient pack");
            }
            struct stat info{};
            if (::fstat(input.fd, &info)) io_error("reading coefficient pack metadata");
            if (info.st_size < 32) fail("truncated coefficient pack");
            inode = info.st_ino;
            mapping.size = info.st_size;
            mapping.address = ::mmap(nullptr, mapping.size, PROT_READ, MAP_PRIVATE, input.fd, 0);
            if (mapping.address == MAP_FAILED) io_error("mapping coefficient pack");
            std::string_view data(static_cast<const char *>(mapping.address), mapping.size);
            // Footer: record count, index checksum, magic.
            if (!data.ends_with("BUCKPK01")) fail("invalid coefficient pack footer");
            const auto count = word(data.substr(data.size()-24), 8);
            if (!count || count > span) fail("invalid coefficient pack count");
            const auto bytes = entry_size * ((count + stride - 1) / stride);
            if (bytes >= data.size()-24) fail("truncated coefficient pack index");
            records = data.substr(0, data.size()-24-bytes);
            index = data.substr(records.size(), bytes);
            if (word(data.substr(data.size()-16), 8) !=
                hash(data.substr(records.size(), bytes+8), hash(base.get_str())))
                fail("coefficient pack index checksum mismatch");
            checked = std::make_unique<std::once_flag[]>(bytes / entry_size);
        }
        std::string_view block(std::size_t i) const {
            const auto start = word(index.substr(i*entry_size+4), 8);
            const auto end = (i+1)*entry_size == index.size() ? records.size()
                              : word(index.substr((i+1)*entry_size+4), 8);
            if (end > records.size() || start >= end || end-start < 8)
                fail("invalid coefficient pack block offset");
            const auto bytes = records.substr(start, end-start-8);
            std::call_once(checked[i], [&] {
                if (hash(bytes) != word(records.substr(end-8), 8))
                    fail("coefficient pack block checksum mismatch");
            });
            return bytes;
        }
        Values lookup(const Integer &d, unsigned key) const {
            std::size_t lo = 0, hi = index.size()/entry_size;
            while (lo < hi) {
                const auto mid = lo + (hi-lo)/2;
                if (word(index.substr(mid*entry_size), 4) <= key) lo = mid+1;
                else hi = mid;
            }
            if (!lo) return {};
            auto block = this->block(lo-1);
            while (!block.empty()) {
                const auto row = next(block);
                if (row.key == key) return decode(d, row.values);
                if (row.key > key) break;
            }
            return {};
        }
    };
    std::map<Integer, std::shared_ptr<Pack>> packs_;

    // Spill private batches during calculation. A single unlinked file holds all
    // batches; only chunk positions per pack remain resident until publication.
    struct Chunk { off_t offset; std::size_t bytes; };
    Descriptor scratch_{-1};
    std::map<Integer, std::vector<Chunk>> chunks_;
    static void write_all(int fd, std::string_view bytes) {
        while (!bytes.empty()) {
            const auto n = ::write(fd, bytes.data(), bytes.size());
            if (n < 0 && errno == EINTR) continue;
            if (n <= 0) io_error("writing coefficient cache");
            bytes.remove_prefix(n);
        }
    }
    void spill() {
        if (pending_.empty()) return;
        if (scratch_.fd < 0) {
            auto name = (root_ / ".pending-XXXXXX").string();
            scratch_.fd = ::mkstemp(name.data());
            if (scratch_.fd < 0) io_error("creating private coefficient buffer");
            if (::unlink(name.c_str())) io_error("unlinking private coefficient buffer");
        }
        auto it = pending_.begin();
        while (it != pending_.end()) {
            const Integer base = it->first - mpz_fdiv_ui(it->first.get_mpz_t(), span);
            const Integer end = base + span;
            std::string bytes;
            do {
                record(bytes, mpz_fdiv_ui(it->first.get_mpz_t(), span), encode(it->first, it->second));
                ++it;
            } while (it != pending_.end() && it->first < end);
            const auto pos = ::lseek(scratch_.fd, 0, SEEK_CUR);
            if (pos < 0) io_error("locating private coefficient buffer");
            write_all(scratch_.fd, bytes);
            chunks_[base].push_back({pos, bytes.size()});
        }
        pending_.clear();
    }
    std::map<unsigned, Values> additions(const Integer &base, const std::vector<Chunk> &chunks) {
        std::map<unsigned, Values> updates;
        for (const auto &chunk : chunks) {
            std::string bytes(chunk.bytes, '\0');
            std::size_t used = 0;
            while (used < bytes.size()) {
                const auto n = ::pread(scratch_.fd, bytes.data()+used, bytes.size()-used, chunk.offset+used);
                if (n < 0 && errno == EINTR) continue;
                if (n <= 0) fail("cannot read private coefficient buffer");
                used += n;
            }
            std::string_view data(bytes);
            while (!data.empty()) {
                const auto row = next(data);
                merge(updates[row.key], decode(base+row.key, row.values));
            }
        }
        return updates;
    }
    static std::string replacement(const Pack &old, const Integer &base, const std::map<unsigned, Values> &updates) {
        std::string_view input;
        std::size_t old_block = 0;
        auto update = updates.begin();
        std::string bytes, index;
        std::size_t count = 0, block_start = 0;
        bool changed = false;
        auto emit = [&](unsigned key, std::string_view values) {
            if (count++ % stride == 0) {
                if (count > 1) append_word(bytes, hash(std::string_view(bytes).substr(block_start)), 8);
                block_start = bytes.size();
                append_word(index, key, 4);
                append_word(index, bytes.size(), 8);
            }
            record(bytes, key, values);
        };
        while (!input.empty() || old_block < old.index.size()/entry_size || update != updates.end()) {
            if (input.empty() && old_block < old.index.size()/entry_size)
                input = old.block(old_block++);
            const auto key = input.empty() ? span : word(input, 4);
            if (update != updates.end() && update->first < key) {
                emit(update->first, encode(base+update->first, update->second));
                changed = true;
                ++update;
                continue;
            }
            const auto row = next(input);
            if (update != updates.end() && update->first == row.key) {
                auto values = decode(base+row.key, row.values);
                changed |= merge(values, update->second);
                emit(row.key, encode(base+row.key, values));
                ++update;
            } else emit(row.key, row.values);
        }
        if (!changed) return {};
        append_word(bytes, hash(std::string_view(bytes).substr(block_start)), 8);
        append_word(index, count, 8);
        const auto digest = hash(index, hash(base.get_str()));
        bytes.append(index);
        append_word(bytes, digest, 8);
        bytes += "BUCKPK01";
        return bytes;
    }
    void publish(const Integer &base, const std::map<unsigned, Values> &updates) {
        const auto file = path(base);
        std::filesystem::create_directories(file.parent_path());
        for (;;) {
            Pack old(file, base);
            const auto bytes = replacement(old, base, updates);
            if (bytes.empty()) return;
            auto name = file.string() + ".XXXXXX";
            Descriptor out(::mkstemp(name.data()));
            if (out.fd < 0) io_error("creating replacement coefficient pack");
            struct Cleanup {
                const std::string &name;
                ~Cleanup() { ::unlink(name.c_str()); }
            } cleanup{name};
            write_all(out.fd, bytes);
            if (::close(std::exchange(out.fd, -1))) io_error("closing coefficient pack");
            Descriptor lock(::open(root_.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC));
            if (lock.fd < 0) io_error("opening coefficient cache directory");
            int status;
            do { status = ::flock(lock.fd, LOCK_EX | LOCK_NB); } while (status && errno == EINTR);
            if (status && errno != EWOULDBLOCK && errno != EAGAIN) io_error("locking coefficient cache");
            if (status) {
                Pack latest(file, base);
                for (const auto &[key, wanted] : updates) {
                    auto values = latest.lookup(base+key, key);
                    if (merge(values, wanted)) fail("cache writer lock busy: " + root_.string());
                }
                return;
            }
            struct stat current{};
            ino_t current_inode = 0;
            if (::stat(file.c_str(), &current) == 0) current_inode = current.st_ino;
            else if (errno != ENOENT) io_error("checking published pack");
            // A competing writer published while we prepared the replacement.
            // Release the lock and rebuild against its immutable snapshot.
            if (current_inode != old.inode) continue;
            if (::rename(name.c_str(), file.c_str())) io_error("publishing coefficient pack");
            return;
        }
    }
  public:
    void open(const std::filesystem::path &root, const std::string &output) {
        if (output != "-") {
            const auto destination = std::filesystem::weakly_canonical(std::filesystem::absolute(output));
            const auto cache = std::filesystem::weakly_canonical(std::filesystem::absolute(root));
            const auto relative = destination.lexically_relative(cache);
            if (!relative.empty() && *relative.begin() != "..") fail("output collides with prefix cache");
        }
        std::filesystem::create_directories(root);
        root_ = root;
    }
    bool get(unsigned weight, const Integer &d, Integer &result) {
        if (root_.empty()) return false;
        const unsigned key = mpz_fdiv_ui(d.get_mpz_t(), span);
        const Integer base = d-key;
        std::shared_ptr<Pack> pack;
        {
            const std::lock_guard guard(mutex_);
            const auto it = pending_.find(d);
            if (it != pending_.end() && it->second[weight-1]) {
                result = *it->second[weight-1];
                return true;
            }
            if (!packs_.contains(base)) {
                if (packs_.size() >= 64) packs_.clear();
                packs_[base] = std::make_shared<Pack>(path(base), base);
            }
            pack = packs_.at(base);
        }
        const auto values = pack->lookup(d, key);
        if (!values[weight-1]) return false;
        result = *values[weight-1];
        return true;
    }
    void put(unsigned weight, const Integer &d, const Integer &value) {
        if (root_.empty()) return;
        const std::lock_guard guard(mutex_);
        pending_[d][weight-1] = value;
        if (pending_.size() >= 4096) spill();
    }
    void finish(bool emitted = true) {
        if (pending_.empty() && chunks_.empty()) return;
        std::cerr << (emitted ? "result emitted" : "result buffered") << ", sorting coefficient caches\n";
        spill();
        for (const auto &[base, chunks] : chunks_)
            publish(base, additions(base, chunks));
    }
};

inline PrefixCache &prefix_cache() {
    static PrefixCache cache;
    return cache;
}
