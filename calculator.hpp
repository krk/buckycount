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

// Atomic file and standard output for all commands.
#include <array>
#include <atomic>
#include <pthread.h>
#include <filesystem>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <unistd.h>

#include <cerrno>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <random>
#include <streambuf>
#include <sys/stat.h>
#include <utility>

namespace {

// Only the signal-wait thread handles termination. Workers register and remove
// files under the same mutex so cleanup sees every live temporary exactly once.
std::mutex temporary_mutex;
std::map<int, std::string> temporary_files;
mode_t creation_mode;

[[noreturn]] void output_error(const std::string& message) { throw std::runtime_error(message); }

class OutputCleanup {
    sigset_t previous_{};
    std::atomic_bool finished_{false};
    std::thread waiter_;

public:
    OutputCleanup() {
        const auto mask = ::umask(0);
        ::umask(mask);
        creation_mode = 0666 & ~mask;
        sigset_t signals;
        sigemptyset(&signals);
        sigaddset(&signals, SIGINT);
        sigaddset(&signals, SIGTERM);
        if (::pthread_sigmask(SIG_BLOCK, &signals, &previous_))
            output_error("cannot block output cleanup signals");
        try {
            waiter_ = std::thread([this, signals] {
                int signal;
                if (::sigwait(&signals, &signal)) ::_exit(1);
                if (finished_) return;
                const std::lock_guard guard(temporary_mutex);
                for (const auto& [directory, name] : temporary_files)
                    ::unlinkat(directory, name.c_str(), 0);
                ::_exit(128 + signal);
            });
        } catch (...) {
            ::pthread_sigmask(SIG_SETMASK, &previous_, nullptr);
            throw;
        }
    }

    ~OutputCleanup() {
        finished_ = true;
        ::pthread_kill(waiter_.native_handle(), SIGTERM);
        waiter_.join();
        ::pthread_sigmask(SIG_SETMASK, &previous_, nullptr);
    }
};

void synchronize(int descriptor, const std::string& failure) {
    while (::fsync(descriptor) != 0) {
        if (errno != EINTR) output_error(failure + ": " + std::strerror(errno));
    }
}

// Buffer writes directly to the exclusively created temporary file.
class OutputBuffer final : public std::streambuf {
    int descriptor_;
    std::array<char, 16384> buffer_{};

    int sync() override {
        const char* next = pbase();
        while (next != pptr()) {
            const auto count = ::write(descriptor_, next, static_cast<std::size_t>(pptr() - next));
            if (count < 0 && errno == EINTR) continue;
            if (count <= 0) return -1;
            next += count;
        }
        setp(buffer_.data(), buffer_.data() + buffer_.size());
        return 0;
    }

    int_type overflow(int_type value) override {
        if (sync() != 0) return traits_type::eof();
        if (!traits_type::eq_int_type(value, traits_type::eof())) {
            *pptr() = traits_type::to_char_type(value);
            pbump(1);
        }
        return traits_type::not_eof(value);
    }

public:
    explicit OutputBuffer(int descriptor) : descriptor_(descriptor) {
        setp(buffer_.data(), buffer_.data() + buffer_.size());
    }
};

// Validate the destination before calculation and publish only a complete result.
class OutputFile {
    std::filesystem::path parent_, name_;
    std::string temporary_;
    int descriptor_ = -1;
    int directory_ = -1;


    // An explicit symlink to a regular file (or a dangling symlink) is
    // replaced as a directory entry. Only actual regular files supply mode bits.
    mode_t destination_mode() const {
        struct stat existing{};
        if (::fstatat(directory_, name_.c_str(), &existing, AT_SYMLINK_NOFOLLOW) == 0) {
            if (S_ISREG(existing.st_mode)) return existing.st_mode & 0777;
            if (!S_ISLNK(existing.st_mode))
                output_error("output is not a regular file: " + (parent_ / name_).string());
            if (::fstatat(directory_, name_.c_str(), &existing, 0) == 0) {
                if (S_ISREG(existing.st_mode)) return creation_mode;
                output_error("output refers to a non-regular file: " + (parent_ / name_).string());
            }
        }
        if (errno != ENOENT)
            output_error("cannot inspect output " + (parent_ / name_).string() + ": " + std::strerror(errno));
        return creation_mode;
    }

public:
    explicit OutputFile(const std::string& output) {
        if (output == "-") return;
        const std::lock_guard guard(temporary_mutex);
        const auto path = std::filesystem::absolute(output);
        parent_ = path.parent_path();
        name_ = path.filename();
        if (name_.empty()) output_error("output requires a filename: " + output);
        directory_ = ::open(parent_.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
        if (directory_ < 0) output_error("cannot open output directory: " + std::string(std::strerror(errno)));
        try {
            destination_mode();

            std::random_device random;
            do {
                temporary_ = ".buckycount-" + std::to_string(random());
                if (temporary_ == name_) continue;
                descriptor_ = ::openat(directory_, temporary_.c_str(),
                                       O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0600);
                if (descriptor_ < 0 && errno != EEXIST)
                    output_error("cannot create output in " + output + ": " + std::strerror(errno));
            } while (descriptor_ < 0);
            temporary_files.emplace(directory_, temporary_);
        } catch (...) {
            if (descriptor_ >= 0) {
                ::close(descriptor_);
                ::unlinkat(directory_, temporary_.c_str(), 0);
            }
            ::close(directory_);
            throw;
        }
    }

    OutputFile(const OutputFile&) = delete;
    OutputFile& operator=(const OutputFile&) = delete;

    ~OutputFile() {
        const std::lock_guard guard(temporary_mutex);
        if (descriptor_ >= 0) ::close(descriptor_);
        if (!temporary_.empty()) ::unlinkat(directory_, temporary_.c_str(), 0);
        temporary_files.erase(directory_);
        if (directory_ >= 0) ::close(directory_);
    }

    void publish(const std::string& output, const std::function<void(std::ostream&)>& write) {
        OutputBuffer buffer(descriptor_);
        std::ostream file(&buffer);
        std::ostream& stream = directory_ < 0 ? std::cout : file;
        write(stream);
        stream.flush();
        if (!stream) output_error("cannot write output " + output);
        if (directory_ < 0) return;
        const mode_t mode = destination_mode();
        if (::fchmod(descriptor_, mode) != 0)
            output_error("cannot set output permissions: " + std::string(std::strerror(errno)));
        synchronize(descriptor_, "cannot synchronize output " + output);
        const int descriptor = std::exchange(descriptor_, -1);
        if (::close(descriptor) != 0)
            output_error("cannot close output " + output + ": " + std::strerror(errno));
        {
            const std::lock_guard guard(temporary_mutex);
            if (::renameat(directory_, temporary_.c_str(), directory_, name_.c_str()) != 0)
                output_error("cannot publish " + output + ": " + std::strerror(errno));
            temporary_files.erase(directory_);
            temporary_.clear();
        }
        synchronize(directory_, "output published but cannot synchronize its directory");

    }
};

} // namespace
