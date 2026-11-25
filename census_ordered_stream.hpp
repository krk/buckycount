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
#include "census_heap_scopes.hpp"

#include "genfullerene.hpp"
#include <condition_variable>
#include <mutex>
#include <sstream>
#include <thread>

namespace gf::fullerene_rank {

// The producer calls emit in census order. The decoder may be called concurrently
// on separate streams. The caller alone owns the destination and writes each
// completed batch in submission order.
template<class Produce, class Write>
void census_write_ordered(std::ostream& out, unsigned workers, std::size_t batch_size,
                          Produce&& produce, Write&& write) {
    CENSUS_HEAP_SCOPE(pipeline);
    need(workers > 0 && batch_size > 0, "ordered stream needs positive worker and batch sizes");
    // Each lane belongs to one decoder thread for its entire lifetime.
    auto write_one = [&](std::ostream& destination, const IdentifierCusp& period, unsigned lane) {
        if constexpr (std::is_invocable_v<Write&, std::ostream&, const IdentifierCusp&, unsigned>)
            write(destination, period, lane);
        else
            write(destination, period);
    };
    if (workers == 1) {
        produce([&](const IdentifierCusp& period) {
            write_one(out, period, 0);
            need(bool(out), "ordered graph output failed");
        });
        return;
    }

    // Two batches per worker can be in flight. The caller fills one additional
    // slot, then waits for the oldest submitted batch before reusing its slot.
    struct Batch {
        std::vector<IdentifierCusp> periods;
        std::string bytes;
        std::exception_ptr failure;
        bool complete = false;
    };
    const auto capacity = std::size_t{2} * workers;
    std::vector<Batch> slots(capacity + 1);
    auto advance = [&](std::size_t i) { return i + 1 == slots.size() ? 0 : i + 1; };

    std::mutex mutex;
    std::condition_variable available, completed;
    std::size_t next_job = 0, jobs = 0;
    bool stopped = false;
    auto run = [&](unsigned lane) {
        CENSUS_HEAP_SCOPE(pipeline);
        for (;;) {
            std::size_t i;
            {
                std::unique_lock lock(mutex);
                available.wait(lock, [&] { return stopped || jobs; });
                if (stopped) return;
                i = next_job;
                next_job = advance(next_job);
                --jobs;
            }
            auto& slot = slots[i];
            try {
                // Starting a fresh stream also resets formatting and error state.
                // Moving the cleared string back in retains its allocated buffer.
                slot.bytes.clear();
                std::ostringstream bytes(std::move(slot.bytes));
                for (const auto& period : slot.periods) write_one(bytes, period, lane);
                need(bool(bytes), "ordered graph decoder failed");
                slot.bytes = std::move(bytes).str();
            } catch (...) {
                slot.failure = std::current_exception();
            }
            {
                std::lock_guard lock(mutex);
                slot.complete = true;
            }
            completed.notify_one();
        }
    };
    auto stop = [&] {
        {
            std::lock_guard lock(mutex);
            stopped = true;
        }
        available.notify_all();
    };

    std::size_t input = 0, output = 0, outstanding = 0;
    auto drain_one = [&] {
        auto& slot = slots[output];
        {
            std::unique_lock lock(mutex);
            completed.wait(lock, [&] { return slot.complete; });
            slot.complete = false;
        }
        if (slot.failure) std::rethrow_exception(slot.failure);
        out.write(slot.bytes.data(), static_cast<std::streamsize>(slot.bytes.size()));
        need(bool(out), "ordered graph output failed");
        slot.periods.clear();
        output = advance(output);
        --outstanding;
    };
    auto submit = [&] {
        if (outstanding == capacity) drain_one();
        ++outstanding;
        input = advance(input);
        {
            std::lock_guard lock(mutex);
            ++jobs;
        }
        available.notify_one();
    };

    std::vector<std::jthread> threads;
    try {
        for (unsigned i = 0; i < workers; ++i) threads.emplace_back(run, i);
        produce([&](const IdentifierCusp& period) {
            auto& batch = slots[input].periods;
            if (batch.empty()) batch.reserve(batch_size);
            batch.push_back(period);
            if (batch.size() == batch_size) submit();
        });
        if (!slots[input].periods.empty()) submit();
        while (outstanding) drain_one();
    } catch (...) {
        stop();
        throw;
    }
    stop();
}

} // namespace gf::fullerene_rank
