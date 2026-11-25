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

// Included inside the arithmetic engine namespace.
#pragma once

// FLINT's polynomial objects and per-worker scratch caches have separate
// lifetimes. A worker can reuse its caches across several products.
class FlintCall {
#if FLINT_USES_TLS
public:
    FlintCall() {
        struct Cleanup { ~Cleanup() { flint_cleanup(); } };
        thread_local Cleanup cleanup;
        static_cast<void>(cleanup);
    }
#else
    // Libraries built without thread-local caches can share mutable scratch.
    inline static std::mutex mutex_;
    std::lock_guard<std::mutex> lock_{mutex_};
#endif
};

class ExactFactorization {
    FlintCall call_;
    fmpz_factor_t factors_;
public:
    explicit ExactFactorization(const Integer& n) {
        fmpz_t input;
        fmpz_init(input);
        fmpz_set_mpz(input,n.get_mpz_t());
        fmpz_factor_init(factors_);
        fmpz_factor(factors_,input);
        fmpz_clear(input);
    }
    ExactFactorization(const ExactFactorization&) = delete;
    ~ExactFactorization() { fmpz_factor_clear(factors_); }
    std::vector<std::pair<Integer,unsigned long>> primes() const {
        std::vector<std::pair<Integer,unsigned long>> result;
        for (slong i=0; i<factors_->num; ++i) {
            if (!fmpz_is_prime(factors_->p+i)) fail("FLINT returned a composite factor");
            Integer prime;
            fmpz_get_mpz(prime.get_mpz_t(), factors_->p+i);
            result.emplace_back(std::move(prime),factors_->exp[i]);
        }
        return result;
    }
};
