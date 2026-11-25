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
// The allocation benchmark injects scope markers at compile time.
#ifndef CENSUS_HEAP_SCOPE
#define CENSUS_HEAP_SCOPE(phase) ((void)0)
#define CENSUS_HEAP_EVAL(phase, ...) (__VA_ARGS__)
#define CENSUS_HEAP_RUN() ((void)0)
#endif
