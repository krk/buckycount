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

// Direct unit-triangle centers in a simple developed shortest-star polygon.
#pragma once
#include "genfullerene.hpp"

namespace gf {
struct StarCenterStats {
    U graphs = 0, centers = 0, boundary_centers = 0, darts = 0, crossings = 0;
    U maximum_crossings = 0;
};
struct StarCenterHash {
    std::size_t operator()(Point p) const {
        return static_cast<std::size_t>((static_cast<U>(p.x) * 0x9e3779b97f4a7c15ULL) ^
                                       (static_cast<U>(p.y) * 0xbf58476d1ce4e5b9ULL));
    }
};
struct StarCenterSeam {
    Point a, b, offset;
    int phase;
    Point image(Point p) const { return rotate(p, phase) + offset; }
};
I star_center_floor(I a, I b) {
    need(b > 0, "nonpositive center divisor");
    I q = a / b, r = a % b;
    return q - (r < 0);
}
I star_center_ceil(I a, I b) { return -star_center_floor(-a, b); }
Patch star_center_triangle(Point p) {
    I r = p.x - 3 * star_center_floor(p.x, 3);
    need((r == 1 || r == 2) && (p.y - r) % 3 == 0, "invalid triangle center");
    I x = (p.x - r) / 3, y = (p.y - r) / 3;
    if (r == 1) return {{{x,y},{x+1,y},{x,y+1}}};
    return {{{x+1,y+1},{x,y+1},{x+1,y}}};
}

} // namespace gf
