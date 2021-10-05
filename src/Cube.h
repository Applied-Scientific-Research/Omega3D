/*
 * Cube.h - Two versions of a simple cube
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <cstdint>
#include <vector>

// two tris per side
const std::vector<float> cube0 = {{
0.0, 0.0, 0.0,
0.0, 0.0, 1.0,
0.0, 1.0, 0.0,
0.0, 1.0, 1.0,
1.0, 0.0, 0.0,
1.0, 0.0, 1.0,
1.0, 1.0, 0.0,
1.0, 1.0, 1.0
}};

const std::vector<uint32_t> cube0idx = {{
0,4,5,
0,5,1,
4,6,7,
4,7,5,
1,5,7,
1,7,3,
6,2,3,
6,3,7,
2,0,1,
2,1,3,
2,6,4,
2,4,0
}};


// 8 tris per side - corners all have 6 tris meet there
// THIS IS INCOMPLETE
const std::vector<float> cube1 = {{
0.0, 0.0, 0.0,
0.0, 0.0, 1.0,
0.0, 1.0, 0.0,
0.0, 1.0, 1.0,
1.0, 0.0, 0.0,
1.0, 0.0, 1.0,
1.0, 1.0, 0.0,
1.0, 1.0, 1.0
}};

const std::vector<uint32_t> cube1idx = {{
0,4,5,
0,5,1,
4,6,7,
4,7,5,
1,5,7,
1,7,3,
6,2,3,
6,3,7,
2,0,1,
2,1,3,
2,6,4,
2,4,0
}};

