/*
 * GeomHelper.h - Generators for basic solid shapes
 *
 * (c)2021 Applied Scientific Research, Inc.
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

#include "ElementPacket.h"

// Create a closed spherical/ovoid object from a icosahedron
ElementPacket<float> generate_ovoid(const float, const float, const float, const float);

// Create a closed rectangular object
ElementPacket<float> generate_cuboid(const float, const float, const float, const float);

// Create a closed torus
ElementPacket<float> generate_torus(const float, const float, const float);

// Create a closed discoid
ElementPacket<float> generate_discoid(const float, const float, const float);

