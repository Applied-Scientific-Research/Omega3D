/*
 * GeomHelper.h - Generators for basic solid shapes
 *
 * (c)2021 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
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

