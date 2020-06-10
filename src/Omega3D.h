/*
 * Omega3D.h - Useful definitions for anywhere in the code
 *
 * (c)2018-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#pragma once

// Use this for indexes into panels or bodies
// using 32-bit because we may have more than 65536 triangles/elements in the system
#include <cstdint>
using Int = uint32_t;
#include <cstdlib>

const size_t Dimensions = 3;

// element type
enum elem_t {
  active   = 1,  // active vorticity
  reactive = 2,  // active once strength is solved
  inert    = 3   // does not affect flow
};

// movement type
enum move_t {
  lagrangian = 1, // moves with local velocity
  bodybound  = 2, // moves with attached body
  fixed      = 3  // does not move
};

