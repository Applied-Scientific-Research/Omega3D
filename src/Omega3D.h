/*
 * Omega3D.h - Useful definitions for anywhere in the code
 *
 * (c)2018 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <cstdlib>

const size_t Dimensions = 3;

// solver order and method
enum solution_t {
  direct_cpu   = 1,
  direct_vc    = 2,
  direct_glsl  = 3,
  treecode_cpu = 4,
  treecode_vc  = 5
};

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

