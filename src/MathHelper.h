/*
 * MathHelper.h - Useful functions to operate on 3-tuples ("vectors", not "std::vector")
 *
 * (c)2018 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <cmath>
#include <array>

// normalize a 3-vector
template <class S>
inline void normalizeVec (std::array<S,3>& in) {
  // note there is no std::rsqrt()
  const S len = 1.0 / std::sqrt(in[0]*in[0] + in[1]*in[1] + in[2]*in[2]);
  in[0] *= len;
  in[1] *= len;
  in[2] *= len;
}

// basic dot product, like Fortran
template <class S>
inline S dot_product (std::array<S,3> const & v1, std::array<S,3> const & v2) {
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

// length of a 3-tuple
template <class S>
inline S length (std::array<S,3> const & vec) {
  return std::sqrt(dot_product<S>(vec, vec));
}

// basic cross product v1 x v2
template <class S>
inline void cross_product (std::array<S,3> const & v1, std::array<S,3> const & v2, std::array<S,3>& result) {
  result[0] = v1[1]*v2[2] - v1[2]*v2[1];
  result[1] = v1[2]*v2[0] - v1[0]*v2[2];
  result[2] = v1[0]*v2[1] - v1[1]*v2[0];
  return;
}

//
// from a single 3D vector, generate a set of orthonormal basis vectors
//   code is from https://graphics.pixar.com/library/OrthonormalB/paper.pdf
//
template <class S>
inline void branchlessONB (const std::array<S,3>& n, std::array<S,3>& b1, std::array<S,3>& b2) {
  const S sign = std::copysign(1.0, n[2]);
  const S a = -1.0 / (sign + n[2]);
  const S b = n[0] * n[1] * a;
  b1[0] = 1.0 + sign * n[0] * n[0] * a;
  b1[1] = sign * b;
  b1[2] = -sign * n[0];
  b2[0] = b;
  b2[1] = sign + n[1] * n[1] * a;
  b2[2] = -n[1];
}

