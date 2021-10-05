/*
 * MathHelper.h - Useful functions to operate on 3-tuples ("vectors", not "std::vector")
 *
 * (c)2018 Applied Scientific Research, Inc.
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

#ifdef _WIN32
#define __restrict__ __restrict
#endif

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <cmath>
#include <array>

// helper functions: exp, recip, sqrt, rsqrt, oor1p5

#ifdef USE_VC
template <class S>
static inline S my_exp(const S _in) {
  return Vc::exp(_in);
}
template <>
inline float my_exp(const float _in) {
  return std::exp(_in);
}
template <>
inline double my_exp(const double _in) {
  return std::exp(_in);
}
#else
template <class S>
static inline S my_exp(const S _in) {
  return std::exp(_in);
}
#endif

#ifdef USE_VC
template <class S>
static inline S my_recip(const S _in) {
  return Vc::reciprocal(_in);
}
template <>
inline float my_recip(const float _in) {
  return 1.0f / _in;
}
template <>
inline double my_recip(const double _in) {
  return 1.0 / _in;
}
#else
template <class S>
static inline S my_recip(const S _in) {
  return S(1.0) / _in;
}
#endif

#ifdef USE_VC
template <class S>
static inline S my_sqrt(const S _in) {
  return Vc::sqrt(_in);
}
template <>
inline float my_sqrt(const float _in) {
  return std::sqrt(_in);
}
template <>
inline double my_sqrt(const double _in) {
  return std::sqrt(_in);
}
#else
template <class S>
static inline S my_sqrt(const S _in) {
  return std::sqrt(_in);
}
#endif

#ifdef USE_VC
template <class S>
static inline S my_rsqrt(const S _in) {
  return Vc::rsqrt(_in);
}
template <>
inline float my_rsqrt(const float _in) {
  return 1.0f / std::sqrt(_in);
}
template <>
inline double my_rsqrt(const double _in) {
  return 1.0 / std::sqrt(_in);
}
#else
template <class S>
static inline S my_rsqrt(const S _in) {
  return S(1.0) / std::sqrt(_in);
}
#endif

#ifdef USE_VC
template <class S>
static inline S oor2p5(const S _in) {
  //return Vc::reciprocal(_in*_in*Vc::sqrt(_in));	// 234 GFlop/s
  return Vc::rsqrt(_in) * Vc::reciprocal(_in*_in);	// 269 GFlop/s
}
template <>
inline float oor2p5(const float _in) {
  return 1.0f / (_in*_in*std::sqrt(_in));
}
template <>
inline double oor2p5(const double _in) {
  return 1.0 / (_in*_in*std::sqrt(_in));
}
#else
template <class S>
static inline S oor2p5(const S _in) {
  return S(1.0) / (_in*_in*std::sqrt(_in));
}
#endif

#ifdef USE_VC
template <class S>
static inline S oor1p5(const S _in) {
  //return Vc::reciprocal(_in*Vc::sqrt(_in));		// 243 GFlop/s
  return Vc::rsqrt(_in) * Vc::reciprocal(_in);		// 302 GFlop/s
}
template <>
inline float oor1p5(const float _in) {
  return 1.0f / (_in*std::sqrt(_in));
}
template <>
inline double oor1p5(const double _in) {
  return 1.0 / (_in*std::sqrt(_in));
}
#else
template <class S>
static inline S oor1p5(const S _in) {
  return S(1.0) / (_in*std::sqrt(_in));
}
#endif

#ifdef USE_VC
template <class S>
static inline S oor0p75(const S _in) {
  const S rsqd = Vc::rsqrt(_in);
  //return rsqd*Vc::sqrt(rsqd);				// 265 GFlop/s
  return rsqd*rsqd*Vc::rsqrt(rsqd);			// 301 GFlop/s
}
template <>
inline float oor0p75(const float _in) {
  const float sqd = std::sqrt(_in);
  return 1.0f / (sqd*std::sqrt(sqd));
}
template <>
inline double oor0p75(const double _in) {
  const double sqd = std::sqrt(_in);
  return 1.0 / (sqd*std::sqrt(sqd));
}
#else
template <class S>
static inline S oor0p75(const S _in) {
  const S sqd = std::sqrt(_in);
  return S(1.0) / (sqd*std::sqrt(sqd));
}
#endif


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

