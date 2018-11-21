#pragma once

#include <cmath>
#include <array>

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

