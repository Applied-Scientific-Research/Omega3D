/*
 * Kernels.h - Non-class inner kernels for influence calculations
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#ifdef _WIN32
#define __restrict__ __restrict
#endif

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <cmath>


// velocity influence functions

// thick-cored particle on thick-cored point, no gradients
template <class S, class A>
static inline void kernel_0_0s (const S sx, const S sy, const S sz,
                                const S sr,
                                const S ssx, const S ssy, const S ssz,
                                const S tx, const S ty, const S tz,
                                const S tr,
                                A* __restrict__ tu, A* __restrict__ tv, A* __restrict__ tw) {
  // 30 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr + tr*tr;
#ifdef USE_VC
  r2 = Vc::reciprocal(r2*Vc::sqrt(r2));
#else
  r2 = 1.0 / (r2*std::sqrt(r2));
#endif
  const S dxxw = dz*ssy - dy*ssz;
  const S dyxw = dx*ssz - dz*ssx;
  const S dzxw = dy*ssx - dx*ssy;
  *tu += r2 * dxxw;
  *tv += r2 * dyxw;
  *tw += r2 * dzxw;
}

// thick-cored particle on singular point, no gradients
template <class S, class A>
static inline void kernel_0v_0p (const S sx, const S sy, const S sz,
                                 const S sr,
                                 const S ssx, const S ssy, const S ssz,
                                 const S tx, const S ty, const S tz,
                                 A* __restrict__ tu, A* __restrict__ tv, A* __restrict__ tw) {
  // 28 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr;
#ifdef USE_VC
  r2 = Vc::reciprocal(r2*Vc::sqrt(r2));
#else
  r2 = 1.0 / (r2*std::sqrt(r2));
#endif
  const S dxxw = dz*ssy - dy*ssz;
  const S dyxw = dx*ssz - dz*ssx;
  const S dzxw = dy*ssx - dx*ssy;
  *tu += r2 * dxxw;
  *tv += r2 * dyxw;
  *tw += r2 * dzxw;
}

template <class S, class A>
static inline void kernel_0_0v (const S* __restrict__ sx, const S __restrict__ sr, const S* __restrict__ ss,
                                const S* __restrict__ tx, const S __restrict__ tr, A* __restrict__ tu) {
  // 30 flops
  const S dx = tx[0] - sx[0];
  const S dy = tx[1] - sx[1];
  const S dz = tx[2] - sx[2];
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr + tr*tr;
  r2 = 1.0 / (r2*std::sqrt(r2));
  const S dxxw = dz*ss[1] - dy*ss[2];
  const S dyxw = dx*ss[2] - dz*ss[0];
  const S dzxw = dy*ss[0] - dx*ss[1];
  tu[0] += r2 * dxxw;
  tu[1] += r2 * dyxw;
  tu[2] += r2 * dzxw;
}

// thick-cored particle on thick-cored point, with gradients
template <class S, class A>
static inline void kernel_0_0sg (const S sx, const S sy, const S sz,
                                 const S sr,
                                 const S ssx, const S ssy, const S ssz,
                                 const S tx, const S ty, const S tz,
                                 const S tr,
                                 A* __restrict__ tu, A* __restrict__ tv, A* __restrict__ tw,
                                 A* __restrict__ tux, A* __restrict__ tvx, A* __restrict__ twx,
                                 A* __restrict__ tuy, A* __restrict__ tvy, A* __restrict__ twy,
                                 A* __restrict__ tuz, A* __restrict__ tvz, A* __restrict__ twz) {
  // 30 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  const S r2 = dx*dx + dy*dy + dz*dz + sr*sr + tr*tr;
#ifdef USE_VC
  const S r3 = Vc::reciprocal(r2*Vc::sqrt(r2));
#else
  const S r3 = 1.0 / (r2*std::sqrt(r2));
#endif
  S dxxw = dz*ssy - dy*ssz;
  S dyxw = dx*ssz - dz*ssx;
  S dzxw = dy*ssx - dx*ssy;
  *tu += r3 * dxxw;
  *tv += r3 * dyxw;
  *tw += r3 * dzxw;

  // accumulate velocity gradients
#ifdef USE_VC
  //const S bbb = -3.f * r3 * Vc::reciprocal(r2);
  const S bbb = S(-3.0) * r3 * Vc::reciprocal(r2);
#else
  const S bbb = -3.0 * r3 / r2;
#endif
  // continuing with grads - this section is 33 flops
  dxxw *= bbb;
  dyxw *= bbb;
  dzxw *= bbb;
  *tux += dx*dxxw;
  *tvx += dx*dyxw + ssz*r3;
  *twx += dx*dzxw - ssy*r3;
  *tuy += dy*dxxw - ssz*r3;
  *tvy += dy*dyxw;
  *twy += dy*dzxw + ssx*r3;
  *tuz += dz*dxxw + ssy*r3;
  *tvz += dz*dyxw - ssx*r3;
  *twz += dz*dzxw;
}

// thick-cored particle on singular point, with gradients
template <class S, class A>
static inline void kernel_0v_0pg (const S sx, const S sy, const S sz,
                                  const S sr,
                                  const S ssx, const S ssy, const S ssz,
                                  const S tx, const S ty, const S tz,
                                  A* __restrict__ tu, A* __restrict__ tv, A* __restrict__ tw,
                                  A* __restrict__ tux, A* __restrict__ tvx, A* __restrict__ twx,
                                  A* __restrict__ tuy, A* __restrict__ tvy, A* __restrict__ twy,
                                  A* __restrict__ tuz, A* __restrict__ tvz, A* __restrict__ twz) {
  // 30 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  const S r2 = dx*dx + dy*dy + dz*dz + sr*sr;
#ifdef USE_VC
  const S r3 = Vc::reciprocal(r2*Vc::sqrt(r2));
#else
  const S r3 = 1.0 / (r2*std::sqrt(r2));
#endif
  S dxxw = dz*ssy - dy*ssz;
  S dyxw = dx*ssz - dz*ssx;
  S dzxw = dy*ssx - dx*ssy;
  *tu += r3 * dxxw;
  *tv += r3 * dyxw;
  *tw += r3 * dzxw;

  // accumulate velocity gradients
#ifdef USE_VC
  //const S bbb = -3.f * r3 * Vc::reciprocal(r2);
  const S bbb = S(-3.0) * r3 * Vc::reciprocal(r2);
#else
  const S bbb = -3.0 * r3 / r2;
#endif
  // continuing with grads - this section is 33 flops
  dxxw *= bbb;
  dyxw *= bbb;
  dzxw *= bbb;
  *tux += dx*dxxw;
  *tvx += dx*dyxw + ssz*r3;
  *twx += dx*dzxw - ssy*r3;
  *tuy += dy*dxxw - ssz*r3;
  *tvy += dy*dyxw;
  *twy += dy*dzxw + ssx*r3;
  *tuz += dz*dxxw + ssy*r3;
  *tvz += dz*dyxw - ssx*r3;
  *twz += dz*dzxw;
}

template <class S, class A>
static inline void kernel_0_0vg (const S* __restrict__ sx, const S __restrict__ sr, const S* __restrict__ ss,
                                 const S* __restrict__ tx, const S __restrict__ tr, A* __restrict__ tu) {
  // 30 flops
  const S dx = tx[0] - sx[0];
  const S dy = tx[1] - sx[1];
  const S dz = tx[2] - sx[2];
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr + tr*tr;
  r2 = 1.0 / (r2*std::sqrt(r2));
  S dxxw = dz*ss[1] - dy*ss[2];
  S dyxw = dx*ss[2] - dz*ss[0];
  S dzxw = dy*ss[0] - dx*ss[1];
  tu[0] += r2 * dxxw;
  tu[1] += r2 * dyxw;
  tu[2] += r2 * dzxw;

  // HACK - you need to figure out what this term is
  const S bbb = r2 / std::sqrt(r2);
  // continuing with grads - this section is 33 flops
  dxxw *= bbb;
  dyxw *= bbb;
  dzxw *= bbb;
  tu[3] += dx*dxxw;
  tu[4] += dx*dyxw + ss[2]*r2;
  tu[5] += dx*dzxw - ss[1]*r2;
  tu[6] += dy*dxxw - ss[2]*r2;
  tu[7] += dy*dyxw;
  tu[8] += dy*dzxw + ss[0]*r2;
  tu[9] += dz*dxxw + ss[1]*r2;
  tu[10] += dz*dyxw - ss[0]*r2;
  tu[11] += dz*dzxw;
}


//
// influence of 3d linear constant-strength *vortex* panel on target point
//   ignoring the 1/4pi factor, which will be multiplied later
// uses four thick-cored integration points on the source side
//   168 flops
//
template <class S, class A>
static inline void kernel_1_0v (const S sx0, const S sy0, const S sz0,
                                const S sx1, const S sy1, const S sz1,
                                const S sx2, const S sy2, const S sz2,
                                const S ssx, const S ssy, const S ssz,
                                const S tx, const S ty, const S tz,
                                A* __restrict__ tu, A* __restrict__ tv, A* __restrict__ tw) {

  // scale the strength by 1/4, to account for the 4 calls below (3 flops)
  const S strx = S(0.25) * ssx;
  const S stry = S(0.25) * ssy;
  const S strz = S(0.25) * ssz;

  // first source point (39 flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (30 flops)
    (void) kernel_0_0s (sx, sy, sz, (S)0.0,
                        strx, stry, strz,
                        tx, ty, tz, (S)0.0,
                        tu, tv, tw);
  }

  // second source point (42)
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0_0s (sx, sy, sz, (S)0.0,
                        strx, stry, strz,
                        tx, ty, tz, (S)0.0,
                        tu, tv, tw);
  }

  // third source point (42)
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0_0s (sx, sy, sz, (S)0.0,
                        strx, stry, strz,
                        tx, ty, tz, (S)0.0,
                        tu, tv, tw);
  }

  // final source point (42)
  {
    const S sx = (sx0 + sx1 + S(4.0)*sx2) / S(6.0);
    const S sy = (sy0 + sy1 + S(4.0)*sy2) / S(6.0);
    const S sz = (sz0 + sz1 + S(4.0)*sz2) / S(6.0);
    (void) kernel_0_0s (sx, sy, sz, (S)0.0,
                        strx, stry, strz,
                        tx, ty, tz, (S)0.0,
                        tu, tv, tw);
  }
}

