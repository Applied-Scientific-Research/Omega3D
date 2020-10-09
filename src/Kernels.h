/*
 * Kernels.h - Non-class inner kernels for influence calculations
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#ifdef _WIN32
#define __restrict__ __restrict
#endif

#include "CoreFunc.h"

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <cmath>


//
// velocity influence functions
//
// here is the naming system:
//   kernel_NS_MT
//     N is the number of dimensions of the source element (0=point, 2=surface)
//     M is the number of dimensions of the target element
//     S is the type of the source element ('v'=vortex, 's'=source, 'vs'=vortex and source)
//     T is the type of the target element
//         first character is 'p' for a singular point, 'b' for a vortex blob
//         second character is 'g' if gradients must be returned
//

// thick-cored particle on thick-cored point, no gradients
template <class S> inline size_t flops_0v_0b () { return 23 + flops_tv_nograds<S>(); }
template <class S, class A>
static inline void kernel_0v_0b (const S sx, const S sy, const S sz,
                                 const S sr,
                                 const S ssx, const S ssy, const S ssz,
                                 const S tx, const S ty, const S tz,
                                 const S tr,
                                 A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {
  // 23+(7|12) flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  const S r2 = core_func<S>(dx*dx + dy*dy + dz*dz, sr, tr);
  const S dxxw = dz*ssy - dy*ssz;
  const S dyxw = dx*ssz - dz*ssx;
  const S dzxw = dy*ssx - dx*ssy;
  *tu += r2 * dxxw;
  *tv += r2 * dyxw;
  *tw += r2 * dzxw;
}

// same, but vortex+source strength
template <class S> inline size_t flops_0vs_0b () { return 29 + flops_tv_nograds<S>(); }
template <class S, class A>
static inline void kernel_0vs_0b (const S sx, const S sy, const S sz,
                                 const S sr,
                                 const S ssx, const S ssy, const S ssz, const S ss,
                                 const S tx, const S ty, const S tz,
                                 const S tr,
                                 A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {
  // 29+(7|12) flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  const S r2 = core_func<S>(dx*dx + dy*dy + dz*dz, sr, tr);
  const S dxxw = dz*ssy - dy*ssz + dx*ss;
  const S dyxw = dx*ssz - dz*ssx + dy*ss;
  const S dzxw = dy*ssx - dx*ssy + dz*ss;
  *tu += r2 * dxxw;
  *tv += r2 * dyxw;
  *tw += r2 * dzxw;
}

// thick-cored particle on singular point, no gradients
template <class S> inline size_t flops_0v_0p () { return 23 + flops_tp_nograds<S>(); }
template <class S, class A>
static inline void kernel_0v_0p (const S sx, const S sy, const S sz,
                                 const S sr,
                                 const S ssx, const S ssy, const S ssz,
                                 const S tx, const S ty, const S tz,
                                 A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {
  // 23+(5|9) flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  const S r2 = core_func<S>(dx*dx + dy*dy + dz*dz, sr);
  const S dxxw = dz*ssy - dy*ssz;
  const S dyxw = dx*ssz - dz*ssx;
  const S dzxw = dy*ssx - dx*ssy;
  *tu += r2 * dxxw;
  *tv += r2 * dyxw;
  *tw += r2 * dzxw;
}

// same, but vortex+source strength
template <class S> inline size_t flops_0vs_0p () { return 29 + flops_tp_nograds<S>(); }
template <class S, class A>
static inline void kernel_0vs_0p (const S sx, const S sy, const S sz,
                                  const S sr,
                                  const S ssx, const S ssy, const S ssz, const S ss,
                                  const S tx, const S ty, const S tz,
                                  A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {
  // 29+(5|9) flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  const S r2 = core_func<S>(dx*dx + dy*dy + dz*dz, sr);
  const S dxxw = dz*ssy - dy*ssz + dx*ss;
  const S dyxw = dx*ssz - dz*ssx + dy*ss;
  const S dzxw = dy*ssx - dx*ssy + dz*ss;
  *tu += r2 * dxxw;
  *tv += r2 * dyxw;
  *tw += r2 * dzxw;
}

// same, but source strength only
template <class S> inline size_t flops_0s_0p () { return 14 + flops_tp_nograds<S>(); }
template <class S, class A>
static inline void kernel_0s_0p (const S sx, const S sy, const S sz,
                                 const S sr,
                                 const S ss,
                                 const S tx, const S ty, const S tz,
                                 A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {
  // 14+(5|9) flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  const S r2 = ss * core_func<S>(dx*dx + dy*dy + dz*dz, sr);
  *tu += r2 * dx;
  *tv += r2 * dy;
  *tw += r2 * dz;
}

// thick-cored particle on thick-cored point, with gradients
//   54+(9|14) flops total
template <class S> inline size_t flops_0v_0bg () { return 54 + flops_tv_grads<S>(); }
template <class S, class A>
static inline void kernel_0v_0bg (const S sx, const S sy, const S sz,
                                  const S sr,
                                  const S ssx, const S ssy, const S ssz,
                                  const S tx, const S ty, const S tz,
                                  const S tr,
                                  A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw,
                                  A* const __restrict__ tux, A* const __restrict__ tvx, A* const __restrict__ twx,
                                  A* const __restrict__ tuy, A* const __restrict__ tvy, A* const __restrict__ twy,
                                  A* const __restrict__ tuz, A* const __restrict__ tvz, A* const __restrict__ twz) {
  // 21 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  S r3, bbb;
  (void) core_func<S>(dx*dx + dy*dy + dz*dz, sr, tr, &r3, &bbb);
  S dxxw = dz*ssy - dy*ssz;
  S dyxw = dx*ssz - dz*ssx;
  S dzxw = dy*ssx - dx*ssy;
  *tu += r3 * dxxw;
  *tv += r3 * dyxw;
  *tw += r3 * dzxw;

  // accumulate velocity gradients
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

// same, but for vortex+source strengths
//   79+(9|14) flops total
template <class S> inline size_t flops_0vs_0bg () { return 79 + flops_tv_grads<S>(); }
template <class S, class A>
static inline void kernel_0vs_0bg (const S sx, const S sy, const S sz,
                                   const S sr,
                                   const S ssx, const S ssy, const S ssz, const S ss,
                                   const S tx, const S ty, const S tz,
                                   const S tr,
                                   A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw,
                                   A* const __restrict__ tux, A* const __restrict__ tvx, A* const __restrict__ twx,
                                   A* const __restrict__ tuy, A* const __restrict__ tvy, A* const __restrict__ twy,
                                   A* const __restrict__ tuz, A* const __restrict__ tvz, A* const __restrict__ twz) {
  // 36 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  S r3, bbb;
  (void) core_func<S>(dx*dx + dy*dy + dz*dz, sr, tr, &r3, &bbb);
  S dxxw = dz*ssy - dy*ssz;
  S dyxw = dx*ssz - dz*ssx;
  S dzxw = dy*ssx - dx*ssy;
  *tu += r3 * (dxxw + dx*ss);
  *tv += r3 * (dyxw + dy*ss);
  *tw += r3 * (dzxw + dz*ss);

  // accumulate velocity gradients
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
  // and the grads due to the source term - another 19 flops
  const S dxs = dx*bbb*ss;
  const S dys = dy*bbb*ss;
  const S dzs = dz*bbb*ss;
  const S dss = ss*r3;
  *tux += dx*dxs + dss;
  *tvx += dx*dys;
  *twx += dx*dzs;
  *tuy += dy*dxs;
  *tvy += dy*dys + dss;
  *twy += dy*dzs;
  *tuz += dz*dxs;
  *tvz += dz*dys;
  *twz += dz*dzs + dss;
}

// thick-cored particle on singular point, with gradients
//   54+(7|11) flops total
template <class S> inline size_t flops_0v_0pg () { return 54 + flops_tp_grads<S>(); }
template <class S, class A>
static inline void kernel_0v_0pg (const S sx, const S sy, const S sz,
                                  const S sr,
                                  const S ssx, const S ssy, const S ssz,
                                  const S tx, const S ty, const S tz,
                                  A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw,
                                  A* const __restrict__ tux, A* const __restrict__ tvx, A* const __restrict__ twx,
                                  A* const __restrict__ tuy, A* const __restrict__ tvy, A* const __restrict__ twy,
                                  A* const __restrict__ tuz, A* const __restrict__ tvz, A* const __restrict__ twz) {
  // 28 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  S r3, bbb;
  (void) core_func<S>(dx*dx + dy*dy + dz*dz, sr, &r3, &bbb);
  S dxxw = dz*ssy - dy*ssz;
  S dyxw = dx*ssz - dz*ssx;
  S dzxw = dy*ssx - dx*ssy;
  *tu += r3 * dxxw;
  *tv += r3 * dyxw;
  *tw += r3 * dzxw;

  // accumulate velocity gradients
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

// same, but for vortex+source strengths
//   79+(7|11) flops total
template <class S> inline size_t flops_0vs_0pg () { return 79 + flops_tp_grads<S>(); }
template <class S, class A>
static inline void kernel_0vs_0pg (const S sx, const S sy, const S sz,
                                   const S sr,
                                   const S ssx, const S ssy, const S ssz, const S ss,
                                   const S tx, const S ty, const S tz,
                                   A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw,
                                   A* const __restrict__ tux, A* const __restrict__ tvx, A* const __restrict__ twx,
                                   A* const __restrict__ tuy, A* const __restrict__ tvy, A* const __restrict__ twy,
                                   A* const __restrict__ tuz, A* const __restrict__ tvz, A* const __restrict__ twz) {
  // 34 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  S r3, bbb;
  (void) core_func<S>(dx*dx + dy*dy + dz*dz, sr, &r3, &bbb);
  S dxxw = dz*ssy - dy*ssz;
  S dyxw = dx*ssz - dz*ssx;
  S dzxw = dy*ssx - dx*ssy;
  *tu += r3 * (dxxw + dx*ss);
  *tv += r3 * (dyxw + dy*ss);
  *tw += r3 * (dzxw + dz*ss);

  // accumulate velocity gradients
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
  // and the grads due to the source term - another 19 flops
  const S dxs = dx*bbb*ss;
  const S dys = dy*bbb*ss;
  const S dzs = dz*bbb*ss;
  const S dss = ss*r3;
  *tux += dx*dxs + dss;
  *tvx += dx*dys;
  *twx += dx*dzs;
  *tuy += dy*dxs;
  *tvy += dy*dys + dss;
  *twy += dy*dzs;
  *tuz += dz*dxs;
  *tvz += dz*dys;
  *twz += dz*dzs + dss;
}


//
// influence of 3d linear constant-strength *vortex* panel on target point
//   ignoring the 1/4pi factor, which will be multiplied later
// uses four singular integration points on the source side
//   140+(20|36) flops
//
template <class S> inline size_t flops_2v_0p () { return 48 + 4*flops_0v_0p<S>(); }
template <class S, class A>
static inline void kernel_2v_0p (const S sx0, const S sy0, const S sz0,
                                 const S sx1, const S sy1, const S sz1,
                                 const S sx2, const S sy2, const S sz2,
                                 const S ssx, const S ssy, const S ssz,
                                 const S tx, const S ty, const S tz,
                                 A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {

  // scale the strength by 1/4, to account for the 4 calls below (3 flops)
  const S strx = S(0.25) * ssx;
  const S stry = S(0.25) * ssy;
  const S strz = S(0.25) * ssz;

  // first source point (32+(5|9) flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (23+(5|9) flops)
    (void) kernel_0v_0p (sx, sy, sz, (S)0.0,
                        strx, stry, strz,
                        tx, ty, tz,
                        tu, tv, tw);
  }

  // second source point (35+(5|9))
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0v_0p (sx, sy, sz, (S)0.0,
                        strx, stry, strz,
                        tx, ty, tz,
                        tu, tv, tw);
  }

  // third source point (35+(5|9))
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0v_0p (sx, sy, sz, (S)0.0,
                        strx, stry, strz,
                        tx, ty, tz,
                        tu, tv, tw);
  }

  // final source point (35+(5|9))
  {
    const S sx = (sx0 + sx1 + S(4.0)*sx2) / S(6.0);
    const S sy = (sy0 + sy1 + S(4.0)*sy2) / S(6.0);
    const S sz = (sz0 + sz1 + S(4.0)*sz2) / S(6.0);
    (void) kernel_0v_0p (sx, sy, sz, (S)0.0,
                        strx, stry, strz,
                        tx, ty, tz,
                        tu, tv, tw);
  }
}

// same for vortex+source terms
//   165+(20|36) flops
template <class S> inline size_t flops_2vs_0p () { return 49 + 4*flops_0vs_0p<S>(); }
template <class S, class A>
static inline void kernel_2vs_0p (const S sx0, const S sy0, const S sz0,
                                  const S sx1, const S sy1, const S sz1,
                                  const S sx2, const S sy2, const S sz2,
                                  const S ssx, const S ssy, const S ssz, const S ss,
                                  const S tx, const S ty, const S tz,
                                  A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {

  // scale the strength by 1/4, to account for the 4 calls below (4 flops)
  const S strx = S(0.25) * ssx;
  const S stry = S(0.25) * ssy;
  const S strz = S(0.25) * ssz;
  const S strs = S(0.25) * ss;

  // first source point (38+(5|9) flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (29+(5|9) flops)
    (void) kernel_0vs_0p (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz,
                          tu, tv, tw);
  }

  // second source point (41+(5|9))
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0vs_0p (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz,
                          tu, tv, tw);
  }

  // third source point (41+(5|9))
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0vs_0p (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz,
                          tu, tv, tw);
  }

  // final source point (41+(5|9))
  {
    const S sx = (sx0 + sx1 + S(4.0)*sx2) / S(6.0);
    const S sy = (sy0 + sy1 + S(4.0)*sy2) / S(6.0);
    const S sz = (sz0 + sz1 + S(4.0)*sz2) / S(6.0);
    (void) kernel_0vs_0p (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz,
                          tu, tv, tw);
  }
}

// same, but for source source terms only
//   102+(20|36) flops
template <class S> inline size_t flops_2s_0p () { return 46 + 4*flops_0s_0p<S>(); }
template <class S, class A>
static inline void kernel_2s_0p (const S sx0, const S sy0, const S sz0,
                                 const S sx1, const S sy1, const S sz1,
                                 const S sx2, const S sy2, const S sz2,
                                 const S ss,
                                 const S tx, const S ty, const S tz,
                                 A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {

  // scale the strength by 1/4, to account for the 4 calls below (1 flop)
  const S strs = S(0.25) * ss;

  // first source point (23+(5|9) flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (14+(5|9) flops)
    (void) kernel_0s_0p (sx, sy, sz, (S)0.0,
                         strs,
                         tx, ty, tz,
                         tu, tv, tw);
  }

  // second source point (26+(5|9))
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0s_0p (sx, sy, sz, (S)0.0,
                         strs,
                         tx, ty, tz,
                         tu, tv, tw);
  }

  // third source point (26+(5|9))
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0s_0p (sx, sy, sz, (S)0.0,
                         strs,
                         tx, ty, tz,
                         tu, tv, tw);
  }

  // final source point (26+(5|9))
  {
    const S sx = (sx0 + sx1 + S(4.0)*sx2) / S(6.0);
    const S sy = (sy0 + sy1 + S(4.0)*sy2) / S(6.0);
    const S sz = (sz0 + sz1 + S(4.0)*sz2) / S(6.0);
    (void) kernel_0s_0p (sx, sy, sz, (S)0.0,
                         strs,
                         tx, ty, tz,
                         tu, tv, tw);
  }
}


//
// influence of 3d linear constant-strength *vortex* panel on thick-cored target vortex
//   ignoring the 1/4pi factor, which will be multiplied later
// uses four singular integration points on the source side
//   140+(28|48) flops
//
template <class S> inline size_t flops_2v_0b () { return 48 + 4*flops_0v_0b<S>(); }
template <class S, class A>
static inline void kernel_2v_0b (const S sx0, const S sy0, const S sz0,
                                 const S sx1, const S sy1, const S sz1,
                                 const S sx2, const S sy2, const S sz2,
                                 const S ssx, const S ssy, const S ssz,
                                 const S tx, const S ty, const S tz, const S tr,
                                 A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {

  // scale the strength by 1/4, to account for the 4 calls below (3 flops)
  const S strx = S(0.25) * ssx;
  const S stry = S(0.25) * ssy;
  const S strz = S(0.25) * ssz;
  //const S targr = S(0.1) * tr;

  // first source point (32+(7|12) flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (23+(7|12) flops)
    (void) kernel_0v_0b (sx, sy, sz, (S)0.0,
                         strx, stry, strz,
                         tx, ty, tz, tr,
                         tu, tv, tw);
  }

  // second source point (35+(7|12))
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0v_0b (sx, sy, sz, (S)0.0,
                         strx, stry, strz,
                         tx, ty, tz, tr,
                         tu, tv, tw);
  }

  // third source point (35+(7|12))
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0v_0b (sx, sy, sz, (S)0.0,
                         strx, stry, strz,
                         tx, ty, tz, tr,
                         tu, tv, tw);
  }

  // final source point (35+(7|12))
  {
    const S sx = (sx0 + sx1 + S(4.0)*sx2) / S(6.0);
    const S sy = (sy0 + sy1 + S(4.0)*sy2) / S(6.0);
    const S sz = (sz0 + sz1 + S(4.0)*sz2) / S(6.0);
    (void) kernel_0v_0b (sx, sy, sz, (S)0.0,
                         strx, stry, strz,
                         tx, ty, tz, tr,
                         tu, tv, tw);
  }
}

// same thing, but for vortex+source strengths
//   165+(28|48) flops
template <class S> inline size_t flops_2vs_0b () { return 49 + 4*flops_0vs_0b<S>(); }
template <class S, class A>
static inline void kernel_2vs_0b (const S sx0, const S sy0, const S sz0,
                                  const S sx1, const S sy1, const S sz1,
                                  const S sx2, const S sy2, const S sz2,
                                  const S ssx, const S ssy, const S ssz, const S ss,
                                  const S tx, const S ty, const S tz, const S tr,
                                  A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {

  // scale the strength by 1/4, to account for the 4 calls below (4 flops)
  const S strx = S(0.25) * ssx;
  const S stry = S(0.25) * ssy;
  const S strz = S(0.25) * ssz;
  const S strs = S(0.25) * ss;
  //const S targr = S(0.1) * tr;

  // first source point (38+(7|12) flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (29+(7|12) flops)
    (void) kernel_0vs_0b (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz, tr,
                          tu, tv, tw);
  }

  // second source point (41+(7|12))
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0vs_0b (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz, tr,
                          tu, tv, tw);
  }

  // third source point (41+(7|12))
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0vs_0b (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz, tr,
                          tu, tv, tw);
  }

  // final source point (41+(7|12))
  {
    const S sx = (sx0 + sx1 + S(4.0)*sx2) / S(6.0);
    const S sy = (sy0 + sy1 + S(4.0)*sy2) / S(6.0);
    const S sz = (sz0 + sz1 + S(4.0)*sz2) / S(6.0);
    (void) kernel_0vs_0b (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz, tr,
                          tu, tv, tw);
  }
}


//
// influence of 3d linear constant-strength *vortex* panel on thick-cored target vortex with gradients
//   ignoring the 1/4pi factor, which will be multiplied later
// uses four singular integration points on the source side
//   264+(36|56) flops
//
template <class S> inline size_t flops_2v_0bg () { return 48 + 4*flops_0v_0bg<S>(); }
template <class S, class A>
static inline void kernel_2v_0bg (const S sx0, const S sy0, const S sz0,
                                  const S sx1, const S sy1, const S sz1,
                                  const S sx2, const S sy2, const S sz2,
                                  const S ssx, const S ssy, const S ssz,
                                  const S tx, const S ty, const S tz, const S tr,
                                  A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw,
                                  A* const __restrict__ tux, A* const __restrict__ tvx, A* const __restrict__ twx,
                                  A* const __restrict__ tuy, A* const __restrict__ tvy, A* const __restrict__ twy,
                                  A* const __restrict__ tuz, A* const __restrict__ tvz, A* const __restrict__ twz) {

  // scale the strength by 1/4, to account for the 4 calls below (3 flops)
  const S strx = S(0.25) * ssx;
  const S stry = S(0.25) * ssy;
  const S strz = S(0.25) * ssz;
  //const S targr = S(0.1) * tr;

  // first source point (63+(9|14) flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (54+(9|14) flops)
    (void) kernel_0v_0bg (sx, sy, sz, (S)0.0,
                          strx, stry, strz,
                          tx, ty, tz, tr,
                          tu, tv, tw,
                          tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // second source point (66+(9|14))
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0v_0bg (sx, sy, sz, (S)0.0,
                          strx, stry, strz,
                          tx, ty, tz, tr,
                          tu, tv, tw,
                          tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // third source point (66+(9|14))
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0v_0bg (sx, sy, sz, (S)0.0,
                          strx, stry, strz,
                          tx, ty, tz, tr,
                          tu, tv, tw,
                          tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // final source point (66+(9|14))
  {
    const S sx = (sx0 + sx1 + S(4.0)*sx2) / S(6.0);
    const S sy = (sy0 + sy1 + S(4.0)*sy2) / S(6.0);
    const S sz = (sz0 + sz1 + S(4.0)*sz2) / S(6.0);
    (void) kernel_0v_0bg (sx, sy, sz, (S)0.0,
                          strx, stry, strz,
                          tx, ty, tz, tr,
                          tu, tv, tw,
                          tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }
}

// same thing for vortex+source strengths
//   365+(36|56) flops
template <class S> inline size_t flops_2vs_0bg () { return 49 + 4*flops_0vs_0bg<S>(); }
template <class S, class A>
static inline void kernel_2vs_0bg (const S sx0, const S sy0, const S sz0,
                                   const S sx1, const S sy1, const S sz1,
                                   const S sx2, const S sy2, const S sz2,
                                   const S ssx, const S ssy, const S ssz, const S ss,
                                   const S tx, const S ty, const S tz, const S tr,
                                   A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw,
                                   A* const __restrict__ tux, A* const __restrict__ tvx, A* const __restrict__ twx,
                                   A* const __restrict__ tuy, A* const __restrict__ tvy, A* const __restrict__ twy,
                                   A* const __restrict__ tuz, A* const __restrict__ tvz, A* const __restrict__ twz) {

  // scale the strength by 1/4, to account for the 4 calls below (4 flops)
  const S strx = S(0.25) * ssx;
  const S stry = S(0.25) * ssy;
  const S strz = S(0.25) * ssz;
  const S strs = S(0.25) * ss;
  //const S targr = S(0.1) * tr;

  // first source point (88+(9|14) flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (79+(9|14) flops)
    (void) kernel_0vs_0bg (sx, sy, sz, (S)0.0,
                           strx, stry, strz, strs,
                           tx, ty, tz, tr,
                           tu, tv, tw,
                           tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // second source point (91+(9|14))
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0vs_0bg (sx, sy, sz, (S)0.0,
                           strx, stry, strz, strs,
                           tx, ty, tz, tr,
                           tu, tv, tw,
                           tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // third source point (91+(9|14))
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0vs_0bg (sx, sy, sz, (S)0.0,
                           strx, stry, strz, strs,
                           tx, ty, tz, tr,
                           tu, tv, tw,
                           tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // final source point (91+(9|14))
  {
    const S sx = (sx0 + sx1 + S(4.0)*sx2) / S(6.0);
    const S sy = (sy0 + sy1 + S(4.0)*sy2) / S(6.0);
    const S sz = (sz0 + sz1 + S(4.0)*sz2) / S(6.0);
    (void) kernel_0vs_0bg (sx, sy, sz, (S)0.0,
                           strx, stry, strz, strs,
                           tx, ty, tz, tr,
                           tu, tv, tw,
                           tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }
}


//
// influence of 3d linear constant-strength *vortex* panel on singular target point with gradients
//   ignoring the 1/4pi factor, which will be multiplied later
// uses four singular integration points on the source side
//   264+(28|44) flops
//
template <class S> inline size_t flops_2v_0pg () { return 48 + 4*flops_0v_0pg<S>(); }
template <class S, class A>
static inline void kernel_2v_0pg (const S sx0, const S sy0, const S sz0,
                                  const S sx1, const S sy1, const S sz1,
                                  const S sx2, const S sy2, const S sz2,
                                  const S ssx, const S ssy, const S ssz,
                                  const S tx, const S ty, const S tz,
                                  A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw,
                                  A* const __restrict__ tux, A* const __restrict__ tvx, A* const __restrict__ twx,
                                  A* const __restrict__ tuy, A* const __restrict__ tvy, A* const __restrict__ twy,
                                  A* const __restrict__ tuz, A* const __restrict__ tvz, A* const __restrict__ twz) {

  // scale the strength by 1/4, to account for the 4 calls below (3 flops)
  const S strx = S(0.25) * ssx;
  const S stry = S(0.25) * ssy;
  const S strz = S(0.25) * ssz;

  // first source point (63+(7|11) flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (54+(7|11) flops)
    (void) kernel_0v_0pg (sx, sy, sz, (S)0.0,
                          strx, stry, strz,
                          tx, ty, tz,
                          tu, tv, tw,
                          tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // second source point (66+(7|11))
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0v_0pg (sx, sy, sz, (S)0.0,
                          strx, stry, strz,
                          tx, ty, tz,
                          tu, tv, tw,
                          tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // third source point (66+(7|11))
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0v_0pg (sx, sy, sz, (S)0.0,
                          strx, stry, strz,
                          tx, ty, tz,
                          tu, tv, tw,
                          tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // final source point (66+(7|11))
  {
    const S sx = (sx0 + sx1 + S(4.0)*sx2) / S(6.0);
    const S sy = (sy0 + sy1 + S(4.0)*sy2) / S(6.0);
    const S sz = (sz0 + sz1 + S(4.0)*sz2) / S(6.0);
    (void) kernel_0v_0pg (sx, sy, sz, (S)0.0,
                          strx, stry, strz,
                          tx, ty, tz,
                          tu, tv, tw,
                          tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }
}

// same thing for vortex+source strengths
//   365+(28|44) flops
template <class S> inline size_t flops_2vs_0pg () { return 49 + 4*flops_0vs_0pg<S>(); }
template <class S, class A>
static inline void kernel_2vs_0pg (const S sx0, const S sy0, const S sz0,
                                   const S sx1, const S sy1, const S sz1,
                                   const S sx2, const S sy2, const S sz2,
                                   const S ssx, const S ssy, const S ssz, const S ss,
                                   const S tx, const S ty, const S tz,
                                   A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw,
                                   A* const __restrict__ tux, A* const __restrict__ tvx, A* const __restrict__ twx,
                                   A* const __restrict__ tuy, A* const __restrict__ tvy, A* const __restrict__ twy,
                                   A* const __restrict__ tuz, A* const __restrict__ tvz, A* const __restrict__ twz) {

  // scale the strength by 1/4, to account for the 4 calls below (4 flops)
  const S strx = S(0.25) * ssx;
  const S stry = S(0.25) * ssy;
  const S strz = S(0.25) * ssz;
  const S strs = S(0.25) * ss;

  // first source point (88+(7|11) flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (79+(7|11) flops)
    (void) kernel_0vs_0pg (sx, sy, sz, (S)0.0,
                           strx, stry, strz, strs,
                           tx, ty, tz,
                           tu, tv, tw,
                           tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // second source point (91+(7|11))
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0vs_0pg (sx, sy, sz, (S)0.0,
                           strx, stry, strz, strs,
                           tx, ty, tz,
                           tu, tv, tw,
                           tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // third source point (91+(7|11))
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0vs_0pg (sx, sy, sz, (S)0.0,
                           strx, stry, strz, strs,
                           tx, ty, tz,
                           tu, tv, tw,
                           tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // final source point (91+(7|11))
  {
    const S sx = (sx0 + sx1 + S(4.0)*sx2) / S(6.0);
    const S sy = (sy0 + sy1 + S(4.0)*sy2) / S(6.0);
    const S sz = (sz0 + sz1 + S(4.0)*sz2) / S(6.0);
    (void) kernel_0vs_0pg (sx, sy, sz, (S)0.0,
                           strx, stry, strz, strs,
                           tx, ty, tz,
                           tu, tv, tw,
                           tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }
}

// helper functions: dist (others in CoreFunc.h)

template <class S>
static inline S my_dist(const S dx, const S dy, const S dz) {
  return my_sqrt(dx*dx + dy*dy + dz*dz);
}

#ifdef USE_VC
template <class S>
static inline bool my_well_sep(const S _dist, const S _size) {
  return Vc::all_of(_dist > _size*S(4.0));
}
template <>
inline bool my_well_sep(const float _dist, const float _size) {
  return (_dist > _size*4.0f);
}
template <>
inline bool my_well_sep(const double _dist, const double _size) {
  return (_dist > _size*4.0);
}
#else
template <class S>
static inline bool my_well_sep(const S _dist, const S _size) {
  return (_dist > _size*S(4.0));
}
#endif

#ifdef USE_VC
template <class S>
static inline int my_simdwide(const S _in) {
  //return Vc::Vector<S>::size();
  return S::size();
}
template <>
inline int my_simdwide(const float _in) {
  return 1;
}
template <>
inline int my_simdwide(const double _in) {
  return 1;
}
#else
template <class S>
static inline int my_simdwide(const S _in) {
  return 1;
}
#endif

// panel-point influence, including subpaneling
//   initial input strength is a sheet strength, all recursive calls are absolute
//   returns flops
template <class S, class A>
int rkernel_2vs_0p (const S sx0, const S sy0, const S sz0,
                    const S sx1, const S sy1, const S sz1,
                    const S sx2, const S sy2, const S sz2,
                    const S ssx, const S ssy, const S ssz, const S ss,
                    const S tx, const S ty, const S tz,
                    const S sa, const int lev, const int maxlev,
                    A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {

  // accumulate flop count
  int flops = 0;

  // convert from sheet strength into absolute strength - only once
  S strx = ssx;
  S stry = ssy;
  S strz = ssz;
  S strs = ss;
  if (lev == 0) {
    strx *= sa;
    stry *= sa;
    strz *= sa;
    strs *= sa;
    flops += 4;
  }

  // compute the sizes and distance
  const S sx = (sx0 + sx1 + sx2) / S(3.0);
  const S sy = (sy0 + sy1 + sy2) / S(3.0);
  const S sz = (sz0 + sz1 + sz2) / S(3.0);
  flops += 9;
  const S trisize = my_sqrt<S>(sa);
  const S dist = my_dist<S>(tx-sx, ty-sy, tz-sz);
  flops += 10;

  // recurse or solve?
  const bool wellseparated = my_well_sep<S>(dist, trisize);
  flops += 1;
  if (wellseparated or lev == maxlev) {

    // run just one influence calculation
    // desingularize, but only by a little bit
    //(void) kernel_0vs_0p (sx, sy, sz, S(0.5)*trisize,
    (void) kernel_0vs_0p (sx, sy, sz, S(0.0),
                          strx, stry, strz, strs,
                          tx, ty, tz,
                          tu, tv, tw);

    flops += (int)flops_0vs_0p<S>();
    flops *= my_simdwide<S>(sx);

  } else {

    // split source into 4 subpanels and run 4 calls

    // scale the strength by 1/4, to account for the 4 calls below
    strx *= S(0.25);
    stry *= S(0.25);
    strz *= S(0.25);
    strs *= S(0.25);
    // and scale the area by 1/4
    const S sca = S(0.25) * sa;
    flops += 5;

    // find the 6 nodes of the source and target triangles
    const S scx[6] = {sx0, S(0.5)*(sx0+sx1), sx1, S(0.5)*(sx0+sx2), S(0.5)*(sx1+sx2), sx2};
    const S scy[6] = {sy0, S(0.5)*(sy0+sy1), sy1, S(0.5)*(sy0+sy2), S(0.5)*(sy1+sy2), sy2};
    const S scz[6] = {sz0, S(0.5)*(sz0+sz1), sz1, S(0.5)*(sz0+sz2), S(0.5)*(sz1+sz2), sz2};
    flops += 18;
    flops *= my_simdwide<S>(sx);

    // the index pointers to the child triangles
    const int id[4][3] = {{0,1,3}, {1,2,4}, {1,4,3}, {3,4,5}};

    for (int i=0; i<4; ++i) {

      // accumulate the influence
      flops += rkernel_2vs_0p (scx[id[i][0]], scy[id[i][0]], scz[id[i][0]],
                               scx[id[i][1]], scy[id[i][1]], scz[id[i][1]],
                               scx[id[i][2]], scy[id[i][2]], scz[id[i][2]],
                               strx, stry, strz, strs,
                               tx, ty, tz,
                               sca, lev+1, maxlev,
                               tu, tv, tw);
    }
  }

  return flops;
}


// panel-point influence with gradients, including subpaneling
//   returns flops
template <class S, class A>
int rkernel_2vs_0pg (const S sx0, const S sy0, const S sz0,
                     const S sx1, const S sy1, const S sz1,
                     const S sx2, const S sy2, const S sz2,
                     const S ssx, const S ssy, const S ssz, const S ss,
                     const S tx, const S ty, const S tz,
                     const S sa, const int lev, const int maxlev,
                     A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw,
                     A* const __restrict__ tux, A* const __restrict__ tvx, A* const __restrict__ twx,
                     A* const __restrict__ tuy, A* const __restrict__ tvy, A* const __restrict__ twy,
                     A* const __restrict__ tuz, A* const __restrict__ tvz, A* const __restrict__ twz) {

  // accumulate flop count
  int flops = 0;

  // convert from sheet strength into absolute strength - only once
  S strx = ssx;
  S stry = ssy;
  S strz = ssz;
  S strs = ss;
  if (lev == 0) {
    strx *= sa;
    stry *= sa;
    strz *= sa;
    strs *= sa;
    flops += 4;
  }

  // compute the sizes and distance
  const S sx = (sx0 + sx1 + sx2) / S(3.0);
  const S sy = (sy0 + sy1 + sy2) / S(3.0);
  const S sz = (sz0 + sz1 + sz2) / S(3.0);
  flops += 9;
  const S trisize = my_sqrt<S>(sa);
  const S dist = my_dist<S>(tx-sx, ty-sy, tz-sz);
  flops += 10;

  // recurse or solve?
  const bool wellseparated = my_well_sep<S>(dist, trisize);
  flops += 1;
  if (wellseparated or lev == maxlev) {

    // run just one influence calculation
    //(void) kernel_0vs_0pg (sx, sy, sz, S(0.5)*trisize,
    (void) kernel_0vs_0pg (sx, sy, sz, S(0.0),
                           strx, stry, strz, strs,
                           tx, ty, tz,
                           tu, tv, tw,
                           tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);

    flops += (int)flops_0vs_0pg<S>();
    flops *= my_simdwide<S>(sx);

  } else {

    // split source into 4 subpanels and run 4 calls

    // scale the strength by 1/4, to account for the 4 calls below
    strx *= S(0.25);
    stry *= S(0.25);
    strz *= S(0.25);
    strs *= S(0.25);
    // and scale the area by 1/4
    const S sca = S(0.25) * sa;
    flops += 5;

    // find the 6 nodes of the source and target triangles
    const S scx[6] = {sx0, S(0.5)*(sx0+sx1), sx1, S(0.5)*(sx0+sx2), S(0.5)*(sx1+sx2), sx2};
    const S scy[6] = {sy0, S(0.5)*(sy0+sy1), sy1, S(0.5)*(sy0+sy2), S(0.5)*(sy1+sy2), sy2};
    const S scz[6] = {sz0, S(0.5)*(sz0+sz1), sz1, S(0.5)*(sz0+sz2), S(0.5)*(sz1+sz2), sz2};
    flops += 18;
    flops *= my_simdwide<S>(sx);

    // the index pointers to the child triangles
    const int id[4][3] = {{0,1,3}, {1,2,4}, {1,4,3}, {3,4,5}};

    for (int i=0; i<4; ++i) {

      // accumulate the influence
      flops += rkernel_2vs_0pg (scx[id[i][0]], scy[id[i][0]], scz[id[i][0]],
                                scx[id[i][1]], scy[id[i][1]], scz[id[i][1]],
                                scx[id[i][2]], scy[id[i][2]], scz[id[i][2]],
                                strx, stry, strz, strs,
                                tx, ty, tz,
                                sca, lev+1, maxlev,
                                tu, tv, tw,
                                tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
    }
  }

  return flops;
}


// panel-panel interaction, allowing subpaneling
//   strengths are assumed to be sheet strengths
//   return value is flops count
template <class S, class A>
int rkernel_2vs_2p (const S sx0, const S sy0, const S sz0,
                    const S sx1, const S sy1, const S sz1,
                    const S sx2, const S sy2, const S sz2,
                    const S ssx, const S ssy, const S ssz, const S ss,
                    const S tx0, const S ty0, const S tz0,
                    const S tx1, const S ty1, const S tz1,
                    const S tx2, const S ty2, const S tz2,
                    const S sa, const S ta, const int lev, const int maxlev,
                    A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {

  // accumulate flop count
  int flops = 0;

  // convert from sheet strength into absolute strength - only once
  S strx = ssx;
  S stry = ssy;
  S strz = ssz;
  S strs = ss;
  if (lev == 0) {
    strx *= sa;
    stry *= sa;
    strz *= sa;
    strs *= sa;
  }

  // compute the sizes and distance
  const S sx = (sx0 + sx1 + sx2) / S(3.0);
  const S sy = (sy0 + sy1 + sy2) / S(3.0);
  const S sz = (sz0 + sz1 + sz2) / S(3.0);
  const S tx = (tx0 + tx1 + tx2) / S(3.0);
  const S ty = (ty0 + ty1 + ty2) / S(3.0);
  const S tz = (tz0 + tz1 + tz2) / S(3.0);
  flops += 18;
  const S trisize = my_sqrt<S>(sa) + my_sqrt<S>(ta);
  const S dist = my_dist<S>(tx-sx, ty-sy, tz-sz);
  flops += 12;

  // recurse or solve?
  const bool wellseparated = my_well_sep<S>(dist, trisize);
  flops += 1;
  if (wellseparated or lev == maxlev) {

    // run just one influence calculation
    //(void) kernel_0vs_0p (sx, sy, sz, S(0.5)*trisize,
    (void) kernel_0vs_0p (sx, sy, sz, S(0.0),
                          strx, stry, strz, strs,
                          tx, ty, tz,
                          tu, tv, tw);

    flops += (int)flops_0vs_0p<S>();
    flops *= my_simdwide<S>(sx);

  } else {

    // split source and target into 4 each and run 16 calls

    // scale the strength by 1/16, to account for reduced strength and reduced target area
    strx *= S(0.0625);
    stry *= S(0.0625);
    strz *= S(0.0625);
    strs *= S(0.0625);
    // scale the area by 1/4
    const S sca = S(0.25) * sa;
    const S tca = S(0.25) * ta;
    flops += 6;

    // find the 6 nodes of the source and target triangles
    const S scx[6] = {sx0, S(0.5)*(sx0+sx1), sx1, S(0.5)*(sx0+sx2), S(0.5)*(sx1+sx2), sx2};
    const S scy[6] = {sy0, S(0.5)*(sy0+sy1), sy1, S(0.5)*(sy0+sy2), S(0.5)*(sy1+sy2), sy2};
    const S scz[6] = {sz0, S(0.5)*(sz0+sz1), sz1, S(0.5)*(sz0+sz2), S(0.5)*(sz1+sz2), sz2};
    const S tcx[6] = {tx0, S(0.5)*(tx0+tx1), tx1, S(0.5)*(tx0+tx2), S(0.5)*(tx1+tx2), tx2};
    const S tcy[6] = {ty0, S(0.5)*(ty0+ty1), ty1, S(0.5)*(ty0+ty2), S(0.5)*(ty1+ty2), ty2};
    const S tcz[6] = {tz0, S(0.5)*(tz0+tz1), tz1, S(0.5)*(tz0+tz2), S(0.5)*(tz1+tz2), tz2};
    flops += 36;
    flops *= my_simdwide<S>(sx);

    // the index pointers to the child triangles
    const int id[4][3] = {{0,1,3}, {1,2,4}, {1,4,3}, {3,4,5}};

    for (int i=0; i<4; ++i) {
      for (int j=0; j<4; ++j) {

        // accumulate the influence
        flops += rkernel_2vs_2p (scx[id[i][0]], scy[id[i][0]], scz[id[i][0]],
                                 scx[id[i][1]], scy[id[i][1]], scz[id[i][1]],
                                 scx[id[i][2]], scy[id[i][2]], scz[id[i][2]],
                                 strx, stry, strz, strs,
                                 tcx[id[j][0]], tcy[id[j][0]], tcz[id[j][0]],
                                 tcx[id[j][1]], tcy[id[j][1]], tcz[id[j][1]],
                                 tcx[id[j][2]], tcy[id[j][2]], tcz[id[j][2]],
                                 sca, tca, lev+1, maxlev,
                                 tu, tv, tw);
      }
    }
  }

  return flops;
}

