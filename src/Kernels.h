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
template <class S, class A>
static inline void kernel_0v_0b (const S sx, const S sy, const S sz,
                                 const S sr,
                                 const S ssx, const S ssy, const S ssz,
                                 const S tx, const S ty, const S tz,
                                 const S tr,
                                 A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {
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

// same, but vortex+source strength
template <class S, class A>
static inline void kernel_0vs_0b (const S sx, const S sy, const S sz,
                                 const S sr,
                                 const S ssx, const S ssy, const S ssz, const S ss,
                                 const S tx, const S ty, const S tz,
                                 const S tr,
                                 A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {
  // 36 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr + tr*tr;
#ifdef USE_VC
  r2 = Vc::reciprocal(r2*Vc::sqrt(r2));
#else
  r2 = 1.0 / (r2*std::sqrt(r2));
#endif
  const S dxxw = dz*ssy - dy*ssz + dx*ss;
  const S dyxw = dx*ssz - dz*ssx + dy*ss;
  const S dzxw = dy*ssx - dx*ssy + dz*ss;
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
                                 A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {
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

// same, but vortex+source strength
template <class S, class A>
static inline void kernel_0vs_0p (const S sx, const S sy, const S sz,
                                  const S sr,
                                  const S ssx, const S ssy, const S ssz, const S ss,
                                  const S tx, const S ty, const S tz,
                                  A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {
  // 34 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr;
#ifdef USE_VC
  r2 = Vc::reciprocal(r2*Vc::sqrt(r2));
#else
  r2 = 1.0 / (r2*std::sqrt(r2));
#endif
  const S dxxw = dz*ssy - dy*ssz + dx*ss;
  const S dyxw = dx*ssz - dz*ssx + dy*ss;
  const S dzxw = dy*ssx - dx*ssy + dz*ss;
  *tu += r2 * dxxw;
  *tv += r2 * dyxw;
  *tw += r2 * dzxw;
}

// same, but source strength only
template <class S, class A>
static inline void kernel_0s_0p (const S sx, const S sy, const S sz,
                                 const S sr,
                                 const S ss,
                                 const S tx, const S ty, const S tz,
                                 A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {
  // 19 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr;
#ifdef USE_VC
  r2 = (S)ss * Vc::reciprocal(r2*Vc::sqrt(r2));
#else
  r2 = ss / (r2*std::sqrt(r2));
#endif
  *tu += r2 * dx;
  *tv += r2 * dy;
  *tw += r2 * dz;
}

// thick-cored particle on thick-cored point, with gradients - 65 flops total
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

// same, but for vortex+source strengths
//   90 flops total
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
  const S r2 = dx*dx + dy*dy + dz*dz + sr*sr + tr*tr;
#ifdef USE_VC
  const S r3 = Vc::reciprocal(r2*Vc::sqrt(r2));
#else
  const S r3 = 1.0 / (r2*std::sqrt(r2));
#endif
  S dxxw = dz*ssy - dy*ssz;
  S dyxw = dx*ssz - dz*ssx;
  S dzxw = dy*ssx - dx*ssy;
  *tu += r3 * (dxxw + dx*ss);
  *tv += r3 * (dyxw + dy*ss);
  *tw += r3 * (dzxw + dz*ss);

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
//   63 flops total
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

// same, but for vortex+source strengths
//   88 flops total
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
  const S r2 = dx*dx + dy*dy + dz*dz + sr*sr;
#ifdef USE_VC
  const S r3 = Vc::reciprocal(r2*Vc::sqrt(r2));
#else
  const S r3 = 1.0 / (r2*std::sqrt(r2));
#endif
  S dxxw = dz*ssy - dy*ssz;
  S dyxw = dx*ssz - dz*ssx;
  S dzxw = dy*ssx - dx*ssy;
  *tu += r3 * (dxxw + dx*ss);
  *tv += r3 * (dyxw + dy*ss);
  *tw += r3 * (dzxw + dz*ss);

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
//   160 flops
//
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

  // first source point (37 flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (28 flops)
    (void) kernel_0v_0p (sx, sy, sz, (S)0.0,
                        strx, stry, strz,
                        tx, ty, tz,
                        tu, tv, tw);
  }

  // second source point (40)
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0v_0p (sx, sy, sz, (S)0.0,
                        strx, stry, strz,
                        tx, ty, tz,
                        tu, tv, tw);
  }

  // third source point (40)
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0v_0p (sx, sy, sz, (S)0.0,
                        strx, stry, strz,
                        tx, ty, tz,
                        tu, tv, tw);
  }

  // final source point (40)
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
//   185 flops
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

  // first source point (43 flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (34 flops)
    (void) kernel_0vs_0p (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz,
                          tu, tv, tw);
  }

  // second source point (46)
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0vs_0p (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz,
                          tu, tv, tw);
  }

  // third source point (46)
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0vs_0p (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz,
                          tu, tv, tw);
  }

  // final source point (46)
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
//   122 flops
template <class S, class A>
static inline void kernel_2s_0p (const S sx0, const S sy0, const S sz0,
                                 const S sx1, const S sy1, const S sz1,
                                 const S sx2, const S sy2, const S sz2,
                                 const S ss,
                                 const S tx, const S ty, const S tz,
                                 A* const __restrict__ tu, A* const __restrict__ tv, A* const __restrict__ tw) {

  // scale the strength by 1/4, to account for the 4 calls below (1 flop)
  const S strs = S(0.25) * ss;

  // first source point (28 flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (19 flops)
    (void) kernel_0s_0p (sx, sy, sz, (S)0.0,
                         strs,
                         tx, ty, tz,
                         tu, tv, tw);
  }

  // second source point (31)
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0s_0p (sx, sy, sz, (S)0.0,
                         strs,
                         tx, ty, tz,
                         tu, tv, tw);
  }

  // third source point (31)
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0s_0p (sx, sy, sz, (S)0.0,
                         strs,
                         tx, ty, tz,
                         tu, tv, tw);
  }

  // final source point (31)
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
//   168 flops
//
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

  // first source point (39 flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (30 flops)
    (void) kernel_0v_0b (sx, sy, sz, (S)0.0,
                         strx, stry, strz,
                         tx, ty, tz, tr,
                         tu, tv, tw);
  }

  // second source point (42)
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0v_0b (sx, sy, sz, (S)0.0,
                         strx, stry, strz,
                         tx, ty, tz, tr,
                         tu, tv, tw);
  }

  // third source point (42)
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0v_0b (sx, sy, sz, (S)0.0,
                         strx, stry, strz,
                         tx, ty, tz, tr,
                         tu, tv, tw);
  }

  // final source point (42)
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
//   193 flops
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

  // first source point (45 flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (36 flops)
    (void) kernel_0vs_0b (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz, tr,
                          tu, tv, tw);
  }

  // second source point (48)
  {
    const S sx = (S(4.0)*sx0 + sx1 + sx2) / S(6.0);
    const S sy = (S(4.0)*sy0 + sy1 + sy2) / S(6.0);
    const S sz = (S(4.0)*sz0 + sz1 + sz2) / S(6.0);
    (void) kernel_0vs_0b (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz, tr,
                          tu, tv, tw);
  }

  // third source point (48)
  {
    const S sx = (sx0 + S(4.0)*sx1 + sx2) / S(6.0);
    const S sy = (sy0 + S(4.0)*sy1 + sy2) / S(6.0);
    const S sz = (sz0 + S(4.0)*sz1 + sz2) / S(6.0);
    (void) kernel_0vs_0b (sx, sy, sz, (S)0.0,
                          strx, stry, strz, strs,
                          tx, ty, tz, tr,
                          tu, tv, tw);
  }

  // final source point (48)
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
//   308 flops
//
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

  // first source point (74 flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (65 flops)
    (void) kernel_0v_0bg (sx, sy, sz, (S)0.0,
                          strx, stry, strz,
                          tx, ty, tz, tr,
                          tu, tv, tw,
                          tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // second source point (77)
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

  // third source point (77)
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

  // final source point (77)
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
//   409 flops
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

  // first source point (99 flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (90 flops)
    (void) kernel_0vs_0bg (sx, sy, sz, (S)0.0,
                           strx, stry, strz, strs,
                           tx, ty, tz, tr,
                           tu, tv, tw,
                           tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // second source point (102)
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

  // third source point (102)
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

  // final source point (102)
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
//   300 flops
//
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

  // first source point (72 flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (63 flops)
    (void) kernel_0v_0pg (sx, sy, sz, (S)0.0,
                          strx, stry, strz,
                          tx, ty, tz,
                          tu, tv, tw,
                          tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // second source point (75)
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

  // third source point (75)
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

  // final source point (75)
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
//   401 flops
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

  // first source point (97 flops)
  {
    // prepare the source points (9 flops)
    const S sx = (sx0 + sx1 + sx2) / S(3.0);
    const S sy = (sy0 + sy1 + sy2) / S(3.0);
    const S sz = (sz0 + sz1 + sz2) / S(3.0);

    // accumulate the influence (88 flops)
    (void) kernel_0vs_0pg (sx, sy, sz, (S)0.0,
                           strx, stry, strz, strs,
                           tx, ty, tz,
                           tu, tv, tw,
                           tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);
  }

  // second source point (100)
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

  // third source point (100)
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

  // final source point (100)
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

// panel-point influence, including subpaneling
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

  // compute the sizes and distance
  const S sx = (sx0 + sx1 + sx2) / S(3.0);
  const S sy = (sy0 + sy1 + sy2) / S(3.0);
  const S sz = (sz0 + sz1 + sz2) / S(3.0);
  flops += 9;
#ifdef USE_VC
  const S trisize = Vc::sqrt(sa);
  const S dist = Vc::sqrt(Vc::pow(tx-sx,2) + Vc::pow(ty-sy,2) + Vc::pow(tz-sz,2));
#else
  const S trisize = std::sqrt(sa);
  const S dist = std::sqrt(std::pow(tx-sx,2) + std::pow(ty-sy,2) + std::pow(tz-sz,2));
#endif
  flops += 10;

  // recurse or solve?
#ifdef USE_VC
  const bool wellseparated = Vc::all_of(dist > trisize*S(4.0));
#else
  const bool wellseparated = dist > trisize*S(4.0);
#endif
  flops += 1;
  if (wellseparated or lev == maxlev) {

    // run just one influence calculation
    // desingularize, but only by a little bit
    //(void) kernel_0vs_0p (sx, sy, sz, S(0.25)*trisize,
    (void) kernel_0vs_0p (sx, sy, sz, S(0.0),
                          ssx, ssy, ssz, ss,
                          tx, ty, tz,
                          tu, tv, tw);

    flops += 34;

#ifdef USE_VC
    flops *= 8;
#endif

  } else {

    // split source into 4 subpanels and run 4 calls

    // scale the strength by 1/4, to account for the 4 calls below
    const S strx = S(0.25) * ssx;
    const S stry = S(0.25) * ssy;
    const S strz = S(0.25) * ssz;
    const S strs = S(0.25) * ss;
    const S sca = S(0.25) * sa;
    flops += 5;

    // find the 6 nodes of the source and target triangles
    const S scx[6] = {sx0, S(0.5)*(sx0+sx1), sx1, S(0.5)*(sx0+sx2), S(0.5)*(sx1+sx2), sx2};
    const S scy[6] = {sy0, S(0.5)*(sy0+sy1), sy1, S(0.5)*(sy0+sy2), S(0.5)*(sy1+sy2), sy2};
    const S scz[6] = {sz0, S(0.5)*(sz0+sz1), sz1, S(0.5)*(sz0+sz2), S(0.5)*(sz1+sz2), sz2};
    flops += 18;
#ifdef USE_VC
    flops *= 8;
#endif

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

  // compute the sizes and distance
  const S sx = (sx0 + sx1 + sx2) / S(3.0);
  const S sy = (sy0 + sy1 + sy2) / S(3.0);
  const S sz = (sz0 + sz1 + sz2) / S(3.0);
  flops += 9;
#ifdef USE_VC
  const S trisize = Vc::sqrt(sa);
  const S dist = Vc::sqrt(Vc::pow(tx-sx,2) + Vc::pow(ty-sy,2) + Vc::pow(tz-sz,2));
#else
  const S trisize = std::sqrt(sa);
  const S dist = std::sqrt(std::pow(tx-sx,2) + std::pow(ty-sy,2) + std::pow(tz-sz,2));
#endif
  flops += 10;

  // recurse or solve?
  flops += 1;
  if (dist > trisize*S(4.0) or lev == maxlev) {

    // run just one influence calculation
    (void) kernel_0vs_0pg (sx, sy, sz, S(0.0),
                           ssx, ssy, ssz, ss,
                           tx, ty, tz,
                           tu, tv, tw,
                           tux, tvx, twx, tuy, tvy, twy, tuz, tvz, twz);

    flops += 100;

#ifdef USE_VC
    flops *= 8;
#endif

  } else {

    // split source into 4 subpanels and run 4 calls

    // scale the strength by 1/4, to account for the 4 calls below
    const S strx = S(0.25) * ssx;
    const S stry = S(0.25) * ssy;
    const S strz = S(0.25) * ssz;
    const S strs = S(0.25) * ss;
    const S sca = S(0.25) * sa;
    flops += 5;

    // find the 6 nodes of the source and target triangles
    const S scx[6] = {sx0, S(0.5)*(sx0+sx1), sx1, S(0.5)*(sx0+sx2), S(0.5)*(sx1+sx2), sx2};
    const S scy[6] = {sy0, S(0.5)*(sy0+sy1), sy1, S(0.5)*(sy0+sy2), S(0.5)*(sy1+sy2), sy2};
    const S scz[6] = {sz0, S(0.5)*(sz0+sz1), sz1, S(0.5)*(sz0+sz2), S(0.5)*(sz1+sz2), sz2};
    flops += 18;
#ifdef USE_VC
    flops *= 8;
#endif

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

  // compute the sizes and distance
  const S sx = (sx0 + sx1 + sx2) / S(3.0);
  const S sy = (sy0 + sy1 + sy2) / S(3.0);
  const S sz = (sz0 + sz1 + sz2) / S(3.0);
  const S tx = (tx0 + tx1 + tx2) / S(3.0);
  const S ty = (ty0 + ty1 + ty2) / S(3.0);
  const S tz = (tz0 + tz1 + tz2) / S(3.0);
  flops += 18;
#ifdef USE_VC
  const S trisize = Vc::sqrt(sa) + Vc::sqrt(ta);
  const S dist = Vc::sqrt(Vc::pow(tx-sx,2) + Vc::pow(ty-sy,2) + Vc::pow(tz-sz,2));
#else
  const S trisize = std::sqrt(sa) + std::sqrt(ta);
  const S dist = std::sqrt(std::pow(tx-sx,2) + std::pow(ty-sy,2) + std::pow(tz-sz,2));
#endif
  flops += 12;

  // recurse or solve?
  flops += 1;
  if (dist > trisize*S(4.0) or lev == maxlev) {

    // run just one influence calculation
    (void) kernel_0vs_0p (sx, sy, sz, S(0.0),
                          ssx, ssy, ssz, ss,
                          tx, ty, tz,
                          tu, tv, tw);

    flops += 34;

#ifdef USE_VC
    flops *= 8;
#endif

  } else {

    // split source and target into 4 each and run 16 calls

    // scale the strength by 1/4, to account for the 4 calls below
    const S strx = S(0.0625) * ssx;
    const S stry = S(0.0625) * ssy;
    const S strz = S(0.0625) * ssz;
    const S strs = S(0.0625) * ss;
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
#ifdef USE_VC
    flops *= 8;
#endif

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

