/*
 * CoreFunc.h - Non-class core function inlines for influence calculations
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#ifdef _WIN32
#define __restrict__ __restrict
#endif

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <cmath>

//#define USE_RM_KERNEL
#define USE_EXPONENTIAL_KERNEL
//#define USE_WL_KERNEL
//#define USE_V2_KERNEL


#ifdef USE_RM_KERNEL
//
// core functions - Rosenhead-Moore
//

// this is always 7 flops
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = distsq + sr*sr + tr*tr;
#ifdef USE_VC
  return Vc::reciprocal(r2*Vc::sqrt(r2));
#else
  return S(1.0) / (r2*std::sqrt(r2));
#endif
}

// this is always 5 flops
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = distsq + sr*sr;
#ifdef USE_VC
  return Vc::reciprocal(r2*Vc::sqrt(r2));
#else
  return S(1.0) / (r2*std::sqrt(r2));
#endif
}

// core functions - Rosenhead-Moore with gradients

// this is always 9 flops
template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S r2 = distsq + sr*sr + tr*tr;
#ifdef USE_VC
  *r3 = Vc::reciprocal(r2*Vc::sqrt(r2));
  *bbb = S(-3.0) * (*r3) * Vc::reciprocal(r2);
#else
  *r3 = S(1.0) / (r2*std::sqrt(r2));
  *bbb = S(-3.0) * (*r3) / r2;
#endif
}

// this is always 7 flops
template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S r2 = distsq + sr*sr;
#ifdef USE_VC
  *r3 = Vc::reciprocal(r2*Vc::sqrt(r2));
  *bbb = S(-3.0) * (*r3) * Vc::reciprocal(r2);
#else
  *r3 = S(1.0) / (r2*std::sqrt(r2));
  *bbb = S(-3.0) * (*r3) / r2;
#endif
}
#endif


#ifdef USE_EXPONENTIAL_KERNEL
//
// core functions - compact exponential
//

// this probably averages out to 12 flops
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
#ifdef USE_VC
  const S dist = Vc::sqrt(distsq);
  const S ood3 = Vc::reciprocal(distsq * dist);
  const S corefac = Vc::reciprocal(sr*sr*sr + tr*tr*tr);
#else
  const S dist = std::sqrt(distsq);
  const S ood3 = S(1.0) / (distsq * dist);
  const S corefac = S(1.0) / (std::pow(sr,3) + std::pow(tr,3));
#endif
  const S reld3 = corefac / ood3;
#ifdef USE_VC
  S returnval = ood3;
  returnval(reld3 < S(16.0)) = ood3 * (S(1.0) - Vc::exp(-reld3));
  returnval(reld3 < S(0.001)) = corefac;
  return returnval;
#else
  if (reld3 > S(16.0)) {
    return ood3;
  } else if (reld3 < S(0.001)) {
    return corefac;
  } else {
    return ood3 * (S(1.0) - std::exp(-reld3));
  }
#endif
}

// this probably averages out to 9 flops
template <class S>
static inline S core_func (const S distsq, const S sr) {
#ifdef USE_VC
  const S dist = Vc::sqrt(distsq);
  const S ood3 = Vc::reciprocal(distsq * dist);
  const S corefac = Vc::reciprocal(sr*sr*sr);
#else
  const S dist = std::sqrt(distsq);
  const S ood3 = S(1.0) / (distsq * dist);
  const S corefac = S(1.0) / std::pow(sr,3);
#endif
  const S reld3 = corefac / ood3;
  // 7 flops to here
#ifdef USE_VC
  S returnval = ood3;
  returnval(reld3 < S(16.0)) = ood3 * (S(1.0) - Vc::exp(-reld3));
  returnval(reld3 < S(0.001)) = corefac;
  //std::cout << std::endl << dist << std::endl << reld3 << std::endl << returnval << std::endl;
  return returnval;
#else
  if (reld3 > S(16.0)) {
    return ood3;
    // 1 flop (comparison)
  } else if (reld3 < S(0.001)) {
    return corefac;
    // 2 flops
  } else {
    return ood3 * (S(1.0) - std::exp(-reld3));
    // 3 flops
  }
#endif
}

// core functions - compact exponential with gradients

// call this one 14 flops average
template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
#ifdef USE_VC
  const S dist = Vc::sqrt(distsq);
  const S corefac = Vc::reciprocal(sr*sr*sr + tr*tr*tr);
#else
  const S dist = std::sqrt(distsq);
  const S corefac = S(1.0) / (std::pow(sr,3) + std::pow(tr,3));
#endif
  const S d3 = distsq * dist;
  const S reld3 = d3 * corefac;
  // 9 flops to here
#ifdef USE_VC
  S myr3, mybbb;
  myr3(reld3 > S(16.0)) = Vc::reciprocal(d3);
  mybbb(reld3 > S(16.0)) = S(-3.0) / (d3 * distsq);
  const S expreld3 = Vc::exp(-reld3);
  myr3(reld3 < S(16.0)) = (S(1.0) - expreld3) / d3;
  mybbb(reld3 < S(16.0)) = S(3.0) * (corefac*expreld3 - myr3) / distsq;
  myr3(reld3 < S(0.001)) = corefac;
  mybbb(reld3 < S(0.001)) = S(-1.5) * dist * corefac * corefac;
  *r3 = myr3;
  *bbb = mybbb;
#else
  if (reld3 > S(16.0)) {
    *r3 = S(1.0) / d3;
    *bbb = S(-3.0) / (d3 * distsq);
    // this is 4 flops and is most likely
  } else if (reld3 < S(0.001)) {
    *r3 = corefac;
    *bbb = S(-1.5) * dist * corefac * corefac;
    // this is 5 flops
  } else {
    const S expreld3 = std::exp(-reld3);
    *r3 = (S(1.0) - expreld3) / d3;
    *bbb = S(3.0) * (corefac*expreld3 - *r3) / distsq;
    // this is 9 flops
  }
#endif
}

// call this one 11 flops average
template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
#ifdef USE_VC
  const S dist = Vc::sqrt(distsq);
  const S corefac = Vc::reciprocal(sr*sr*sr);
#else
  const S dist = std::sqrt(distsq);
  const S corefac = S(1.0) / std::pow(sr,3);
#endif
  const S d3 = distsq * dist;
  const S reld3 = d3 * corefac;
  // 6 flops to here
#ifdef USE_VC
  S myr3, mybbb;
  myr3(reld3 > S(16.0)) = Vc::reciprocal(d3);
  mybbb(reld3 > S(16.0)) = S(-3.0) / (d3 * distsq);
  const S expreld3 = Vc::exp(-reld3);
  myr3(reld3 < S(16.0)) = (S(1.0) - expreld3) / d3;
  mybbb(reld3 < S(16.0)) = S(3.0) * (corefac*expreld3 - myr3) / distsq;
  myr3(reld3 < S(0.001)) = corefac;
  mybbb(reld3 < S(0.001)) = S(-1.5) * dist * corefac * corefac;
  *r3 = myr3;
  *bbb = mybbb;
#else
  if (reld3 > S(16.0)) {
    *r3 = S(1.0) / d3;
    *bbb = S(-3.0) / (d3 * distsq);
    // this is 4 flops and is most likely
  } else if (reld3 < S(0.001)) {
    *r3 = corefac;
    *bbb = S(-1.5) * dist * corefac * corefac;
    // this is 5 flops
  } else {
    const S expreld3 = std::exp(-reld3);
    *r3 = (S(1.0) - expreld3) / d3;
    *bbb = S(3.0) * (corefac*expreld3 - *r3) / distsq;
    // this is 9 flops
  }
#endif
}
#endif


#ifdef USE_WL_KERNEL
//
// core functions - Winckelmans–Leonard
//

// this is always 10 flops
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = sr*sr + tr*tr;
  const S d2 = distsq + r2;
#ifdef USE_VC
  return (distsq + S(2.5)*r2) / (d2*d2*Vc::sqrt(d2));
#else
  return (distsq + S(2.5)*r2) / (d2*d2*std::sqrt(d2));
#endif
}

// this is always 8 flops
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = sr*sr;
  const S d2 = distsq + r2;
#ifdef USE_VC
  return (distsq + S(2.5)*r2) / (d2*d2*Vc::sqrt(d2));
#else
  return (distsq + S(2.5)*r2) / (d2*d2*std::sqrt(d2));
#endif
}

template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
}

template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
}
#endif

