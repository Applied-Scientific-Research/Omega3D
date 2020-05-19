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
//#define USE_EXPONENTIAL_KERNEL
//#define USE_WL_KERNEL
#define USE_V2_KERNEL
//#define USE_V3_KERNEL	// not programmed

// helper functions: recip, oor1p5
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


#ifdef USE_RM_KERNEL
//
// Rosenhead-Moore velocity-only
//
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = distsq + sr*sr + tr*tr;
  return oor1p5(r2);
}
template <class S> inline size_t flops_tv_nograds () { return 7; }

template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = distsq + sr*sr;
  return oor1p5(r2);
}
template <class S> inline size_t flops_tp_nograds () { return 5; }

//
// Rosenhead-Moore with gradients
//
template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S r2 = distsq + sr*sr + tr*tr;
  *r3 = oor1p5(r2);
  *bbb = S(-3.0) * (*r3) * my_recip(r2);
}
template <class S> inline size_t flops_tv_grads () { return 9; }

template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S r2 = distsq + sr*sr;
  *r3 = oor1p5(r2);
  *bbb = S(-3.0) * (*r3) * my_recip(r2);
}
template <class S> inline size_t flops_tp_grads () { return 7; }
#endif


#ifdef USE_EXPONENTIAL_KERNEL
//
// exponential core - velocity only
//
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
template <class S> inline size_t flops_tp_nograds () { return 9; }

// non-singular targets
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
template <class S> inline size_t flops_tv_nograds () { return 12; }

//
// exponential core - with gradients
//
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
template <class S> inline size_t flops_tp_grads () { return 11; }

// non-singular targets
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
template <class S> inline size_t flops_tv_grads () { return 14; }
#endif


#ifdef USE_WL_KERNEL
//
// Winckelmans–Leonard - velocity only
//
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = sr*sr + tr*tr;
  const S d2 = distsq + r2;
  return (distsq + S(2.5)*r2) * oor2p5<S>(d2);
}
template <class S> inline size_t flops_tv_nograds () { return 10; }

template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = sr*sr;
  const S d2 = distsq + r2;
  return (distsq + S(2.5)*r2) * oor2p5<S>(d2);
}
template <class S> inline size_t flops_tp_nograds () { return 8; }

template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S r2 = sr*sr + tr*tr;
  const S d2 = distsq + r2;
  const S d2top = distsq + S(2.5)*r2;
  const S dn5 = oor2p5<S>(d2);
  *r3 = d2top * dn5;
  *bbb = S(2.0)*dn5 - S(5.0)*d2top*dn5/d2;
}
template <class S> inline size_t flops_tv_grads () { return 16; }

template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S r2 = sr*sr;
  const S d2 = distsq + r2;
  const S d2top = distsq + S(2.5)*r2;
  const S dn5 = oor2p5<S>(d2);
  *r3 = d2top * dn5;
  *bbb = S(2.0)*dn5 - S(5.0)*d2top*dn5/d2;
}
template <class S> inline size_t flops_tp_grads () { return 14; }
#endif


#ifdef USE_V2_KERNEL
//
// Vatistas n=2 - velocity only
//
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S s2 = sr*sr;
  const S t2 = tr*tr;
  const S denom = distsq*distsq + s2*s2 + t2*t2;
  return oor0p75<S>(denom);
}
template <class S> inline size_t flops_tv_nograds () { return 11; }

template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S s2 = sr*sr;
  const S denom = distsq*distsq + s2*s2;
  return oor0p75<S>(denom);
}
template <class S> inline size_t flops_tp_nograds () { return 8; }

// core functions - Vatistas n=2 with gradients
template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S s2 = sr*sr;
  const S t2 = tr*tr;
  const S denom = distsq*distsq + s2*s2 + t2*t2;
  *r3 = oor0p75<S>(denom);
  *bbb = S(-3.0) * distsq / denom;
}
template <class S> inline size_t flops_tv_grads () { return 13; }

template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S s2 = sr*sr;
  const S denom = distsq*distsq + s2*s2;
  *r3 = oor0p75<S>(denom);
  *bbb = S(-3.0) * distsq / denom;
}
template <class S> inline size_t flops_tp_grads () { return 10; }
#endif


#ifdef USE_V3_KERNEL
//
// Vatistas n=3 - velocity only
//
#endif

