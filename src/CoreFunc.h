/*
 * CoreFunc.h - Non-class core function inlines for influence calculations
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
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
#define USE_WL_KERNEL
//#define USE_V2_KERNEL
//#define USE_V3_KERNEL	// not programmed

// helper functions: sqrt, recip, oor1p5
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
// a helper conditional
#ifdef USE_VC
template <class S>
static inline S exp_cond (const S ood3, const S corefac, const S reld3) {
  S returnval = ood3;
  returnval(reld3 < S(16.0)) = ood3 * (S(1.0) - Vc::exp(-reld3));
  returnval(reld3 < S(0.001)) = corefac;
  return returnval;
}
template <>
inline float exp_cond (const float ood3, const float corefac, const float reld3) {
  if (reld3 > 16.0f) {
    return ood3;
  } else if (reld3 < 0.001f) {
    return corefac;
  } else {
    return ood3 * (1.0f - std::exp(-reld3));
  }
}
template <>
inline double exp_cond (const double ood3, const double corefac, const double reld3) {
  if (reld3 > 16.0) {
    return ood3;
  } else if (reld3 < 0.001) {
    return corefac;
  } else {
    return ood3 * (1.0 - std::exp(-reld3));
  }
}
#else
template <class S>
static inline S exp_cond (const S ood3, const S corefac, const S reld3) {
  if (reld3 > 16.0f) {
    return ood3;
    // 1 flop (comparison)
  } else if (reld3 < 0.001f) {
    return corefac;
    // 2 flops
  } else {
    return ood3 * (1.0f - std::exp(-reld3));
    // 3 flops
  }
}
#endif
// and another one for bbb
#ifdef USE_VC
template <class S>
static inline S exp_bbb (const S r3, const S corefac, const S reld3, const S dist, const S distsq) {
  S mybbb;
  mybbb(reld3 > S(16.0)) = S(-3.0) * r3 / distsq;
  const S expreld3 = Vc::exp(-reld3);
  mybbb(reld3 < S(16.0)) = S(3.0) * (corefac*expreld3 - r3) / distsq;
  mybbb(reld3 < S(0.001)) = S(-1.5) * dist * r3 * r3;
  return mybbb;
}
template <>
inline float exp_bbb (const float r3, const float corefac, const float reld3, const float dist, const float distsq) {
  if (reld3 > 16.0f) {
    return -3.0f * r3 / distsq;
  } else if (reld3 < 0.001f) {
    return -1.5f * dist * r3 * r3;
  } else {
    const float expreld3 = std::exp(-reld3);
    return 3.0f * (corefac*expreld3 - r3) / distsq;
  }
}
template <>
inline double exp_bbb (const double r3, const double corefac, const double reld3, const double dist, const double distsq) {
  if (reld3 > 16.0) {
    return -3.0 * r3 / distsq;
  } else if (reld3 < 0.001) {
    return -1.5 * dist * r3 * r3;
  } else {
    const double expreld3 = std::exp(-reld3);
    return 3.0 * (corefac*expreld3 - r3) / distsq;
  }
}
#else
template <class S>
static inline S exp_bbb (const S r3, const S corefac, const S reld3, const S dist, const S distsq) {
  if (reld3 > 16.0f) {
    return -3.0f * r3 / distsq;
    // this is 2 flops and is most likely
  } else if (reld3 < 0.001f) {
    return -1.5f * dist * r3 * r3;
    // this is 3 flops
  } else {
    const S expreld3 = std::exp(-reld3);
    return 3.0f * (corefac*expreld3 - r3) / distsq;
    // this is 5 flops
  }
}
#endif
//
// exponential core - velocity only
//
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S dist = my_sqrt(distsq);
  const S ood3 = my_recip(distsq*dist);
  const S corefac = my_recip(sr*sr*sr);
  const S reld3 = corefac / ood3;
  // 7 flops to here
  return exp_cond(ood3, corefac, reld3);
}
template <class S> inline size_t flops_tp_nograds () { return 9; }

// non-singular targets
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S dist = my_sqrt(distsq);
  const S ood3 = my_recip(distsq*dist);
  const S corefac = my_recip(sr*sr*sr + tr*tr*tr);
  const S reld3 = corefac / ood3;
  return exp_cond(ood3, corefac, reld3);
}
template <class S> inline size_t flops_tv_nograds () { return 12; }

//
// exponential core - with gradients
//
template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S dist = my_sqrt(distsq);
  const S corefac = my_recip(sr*sr*sr);
  const S d3 = distsq * dist;
  const S reld3 = d3 * corefac;
  // 6 flops to here
  const S ood3 = my_recip(d3);
  *r3 = exp_cond(ood3, corefac, reld3);
  *bbb = exp_bbb(*r3, corefac, reld3, dist, distsq);
}
template <class S> inline size_t flops_tp_grads () { return 11; }

// non-singular targets
template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S dist = my_sqrt(distsq);
  const S corefac = my_recip(sr*sr*sr + tr*tr*tr);
  const S d3 = distsq * dist;
  const S reld3 = d3 * corefac;
  // 9 flops to here
  const S ood3 = my_recip(d3);
  *r3 = exp_cond(ood3, corefac, reld3);
  *bbb = exp_bbb(*r3, corefac, reld3, dist, distsq);
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

//
// Winckelmans–Leonard - with gradients
//
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

//
// Vatistas n=2 - with gradients
//
template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S s2 = sr*sr;
  const S t2 = tr*tr;
  const S denom = distsq*distsq + s2*s2 + t2*t2;
  *r3 = oor0p75<S>(denom);
  // this does not seem right
  //*bbb = S(-3.0) * distsq / denom;
  // this matches the other kernels better
  const S denomg = distsq*distsq*my_sqrt(distsq) + s2*s2 + t2*t2;
  *bbb = S(-3.0) / denom;
  // but still returns elongation half of what it should be
  // gnuplot this:
  // set logscale y
  // plot [0:10][1e-5:1000] 3*(x*x+1)**-2.5, -2*(x*x+1)**-2.5+5*(x*x+2.5)*(x*x+1)**-3.5, 3*x**-5, 3/(1+x**5)
}
template <class S> inline size_t flops_tv_grads () { return 13; }

template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S s2 = sr*sr;
  const S denom = distsq*distsq + s2*s2;
  *r3 = oor0p75<S>(denom);
  // see above
  //*bbb = S(-3.0) * distsq / denom;
  const S denomg = distsq*distsq*my_sqrt(distsq) + s2*s2;
  *bbb = S(-3.0) / denom;
}
template <class S> inline size_t flops_tp_grads () { return 10; }
#endif


#ifdef USE_V3_KERNEL
//
// Vatistas n=3 - velocity only
//
#endif

