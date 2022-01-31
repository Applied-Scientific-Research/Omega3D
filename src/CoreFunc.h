/*
 * CoreFunc.h - Non-class core function inlines for influence calculations
 *
 * (c)2020,2 Applied Scientific Research, Inc.
 *           Mark J Stock <markjstock@gmail.com>
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

#include "MathHelper.h"

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <cmath>

//#define USE_RM_KERNEL
//#define USE_EXPONENTIAL_KERNEL
//#define USE_WL_KERNEL
#define USE_V2_KERNEL
//#define USE_V3_KERNEL	// not programmed



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

// and the one for singular targets
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

// and the one for singular targets
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
// Winckelmans-Leonard - velocity only
//
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = sr*sr + tr*tr;
  const S d2 = distsq + r2;
  return (distsq + S(2.5)*r2) * oor2p5<S>(d2);
}
template <class S> inline size_t flops_tv_nograds () { return 10; }

// and the one for singular targets
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = sr*sr;
  const S d2 = distsq + r2;
  return (distsq + S(2.5)*r2) * oor2p5<S>(d2);
}
template <class S> inline size_t flops_tp_nograds () { return 8; }

//
// Winckelmans-Leonard - with gradients
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

// and the one for singular targets
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

// and for singular targets
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
  *bbb = S(-3.0) * distsq * (*r3) * my_recip(denom);
}
template <class S> inline size_t flops_tv_grads () { return 13; }

// and the one for singular targets
template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r3, S* const __restrict__ bbb) {
  const S s2 = sr*sr;
  const S denom = distsq*distsq + s2*s2;
  *r3 = oor0p75<S>(denom);
  *bbb = S(-3.0) * distsq * (*r3) * my_recip(denom);
}
template <class S> inline size_t flops_tp_grads () { return 10; }
#endif


#ifdef USE_V3_KERNEL
//
// Vatistas n=3 - velocity only
//
#endif

