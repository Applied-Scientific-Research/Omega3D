/*
 * Core.h - Useful values for diffusion for various core functions
 *
 * (c)2018 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <cmath>

enum CoreType { gaussian, compactg };

//
// non-class templated functions (on real type)
//

// this is the coefficient in the core function to guarantee unit-volume
template <class RT>
RT get_core_coefficient(const CoreType thiscore) {
  return (thiscore == CoreType::gaussian) ? (1.0/(M_PI*std::sqrt(M_PI))) : (0.23873241);
}

// this is the per-axis second moment of a unit-volume core
template <class RT>
RT get_core_second_mom(const CoreType thiscore) {
  return (thiscore == CoreType::gaussian) ? (0.5) : (0.300915092);
}

// this is the per-axis fourth moment of a unit-volume core
template <class RT>
RT get_core_fourth_mom(const CoreType thiscore) {
  return (thiscore == CoreType::gaussian) ? (0.75) : (0.238127865);
}

// this is the diffusion time from a delta function to the core function
template <class RT>
RT get_core_diffusion_time(const CoreType thiscore) {
  return (thiscore == CoreType::gaussian) ? (1.0) : (0.601830184);
}
