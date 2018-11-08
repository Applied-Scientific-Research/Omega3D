#pragma once

#include "Kernels.h"
#include "Points.h"
#include "Panels.h"

#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#define _USE_MATH_DEFINES
#include <cmath>


template <class S, class A>
void points_affect_points (Points<S> const& src, Points<S>& targ) {
  std::cout << "    0_0 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // get references to use locally
  const std::vector<S>& sx = src.get_pos();
  const std::vector<S>& sr = src.get_rad();
  const std::vector<S>& ss = src.get_str();
  const std::vector<S>& tx = targ.get_pos();
  const std::vector<S>& tr = targ.get_rad();
  std::vector<S>& tu = targ.get_vel();
  auto& opttug = targ.get_velgrad();

  // here is where we can dispatch on solver type, grads-or-not, core function, etc.

  // dispatch on presence of val grads
  if (opttug) {
    //std::vector<S>& tug = targ.get_vel();
    // get the pointer from the optional
    auto tug = *opttug;

    // velocity+grads kernel
    for (size_t i=0; i<targ.getn(); ++i) {
      std::array<A,12> accum = {0.0};
      for (size_t j=0; j<src.getn(); ++j) {
        kernel_0_0g<S,A>(&sx[3*j], sr[j], &ss[3*j],
                         &tx[3*i], tr[i], accum.data());
      }
      tu[3*i+0] += accum[0];
      tu[3*i+1] += accum[1];
      tu[3*i+2] += accum[2];
      tug[9*i+0] += accum[3];
      tug[9*i+1] += accum[4];
      tug[9*i+2] += accum[5];
      tug[9*i+3] += accum[6];
      tug[9*i+4] += accum[7];
      tug[9*i+5] += accum[8];
      tug[9*i+6] += accum[9];
      tug[9*i+7] += accum[10];
      tug[9*i+8] += accum[11];
    }

  } else {
    // velocity-only kernel
    for (size_t i=0; i<targ.getn(); ++i) {
      std::array<A,3> accum = {0.0};
      for (size_t j=0; j<src.getn(); ++j) {
        kernel_0_0<S,A>(&sx[3*j], sr[j], &ss[3*j],
                        &tx[3*i], tr[i], accum.data());
      }
      tu[3*i+0] += accum[0];
      tu[3*i+1] += accum[1];
      tu[3*i+2] += accum[2];
    }
  }
}


template <class S, class A>
void panels_affect_points (Panels<S> const& src, Points<S>& targ) {
  std::cout << "    1_0 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // get references to use locally
  const std::vector<S>& sx = src.get_pos();
  //const std::vector<S>& sr = src.get_rad();
  const std::vector<uint16_t>& si = src.get_idx();
  const std::vector<S>& ss = src.get_str();
  const std::vector<S>& tx = targ.get_pos();
  //const std::vector<S>& tr = targ.get_rad();
  std::vector<S>& tu = targ.get_vel();

  for (size_t i=0; i<targ.getn(); ++i) {
    std::array<A,2> accum = {0.0};
    for (size_t j=0; j<src.getn(); ++j) {
      kernel_1_0<S,A>(&sx[2*si[2*j]], &sx[2*si[2*j+1]], ss[j],
                      &tx[2*i], accum.data());
    }
    tu[2*i]   += accum[0];
    tu[2*i+1] += accum[1];
  }
}


template <class S, class A>
void points_affect_panels (Points<S> const& src, Panels<S>& targ) {
  std::cout << "    0_1 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // get references to use locally
  const std::vector<S>& sx = src.get_pos();
  const std::vector<S>& ss = src.get_str();
  const std::vector<S>& tx = targ.get_pos();
  const std::vector<uint16_t>& ti = targ.get_idx();
  std::vector<S>& tu = targ.get_vel();

  for (size_t i=0; i<targ.getn(); ++i) {
    std::array<A,2> accum = {0.0};
    for (size_t j=0; j<src.getn(); ++j) {
      // note that this is the same kernel as panels_affect_points!
      kernel_1_0<S,A>(&tx[2*ti[2*i]], &tx[2*ti[2*i+1]], ss[j],
                      &sx[2*j], accum.data());
    }
    // we use it backwards, so the resulting velocities are negative
    tu[2*i]   -= accum[0];
    tu[2*i+1] -= accum[1];
  }
}

template <class S, class A>
void panels_affect_panels (Panels<S> const& src, Panels<S>& targ) {
  std::cout << "    1_1 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
  // not sure how to do this - find field points of one and apply a function above?
}


// helper struct for dispatching through a variant
template <class A>
struct InfluenceVisitor {
  // source collection, target collection
  void operator()(Points<float> const& src, Points<float>& targ) { points_affect_points<float,A>(src, targ); } 
  void operator()(Panels<float> const& src, Points<float>& targ) { panels_affect_points<float,A>(src, targ); } 
  void operator()(Points<float> const& src, Panels<float>& targ) { points_affect_panels<float,A>(src, targ); } 
  void operator()(Panels<float> const& src, Panels<float>& targ) { panels_affect_panels<float,A>(src, targ); } 
};

