#pragma once

#include "Kernels.h"
#include "Points.h"
#include "Panels.h"

#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>


template <class S, class A>
void points_affect_points (Points<S> const& src, Points<S>& targ) {
  std::cout << "    0_0 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  // get references to use locally
  const std::array<std::vector<S>,Dimensions>& sx = src.get_pos();
  const std::vector<S>&                        sr = src.get_rad();
  const std::vector<S>&                        ss = src.get_str();
  const std::array<std::vector<S>,Dimensions>& tx = targ.get_pos();
  const std::vector<S>&                        tr = targ.get_rad();
  std::array<std::vector<S>,Dimensions>&       tu = targ.get_vel();
  std::optional<std::array<std::vector<S>,9>>& opttug = targ.get_velgrad();
  float flops = (float)targ.getn();

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
        //kernel_0_0vg<S,A>(&sx[3*j], sr[j], &ss[3*j],
        //                  &tx[3*i], tr[i], accum.data());
        kernel_0_0sg<S,A>(sx[0][j], sx[1][j], sx[2][j], sr[j], &ss[3*j],
                          tx[0][i], tx[1][i], tx[2][i], tr[i], accum.data());
      }
      tu[0][i] += accum[0];
      tu[1][i] += accum[1];
      tu[2][i] += accum[2];
      tug[0][i] += accum[3];
      tug[1][i] += accum[4];
      tug[2][i] += accum[5];
      tug[3][i] += accum[6];
      tug[4][i] += accum[7];
      tug[5][i] += accum[8];
      tug[6][i] += accum[9];
      tug[7][i] += accum[10];
      tug[8][i] += accum[11];
    }
    flops *= 12.0 + 65.0*(float)src.getn();

  } else {
    // velocity-only kernel
    for (size_t i=0; i<targ.getn(); ++i) {
      std::array<A,3> accum = {0.0};
      for (size_t j=0; j<src.getn(); ++j) {
        //kernel_0_0v<S,A>(&sx[3*j], sr[j], &ss[3*j],
        //                 &tx[3*i], tr[i], accum.data());
        kernel_0_0s<S,A>(sx[0][j], sx[1][j], sx[2][j], sr[j], &ss[3*j],
                         tx[0][i], tx[1][i], tx[2][i], tr[i], accum.data());
      }
      tu[0][i] += accum[0];
      tu[1][i] += accum[1];
      tu[2][i] += accum[2];
    }
    flops *= 3.0 + 30.0*(float)src.getn();
  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "    time " << (float)elapsed_seconds.count() << " at " << (1.e-9*flops/elapsed_seconds.count()) << " GFlop/s" << std::endl;
}


template <class S, class A>
void panels_affect_points (Panels<S> const& src, Points<S>& targ) {
  std::cout << "    1_0 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // get references to use locally
  const std::array<std::vector<S>,Dimensions>& sx = src.get_pos();
  //const std::vector<S>&                      sr = src.get_rad();
  const std::vector<uint16_t>&                 si = src.get_idx();
  const std::vector<S>&                        ss = src.get_str();
  const std::array<std::vector<S>,Dimensions>& tx = targ.get_pos();
  //const std::vector<S>&                      tr = targ.get_rad();
  std::array<std::vector<S>,Dimensions>&       tu = targ.get_vel();

  for (size_t i=0; i<targ.getn(); ++i) {
    std::array<A,3> accum = {0.0};
    for (size_t j=0; j<src.getn(); ++j) {
      const size_t jp0 = si[2*j];
      const size_t jp1 = si[2*j+1];
      //kernel_1_0v<S,A>(&sx[2*si[2*j]], &sx[2*si[2*j+1]], ss[j],
      //                &tx[2*i], accum.data());
      kernel_1_0s<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                       sx[0][jp1], sx[1][jp1], sx[2][jp1],
                       ss[j],
                       tx[0][i], tx[1][i], tx[2][i],
                       accum.data());
    }
    tu[0][i] += accum[0];
    tu[1][i] += accum[1];
    tu[2][i] += accum[2];
  }
}


template <class S, class A>
void points_affect_panels (Points<S> const& src, Panels<S>& targ) {
  std::cout << "    0_1 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // get references to use locally
  const std::array<std::vector<S>,Dimensions>& sx = src.get_pos();
  const std::vector<S>&                        ss = src.get_str();
  const std::array<std::vector<S>,Dimensions>& tx = targ.get_pos();
  const std::vector<uint16_t>&                 ti = targ.get_idx();
  std::array<std::vector<S>,Dimensions>&       tu = targ.get_vel();

  for (size_t i=0; i<targ.getn(); ++i) {
    std::array<A,3> accum = {0.0};
    const size_t ip0 = ti[2*i];
    const size_t ip1 = ti[2*i+1];
    for (size_t j=0; j<src.getn(); ++j) {
      // note that this is the same kernel as panels_affect_points!
      //kernel_1_0v<S,A>(&tx[2*ti[2*i]], &tx[2*ti[2*i+1]], ss[j],
      //                 &sx[2*j], accum.data());
      kernel_1_0s<S,A>(tx[0][ip0], tx[1][ip0], tx[2][ip0],
                       tx[0][ip1], tx[1][ip1], tx[2][ip1],
                       ss[j],
                       sx[0][j], sx[1][j], sx[2][j],
                       accum.data());
    }
    // we use it backwards, so the resulting velocities are negative
    tu[0][i] -= accum[0];
    tu[1][i] -= accum[1];
    tu[2][i] -= accum[2];
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

