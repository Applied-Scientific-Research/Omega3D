/*
 * RHS.h - Non-class velocity-to-right-hand-side calculations
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega3D.h"
#include "VectorHelper.h"
#include "Points.h"
#include "Surfaces.h"

#include <iostream>
#include <vector>
#include <cassert>


template <class S>
std::vector<S> vels_to_rhs_points (Points<S> const& targ) {
  std::cout << "    NOT converting vels to RHS vector for " << targ.to_string() << std::endl;

  //auto start = std::chrono::system_clock::now();
  //float flops = 0.0;

  // size the return vector
  size_t ntarg  = targ.get_n();
  std::vector<S> rhs;
  rhs.resize(ntarg);

  //auto end = std::chrono::system_clock::now();
  //std::chrono::duration<double> elapsed_seconds = end-start;
  //const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  //printf("    points_on_points_coeff: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);

  return rhs;
}

template <class S>
std::vector<S> vels_to_rhs_panels (Surfaces<S> const& targ) {
  std::cout << "    convert vels to RHS vector for" << targ.to_string() << std::endl;

  // pull references to the element arrays
  //const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  //const std::vector<Int>&                 ti = targ.get_idx();
  const std::array<Vector<S>,Dimensions>& x1 = targ.get_x1();
  const std::array<Vector<S>,Dimensions>& x2 = targ.get_x2();
  const std::array<Vector<S>,Dimensions>& norm = targ.get_norm();
  //const Vector<S>&                        area = targ.get_area();
  const std::array<Vector<S>,Dimensions>& tu = targ.get_vel();
  //const Vector<S>&                      tb = targ.get_norm_bcs();

  //std::cout << "size of x1 is " << x1.size() << " by " << x1[0].size() << std::endl;
  //std::cout << "size of x2 is " << x2.size() << " by " << x2[0].size() << std::endl;
  //std::cout << "size of norm is " << norm.size() << " by " << norm[0].size() << std::endl;
  //std::cout << "size of tu is " << tu.size() << " by " << tu[0].size() << std::endl;
  //std::cout << "size of tb is " << tb.size() << " by " << tb[0].size() << std::endl;

  // check for data
  assert(x1[0].size() == x2[0].size() && "Array size mismatch");
  assert(x1[0].size() == norm[0].size() && "Array size mismatch");
  assert(x1[0].size() == tu[0].size() && "Array size mismatch");
  //assert(tb.size() > 0 && "Target BC vector empty");
  //assert(x1[0].size() == tb[0].size() && "Array size mismatch");

  // find array sizes
  const size_t ntarg = targ.get_npanels();
  const size_t nunk  = targ.num_unknowns_per_panel();

  // prepare the rhs vector
  std::vector<S> rhs;
  rhs.resize(ntarg*nunk);

  // convert velocity and boundary condition to RHS values
  // ONLY include the influence of the source boundary condition here, apply vortex before shedding

  if (nunk == 1) {
    // normal-only
    for (size_t i=0; i<ntarg; ++i) {
      rhs[i] = -(tu[0][i]*norm[0][i] + tu[1][i]*norm[1][i] + tu[2][i]*norm[2][i]);
    }
  } else if (nunk == 2) {
    // dot product of x1 tangent with local velocity, applying normalization
    for (size_t i=0; i<ntarg; ++i) {
      rhs[2*i+0] = -(tu[0][i]*x1[0][i] + tu[1][i]*x1[1][i] + tu[2][i]*x1[2][i]);
      rhs[2*i+1] = -(tu[0][i]*x2[0][i] + tu[1][i]*x2[1][i] + tu[2][i]*x2[2][i]);
    }
  } else if (nunk == 3) {
    // two tangentials then the normal
    for (size_t i=0; i<ntarg; ++i) {
      rhs[3*i+0] = -(tu[0][i]*x1[0][i] + tu[1][i]*x1[1][i] + tu[2][i]*x1[2][i]);
      rhs[3*i+1] = -(tu[0][i]*x2[0][i] + tu[1][i]*x2[1][i] + tu[2][i]*x2[2][i]);
      rhs[3*i+2] = -(tu[0][i]*norm[0][i] + tu[1][i]*norm[1][i] + tu[2][i]*norm[2][i]);
    }
  }

  // DO NOT include the influence of (at least tangential) boundary condition here, do it before shedding
  //   nunk=1 is normal-only
  //   nunk=2 is tangential-only (x1, x2)
  //   nunk=3 is tangential and normal (x1, x2, norm)
  //for (size_t i=0; i<ntarg; ++i) {
    //for (size_t j=0; j<nunk; ++j) {
      //rhs[nunk*i+j] -= tb[j][i];
    //}
  //}

  return rhs;
}


// helper struct for dispatching through a variant
struct RHSVisitor {
  // source collection, target collection
  std::vector<float> operator()(Points<float> const& targ)   { return vels_to_rhs_points<float>(targ); } 
  std::vector<float> operator()(Surfaces<float> const& targ) { return vels_to_rhs_panels<float>(targ); } 
};

