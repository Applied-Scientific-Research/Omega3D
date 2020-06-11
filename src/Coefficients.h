/*
 * Coefficients.h - Non-class influence coefficients calculations
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#define RECURSIVE_LEVELS 3

#include "Omega3D.h"
#include "VectorHelper.h"
#include "Kernels.h"
#include "Points.h"
#include "Surfaces.h"

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <algorithm>	// for std::transform
#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#include <chrono>
#include <cmath>	// for M_PI


template <class S>
Vector<S> points_on_points_coeff (Points<S> const& src, Points<S>& targ) {
  std::cout << "    0_0 compute coefficients of" << src.to_string() << " on" << targ.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  assert(false && "ERROR: points_on_points_coeff is not implemented");

  // when we need this, copy it from Influence.h
  Vector<S> coeffs;
  float flops = 0.0;

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    points_on_points_coeff: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);

  return coeffs;
}

template <class S>
Vector<S> panels_on_points_coeff (Surfaces<S> const& src, Points<S>& targ) {
  std::cout << "    1_0 compute coefficients of" << src.to_string() << " on" << targ.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  assert(false && "ERROR: panels_on_points_coeff is not implemented");

  Vector<S> coeffs;
  float flops = 0.0;

/*
  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  //const Vector<S>&                      sr = src.get_rad();
  const std::vector<uint16_t>&            si = src.get_idx();
  const std::array<Vector<S>,Dimensions>& ss = src.get_str();
  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  //const Vector<S>&                      tr = targ.get_rad();
  std::array<Vector<S>,Dimensions>&       tu = targ.get_vel();

  #pragma omp parallel for
  for (size_t i=0; i<targ.get_n(); ++i) {
    //std::array<A,3> accum = {0.0};
    A accumu = 0.0;
    A accumv = 0.0;
    A accumw = 0.0;
    for (size_t j=0; j<src.get_n(); ++j) {
      const size_t jp0 = si[2*j];
      const size_t jp1 = si[2*j+1];
      kernel_2v_0p<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                        sx[0][jp1], sx[1][jp1], sx[2][jp1],
                        ss[0][j], ss[1][j], ss[2][j],
                        tx[0][i], tx[1][i], tx[2][i],
                        &accumu, &accumv, &accumw);
    }
    //tu[0][i] += accum[0];
    //tu[1][i] += accum[1];
    //tu[2][i] += accum[2];
    tu[0][i] += accumu;
    tu[1][i] += accumv;
    tu[2][i] += accumw;
  }
*/

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    panels_on_points_coeff: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);

  return coeffs;
}


template <class S>
Vector<S> points_on_panels_coeff (Points<S> const& src, Surfaces<S>& targ) {
  std::cout << "    0_1 compute coefficients of" << src.to_string() << " on" << targ.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  assert(false && "ERROR: points_on_panels_coeff is not implemented");

  Vector<S> coeffs;
  float flops = 0.0;

/*
  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const std::array<Vector<S>,Dimensions>& ss = src.get_str();
  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  const std::vector<uint16_t>&            ti = targ.get_idx();
  std::array<Vector<S>,Dimensions>&       tu = targ.get_vel();

  #pragma omp parallel for
  for (size_t i=0; i<targ.get_n(); ++i) {
    //std::array<A,3> accum = {0.0};
    A accumu = 0.0;
    A accumv = 0.0;
    A accumw = 0.0;
    const size_t ip0 = ti[2*i];
    const size_t ip1 = ti[2*i+1];
    for (size_t j=0; j<src.get_n(); ++j) {
      // note that this is the same kernel as panels_on_points_coeff!
      kernel_2v_0p<S,A>(tx[0][ip0], tx[1][ip0], tx[2][ip0],
                        tx[0][ip1], tx[1][ip1], tx[2][ip1],
                        ss[0][j], ss[1][j], ss[2][j],
                        sx[0][j], sx[1][j], sx[2][j],
                        &accumu, &accumv, &accumw);
    }
    // we use it backwards, so the resulting velocities are negative
    //tu[0][i] -= accum[0];
    //tu[1][i] -= accum[1];
    //tu[2][i] -= accum[2];
    tu[0][i] -= accumu;
    tu[1][i] -= accumv;
    tu[2][i] -= accumw;
  }
*/

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    points_on_panels_coeff: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);

  return coeffs;
}

template <class S>
Vector<S> panels_on_panels_coeff (Surfaces<S> const& src, Surfaces<S>& targ) {
  std::cout << "    2_2 compute coefficients of" << src.to_string() << " on" << targ.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  // how large of a problem do we have?
  const size_t nsrc  = src.get_npanels();
  const size_t ntarg = targ.get_npanels();

  // and how many rows and cols does that mean?
  const bool src_has_src  = src.src_is_unknown();
  const bool targ_has_src = targ.src_is_unknown();
  const size_t snunk = src.num_unknowns_per_panel();
  const size_t tnunk = targ.num_unknowns_per_panel();
  const size_t oldncols =  nsrc * snunk;
  const size_t oldnrows = ntarg * tnunk;
  const size_t nunk = targ.num_unknowns_per_panel();
  assert(targ.num_unknowns_per_panel() == src.num_unknowns_per_panel() && "nunk are not the same");

  // pull references to the element arrays
  const std::array<Vector<S>,Dimensions>&  sx = src.get_pos();
  const std::vector<Int>&                  si = src.get_idx();
  const std::array<Vector<S>,Dimensions>& sb1 = src.get_x1();
  const std::array<Vector<S>,Dimensions>& sb2 = src.get_x2();
  const Vector<S>&                         sa = src.get_area();

  const std::array<Vector<S>,Dimensions>&  tx = targ.get_pos();
  const std::vector<Int>&                  ti = targ.get_idx();
  const std::array<Vector<S>,Dimensions>& tb1 = targ.get_x1();
  const std::array<Vector<S>,Dimensions>& tb2 = targ.get_x2();
  const std::array<Vector<S>,Dimensions>&  tn = targ.get_norm();
  const Vector<S>&                         ta = targ.get_area();

#ifdef USE_VC
  // define vector types for Vc (still only S==A supported here)
  typedef Vc::Vector<S> StoreVec;
#endif

  // allocate space for the output array
  Vector<S> coeffs;
  coeffs.resize(oldncols*oldnrows);

  // use floats to prevent overruns
  float flops = 0.0;

  // run a panels-on-points algorithm - THIS CAN BE MORE EFFICIENT
  #pragma omp parallel for
  for (int32_t j=0; j<(int32_t)nsrc; j++) {
    // we are looping over sources first, so we're computing the next nunk *columns* in the A matrix

    // store separate pointers for each of the nunk columns
    std::vector<size_t> jptr;
    for (size_t i=0; i<snunk; ++i) jptr[i] = (j*snunk + i) * (ntarg*tnunk);

    // source triangular panel stays the same
    const Int sfirst  = si[3*j];
    const Int ssecond = si[3*j+1];
    const Int sthird  = si[3*j+2];

#ifdef USE_VC
    const StoreVec sx0 = sx[0][sfirst];
    const StoreVec sy0 = sx[1][sfirst];
    const StoreVec sz0 = sx[2][sfirst];
    const StoreVec sx1 = sx[0][ssecond];
    const StoreVec sy1 = sx[1][ssecond];
    const StoreVec sz1 = sx[2][ssecond];
    const StoreVec sx2 = sx[0][sthird];
    const StoreVec sy2 = sx[1][sthird];
    const StoreVec sz2 = sx[2][sthird];
    const S sarea = sa[j];

    const size_t ntargvec = 1 + (ntarg-1) / StoreVec::size();

    for (size_t i=0; i<ntargvec; i++) {

      // fill a 4- or 8-wide vector with the target coordinates
      StoreVec tx0, ty0, tz0, tx1, ty1, tz1, tx2, ty2, tz2, tav;
      for (size_t ii=0; ii<StoreVec::size() && i*StoreVec::size()+ii<ntarg; ++ii) {
        const size_t idx = i*StoreVec::size() + ii;
        const Int tfirst  = ti[3*idx];
        const Int tsecond = ti[3*idx+1];
        const Int tthird  = ti[3*idx+2];
        tx0[ii] = tx[0][tfirst];
        ty0[ii] = tx[1][tfirst];
        tz0[ii] = tx[2][tfirst];
        tx1[ii] = tx[0][tsecond];
        ty1[ii] = tx[1][tsecond];
        tz1[ii] = tx[2][tsecond];
        tx2[ii] = tx[0][tthird];
        ty2[ii] = tx[1][tthird];
        tz2[ii] = tx[2][tthird];
        tav[ii] = ta[idx];
      }

      // influence of vortex panel j with unit circulation on center of panel i
      StoreVec resultu, resultv, resultw;
      resultu = 0.0; resultv = 0.0; resultw = 0.0;

      // first, strength along panel-centric x1 vector
      flops += rkernel_2vs_2p<StoreVec,StoreVec>(sx0, sy0, sz0,
                                      sx1, sy1, sz1,
                                      sx2, sy2, sz2,
                                      StoreVec(sb1[0][j]), StoreVec(sb1[1][j]), StoreVec(sb1[2][j]),
                                      StoreVec(0.0),
                                      tx0, ty0, tz0,
                                      tx1, ty1, tz1,
                                      tx2, ty2, tz2,
                                      StoreVec(sarea), tav, 0, RECURSIVE_LEVELS,
                                      &resultu, &resultv, &resultw);

      // dot product with tangent vector, applying normalization here
      // spread the results from a vector register back to the primary array
      for (size_t ii=0; ii<StoreVec::size() && i*StoreVec::size()+ii<ntarg; ++ii) {
        const size_t idx = i*StoreVec::size() + ii;
        coeffs[jptr[0]++] = resultu[ii]*tb1[0][idx] + resultv[ii]*tb1[1][idx] + resultw[ii]*tb1[2][idx];

        // recompute for other target vector
        coeffs[jptr[0]++] = resultu[ii]*tb2[0][idx] + resultv[ii]*tb2[1][idx] + resultw[ii]*tb2[2][idx];

        // recompute for normal
        if (targ_has_src) coeffs[jptr[0]++] = resultu[ii]*tn[0][idx] + resultv[ii]*tn[1][idx] + resultw[ii]*tn[2][idx];
      }

      // now, along x2 direction (another 168+20 flops)
      resultu = 0.0; resultv = 0.0; resultw = 0.0;
      flops += rkernel_2vs_2p<StoreVec,StoreVec>(sx0, sy0, sz0,
                                      sx1, sy1, sz1,
                                      sx2, sy2, sz2,
                                      StoreVec(sb2[0][j]), StoreVec(sb2[1][j]), StoreVec(sb2[2][j]),
                                      StoreVec(0.0),
                                      tx0, ty0, tz0,
                                      tx1, ty1, tz1,
                                      tx2, ty2, tz2,
                                      StoreVec(sarea), tav, 0, RECURSIVE_LEVELS,
                                      &resultu, &resultv, &resultw);
      for (size_t ii=0; ii<StoreVec::size() && i*StoreVec::size()+ii<ntarg; ++ii) {
        const size_t idx = i*StoreVec::size() + ii;
        coeffs[jptr[1]++] = resultu[ii]*tb1[0][idx] + resultv[ii]*tb1[1][idx] + resultw[ii]*tb1[2][idx];
        coeffs[jptr[1]++] = resultu[ii]*tb2[0][idx] + resultv[ii]*tb2[1][idx] + resultw[ii]*tb2[2][idx];
        if (targ_has_src) coeffs[jptr[1]++] = resultu[ii]*tn[0][idx] + resultv[ii]*tn[1][idx] + resultw[ii]*tn[2][idx];
      }

      // finally the influence of a unit-strength source sheet
      if (src_has_src) {
        resultu = 0.0; resultv = 0.0; resultw = 0.0;
        flops += rkernel_2vs_2p<StoreVec,StoreVec>(sx0, sy0, sz0,
                                      sx1, sy1, sz1,
                                      sx2, sy2, sz2,
                                      StoreVec(0.0), StoreVec(0.0), StoreVec(0.0),
                                      StoreVec(1.0),
                                      tx0, ty0, tz0,
                                      tx1, ty1, tz1,
                                      tx2, ty2, tz2,
                                      StoreVec(sarea), tav, 0, RECURSIVE_LEVELS,
                                      &resultu, &resultv, &resultw);
        for (size_t ii=0; ii<StoreVec::size() && i*StoreVec::size()+ii<ntarg; ++ii) {
          const size_t idx = i*StoreVec::size() + ii;
          coeffs[jptr[2]++] = resultu[ii]*tb1[0][idx] + resultv[ii]*tb1[1][idx] + resultw[ii]*tb1[2][idx];
          coeffs[jptr[2]++] = resultu[ii]*tb2[0][idx] + resultv[ii]*tb2[1][idx] + resultw[ii]*tb2[2][idx];
          if (targ_has_src) coeffs[jptr[2]++] = resultu[ii]*tn[0][idx] + resultv[ii]*tn[1][idx] + resultw[ii]*tn[2][idx];
        }
      }
    }
#else	// no Vc

    const S sx0 = sx[0][sfirst];
    const S sy0 = sx[1][sfirst];
    const S sz0 = sx[2][sfirst];
    const S sx1 = sx[0][ssecond];
    const S sy1 = sx[1][ssecond];
    const S sz1 = sx[2][ssecond];
    const S sx2 = sx[0][sthird];
    const S sy2 = sx[1][sthird];
    const S sz2 = sx[2][sthird];
    const S sarea = sa[j];

    for (size_t i=0; i<ntarg; i++) {
      // 382 flops for each of these, assuming vortex-only interactions

      const Int tfirst  = ti[3*i];
      const Int tsecond = ti[3*i+1];
      const Int tthird  = ti[3*i+2];

      // influence of panel j with unit vortex sheet strength on entirety of panel i
      S resultu, resultv, resultw;
      resultu = 0.0; resultv = 0.0; resultw = 0.0;

      // first, strength along panel-centric x1 vector
      flops += rkernel_2vs_2p<S,S> (sx0, sy0, sz0,
                                    sx1, sy1, sz1,
                                    sx2, sy2, sz2,
                                    sb1[0][j], sb1[1][j], sb1[2][j], S(0.0),
                                    tx[0][tfirst], tx[1][tfirst], tx[2][tfirst],
                                    tx[0][tsecond], tx[1][tsecond], tx[2][tsecond],
                                    tx[0][tthird], tx[1][tthird], tx[2][tthird],
                                    sarea, ta[i], 0, RECURSIVE_LEVELS,
                                    &resultu, &resultv, &resultw);

      // dot product with tangent vector, applying normalization here
      coeffs[jptr[0]++] = resultu*tb1[0][i] + resultv*tb1[1][i] + resultw*tb1[2][i];

      // recompute for other target vector
      coeffs[jptr[0]++] = resultu*tb2[0][i] + resultv*tb2[1][i] + resultw*tb2[2][i];

      // recompute for normal
      if (targ_has_src) coeffs[jptr[0]++] = resultu*tn[0][i] + resultv*tn[1][i] + resultw*tn[2][i];


      // now, along x2 direction (another 168+20 flops)
      resultu = 0.0; resultv = 0.0; resultw = 0.0;
      flops += rkernel_2vs_2p<S,S> (sx0, sy0, sz0,
                                    sx1, sy1, sz1,
                                    sx2, sy2, sz2,
                                    sb2[0][j], sb2[1][j], sb2[2][j], S(0.0),
                                    tx[0][tfirst], tx[1][tfirst], tx[2][tfirst],
                                    tx[0][tsecond], tx[1][tsecond], tx[2][tsecond],
                                    tx[0][tthird], tx[1][tthird], tx[2][tthird],
                                    sarea, ta[i], 0, RECURSIVE_LEVELS,
                                    &resultu, &resultv, &resultw);
      coeffs[jptr[1]++] = resultu*tb1[0][i] + resultv*tb1[1][i] + resultw*tb1[2][i];
      coeffs[jptr[1]++] = resultu*tb2[0][i] + resultv*tb2[1][i] + resultw*tb2[2][i];
      if (targ_has_src) coeffs[jptr[1]++] = resultu*tn[0][i] + resultv*tn[1][i] + resultw*tn[2][i];


      // finally the influence of a unit-strength source sheet
      if (src_has_src) {
        resultu = 0.0; resultv = 0.0; resultw = 0.0;

        flops += rkernel_2vs_2p<S,S> (sx0, sy0, sz0,
                                      sx1, sy1, sz1,
                                      sx2, sy2, sz2,
                                      (S)0.0, (S)0.0, (S)0.0, (S)1.0,
                                      tx[0][tfirst], tx[1][tfirst], tx[2][tfirst],
                                      tx[0][tsecond], tx[1][tsecond], tx[2][tsecond],
                                      tx[0][tthird], tx[1][tthird], tx[2][tthird],
                                      sarea, ta[i], 0, RECURSIVE_LEVELS,
                                      &resultu, &resultv, &resultw);

        coeffs[jptr[2]++] = resultu*tb1[0][i] + resultv*tb1[1][i] + resultw*tb1[2][i];
        coeffs[jptr[2]++] = resultu*tb2[0][i] + resultv*tb2[1][i] + resultw*tb2[2][i];
        if (targ_has_src) coeffs[jptr[2]++] = resultu*tn[0][i] + resultv*tn[1][i] + resultw*tn[2][i];
      }
    }
#endif

    // special case: self-influence
    if (&src == &targ) {
      // find the diagonal components
      size_t dptr = j*ntarg*nunk*nunk + j*nunk;
      // and set to 0 or pi
      coeffs[dptr]   = 0.0;
      coeffs[dptr+1] = 2.0*M_PI;
      if (targ_has_src) coeffs[dptr+2] = 0.0;

      // next column
      dptr += ntarg*nunk;
      coeffs[dptr]   = -2.0*M_PI;
      coeffs[dptr+1] = 0.0;
      if (targ_has_src) coeffs[dptr+2] = 0.0;

      // last (source) column
      if (src_has_src) {
        dptr += ntarg*nunk;
        coeffs[dptr]   = 0.0;
        coeffs[dptr+1] = 0.0;
        if (targ_has_src) coeffs[dptr+2] = 2.0*M_PI;
      }
    }
  }

#ifdef USE_VC
  flops += (float)nsrc*(float)ntarg*382.0;
#endif

  // scale all influences by the constant
  const S fac = 1.0 / (4.0 * M_PI);
  std::transform(coeffs.begin(), coeffs.end(), coeffs.begin(),
                 [fac](S elem) { return elem * fac; });
  flops += 2.0 + (float)coeffs.size();

  // debug print the top-left and bottom-right corners
  if (false) {
    const size_t nrows = nunk*ntarg;
    const size_t ncols = nunk*nsrc;
    std::cout << "Influence matrix is " << nrows << " by " << ncols << std::endl;
    std::cout << "Top-left corner of influence matrix:" << std::endl;
    for (size_t i=0; i<6; ++i) {
      for (size_t j=0; j<6; ++j) {
        std::cout << " \t" << coeffs[nrows*j+i];
      }
      std::cout << std::endl;
    }
    std::cout << "Top-right corner of influence matrix:" << std::endl;
    for (size_t i=0; i<6; ++i) {
      for (size_t j=ncols-6; j<ncols; ++j) {
        std::cout << " \t" << coeffs[nrows*j+i];
      }
      std::cout << std::endl;
    }
    std::cout << "Bottom-right corner of influence matrix:" << std::endl;
    for (size_t i=nrows-6; i<nrows; ++i) {
      for (size_t j=ncols-6; j<ncols; ++j) {
        std::cout << " \t" << coeffs[nrows*j+i];
      }
      std::cout << std::endl;
    }
  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    matrix block:\t[%.4f] cpu seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);

  // SKIP THE REST IF THERES NO AUGMENTATION - TESTING ONLY
  return coeffs;


  //
  // now we augment this "matrix" with an optional new row and column
  //

  const size_t nrows = nunk*ntarg + (targ.is_augmented() ? 3 : 0);
  const size_t ncols = nunk*nsrc  + ( src.is_augmented() ? 3 : 0);
  std::cout << "    augmenting the " << ntarg << " x " << nsrc << " block to " << nrows << " x " << ncols << std::endl;

  // make a new 1-D vector to contain the coefficients
  Vector<S> augcoeff;
  augcoeff.resize(nrows*ncols);

  // first, copy the old vector into the new and generate the new bottom row

  // original vector is column-major (src/column index changes slowest), so let's keep that
  for (size_t j=0; j<nsrc; j++) {

    // iterators in each vector
    auto c_iter = coeffs.begin() + j*ntarg;
    auto a_iter = augcoeff.begin() + j*nrows;

    // copy the next nsrc numbers into the new vector
    std::copy(c_iter, c_iter+ntarg, a_iter);

    // and add the bottom value to this column
    if (src.is_augmented()) {
      // always include the panel lengths of the source body
      a_iter += nsrc;
      // then write the last value in this column - the length of this panel
      const Int tfirst  = ti[2*j];
      const Int tsecond = ti[2*j+1];
      // target panel vector
      const S panelx = tx[0][tsecond] - tx[0][tfirst];
      const S panely = tx[1][tsecond] - tx[1][tfirst];
      const S panell = std::sqrt(panelx*panelx + panely*panely);
      // coefficient in matrix is the panel length
      *a_iter = panell;
    }
  }

  // no longer need coeffs
  coeffs.clear();

  // then add the last column, if necessary, if target rotates
  if (targ.is_augmented()) {
    auto a_iter = augcoeff.begin() + nsrc*nrows;

    // this is the velocity influence from the source body with unit rotational rate on these target panels
    std::array<Vector<S>,Dimensions> const& vel = targ.get_vel();

    // fill in the entries - the influence of the source body's panels, when vort and source terms are set
    //   such that the integration results in the flow imposed by the body's volume of vorticity,
    //   on this target panel.
    for (size_t i=0; i<ntarg; ++i) {

      // find target point - just above the panel
      // yes, I know I am calling these "source"
      const Int tfirst  = ti[2*i];
      const Int tsecond = ti[2*i+1];
      const S panelx = tx[0][tsecond] - tx[0][tfirst];
      const S panely = tx[1][tsecond] - tx[1][tfirst];
      const S panell = std::sqrt(panelx*panelx + panely*panely);
      //std::cout << "targ panel " << i << " at " << tx[0][tfirst] << " " << tx[1][tfirst] << std::endl;

      // and find the component of that velocity along the "target" panel
      *a_iter = (vel[0][i]*panelx + vel[1][i]*panely) / panell;
      ++a_iter;
    }

    // finally, the bottom corner is a 3x3 of the circulation at unit rotation of the body
    if (&src == &targ) {
      *a_iter = 2.0 * src.get_vol();
    } else {
      *a_iter = 0.0;
    }
  } else {
    // if there's no body pointer, then what?
  }

  // debug print the top-left and bottom-right corners
  if (true) {
    std::cout << "Influence matrix is " << nrows << " by " << ncols << std::endl;
    std::cout << "Top-left corner of influence matrix:" << std::endl;
    for (size_t i=0; i<6; ++i) {
      for (size_t j=0; j<6; ++j) {
        std::cout << " \t" << augcoeff[nrows*j+i];
      }
      std::cout << std::endl;
    }
    std::cout << "Top-right corner of influence matrix:" << std::endl;
    for (size_t i=0; i<6; ++i) {
      for (size_t j=ncols-6; j<ncols; ++j) {
        std::cout << " \t" << augcoeff[nrows*j+i];
      }
      std::cout << std::endl;
    }
    std::cout << "Bottom-right corner of influence matrix:" << std::endl;
    for (size_t i=nrows-6; i<nrows; ++i) {
      for (size_t j=ncols-6; j<ncols; ++j) {
        std::cout << " \t" << augcoeff[nrows*j+i];
      }
      std::cout << std::endl;
    }
  }

  return augcoeff;
}


// helper struct for dispatching through a variant
struct CoefficientVisitor {
  // source collection, target collection
  Vector<float> operator()(Points<float> const& src,   Points<float>& targ)   { return points_on_points_coeff<float>(src, targ); } 
  Vector<float> operator()(Surfaces<float> const& src, Points<float>& targ)   { return panels_on_points_coeff<float>(src, targ); } 
  Vector<float> operator()(Points<float> const& src,   Surfaces<float>& targ) { return points_on_panels_coeff<float>(src, targ); } 
  Vector<float> operator()(Surfaces<float> const& src, Surfaces<float>& targ) { return panels_on_panels_coeff<float>(src, targ); } 
};

