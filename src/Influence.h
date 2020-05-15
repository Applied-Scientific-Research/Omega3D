/*
 * Influence.h - Non-class influence calculations
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

#ifdef EXTERNAL_VEL_SOLVE
extern "C" float external_vel_solver_f_(int*, const float*, const float*, const float*,
                                              const float*, const float*, const float*,
                                              const float*,
                                        int*, const float*, const float*, const float*,
                                              float*, float*, float*, float*,
                                              float*, float*, float*, float*,
                                              float*, float*, float*, float*);
extern "C" float external_vel_solver_d_(int*, const double*, const double*, const double*,
                                              const double*, const double*, const double*,
                                              const double*,
                                        int*, const double*, const double*, const double*,
                                              double*, double*, double*, double*,
                                              double*, double*, double*, double*,
                                              double*, double*, double*, double*);
#endif

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#include <chrono>
#include <thread>
#include <cmath>
#include <cassert>


//
// Vc and x86 versions of Points/Particles affecting Points/Particles
//
template <class S, class A>
void points_affect_points (Points<S> const& src, Points<S>& targ) {
  auto start = std::chrono::system_clock::now();
  float flops = (float)targ.get_n();


  // get references to use locally
  const std::array<Vector<S>,Dimensions>&     sx = src.get_pos();
  const Vector<S>&                            sr = src.get_rad();
  const std::array<Vector<S>,Dimensions>&     ss = src.get_str();

  const std::array<Vector<S>,Dimensions>&     tx = targ.get_pos();
  std::array<Vector<S>,Dimensions>&           tu = targ.get_vel();
  std::optional<std::array<Vector<S>,9>>& opttug = targ.get_velgrad();

#ifdef EXTERNAL_VEL_SOLVE
  std::cout << "    external influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
  assert(opttug && "Targets do not have velocity gradients in points_affect_points");
  if (opttug) {
    std::array<Vector<S>,9>& tug = *opttug;
    int ns = src.get_n();
    int nt = targ.get_n();
    flops = external_vel_solver_f_(&ns, sx[0].data(), sx[1].data(), sx[2].data(),
                                        ss[0].data(), ss[1].data(), ss[2].data(), sr.data(), 
                                   &nt, tx[0].data(), tx[1].data(), tx[2].data(),
                                        tu[0].data(), tu[1].data(), tu[2].data(),
                                        tug[0].data(), tug[1].data(), tug[2].data(),
                                        tug[3].data(), tug[4].data(), tug[5].data(),
                                        tug[6].data(), tug[7].data(), tug[8].data());
  }
#else  // no external fast solve, perform O(N^2) calculations here

#ifdef USE_OGL_COMPUTE
  // atomic must be ready to accept computation
  if (not targ.is_compute_still_working()) {
    //std::cout << "starting work" << std::endl << std::flush;

    auto gcs = targ.get_gcs();

    // fill in gcs's idea of sources
    // the positions
    gcs->hsx.resize(4*src.get_n());
    for (size_t i=0; i<src.get_n(); ++i) {
      gcs->hsx[4*i+0] = sx[0][i];
      gcs->hsx[4*i+1] = sx[1][i];
      gcs->hsx[4*i+2] = sx[2][i];
      gcs->hsx[4*i+3] = sr[i];
    }
    // the strengths
    gcs->hss.resize(4*src.get_n());
    for (size_t i=0; i<src.get_n(); ++i) {
      gcs->hss[4*i+0] = ss[0][i];
      gcs->hss[4*i+1] = ss[1][i];
      gcs->hss[4*i+2] = ss[2][i];
      gcs->hss[4*i+3] = 0.0;
    }

    // fill in gcs's idea of targets
    gcs->htx.resize(4*targ.get_n());
    for (size_t i=0; i<targ.get_n(); ++i) {
      gcs->htx[4*i+0] = tx[0][i];
      gcs->htx[4*i+1] = tx[1][i];
      gcs->htx[4*i+2] = tx[2][i];
      gcs->htx[4*i+3] = 0.0;
    }
    if (not targ.is_inert()) {
      const Vector<S>& tr = targ.get_rad();
      for (size_t i=0; i<targ.get_n(); ++i) gcs->htx[4*i+3] = tr[i];
    }

    // fill in gcs's space for results
    gcs->hr1.resize(4*targ.get_n());
    gcs->hr2.resize(4*targ.get_n());
    gcs->hr3.resize(4*targ.get_n());

    // then tell the graphics thread to begin computing
    targ.trigger_compute();

    // then hold here until its done
    while (targ.is_compute_still_working()) {
      // check every millisecond
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    // retrieve the results from the gcs vector
    for (size_t i=0; i<targ.get_n(); ++i) {
      tu[0][i] = gcs->hr1[4*i+0];
      tu[1][i] = gcs->hr1[4*i+1];
      tu[2][i] = gcs->hr1[4*i+2];
    }
    if (opttug) {
      std::array<Vector<S>,9>& tug = *opttug;
      for (size_t i=0; i<targ.get_n(); ++i) {
        tug[0][i] = gcs->hr1[4*i+3];
        tug[1][i] = gcs->hr2[4*i+0];
        tug[2][i] = gcs->hr2[4*i+1];
        tug[3][i] = gcs->hr2[4*i+2];
        tug[4][i] = gcs->hr2[4*i+3];
        tug[5][i] = gcs->hr3[4*i+0];
        tug[6][i] = gcs->hr3[4*i+1];
        tug[7][i] = gcs->hr3[4*i+2];
        tug[8][i] = gcs->hr3[4*i+3];
      }
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    float flops = (float)targ.get_n() * (12.0 + (float)flops_0v_0bg<S,A>()*(float)src.get_n());
    const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
    printf("    ptptvelgrad shader: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);

    return;
  }
#endif

#ifdef USE_VC
  // define vector types for Vc
  typedef Vc::Vector<S> StoreVec;
  typedef Vc::SimdArray<A, Vc::Vector<S>::size()> AccumVec;

  // create float_v versions of the source vectors
  const Vc::Memory<StoreVec> sxv  = stdvec_to_vcvec<S>(sx[0], 0.0);
  const Vc::Memory<StoreVec> syv  = stdvec_to_vcvec<S>(sx[1], 0.0);
  const Vc::Memory<StoreVec> szv  = stdvec_to_vcvec<S>(sx[2], 0.0);
  const Vc::Memory<StoreVec> srv  = stdvec_to_vcvec<S>(sr,    1.0);
  const Vc::Memory<StoreVec> ssxv = stdvec_to_vcvec<S>(ss[0], 0.0);
  const Vc::Memory<StoreVec> ssyv = stdvec_to_vcvec<S>(ss[1], 0.0);
  const Vc::Memory<StoreVec> sszv = stdvec_to_vcvec<S>(ss[2], 0.0);
#endif

  // We need 8 different loops here, for the options:
  //   target radii or no target radii
  //   grads or no grads
  //   Vc or no Vc

  //
  // targets are field points, with no core radius ===============================================
  //
  if (targ.is_inert()) {

  // dispatch on presence of val grads
  if (opttug) {
    std::cout << "    0v_0pg compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
    // get the pointer from the optional
    std::array<Vector<S>,9>& tug = *opttug;

    // velocity+grads kernel
    #pragma omp parallel for
    for (int32_t i=0; i<(int32_t)targ.get_n(); ++i) {
#ifdef USE_VC
      const StoreVec txv(tx[0][i]);
      const StoreVec tyv(tx[1][i]);
      const StoreVec tzv(tx[2][i]);
      AccumVec accumu(0.0);
      AccumVec accumv(0.0);
      AccumVec accumw(0.0);
      AccumVec accumux = 0.0;
      AccumVec accumvx = 0.0;
      AccumVec accumwx = 0.0;
      AccumVec accumuy = 0.0;
      AccumVec accumvy = 0.0;
      AccumVec accumwy = 0.0;
      AccumVec accumuz = 0.0;
      AccumVec accumvz = 0.0;
      AccumVec accumwz = 0.0;
      // loop over source particles
      for (size_t j=0; j<sxv.vectorsCount(); ++j) {
        // NOTE: .vectorAt(i) gets the vector at scalar position i
        //       .vector(i) gets the i'th vector!!!
        kernel_0v_0pg<StoreVec,AccumVec>(
                           sxv.vector(j), syv.vector(j), szv.vector(j), srv.vector(j),
                           ssxv.vector(j), ssyv.vector(j), sszv.vector(j),
                           txv, tyv, tzv,
                           &accumu,  &accumv,  &accumw,
                           &accumux, &accumvx, &accumwx,
                           &accumuy, &accumvy, &accumwy,
                           &accumuz, &accumvz, &accumwz);
      }
      tu[0][i] += accumu.sum();
      tu[1][i] += accumv.sum();
      tu[2][i] += accumw.sum();
      tug[0][i] += accumux.sum();
      tug[1][i] += accumvx.sum();
      tug[2][i] += accumwx.sum();
      tug[3][i] += accumuy.sum();
      tug[4][i] += accumvy.sum();
      tug[5][i] += accumwy.sum();
      tug[6][i] += accumuz.sum();
      tug[7][i] += accumvz.sum();
      tug[8][i] += accumwz.sum();
#else  // no Vc
      A accumu = 0.0;
      A accumv = 0.0;
      A accumw = 0.0;
      A accumux = 0.0;
      A accumvx = 0.0;
      A accumwx = 0.0;
      A accumuy = 0.0;
      A accumvy = 0.0;
      A accumwy = 0.0;
      A accumuz = 0.0;
      A accumvz = 0.0;
      A accumwz = 0.0;
      // loop over source particles
      for (size_t j=0; j<src.get_n(); ++j) {
        kernel_0v_0pg<S,A>(sx[0][j], sx[1][j], sx[2][j], sr[j],
                           ss[0][j], ss[1][j], ss[2][j],
                           tx[0][i], tx[1][i], tx[2][i],
                           &accumu,  &accumv,  &accumw,
                           &accumux, &accumvx, &accumwx,
                           &accumuy, &accumvy, &accumwy,
                           &accumuz, &accumvz, &accumwz);
      }
      tu[0][i] += accumu;
      tu[1][i] += accumv;
      tu[2][i] += accumw;
      tug[0][i] += accumux;
      tug[1][i] += accumvx;
      tug[2][i] += accumwx;
      tug[3][i] += accumuy;
      tug[4][i] += accumvy;
      tug[5][i] += accumwy;
      tug[6][i] += accumuz;
      tug[7][i] += accumvz;
      tug[8][i] += accumwz;
#endif // no Vc
    }
    flops *= 12.0 + (float)flops_0v_0pg<S,A>() * (float)src.get_n();

  } else {

    // velocity-only kernel
    std::cout << "    0v_0p compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
    #pragma omp parallel for
    for (size_t i=0; i<targ.get_n(); ++i) {
#ifdef USE_VC
      const StoreVec txv = tx[0][i];
      const StoreVec tyv = tx[1][i];
      const StoreVec tzv = tx[2][i];
      AccumVec accumu = 0.0;
      AccumVec accumv = 0.0;
      AccumVec accumw = 0.0;
      for (size_t j=0; j<sxv.vectorsCount(); ++j) {
        kernel_0v_0p<StoreVec,AccumVec>(
                          sxv.vector(j), syv.vector(j), szv.vector(j), srv.vector(j),
                          ssxv.vector(j), ssyv.vector(j), sszv.vector(j),
                          txv, tyv, tzv,
                          &accumu, &accumv, &accumw);
      }
      tu[0][i] += accumu.sum();
      tu[1][i] += accumv.sum();
      tu[2][i] += accumw.sum();
#else  // no Vc
      A accumu = 0.0;
      A accumv = 0.0;
      A accumw = 0.0;
      for (size_t j=0; j<src.get_n(); ++j) {
        kernel_0v_0p<S,A>(sx[0][j], sx[1][j], sx[2][j], sr[j],
                          ss[0][j], ss[1][j], ss[2][j],
                          tx[0][i], tx[1][i], tx[2][i],
                          &accumu, &accumv, &accumw);
      }
      tu[0][i] += accumu;
      tu[1][i] += accumv;
      tu[2][i] += accumw;
#endif // no Vc
    }
    flops *= 3.0 + (float)flops_0v_0p<S,A>() * (float)src.get_n();
  }

  //
  // targets are particles, with a core radius ===================================================
  //
  } else {

    const Vector<S>&                            tr = targ.get_rad();

  // dispatch on presence of val grads
  if (opttug) {
    std::cout << "    0v_0vg compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
    // get the pointer from the optional
    std::array<Vector<S>,9>& tug = *opttug;

    // velocity+grads kernel
    #pragma omp parallel for
    for (int32_t i=0; i<(int32_t)targ.get_n(); ++i) {
#ifdef USE_VC
      const StoreVec txv(tx[0][i]);
      const StoreVec tyv(tx[1][i]);
      const StoreVec tzv(tx[2][i]);
      const StoreVec trv(tr[i]);
      AccumVec accumu(0.0);
      AccumVec accumv(0.0);
      AccumVec accumw(0.0);
      AccumVec accumux = 0.0;
      AccumVec accumvx = 0.0;
      AccumVec accumwx = 0.0;
      AccumVec accumuy = 0.0;
      AccumVec accumvy = 0.0;
      AccumVec accumwy = 0.0;
      AccumVec accumuz = 0.0;
      AccumVec accumvz = 0.0;
      AccumVec accumwz = 0.0;
      // loop over source particles
      for (size_t j=0; j<sxv.vectorsCount(); ++j) {
        // NOTE: .vectorAt(i) gets the vector at scalar position i
        //       .vector(i) gets the i'th vector!!!
        kernel_0v_0bg<StoreVec,AccumVec>(
                          sxv.vector(j), syv.vector(j), szv.vector(j), srv.vector(j),
                          ssxv.vector(j), ssyv.vector(j), sszv.vector(j),
                          txv, tyv, tzv, trv,
                          &accumu,  &accumv,  &accumw,
                          &accumux, &accumvx, &accumwx,
                          &accumuy, &accumvy, &accumwy,
                          &accumuz, &accumvz, &accumwz);
      }
      tu[0][i] += accumu.sum();
      tu[1][i] += accumv.sum();
      tu[2][i] += accumw.sum();
      tug[0][i] += accumux.sum();
      tug[1][i] += accumvx.sum();
      tug[2][i] += accumwx.sum();
      tug[3][i] += accumuy.sum();
      tug[4][i] += accumvy.sum();
      tug[5][i] += accumwy.sum();
      tug[6][i] += accumuz.sum();
      tug[7][i] += accumvz.sum();
      tug[8][i] += accumwz.sum();
#else  // no Vc
      A accumu = 0.0;
      A accumv = 0.0;
      A accumw = 0.0;
      A accumux = 0.0;
      A accumvx = 0.0;
      A accumwx = 0.0;
      A accumuy = 0.0;
      A accumvy = 0.0;
      A accumwy = 0.0;
      A accumuz = 0.0;
      A accumvz = 0.0;
      A accumwz = 0.0;
      // loop over source particles
      for (size_t j=0; j<src.get_n(); ++j) {
        kernel_0v_0bg<S,A>(sx[0][j], sx[1][j], sx[2][j], sr[j],
                          ss[0][j], ss[1][j], ss[2][j],
                          tx[0][i], tx[1][i], tx[2][i], tr[i],
                          &accumu,  &accumv,  &accumw,
                          &accumux, &accumvx, &accumwx,
                          &accumuy, &accumvy, &accumwy,
                          &accumuz, &accumvz, &accumwz);
      }
      tu[0][i] += accumu;
      tu[1][i] += accumv;
      tu[2][i] += accumw;
      tug[0][i] += accumux;
      tug[1][i] += accumvx;
      tug[2][i] += accumwx;
      tug[3][i] += accumuy;
      tug[4][i] += accumvy;
      tug[5][i] += accumwy;
      tug[6][i] += accumuz;
      tug[7][i] += accumvz;
      tug[8][i] += accumwz;
#endif // no Vc
    }
    flops *= 12.0 + (float)flops_0v_0bg<S,A>() * (float)src.get_n();

  } else {

    // velocity-only kernel
    std::cout << "    0v_0v compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
    #pragma omp parallel for
    for (int32_t i=0; i<(int32_t)targ.get_n(); ++i) {
#ifdef USE_VC
      const StoreVec txv = tx[0][i];
      const StoreVec tyv = tx[1][i];
      const StoreVec tzv = tx[2][i];
      const StoreVec trv = tr[i];
      AccumVec accumu = 0.0;
      AccumVec accumv = 0.0;
      AccumVec accumw = 0.0;
      for (size_t j=0; j<sxv.vectorsCount(); ++j) {
        kernel_0v_0b<StoreVec,AccumVec>(
                         sxv.vector(j), syv.vector(j), szv.vector(j), srv.vector(j),
                         ssxv.vector(j), ssyv.vector(j), sszv.vector(j),
                         txv, tyv, tzv, trv,
                         &accumu, &accumv, &accumw);
      }
      tu[0][i] += accumu.sum();
      tu[1][i] += accumv.sum();
      tu[2][i] += accumw.sum();
#else  // no Vc
      A accumu = 0.0;
      A accumv = 0.0;
      A accumw = 0.0;
      for (size_t j=0; j<src.get_n(); ++j) {
        kernel_0v_0b<S,A>(sx[0][j], sx[1][j], sx[2][j], sr[j],
                         ss[0][j], ss[1][j], ss[2][j],
                         tx[0][i], tx[1][i], tx[2][i], tr[i],
                         &accumu, &accumv, &accumw);
      }
      tu[0][i] += accumu;
      tu[1][i] += accumv;
      tu[2][i] += accumw;
#endif // no Vc
    }
    flops *= 3.0 + (float)flops_0v_0b<S,A>() * (float)src.get_n();
  }

  //
  // end conditional over whether targets are field points (with no core radius)
  //
  }
#endif // no external fast solve

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    points_affect_points: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


//
// Vc and x86 versions of Panels/Surfaces affecting Points/Particles
//
template <class S, class A>
void panels_affect_points (Surfaces<S> const& src, Points<S>& targ) {
  auto start = std::chrono::system_clock::now();

  // get references to use locally
  const std::array<Vector<S>,Dimensions>&     sx = src.get_pos();
  const std::vector<Int>&                     si = src.get_idx();
  const std::array<Vector<S>,Dimensions>&     ss = src.get_str();
  const bool                              havess = src.have_src_str();
  const Vector<S>&                           sss = src.get_src_str();
  const Vector<S>&                            sa = src.get_area();
  const std::array<Vector<S>,Dimensions>&     tx = targ.get_pos();
  std::array<Vector<S>,Dimensions>&           tu = targ.get_vel();
  std::optional<std::array<Vector<S>,9>>& opttug = targ.get_velgrad();

  const int32_t ntarg = targ.get_n();
  float flops = 0.0;

#ifdef USE_VC
  // define vector types for Vc (still only S==A supported here)
  typedef Vc::Vector<S> StoreVec;
  typedef Vc::SimdArray<A, Vc::Vector<S>::size()> AccumVec;

  // prepare the source panels for vectorization - first the strengths
  const Vc::Memory<StoreVec> sav  = stdvec_to_vcvec<S>(sa, 0.0);
  Vc::Memory<StoreVec> ssxv = stdvec_to_vcvec<S>(ss[0], 0.0);
  Vc::Memory<StoreVec> ssyv = stdvec_to_vcvec<S>(ss[1], 0.0);
  Vc::Memory<StoreVec> sszv = stdvec_to_vcvec<S>(ss[2], 0.0);
  const Vc::Memory<StoreVec> sssv = stdvec_to_vcvec<S>(sss, 0.0);
  for (size_t j=0; j<src.get_npanels(); ++j) {
    ssxv[j] /= sa[j];
    ssyv[j] /= sa[j];
    sszv[j] /= sa[j];
  }
  // then the triangle nodes
  //Vector<S> sx0, sy0, sz0, sx1, sy1, sz1, sx2, sy2, sz2;
  Vector<S> sx0(src.get_npanels());
  Vector<S> sy0(src.get_npanels());
  Vector<S> sz0(src.get_npanels());
  Vector<S> sx1(src.get_npanels());
  Vector<S> sy1(src.get_npanels());
  Vector<S> sz1(src.get_npanels());
  Vector<S> sx2(src.get_npanels());
  Vector<S> sy2(src.get_npanels());
  Vector<S> sz2(src.get_npanels());
  for (size_t j=0; j<src.get_npanels(); ++j) {
    const size_t jp0 = si[3*j];
    const size_t jp1 = si[3*j+1];
    const size_t jp2 = si[3*j+2];
    sx0[j] = sx[0][jp0];
    sy0[j] = sx[1][jp0];
    sz0[j] = sx[2][jp0];
    sx1[j] = sx[0][jp1];
    sy1[j] = sx[1][jp1];
    sz1[j] = sx[2][jp1];
    sx2[j] = sx[0][jp2];
    sy2[j] = sx[1][jp2];
    sz2[j] = sx[2][jp2];
  }
  const Vc::Memory<StoreVec> sx0v = stdvec_to_vcvec<S>(sx0, -1.0);
  const Vc::Memory<StoreVec> sy0v = stdvec_to_vcvec<S>(sy0, -1.0);
  const Vc::Memory<StoreVec> sz0v = stdvec_to_vcvec<S>(sz0, 9.0);
  const Vc::Memory<StoreVec> sx1v = stdvec_to_vcvec<S>(sx1, 0.0);
  const Vc::Memory<StoreVec> sy1v = stdvec_to_vcvec<S>(sy1, 1.0);
  const Vc::Memory<StoreVec> sz1v = stdvec_to_vcvec<S>(sz1, 9.0);
  const Vc::Memory<StoreVec> sx2v = stdvec_to_vcvec<S>(sx2, 1.0);
  const Vc::Memory<StoreVec> sy2v = stdvec_to_vcvec<S>(sy2, -1.0);
  const Vc::Memory<StoreVec> sz2v = stdvec_to_vcvec<S>(sz2, 9.0);
#endif

  // We need 8 different loops here, for the options:
  //   target radii or no target radii
  //   grads or no grads
  //   Vc or no Vc

  //
  // targets are field points, with no core radius ===============================================
  //
  if (targ.is_inert()) {

    if (opttug) { // velocity-and-grads kernel -------------------------------------------------
      if (havess) {
        std::cout << "    2vs_0pg compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
      } else {
        std::cout << "    2v_0pg compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
      }

      // get the pointer from the optional
      std::array<Vector<S>,9>& tug = *opttug;

      #ifdef USE_VC

        #pragma omp parallel for
        for (int32_t i=0; i<ntarg; ++i) {
          const StoreVec txv(tx[0][i]);
          const StoreVec tyv(tx[1][i]);
          const StoreVec tzv(tx[2][i]);
          AccumVec accumu(0.0);
          AccumVec accumv(0.0);
          AccumVec accumw(0.0);
          AccumVec accumux(0.0);
          AccumVec accumvx(0.0);
          AccumVec accumwx(0.0);
          AccumVec accumuy(0.0);
          AccumVec accumvy(0.0);
          AccumVec accumwy(0.0);
          AccumVec accumuz(0.0);
          AccumVec accumvz(0.0);
          AccumVec accumwz(0.0);
          if (havess) {
            for (size_t j=0; j<sx0v.vectorsCount(); ++j) {
              // NOTE: .vectorAt(i) gets the vector at scalar position i
              //       .vector(i) gets the i'th vector!!!
              flops += rkernel_2vs_0pg<StoreVec,AccumVec>(sx0v.vector(j), sy0v.vector(j), sz0v.vector(j),
                                   sx1v.vector(j), sy1v.vector(j), sz1v.vector(j),
                                   sx2v.vector(j), sy2v.vector(j), sz2v.vector(j),
                                   ssxv.vector(j), ssyv.vector(j), sszv.vector(j),
                                   sssv.vector(j),
                                   txv, tyv, tzv,
                                   sav.vector(j), 0, RECURSIVE_LEVELS,
                                   &accumu, &accumv, &accumw,
                                   &accumux, &accumvx, &accumwx,
                                   &accumuy, &accumvy, &accumwy,
                                   &accumuz, &accumvz, &accumwz);
            }
          } else {
            for (size_t j=0; j<sx0v.vectorsCount(); ++j) {
              flops += rkernel_2vs_0pg<StoreVec,AccumVec>(sx0v.vector(j), sy0v.vector(j), sz0v.vector(j),
                                   sx1v.vector(j), sy1v.vector(j), sz1v.vector(j),
                                   sx2v.vector(j), sy2v.vector(j), sz2v.vector(j),
                                   ssxv.vector(j), ssyv.vector(j), sszv.vector(j),
                                   StoreVec(0.0),
                                   txv, tyv, tzv,
                                   sav.vector(j), 0, RECURSIVE_LEVELS,
                                   &accumu, &accumv, &accumw,
                                   &accumux, &accumvx, &accumwx,
                                   &accumuy, &accumvy, &accumwy,
                                   &accumuz, &accumvz, &accumwz);
            }
          }
          tu[0][i] += accumu.sum();
          tu[1][i] += accumv.sum();
          tu[2][i] += accumw.sum();
          tug[0][i] += accumux.sum();
          tug[1][i] += accumvx.sum();
          tug[2][i] += accumwx.sum();
          tug[3][i] += accumuy.sum();
          tug[4][i] += accumvy.sum();
          tug[5][i] += accumwy.sum();
          tug[6][i] += accumuz.sum();
          tug[7][i] += accumvz.sum();
          tug[8][i] += accumwz.sum();
          flops += 12.0;
        }

      #else  // no Vc

        #pragma omp parallel for
        for (int32_t i=0; i<ntarg; ++i) {
          A accumu = 0.0;
          A accumv = 0.0;
          A accumw = 0.0;
          A accumux = 0.0;
          A accumvx = 0.0;
          A accumwx = 0.0;
          A accumuy = 0.0;
          A accumvy = 0.0;
          A accumwy = 0.0;
          A accumuz = 0.0;
          A accumvz = 0.0;
          A accumwz = 0.0;
          if (havess) {
            for (size_t j=0; j<src.get_npanels(); ++j) {
              const size_t jp0 = si[3*j];
              const size_t jp1 = si[3*j+1];
              const size_t jp2 = si[3*j+2];
              flops += rkernel_2vs_0pg<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                                   sx[0][jp1], sx[1][jp1], sx[2][jp1],
                                   sx[0][jp2], sx[1][jp2], sx[2][jp2],
                                   ss[0][j]/sa[j], ss[1][j]/sa[j], ss[2][j]/sa[j], sss[j],
                                   tx[0][i], tx[1][i], tx[2][i],
                                   sa[j], 0, RECURSIVE_LEVELS,
                                   &accumu, &accumv, &accumw,
                                   &accumux, &accumvx, &accumwx,
                                   &accumuy, &accumvy, &accumwy,
                                   &accumuz, &accumvz, &accumwz);
            }
          } else {
            for (size_t j=0; j<src.get_npanels(); ++j) {
              const size_t jp0 = si[3*j];
              const size_t jp1 = si[3*j+1];
              const size_t jp2 = si[3*j+2];
              flops += rkernel_2vs_0pg<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                                   sx[0][jp1], sx[1][jp1], sx[2][jp1],
                                   sx[0][jp2], sx[1][jp2], sx[2][jp2],
                                   ss[0][j]/sa[j], ss[1][j]/sa[j], ss[2][j]/sa[j], S(0.0),
                                   tx[0][i], tx[1][i], tx[2][i],
                                   sa[j], 0, RECURSIVE_LEVELS,
                                   &accumu, &accumv, &accumw,
                                   &accumux, &accumvx, &accumwx,
                                   &accumuy, &accumvy, &accumwy,
                                   &accumuz, &accumvz, &accumwz);
            }
          }
          tu[0][i] += accumu;
          tu[1][i] += accumv;
          tu[2][i] += accumw;
          tug[0][i] += accumux;
          tug[1][i] += accumvx;
          tug[2][i] += accumwx;
          tug[3][i] += accumuy;
          tug[4][i] += accumvy;
          tug[5][i] += accumwy;
          tug[6][i] += accumuz;
          tug[7][i] += accumvz;
          tug[8][i] += accumwz;
          flops += 12.0;
        }
      #endif // no Vc

    } else { // velocity-only kernel -----------------------------------------------------------
      if (havess) {
        std::cout << "    2vs_0p compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
      } else {
        std::cout << "    2v_0p compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
      }

      #ifdef USE_VC
        #pragma omp parallel for
        for (int32_t i=0; i<ntarg; ++i) {
          const StoreVec txv(tx[0][i]);
          const StoreVec tyv(tx[1][i]);
          const StoreVec tzv(tx[2][i]);
          AccumVec accumu(0.0);
          AccumVec accumv(0.0);
          AccumVec accumw(0.0);

          if (havess) {
            for (size_t j=0; j<sx0v.vectorsCount(); ++j) {
              // NOTE: .vectorAt(i) gets the vector at scalar position i
              //       .vector(i) gets the i'th vector!!!
              flops += rkernel_2vs_0p<StoreVec,AccumVec>(sx0v.vector(j), sy0v.vector(j), sz0v.vector(j),
                                               sx1v.vector(j), sy1v.vector(j), sz1v.vector(j),
                                               sx2v.vector(j), sy2v.vector(j), sz2v.vector(j),
                                               ssxv.vector(j), ssyv.vector(j), sszv.vector(j),
                                               sssv.vector(j),
                                               txv, tyv, tzv,
                                               sav.vector(j), 0, RECURSIVE_LEVELS,
                                               &accumu, &accumv, &accumw);
            }
          } else {
            for (size_t j=0; j<sx0v.vectorsCount(); ++j) {
              flops += rkernel_2vs_0p<StoreVec,AccumVec>(sx0v.vector(j), sy0v.vector(j), sz0v.vector(j),
                                               sx1v.vector(j), sy1v.vector(j), sz1v.vector(j),
                                               sx2v.vector(j), sy2v.vector(j), sz2v.vector(j),
                                               ssxv.vector(j), ssyv.vector(j), sszv.vector(j),
                                               StoreVec(0.0),
                                               txv, tyv, tzv,
                                               sav.vector(j), 0, RECURSIVE_LEVELS,
                                               &accumu, &accumv, &accumw);
            }
          }

          tu[0][i] += accumu.sum();
          tu[1][i] += accumv.sum();
          tu[2][i] += accumw.sum();
          flops += 3.0;
        }

      #else  // no Vc
        #pragma omp parallel for
        for (int32_t i=0; i<ntarg; ++i) {
          A accumu = 0.0;
          A accumv = 0.0;
          A accumw = 0.0;
          if (havess) {
            for (size_t j=0; j<src.get_npanels(); ++j) {
              const size_t jp0 = si[3*j];
              const size_t jp1 = si[3*j+1];
              const size_t jp2 = si[3*j+2];
              flops += rkernel_2vs_0p<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                                  sx[0][jp1], sx[1][jp1], sx[2][jp1],
                                  sx[0][jp2], sx[1][jp2], sx[2][jp2],
                                  ss[0][j]/sa[j], ss[1][j]/sa[j], ss[2][j]/sa[j], sss[j],
                                  tx[0][i], tx[1][i], tx[2][i],
                                  sa[j], 0, RECURSIVE_LEVELS,
                                  &accumu, &accumv, &accumw);
            }
          } else {
            for (size_t j=0; j<src.get_npanels(); ++j) {
              const size_t jp0 = si[3*j];
              const size_t jp1 = si[3*j+1];
              const size_t jp2 = si[3*j+2];
              flops += rkernel_2vs_0p<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                                  sx[0][jp1], sx[1][jp1], sx[2][jp1],
                                  sx[0][jp2], sx[1][jp2], sx[2][jp2],
                                  ss[0][j]/sa[j], ss[1][j]/sa[j], ss[2][j]/sa[j], S(0.0),
                                  tx[0][i], tx[1][i], tx[2][i],
                                  sa[j], 0, RECURSIVE_LEVELS,
                                  &accumu, &accumv, &accumw);
            }
          }
          tu[0][i] += accumu;
          tu[1][i] += accumv;
          tu[2][i] += accumw;
          flops += 3.0;
        }
      #endif // no Vc
    }

  //
  // targets are particles, with a core radius ===================================================
  //
  } else {

    // get the core radius
    //const Vector<S>&                            tr = targ.get_rad();

    if (opttug) { // velocity-and-grads kernel -------------------------------------------------
      if (havess) {
        std::cout << "    2vs_0bg compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
      } else {
        std::cout << "    2v_0bg compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
      }

      // get the pointer from the optional
      std::array<Vector<S>,9>& tug = *opttug;

      #ifdef USE_VC
        #pragma omp parallel for
        for (int32_t i=0; i<ntarg; ++i) {
          const StoreVec txv(tx[0][i]);
          const StoreVec tyv(tx[1][i]);
          const StoreVec tzv(tx[2][i]);
          //const StoreVec trv(tr[i]);
          AccumVec accumu(0.0);
          AccumVec accumv(0.0);
          AccumVec accumw(0.0);
          AccumVec accumux(0.0);
          AccumVec accumvx(0.0);
          AccumVec accumwx(0.0);
          AccumVec accumuy(0.0);
          AccumVec accumvy(0.0);
          AccumVec accumwy(0.0);
          AccumVec accumuz(0.0);
          AccumVec accumvz(0.0);
          AccumVec accumwz(0.0);

          if (havess) {
            for (size_t j=0; j<sx0v.vectorsCount(); ++j) {
              // NOTE: .vectorAt(i) gets the vector at scalar position i
              //       .vector(i) gets the i'th vector!!!
              flops += rkernel_2vs_0pg<StoreVec,AccumVec>(sx0v.vector(j), sy0v.vector(j), sz0v.vector(j),
                                                sx1v.vector(j), sy1v.vector(j), sz1v.vector(j),
                                                sx2v.vector(j), sy2v.vector(j), sz2v.vector(j),
                                                ssxv.vector(j), ssyv.vector(j), sszv.vector(j),
                                                sssv.vector(j),
                                                txv, tyv, tzv,
                                                sav.vector(j), 0, RECURSIVE_LEVELS,
                                                &accumu, &accumv, &accumw,
                                                &accumux, &accumvx, &accumwx,
                                                &accumuy, &accumvy, &accumwy,
                                                &accumuz, &accumvz, &accumwz);
            }
          } else {
            for (size_t j=0; j<sx0v.vectorsCount(); ++j) {
              flops += rkernel_2vs_0pg<StoreVec,AccumVec>(sx0v.vector(j), sy0v.vector(j), sz0v.vector(j),
                                                sx1v.vector(j), sy1v.vector(j), sz1v.vector(j),
                                                sx2v.vector(j), sy2v.vector(j), sz2v.vector(j),
                                                ssxv.vector(j), ssyv.vector(j), sszv.vector(j),
                                                StoreVec(0.0),
                                                txv, tyv, tzv,
                                                sav.vector(j), 0, RECURSIVE_LEVELS,
                                                &accumu, &accumv, &accumw,
                                                &accumux, &accumvx, &accumwx,
                                                &accumuy, &accumvy, &accumwy,
                                                &accumuz, &accumvz, &accumwz);
            }
          }

          tu[0][i] += accumu.sum();
          tu[1][i] += accumv.sum();
          tu[2][i] += accumw.sum();
          tug[0][i] += accumux.sum();
          tug[1][i] += accumvx.sum();
          tug[2][i] += accumwx.sum();
          tug[3][i] += accumuy.sum();
          tug[4][i] += accumvy.sum();
          tug[5][i] += accumwy.sum();
          tug[6][i] += accumuz.sum();
          tug[7][i] += accumvz.sum();
          tug[8][i] += accumwz.sum();
        }
        flops += 12.0*(float)ntarg;

      #else  // no Vc
        #pragma omp parallel for
        for (int32_t i=0; i<ntarg; ++i) {
          A accumu = 0.0;
          A accumv = 0.0;
          A accumw = 0.0;
          A accumux = 0.0;
          A accumvx = 0.0;
          A accumwx = 0.0;
          A accumuy = 0.0;
          A accumvy = 0.0;
          A accumwy = 0.0;
          A accumuz = 0.0;
          A accumvz = 0.0;
          A accumwz = 0.0;
          if (havess) {
            for (size_t j=0; j<src.get_npanels(); ++j) {
              const size_t jp0 = si[3*j];
              const size_t jp1 = si[3*j+1];
              const size_t jp2 = si[3*j+2];
              flops += rkernel_2vs_0pg<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                                   sx[0][jp1], sx[1][jp1], sx[2][jp1],
                                   sx[0][jp2], sx[1][jp2], sx[2][jp2],
                                   ss[0][j]/sa[j], ss[1][j]/sa[j], ss[2][j]/sa[j], sss[j],
                                   tx[0][i], tx[1][i], tx[2][i],
                                   sa[j], 0, RECURSIVE_LEVELS,
                                   &accumu, &accumv, &accumw,
                                   &accumux, &accumvx, &accumwx,
                                   &accumuy, &accumvy, &accumwy,
                                   &accumuz, &accumvz, &accumwz);
            }
          } else {
            for (size_t j=0; j<src.get_npanels(); ++j) {
              const size_t jp0 = si[3*j];
              const size_t jp1 = si[3*j+1];
              const size_t jp2 = si[3*j+2];
              flops += rkernel_2vs_0pg<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                                   sx[0][jp1], sx[1][jp1], sx[2][jp1],
                                   sx[0][jp2], sx[1][jp2], sx[2][jp2],
                                   ss[0][j]/sa[j], ss[1][j]/sa[j], ss[2][j]/sa[j], S(0.0),
                                   tx[0][i], tx[1][i], tx[2][i],
                                   sa[j], 0, RECURSIVE_LEVELS,
                                   &accumu, &accumv, &accumw,
                                   &accumux, &accumvx, &accumwx,
                                   &accumuy, &accumvy, &accumwy,
                                   &accumuz, &accumvz, &accumwz);
            }
          }
          tu[0][i] += accumu;
          tu[1][i] += accumv;
          tu[2][i] += accumw;
          tug[0][i] += accumux;
          tug[1][i] += accumvx;
          tug[2][i] += accumwx;
          tug[3][i] += accumuy;
          tug[4][i] += accumvy;
          tug[5][i] += accumwy;
          tug[6][i] += accumuz;
          tug[7][i] += accumvz;
          tug[8][i] += accumwz;
          flops += 12;
        }
      #endif // no Vc


    } else { // velocity-only kernel -----------------------------------------------------------
      if (havess) {
        std::cout << "    2vs_0b compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
      } else {
        std::cout << "    2v_0b compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
      }

      #ifdef USE_VC
        #pragma omp parallel for
        for (int32_t i=0; i<ntarg; ++i) {
          const StoreVec txv(tx[0][i]);
          const StoreVec tyv(tx[1][i]);
          const StoreVec tzv(tx[2][i]);
          //const StoreVec trv(tr[i]);
          AccumVec accumu(0.0);
          AccumVec accumv(0.0);
          AccumVec accumw(0.0);

          if (havess) {
            for (size_t j=0; j<sx0v.vectorsCount(); ++j) {
              // NOTE: .vectorAt(i) gets the vector at scalar position i
              //       .vector(i) gets the i'th vector!!!
              flops += rkernel_2vs_0p<StoreVec,AccumVec>(sx0v.vector(j), sy0v.vector(j), sz0v.vector(j),
                                               sx1v.vector(j), sy1v.vector(j), sz1v.vector(j),
                                               sx2v.vector(j), sy2v.vector(j), sz2v.vector(j),
                                               ssxv.vector(j), ssyv.vector(j), sszv.vector(j),
                                               sssv.vector(j),
                                               txv, tyv, tzv,
                                               sav.vector(j), 0, RECURSIVE_LEVELS,
                                               &accumu, &accumv, &accumw);
            }
          } else {
            for (size_t j=0; j<sx0v.vectorsCount(); ++j) {
              flops += rkernel_2vs_0p<StoreVec,AccumVec>(sx0v.vector(j), sy0v.vector(j), sz0v.vector(j),
                                               sx1v.vector(j), sy1v.vector(j), sz1v.vector(j),
                                               sx2v.vector(j), sy2v.vector(j), sz2v.vector(j),
                                               ssxv.vector(j), ssyv.vector(j), sszv.vector(j),
                                               StoreVec(0.0),
                                               txv, tyv, tzv,
                                               sav.vector(j), 0, RECURSIVE_LEVELS,
                                               &accumu, &accumv, &accumw);
            }
          }

          tu[0][i] += accumu.sum();
          tu[1][i] += accumv.sum();
          tu[2][i] += accumw.sum();
        }
        flops += 3.0*(float)ntarg;

      #else  // no Vc
        #pragma omp parallel for
        for (int32_t i=0; i<ntarg; ++i) {
          A accumu = 0.0;
          A accumv = 0.0;
          A accumw = 0.0;
          if (havess) {
            for (size_t j=0; j<src.get_npanels(); ++j) {
              const size_t jp0 = si[3*j];
              const size_t jp1 = si[3*j+1];
              const size_t jp2 = si[3*j+2];
              flops += rkernel_2vs_0p<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                                  sx[0][jp1], sx[1][jp1], sx[2][jp1],
                                  sx[0][jp2], sx[1][jp2], sx[2][jp2],
                                  ss[0][j]/sa[j], ss[1][j]/sa[j], ss[2][j]/sa[j], sss[j],
                                  tx[0][i], tx[1][i], tx[2][i],
                                  sa[j], 0, RECURSIVE_LEVELS,
                                  &accumu, &accumv, &accumw);
            }
          } else {
            for (size_t j=0; j<src.get_npanels(); ++j) {
              const size_t jp0 = si[3*j];
              const size_t jp1 = si[3*j+1];
              const size_t jp2 = si[3*j+2];
              flops += rkernel_2vs_0p<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                                  sx[0][jp1], sx[1][jp1], sx[2][jp1],
                                  sx[0][jp2], sx[1][jp2], sx[2][jp2],
                                  ss[0][j]/sa[j], ss[1][j]/sa[j], ss[2][j]/sa[j], S(0.0),
                                  tx[0][i], tx[1][i], tx[2][i],
                                  sa[j], 0, RECURSIVE_LEVELS,
                                  &accumu, &accumv, &accumw);
            }
          }
          tu[0][i] += accumu;
          tu[1][i] += accumv;
          tu[2][i] += accumw;
          flops += 3;
        }
      #endif // no Vc
    }
  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    panels_affect_points: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


//
// Vc and x86 versions of Points/Particles affecting Panels/Surfaces
// Should never need grads here - we only use it to calculate RHS
// And sources are never inert points, always active particles
//
template <class S, class A>
void points_affect_panels (Points<S> const& src, Surfaces<S>& targ) {
  std::cout << "    0v_2p compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const std::array<Vector<S>,Dimensions>& ss = src.get_str();

  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  const std::vector<Int>&                 ti = targ.get_idx();
  std::array<Vector<S>,Dimensions>&       tu = targ.get_vel();
  const Vector<S>&                        ta = targ.get_area();

  float flops = 0.0;

#ifdef USE_VC
  // define vector types for Vc (still only S==A supported here)
  typedef Vc::Vector<S> StoreVec;
  typedef Vc::SimdArray<A, Vc::Vector<S>::size()> AccumVec;

  // process source particles into Vc-ready memory format
  const Vc::Memory<StoreVec> sxv  = stdvec_to_vcvec<S>(sx[0], 0.0);
  const Vc::Memory<StoreVec> syv  = stdvec_to_vcvec<S>(sx[1], 0.0);
  const Vc::Memory<StoreVec> szv  = stdvec_to_vcvec<S>(sx[2], 0.0);
  //const Vc::Memory<StoreVec> srv  = stdvec_to_vcvec<S>(sr,    1.0);
  const Vc::Memory<StoreVec> ssxv = stdvec_to_vcvec<S>(ss[0], 0.0);
  const Vc::Memory<StoreVec> ssyv = stdvec_to_vcvec<S>(ss[1], 0.0);
  const Vc::Memory<StoreVec> sszv = stdvec_to_vcvec<S>(ss[2], 0.0);

  #pragma omp parallel for
  for (int32_t i=0; i<(int32_t)targ.get_npanels(); ++i) {

    // prepare vector registers for target accumulators
    const size_t ip0 = ti[3*i];
    const size_t ip1 = ti[3*i+1];
    const size_t ip2 = ti[3*i+2];
    const StoreVec tx0 = tx[0][ip0];
    const StoreVec ty0 = tx[1][ip0];
    const StoreVec tz0 = tx[2][ip0];
    const StoreVec tx1 = tx[0][ip1];
    const StoreVec ty1 = tx[1][ip1];
    const StoreVec tz1 = tx[2][ip1];
    const StoreVec tx2 = tx[0][ip2];
    const StoreVec ty2 = tx[1][ip2];
    const StoreVec tz2 = tx[2][ip2];
    const StoreVec tav = ta[i];

    AccumVec accumu = 0.0;
    AccumVec accumv = 0.0;
    AccumVec accumw = 0.0;

    //const size_t nsrc = src.get_n();
    //const size_t nsrcvec = 1 + (nsrc-1) / StoreVec::size();
    //for (size_t j=0; j<nsrcvec; j++) {

    for (size_t j=0; j<sxv.vectorsCount(); ++j) {
      // NOTE: .vectorAt(i) gets the vector at scalar position i
      //       .vector(i) gets the i'th vector!!!
      flops += rkernel_2vs_0p<StoreVec,AccumVec>(tx0, ty0, tz0,
                                      tx1, ty1, tz1,
                                      tx2, ty2, tz2,
                                      ssxv.vector(j)/tav, ssyv.vector(j)/tav, sszv.vector(j)/tav,
                                      StoreVec(0.0),
                                      sxv.vector(j), syv.vector(j), szv.vector(j),
                                      tav, 0, RECURSIVE_LEVELS,
                                      &accumu, &accumv, &accumw);
    }
    tu[0][i] -= accumu.sum();
    tu[1][i] -= accumv.sum();
    tu[2][i] -= accumw.sum();
  }
  flops += 3.0 * (float)targ.get_npanels();

#else  // no Vc
  #pragma omp parallel for
  for (int32_t i=0; i<(int32_t)targ.get_npanels(); ++i) {
    A accumu = 0.0;
    A accumv = 0.0;
    A accumw = 0.0;
    const size_t ip0 = ti[3*i];
    const size_t ip1 = ti[3*i+1];
    const size_t ip2 = ti[3*i+2];
    for (size_t j=0; j<src.get_n(); ++j) {
      // note that this is the same kernel as panels_affect_points!
      flops += rkernel_2vs_0p<S,A>(tx[0][ip0], tx[1][ip0], tx[2][ip0],
                                   tx[0][ip1], tx[1][ip1], tx[2][ip1],
                                   tx[0][ip2], tx[1][ip2], tx[2][ip2],
                                   ss[0][j]/ta[i], ss[1][j]/ta[i], ss[2][j]/ta[i],
                                   S(0.0),
                                   sx[0][j], sx[1][j], sx[2][j],
                                   ta[i], 0, RECURSIVE_LEVELS,
                                   &accumu, &accumv, &accumw);
    }
    // we use it backwards, so the resulting velocities are negative
    tu[0][i] -= accumu;
    tu[1][i] -= accumv;
    tu[2][i] -= accumw;
    flops += 3;
  }
#endif // no Vc

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    points_affect_panels: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


template <class S, class A>
void panels_affect_panels (Surfaces<S> const& src, Surfaces<S>& targ) {
  std::cout << "    2_2 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // run panels_affect_points instead

  // generate temporary colocation points as Points - is this inefficient?
  std::vector<S> xysr = targ.represent_as_particles(0.001, 0.001);
  Points<float> temppts(xysr, active, lagrangian, nullptr);

  // run the calculation
  panels_affect_points<S,A>(src, temppts);

  // and copy the velocities to the real target
  std::array<Vector<S>,Dimensions>& fromvel = temppts.get_vel();
  std::array<Vector<S>,Dimensions>& tovel   = targ.get_vel();
  for (size_t i=0; i<Dimensions; ++i) {
    std::copy(fromvel[i].begin(), fromvel[i].end(), tovel[i].begin());
  }
}


//
// helper struct for dispatching through a variant
//
template <class A>
struct InfluenceVisitor {
  // source collection, target collection
  void operator()(Points<float> const& src,   Points<float>& targ)   { points_affect_points<float,A>(src, targ); } 
  void operator()(Surfaces<float> const& src, Points<float>& targ)   { panels_affect_points<float,A>(src, targ); } 
  void operator()(Points<float> const& src,   Surfaces<float>& targ) { points_affect_panels<float,A>(src, targ); } 
  void operator()(Surfaces<float> const& src, Surfaces<float>& targ) { panels_affect_panels<float,A>(src, targ); } 
};

