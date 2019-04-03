/*
 * Influence.h - Non-class influence calculations
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega3D.h"
#include "VectorHelper.h"
#include "Kernels.h"
#include "Points.h"
#include "Surfaces.h"

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#include <chrono>
#include <cmath>


//
// Vc and x86 versions of Points/Particles affecting Points/Particles
//
template <class S, class A>
void points_affect_points (Points<S> const& src, Points<S>& targ) {
  auto start = std::chrono::system_clock::now();

  // is this where we dispatch the OpenGL compute shader?

  // get references to use locally
  const std::array<Vector<S>,Dimensions>&     sx = src.get_pos();
  const Vector<S>&                            sr = src.get_rad();
  const std::array<Vector<S>,Dimensions>&     ss = src.get_str();

  const std::array<Vector<S>,Dimensions>&     tx = targ.get_pos();
  std::array<Vector<S>,Dimensions>&           tu = targ.get_vel();
  std::optional<std::array<Vector<S>,9>>& opttug = targ.get_velgrad();

#ifdef USE_VC
  // create float_v versions of the source vectors
  const Vc::Memory<Vc::Vector<S>> sxv  = stdvec_to_vcvec<S>(sx[0], 0.0);
  const Vc::Memory<Vc::Vector<S>> syv  = stdvec_to_vcvec<S>(sx[1], 0.0);
  const Vc::Memory<Vc::Vector<S>> szv  = stdvec_to_vcvec<S>(sx[2], 0.0);
  const Vc::Memory<Vc::Vector<S>> srv  = stdvec_to_vcvec<S>(sr,    1.0);
  const Vc::Memory<Vc::Vector<S>> ssxv = stdvec_to_vcvec<S>(ss[0], 0.0);
  const Vc::Memory<Vc::Vector<S>> ssyv = stdvec_to_vcvec<S>(ss[1], 0.0);
  const Vc::Memory<Vc::Vector<S>> sszv = stdvec_to_vcvec<S>(ss[2], 0.0);
#endif

  float flops = (float)targ.get_n();

  // We need 8 different loops here, for the options:
  //   target radii or no target radii
  //   grads or no grads
  //   Vc or no Vc

  //
  // targets are field points, with no core radius
  //
  if (targ.is_inert()) {

  // dispatch on presence of val grads
  if (opttug) {
    std::cout << "    0v_0pg compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
    // get the pointer from the optional
    std::array<Vector<S>,9>& tug = *opttug;

    // velocity+grads kernel
    #pragma omp parallel for
    for (size_t i=0; i<targ.get_n(); ++i) {
#ifdef USE_VC
      const Vc::Vector<S> txv(tx[0][i]);
      const Vc::Vector<S> tyv(tx[1][i]);
      const Vc::Vector<S> tzv(tx[2][i]);
      // care must be taken if S != A, because these vectors must have the same length
      Vc::Vector<A> accumu(0.0);
      Vc::Vector<A> accumv(0.0);
      Vc::Vector<A> accumw(0.0);
      Vc::Vector<A> accumux = 0.0;
      Vc::Vector<A> accumvx = 0.0;
      Vc::Vector<A> accumwx = 0.0;
      Vc::Vector<A> accumuy = 0.0;
      Vc::Vector<A> accumvy = 0.0;
      Vc::Vector<A> accumwy = 0.0;
      Vc::Vector<A> accumuz = 0.0;
      Vc::Vector<A> accumvz = 0.0;
      Vc::Vector<A> accumwz = 0.0;
      // loop over source particles
      for (size_t j=0; j<sxv.vectorsCount(); ++j) {
        // NOTE: .vectorAt(i) gets the vector at scalar position i
        //       .vector(i) gets the i'th vector!!!
        kernel_0v_0pg<Vc::Vector<S>,Vc::Vector<A>>(
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
    flops *= 12.0 + 63.0*(float)src.get_n();

  } else {

    // velocity-only kernel
    std::cout << "    0v_0p compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
    #pragma omp parallel for
    for (size_t i=0; i<targ.get_n(); ++i) {
#ifdef USE_VC
      const Vc::Vector<S> txv = tx[0][i];
      const Vc::Vector<S> tyv = tx[1][i];
      const Vc::Vector<S> tzv = tx[2][i];
      // care must be taken if S != A, because these vectors must have the same length
      Vc::Vector<A> accumu = 0.0;
      Vc::Vector<A> accumv = 0.0;
      Vc::Vector<A> accumw = 0.0;
      for (size_t j=0; j<sxv.vectorsCount(); ++j) {
        kernel_0v_0p<Vc::Vector<S>,Vc::Vector<A>>(
                          sxv[j], syv[j], szv[j], srv[j],
                          ssxv[j], ssyv[j], sszv[j],
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
    flops *= 3.0 + 28.0*(float)src.get_n();
  }

  //
  // targets are particles, with a core radius
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
    for (size_t i=0; i<targ.get_n(); ++i) {
#ifdef USE_VC
      const Vc::Vector<S> txv(tx[0][i]);
      const Vc::Vector<S> tyv(tx[1][i]);
      const Vc::Vector<S> tzv(tx[2][i]);
      const Vc::Vector<S> trv(tr[i]);
      // care must be taken if S != A, because these vectors must have the same length
      Vc::Vector<A> accumu(0.0);
      Vc::Vector<A> accumv(0.0);
      Vc::Vector<A> accumw(0.0);
      Vc::Vector<A> accumux = 0.0;
      Vc::Vector<A> accumvx = 0.0;
      Vc::Vector<A> accumwx = 0.0;
      Vc::Vector<A> accumuy = 0.0;
      Vc::Vector<A> accumvy = 0.0;
      Vc::Vector<A> accumwy = 0.0;
      Vc::Vector<A> accumuz = 0.0;
      Vc::Vector<A> accumvz = 0.0;
      Vc::Vector<A> accumwz = 0.0;
      // loop over source particles
      for (size_t j=0; j<sxv.vectorsCount(); ++j) {
        // NOTE: .vectorAt(i) gets the vector at scalar position i
        //       .vector(i) gets the i'th vector!!!
        kernel_0_0sg<Vc::Vector<S>,Vc::Vector<A>>(
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
        kernel_0_0sg<S,A>(sx[0][j], sx[1][j], sx[2][j], sr[j],
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
    flops *= 12.0 + 65.0*(float)src.get_n();

  } else {

    // velocity-only kernel
    std::cout << "    0v_0v compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
    #pragma omp parallel for
    for (size_t i=0; i<targ.get_n(); ++i) {
#ifdef USE_VC
      const Vc::Vector<S> txv = tx[0][i];
      const Vc::Vector<S> tyv = tx[1][i];
      const Vc::Vector<S> tzv = tx[2][i];
      const Vc::Vector<S> trv = tr[i];
      // care must be taken if S != A, because these vectors must have the same length
      Vc::Vector<A> accumu = 0.0;
      Vc::Vector<A> accumv = 0.0;
      Vc::Vector<A> accumw = 0.0;
      for (size_t j=0; j<sxv.vectorsCount(); ++j) {
        kernel_0_0s<Vc::Vector<S>,Vc::Vector<A>>(
                         sxv[j], syv[j], szv[j], srv[j],
                         ssxv[j], ssyv[j], sszv[j],
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
        kernel_0_0s<S,A>(sx[0][j], sx[1][j], sx[2][j], sr[j],
                         ss[0][j], ss[1][j], ss[2][j],
                         tx[0][i], tx[1][i], tx[2][i], tr[i],
                         &accumu, &accumv, &accumw);
      }
      tu[0][i] += accumu;
      tu[1][i] += accumv;
      tu[2][i] += accumw;
#endif // no Vc
    }
    flops *= 3.0 + 30.0*(float)src.get_n();
  }

  //
  // end conditional over whether targets are field points (with no core radius)
  //
  }

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
  std::cout << "    1_0 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
  //auto start = std::chrono::system_clock::now();

  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  //const Vector<S>&                      sr = src.get_rad();
  const std::vector<Int>&                 si = src.get_idx();
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
      //kernel_1_0v<S,A>(&sx[2*si[2*j]], &sx[2*si[2*j+1]], ss[j],
      //                &tx[2*i], accum.data());
      kernel_1_0s<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                       sx[0][jp1], sx[1][jp1], sx[2][jp1],
                       ss[0][j], ss[1][j], ss[2][j],
                       tx[0][i], tx[1][i], tx[2][i],
                       //accum.data());
                       &accumu, &accumv, &accumw);
    }
    //tu[0][i] += accum[0];
    //tu[1][i] += accum[1];
    //tu[2][i] += accum[2];
    tu[0][i] += accumu;
    tu[1][i] += accumv;
    tu[2][i] += accumw;
  }
}


//
// Vc and x86 versions of Points/Particles affecting Panels/Surfaces
//
template <class S, class A>
void points_affect_panels (Points<S> const& src, Surfaces<S>& targ) {
  std::cout << "    0_1 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const std::array<Vector<S>,Dimensions>& ss = src.get_str();
  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  const std::vector<Int>&                 ti = targ.get_idx();
  std::array<Vector<S>,Dimensions>&       tu = targ.get_vel();

  float flops = (float)targ.get_n();

  #pragma omp parallel for
  for (size_t i=0; i<targ.get_n(); ++i) {
    //std::array<A,3> accum = {0.0};
    A accumu = 0.0;
    A accumv = 0.0;
    A accumw = 0.0;
    const size_t ip0 = ti[2*i];
    const size_t ip1 = ti[2*i+1];
    for (size_t j=0; j<src.get_n(); ++j) {
      // note that this is the same kernel as panels_affect_points!
      //kernel_1_0v<S,A>(&tx[2*ti[2*i]], &tx[2*ti[2*i+1]], ss[j],
      //                 &sx[2*j], accum.data());
      kernel_1_0s<S,A>(tx[0][ip0], tx[1][ip0], tx[2][ip0],
                       tx[0][ip1], tx[1][ip1], tx[2][ip1],
                       ss[0][j], ss[1][j], ss[2][j],
                       sx[0][j], sx[1][j], sx[2][j],
                       //accum.data());
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
  flops *= 11.0 + 0.0*(float)src.get_n();

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    points_affect_panels: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


template <class S, class A>
void panels_affect_panels (Surfaces<S> const& src, Surfaces<S>& targ) {
  std::cout << "    1_1 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
  // not sure how to do this - find field points of one and apply a function above?
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

