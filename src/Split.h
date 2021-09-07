/*
 * Split.h - split elongated particles according to local vorticity
 *
 * (c)2018 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Core.h"
#include "VectorHelper.h"
#include "MathHelper.h"
#include "nanoflann/nanoflann.hpp"

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>


//
// Find elongated particles and split them, maintaining 0,1 moments and approximating
// new second moment and new elongation
//
// templated on storage class S (typically float or double)
//
template <class S>
size_t split_elongated(Vector<S>& x, Vector<S>& y, Vector<S>& z,
                       Vector<S>& r, Vector<S>& elong,
                       Vector<S>& sx, Vector<S>& sy, Vector<S>& sz,
                       const CoreType corefunc,
                       const S particle_overlap,
                       const S threshold) {

  // start timer
  auto start = std::chrono::system_clock::now();

  // make sure all vector sizes are identical
  assert(x.size()==y.size() && "Input array sizes do not match");
  assert(x.size()==z.size() && "Input array sizes do not match");
  assert(x.size()==r.size() && "Input array sizes do not match");
  assert(x.size()==elong.size() && "Input array sizes do not match");
  assert(x.size()==sx.size() && "Input array sizes do not match");
  assert(x.size()==sy.size() && "Input array sizes do not match");
  assert(x.size()==sz.size() && "Input array sizes do not match");
  const size_t n = x.size();

  std::cout << "  Splitting elongated particles with n " << n << std::endl;

  // convert particle positions into something nanoflann can understand
  Eigen::Matrix<S, Eigen::Dynamic, 3> xp;
  xp.resize(n,3);
  xp.col(0) = Eigen::Map<Eigen::Matrix<S, Eigen::Dynamic, 1> >(x.data(), n);
  xp.col(1) = Eigen::Map<Eigen::Matrix<S, Eigen::Dynamic, 1> >(y.data(), n);
  xp.col(2) = Eigen::Map<Eigen::Matrix<S, Eigen::Dynamic, 1> >(z.data(), n);
  
  // generate the searchable data structure
  typedef typename Eigen::Matrix<S, Eigen::Dynamic, Dimensions> EigenMatType;
  typedef typename EigenMatType::Index EigenIndexType;
  typedef nanoflann::KDTreeEigenMatrixAdaptor< EigenMatType >  my_kd_tree_t;
  my_kd_tree_t mat_index(Dimensions, std::cref(xp));
  mat_index.index->buildIndex();

  std::vector<std::pair<EigenIndexType,S> > ret_matches;
  ret_matches.reserve(48);
  nanoflann::SearchParams params;
  params.sorted = true;

  //
  // split elongated particles, add new ones to end of list
  //
  
  // prepare vectors for new particles
  size_t num_split = 0;
  Vector<S> newx; newx.reserve(n/40);
  Vector<S> newy; newy.reserve(n/40);
  Vector<S> newz; newz.reserve(n/40);
  Vector<S> newr; newr.reserve(n/40);
  Vector<S> newe; newe.reserve(n/40);
  Vector<S> newsx; newsx.reserve(n/40);
  Vector<S> newsy; newsy.reserve(n/40);
  Vector<S> newsz; newsz.reserve(n/40);

  // track the largest elongation
  S maxElong = 0.0;
  size_t iElong = 0;

  // do we iterate over this algorithm?
  //for (size_t iters=0; iters<2; ++iters) {

  // now, for every elongated particle, split it into two particles
  for (size_t i=0; i<n; ++i) {

    // check vs. running max
    if (elong[i] > maxElong) {
      maxElong = elong[i];
      iElong = i;
    }

    if (elong[i] > threshold) {

    num_split++;

    // find elongation axis
    std::array<S,3> axis = {sx[i], sy[i], sz[i]};
    normalizeVec(axis);

    // find half-splitting distance
    const S nom_sep = r[i] / particle_overlap;
    const S halfd = 0.25 * elong[i] * nom_sep;

    // determine search radius
    const S search_rad = 3.0 * nom_sep;
    const S distsq_thresh = std::pow(search_rad, 2);

    // tree-based search with nanoflann
    const S query_pt[3] = { x[i], y[i], z[i] };
    (void) mat_index.index->radiusSearch(query_pt, distsq_thresh, ret_matches, params);
    //const size_t nMatches = mat_index.index->radiusSearch(query_pt, distsq_thresh, ret_matches, params);

    // search the new particle list, also
    //LATER
    // my old code says "no need to search np list, as the results won't change the calculation for local curvature"

    //std::cout << "    splitting " << x[i] << " " << y[i] << " " << z[i] << " with str " << sx[i] << " " << sy[i] << " " << sz[i] << " and elong " << elong[i] << std::endl;

    // prepare the accumulators
    S dVol = 0.0;
    std::array<S,3> curv = {0.0, 0.0, 0.0};

    // loop through all matches
    for (size_t j=0; j<ret_matches.size(); ++j) {
      const size_t iother = ret_matches[j].first;
      std::array<S,3> dx = {x[i]-x[iother], y[i]-y[iother], z[i]-z[iother]};
      std::array<S,3> s_other = {sx[iother], sy[iother], sz[iother]};

      // note that distance returned from radiusSearch is already squared
      //const S dist = std::sqrt(ret_matches[j].second);

      // compute this component
      const S factor = std::sqrt( ret_matches[j].second / (r[iother]*r[i])) * particle_overlap;
      const S fac2   = get_core_coefficient<S>(corefunc) * std::exp(-(std::pow(factor,3)));
      const S tStr   = length<S>(s_other);
      const S wgt    = dot_product<S>(axis, s_other);
      //std::cout << "      inear " << iother << " has " << tStr << " " << wgt << std::endl;

      // only add this influence if the strength is non-zero!
      if (tStr > 1.e-15 and wgt > 0.0) {
        // find the axis of the neighbor strength vector
        normalizeVec(s_other);

        const S dsx = axis[0] - sx[iother]/tStr;
        const S dsy = axis[1] - sy[iother]/tStr;
        const S dsz = axis[2] - sz[iother]/tStr;

        const S thisfac = wgt * fac2 * tStr;
        curv[0] += thisfac * (dx[1]*dsz - dx[2]*dsy);
        curv[1] += thisfac * (dx[2]*dsx - dx[0]*dsz);
        curv[2] += thisfac * (dx[0]*dsy - dx[1]*dsx);
        dVol += thisfac;
      }
    }

    // make a new child particle
    newx.resize(num_split);
    newy.resize(num_split);
    newz.resize(num_split);
    newr.resize(num_split);
    newe.resize(num_split);
    newsx.resize(num_split);
    newsy.resize(num_split);
    newsz.resize(num_split);

    newx[num_split-1] = x[i] - halfd*axis[0];
    newy[num_split-1] = y[i] - halfd*axis[1];
    newz[num_split-1] = z[i] - halfd*axis[2];
    newr[num_split-1] = r[i];
    newe[num_split-1] = 0.5 * elong[i];
    newsx[num_split-1] = 0.5 * sx[i];
    newsy[num_split-1] = 0.5 * sy[i];
    newsz[num_split-1] = 0.5 * sz[i];

    // reset the split particle
    x[i] += halfd*axis[0];
    y[i] += halfd*axis[1];
    z[i] += halfd*axis[2];
    elong[i] *= 0.5;
    sx[i] *= 0.5;
    sy[i] *= 0.5;
    sz[i] *= 0.5;

    // only include curvature if dVol is non-zero
    if (dVol > 1e-15) {

      // scale curvature by weight
      curv[0] /= dVol;
      curv[1] /= dVol;
      curv[2] /= dVol;

      // this is a correction factor (don't ask)
      const S factor = particle_overlap * 0.398967 / r[i];
      curv[0] *= factor;
      curv[1] *= factor;
      curv[2] *= factor;

      // calculate extra curvature component
      const S dsx = sy[i]*curv[2] - sz[i]*curv[1];
      const S dsy = sz[i]*curv[0] - sx[i]*curv[2];
      const S dsz = sx[i]*curv[1] - sy[i]*curv[0];

      // add it to the two particles
      // first, new child
      newsx[num_split-1] += dsx;
      newsy[num_split-1] += dsy;
      newsz[num_split-1] += dsz;

      // then, new original
      sx[i] -= dsx;
      sy[i] -= dsy;
      sz[i] -= dsz;
    }

    //std::cout << "      into " << newx[num_split-1] << " " << newy[num_split-1] << " " << newz[num_split-1] << " str " << newsx[num_split-1] << " " << newsy[num_split-1] << " " << newsz[num_split-1] << " and elong " << newe[num_split-1] << std::endl;
    //std::cout << "      and " << x[i] << " " << y[i] << " " << z[i] << " str " << sx[i] << " " << sy[i] << " " << sz[i] << " and elong " << elong[i] << std::endl;

    } // end if elong>thresh
  } // end for

  // merge the original and new arrays
  if (num_split > 0) {

    // now march through the arrays and concatenate them
    x.insert(x.end(), newx.begin(), newx.end());
    y.insert(y.end(), newy.begin(), newy.end());
    z.insert(z.end(), newz.begin(), newz.end());
    r.insert(r.end(), newr.begin(), newr.end());
    elong.insert(elong.end(), newe.begin(), newe.end());
    sx.insert(sx.end(), newsx.begin(), newsx.end());
    sy.insert(sy.end(), newsy.begin(), newsy.end());
    sz.insert(sz.end(), newsz.begin(), newsz.end());

    std::cout << "    split added " << num_split << " particles" << std::endl;
  }
  std::cout << "    max elong was " << maxElong << " at particle " << iElong << std::endl;

  // finish timer and report
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    split time:\t[%.4f] seconds\n", (float)elapsed_seconds.count());

  return num_split;
}


//
// split any collections of Points due to stretching
//
template <class S>
void split_operation(std::vector<Collection>& _vort,
                     const CoreType _corefunc,
                     const S _overlap,
                     const S _threshold) {

  for (auto &coll: _vort) {

    // but only check particles ("Points")
    if (std::holds_alternative<Points<S>>(coll)) {

      Points<S>& pts = std::get<Points<S>>(coll);
      //std::cout << "    check split for " << pts.get_n() << " particles" << std::endl;
      //std::cout << std::endl;

      // none of these are passed as const, because both may be extended with new particles
      std::array<Vector<S>,Dimensions>&       x = pts.get_pos();
      Vector<S>&                              r = pts.get_rad();
      Vector<S>&                              elong = pts.get_elong();
      std::array<Vector<S>, numStrenPerNode>& s = pts.get_str();

      // last two arguments are: relative distance, allow variable core radii
      (void)split_elongated<S>(x[0], x[1], x[2], r, elong, s[0], s[1], s[2],
                               _corefunc,
                               _overlap,
                               _threshold);

      // we probably have a different number of particles now, resize the u, ug, elong arrays
      pts.resize(r.size());
    }
  }
}
