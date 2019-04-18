/*
 * Merge.h - library code for a three-dimensional particle merging scheme
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega3D.h"
#include "VectorHelper.h"
#include "nanoflann.hpp"

#include <Eigen/Dense>

#include <array>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>


//
// Find close particles and merge them, maintaining 0,1 moments and approximating
// new second moment
//
// templated on storage class S (typically float or double)
//
template <class S>
size_t merge_close_particles(std::array<Vector<S>,3>& pos,
                             std::array<Vector<S>,3>& str,
                             Vector<S>&               rad,
                             const S                  particle_overlap,
                             const S                  threshold,
                             const bool               adapt_radii) {

  // make sure all vector sizes are identical
  assert(pos[0].size()==pos[1].size() && "Input array sizes do not match");
  assert(pos[0].size()==pos[2].size() && "Input array sizes do not match");
  assert(pos[0].size()==str[0].size() && "Input array sizes do not match");
  assert(pos[0].size()==str[1].size() && "Input array sizes do not match");
  assert(pos[0].size()==str[2].size() && "Input array sizes do not match");
  assert(pos[0].size()==rad.size() && "Input array sizes do not match");
  const size_t n = rad.size();

  std::cout << "  Merging close particles with n " << n << std::endl;

  // start timer
  auto start = std::chrono::system_clock::now();

  // reference or generate the local set of vectors
  Vector<S>& x = pos[0];
  Vector<S>& y = pos[1];
  Vector<S>& z = pos[2];
  Vector<S>& r = rad;
  Vector<S>& sx = str[0];
  Vector<S>& sy = str[1];
  Vector<S>& sz = str[2];

  // convert particle positions into something nanoflann can understand
  Eigen::Matrix<S, Eigen::Dynamic, Dimensions> xp;
  xp.resize(n,Dimensions);
  xp.col(0) = Eigen::Map<Eigen::Matrix<S, Eigen::Dynamic, 1> >(x.data(), n);
  xp.col(1) = Eigen::Map<Eigen::Matrix<S, Eigen::Dynamic, 1> >(y.data(), n);
  xp.col(2) = Eigen::Map<Eigen::Matrix<S, Eigen::Dynamic, 1> >(z.data(), n);
  
  // generate the searchable data structure
  typedef nanoflann::KDTreeEigenMatrixAdaptor< Eigen::Matrix<S, Eigen::Dynamic, Dimensions> >  my_kd_tree_t;
  my_kd_tree_t mat_index(xp, 20);
  mat_index.index->buildIndex();
  std::vector<std::pair<long int,S> > ret_matches;
  ret_matches.reserve(48);
  nanoflann::SearchParams params;
  params.sorted = true;

  //
  // merge co-located particles with identical radii
  //
  
  // prepare a vector of particles to remove
  std::vector<bool> erase_me;
  erase_me.resize(n);
  std::fill(erase_me.begin(), erase_me.end(), false);

  // now, for every particle, search for a co-located and identical-radius particle!
  for (size_t i=0; i<n; ++i) {

    // nominal separation for this particle
    const S nom_sep = r[i] / particle_overlap;
    const S search_rad = nom_sep * threshold;
    const S distsq_thresh = std::pow(search_rad, 2);

    // tree-based search with nanoflann
    const S query_pt[Dimensions] = { x[i], y[i], z[i] };
    const size_t nMatches = mat_index.index->radiusSearch(query_pt, distsq_thresh, ret_matches, params);

    // one match should be self, but we don't know which one (really?)
    if (nMatches > 1 and false)  {
      std::cout << "part " << i << std::endl;
      for (size_t j=0; j<ret_matches.size(); ++j) {
        std::cout << "\t" << j << "\t" << std::sqrt(ret_matches[j].second)/r[i] << std::endl;
      }
    }

    // if there are more than one, check the radii
    if (nMatches > 1) {
      for (size_t j=0; j<ret_matches.size(); ++j) {
        const size_t iother = ret_matches[j].first;
        if (i != iother) {
          // make sure distance is also less than target particle's threshold
          // note that distance returned from radiusSearch is already squared
          if (std::sqrt(ret_matches[j].second) < threshold*r[iother]/particle_overlap) {
            if (erase_me[i] or erase_me[iother]) {
              // we've already account for this one
            } else {
              //std::cout << "  particles " << i << " and " << iother << " will merge" << std::endl;
              //std::cout << "    first at " << x[i] << " " << y[i] << " " << z[i] << " with str " << sx[i] << " and rad " << r[i] << std::endl;
              //std::cout << "    other at " << x[iother] << " " << y[iother] << " " << z[iother] << " with str " << sx[iother] << " and rad " << r[iother] << std::endl;
              const S si = std::sqrt(sx[i]*sx[i] + sy[i]*sy[i] + sz[i]*sz[i]) + std::numeric_limits<S>::epsilon();
              const S so = std::sqrt(sx[iother]*sx[iother] + sy[iother]*sy[iother] + sz[iother]*sz[iother]) + std::numeric_limits<S>::epsilon();
              const S strength_mag = si + so;
              // find center of strength
              const S newx = (x[i]*si + x[iother]*so) / strength_mag;
              const S newy = (y[i]*si + y[iother]*so) / strength_mag;
              const S newz = (z[i]*si + z[iother]*so) / strength_mag;
              // move strengths to particle i
              x[i] = newx;
              y[i] = newy;
              z[i] = newz;
              if (adapt_radii) {
                r[i] = std::cbrt((si*std::pow(r[i],3) + so*std::pow(r[iother],3))/strength_mag);
              }
              sx[i] = sx[i] + sx[iother];
              sy[i] = sy[i] + sy[iother];
              sz[i] = sz[i] + sz[iother];
              //std::cout << "    result   " << x[i] << " " << y[i] << " with str " << sx[i] << " " << sy[i] << " " << sz[i] << " and rad " << r[i] << std::endl;
              // flag other particle for deletion
              erase_me[iother] = true;
            }
          }
        }
      }
    }
  }

  // how many to be erased?
  const size_t num_removed = std::count (erase_me.begin(), erase_me.end(), true);

  if (num_removed > 0) {

    // now march through the arrays and compress them
    size_t copyto = 0;
    for (size_t i=0; i<n; ++i) {
      if (not erase_me[i]) {
        if (i != copyto) {
          x[copyto] = x[i];
          y[copyto] = y[i];
          z[copyto] = z[i];
          r[copyto] = r[i];
          sx[copyto] = sx[i];
          sy[copyto] = sy[i];
          sz[copyto] = sz[i];
        }
        copyto++;
      }
    }
    // reset n and resize the arrays
    const size_t new_n = copyto;
    x.resize(new_n);
    y.resize(new_n);
    z.resize(new_n);
    r.resize(new_n);
    sx.resize(new_n);
    sy.resize(new_n);
    sz.resize(new_n);

    std::cout << "    merge removed " << num_removed << " particles" << std::endl;
  }

  // finish timer and report
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    merging time:\t[%.4f] seconds\n", (float)elapsed_seconds.count());

  return num_removed;
}


//
// run some number of merge ops on each collection of Point vorts
//
template <class S>
void merge_operation(std::vector<Collection>& _vort,
                     const S                  _overlap,
                     const S                  _thresh,
                     const bool               _isadapt) {

  const size_t maxiters = 1;

  for (auto &coll : _vort) {

    // if inert, no need to merge (keep all tracer particles)
    if (std::visit([=](auto& elem) { return elem.is_inert(); }, coll)) continue;

    // run this step if the collection is Points
    if (std::holds_alternative<Points<S>>(coll)) {

      Points<S>& pts = std::get<Points<S>>(coll);

      //std::cout << "    merging among " << pts.get_n() << " particles" << std::endl;

      // perform possibly multiple iterations
      for (size_t iter=0; iter<maxiters; ++iter) {

        // last two arguments are: relative distance, allow variable core radii
        (void) merge_close_particles(pts.get_pos(),
                                     pts.get_str(),
                                     pts.get_rad(),
                                     _overlap,
                                     _thresh,
                                     _isadapt);

        // resize the rest of the arrays
        pts.resize(pts.get_rad().size());
      }
    }
  }
}
