/*
 * PSE.h - Library code for a three-dimensional Particle-Strength-Exchange scheme
 *
 * (c)2018-21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
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

#include "Core.h"
#include "VectorHelper.h"

#include "nanoflann/nanoflann.hpp"
#include "json/json.hpp"

#include <Eigen/Dense>

#include <algorithm>	// for std::max
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>


//
// Class to hold PSE parameters and temporaries
//
// templatized on storage type ST and compute type CT
//
template <class ST, class CT>
class PSE {
public:
  PSE();

  void set_relative(const bool _in) { thresholds_are_relative = _in; }
  void set_ignore(const float _in) { ignore_thresh = _in; }
  void set_volumes(const bool _in) { use_volumes = _in; }

  const bool get_relative() const { return thresholds_are_relative; }
  const float get_ignore() const { return ignore_thresh; }
  const bool get_volumes() const { return use_volumes; }

  // all-to-all diffuse; can change array sizes
  void diffuse_all(std::array<Vector<ST>,3>&,
                   std::array<Vector<ST>,3>&,
                   Vector<ST>&,
                   const ST,
                   const CoreType,
                   const ST);

  void from_json(const nlohmann::json);
  void add_to_json(nlohmann::json&) const;

protected:
  // search for new target location
  std::array<ST,3> fill_neighborhood_search(const int32_t,
                                            const Vector<ST>&,
                                            const Vector<ST>&,
                                            const Vector<ST>&,
                                            const std::vector<int32_t>&,
                                            const ST);

  // add new particles to boundary
  size_t add_new_boundary_parts(Vector<ST>&,
                                Vector<ST>&,
                                Vector<ST>&,
                                Vector<ST>&,
                                Vector<ST>&,
                                Vector<ST>&,
                                Vector<ST>&,
                                const ST);

  // find particle volumes
  //void find_particle_volumes(Vector<ST>&,
  //                           Vector<ST>&, Vector<ST>&,
  //                           Vector<ST>&,
  //                           Vector<ST>&,
  //                           Vector<ST>&,
  //                           const CoreType,
  //                           const ST);

private:

  // useful constants
  static const int32_t minNearby = 14;
  static const int32_t max_near = 256;

  // new point insertion sites (normalized to nominal separation and centered around origin)
  size_t num_sites;
  std::vector<ST> xsite,ysite,zsite;
  void initialize_sites();

  // do not perform PSE if source particle strength is less than
  //   this fraction of max particle strength
  ST ignore_thresh = 1.e-4;
  // are thresholds absolute or relative to strongest particle?
  bool thresholds_are_relative = true;

  // calculate and use particle volumes? (increases accuracy)
  bool use_volumes = true;

  // use nanoflann for nearest-neighbor searching? false uses direct search
  const bool use_tree = true;
};

// primary constructor
template <class ST, class CT>
PSE<ST,CT>::PSE() {
  initialize_sites();
}

//
// initialize the particle placement sites array
//
template <class ST, class CT>
void PSE<ST,CT>::initialize_sites() {

  const std::vector<ST>& ico = ico2;

  // for all PSE - one layer of test positions
  if (true) {
    // pull data from a multi-level refinement of an icosahedron
    num_sites = ico.size() / Dimensions;
    xsite.resize(num_sites);
    ysite.resize(num_sites);
    zsite.resize(num_sites);
    for (size_t i=0; i<num_sites; ++i) {
      xsite[i] = 2.0 * ico[3*i+0];
      ysite[i] = 2.0 * ico[3*i+1];
      zsite[i] = 2.0 * ico[3*i+2];
    }
  }
}

//
// Using a strength heuristic, add new particles to the boundary of the
// particle distribution; this is very similar to what VRM does
//
// templated on storage class ST (typically float or double)
//
template <class ST, class CT>
size_t PSE<ST,CT>::add_new_boundary_parts(Vector<ST>& x,
                                          Vector<ST>& y,
                                          Vector<ST>& z,
                                          Vector<ST>& r,
                                          Vector<ST>& sx,
                                          Vector<ST>& sy,
                                          Vector<ST>& sz,
                                          const ST particle_overlap) {

  // start timer
  auto start = std::chrono::system_clock::now();

  // make sure all vector sizes are identical
  assert(x.size()==y.size());
  assert(x.size()==z.size());
  assert(x.size()==r.size());
  assert(x.size()==sx.size());
  assert(x.size()==sy.size());
  assert(x.size()==sz.size());
  size_t n = x.size();

  std::cout << "  Adding buffer particles with n " << n << std::endl;

  // convert particle positions into something nanoflann can understand
  Eigen::Matrix<ST, Eigen::Dynamic, Dimensions> xp;
  xp.resize(n,Dimensions);
  xp.col(0) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(x.data(), n);
  xp.col(1) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(y.data(), n);
  xp.col(2) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(z.data(), n);
  
  // generate the searchable data structure
  typedef typename Eigen::Matrix<ST, Eigen::Dynamic, 3> EigenMatType;
  typedef typename EigenMatType::Index EigenIndexType;
  typedef nanoflann::KDTreeEigenMatrixAdaptor< EigenMatType >  my_kd_tree_t;
  my_kd_tree_t mat_index(Dimensions, std::cref(xp));
  if (use_tree) mat_index.index->buildIndex();

  std::vector<std::pair<EigenIndexType,ST> > ret_matches;
  ret_matches.reserve(max_near);
  nanoflann::SearchParams params;
  params.sorted = true;

  // what is maximum strength of all particles?
  ST maxAbsStr = -1.0;
  for (size_t i=0; i<n; ++i) {
    maxAbsStr = std::max(maxAbsStr, sx[i]*sx[i] + sy[i]*sy[i] + sz[i]*sz[i]);
  }
  maxAbsStr = std::sqrt(maxAbsStr);
  std::cout << "    maxAbsStr " << maxAbsStr << std::endl;

  //
  // check each strong-enough particle for an appropriate number of neighbors
  //
  const int32_t initial_n = n;
  for (int32_t i=0; i<initial_n; ++i) {

    const ST thisstr = std::sqrt(sx[i]*sx[i] + sy[i]*sy[i] + sz[i]*sz[i]);

    // if current particle strength is very small, skip out
    if ((thresholds_are_relative && (thisstr < maxAbsStr * ignore_thresh)) or
        (!thresholds_are_relative && (thisstr < ignore_thresh))) continue;

    // nominal separation for this particle (insertion distance)
    const ST nom_sep = r[i] / particle_overlap;

    // what is search radius? (vdelta)
    //const ST search_rad = nom_sep * 1.6;
    const ST search_rad = r[i];

    // initialize vector of indexes of nearest particles
    std::vector<int32_t> inear;

    // switch on search method
    if (use_tree) {
      // tree-based search with nanoflann
      const ST distsq_thresh = std::pow(search_rad, 2);
      const ST query_pt[3] = { x[i], y[i], z[i] };
      (void) mat_index.index->radiusSearch(query_pt, distsq_thresh, ret_matches, params);
      //if (ret_matches.size() > 20) std::cout << "part " << i << " at " << x[i] << " " << y[i] << " " << z[i] << " has " << ret_matches.size() << " matches" << std::endl;

      // copy the indexes into my vector
      for (size_t j=0; j<ret_matches.size(); ++j) inear.push_back((int32_t)ret_matches[j].first);

      // now direct search over all newer particles
      for (size_t j=initial_n; j<n; ++j) {
        const ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2) + std::pow(z[i]-z[j], 2);
        if (distsq < distsq_thresh) inear.push_back((int32_t)j);
      }

    } else {
      // direct search: look for all neighboring particles, include newly-created ones
      //   ideally this would be a tree search
      const ST distsq_thresh = std::pow(search_rad, 2);
      for (size_t j=0; j<n; ++j) {
        const ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2) + std::pow(z[i]-z[j], 2);
        if (distsq < distsq_thresh) inear.push_back((int32_t)j);
      }
    }

    //std::cout << "    found " << inear.size() << " within " << search_rad << std::endl;

    // if there are less than, say, 6, we should just add some now
    while (inear.size() < minNearby) {
      auto newpt = fill_neighborhood_search(i, x, y, z, inear, nom_sep);
      // add the index to the near list
      inear.push_back(n);
      // add it to the master list
      x.push_back(newpt[0]);
      y.push_back(newpt[1]);
      z.push_back(newpt[2]);
      r.push_back(r[i]);
      sx.push_back(0.0);
      sy.push_back(0.0);
      sz.push_back(0.0);
      n++;
    }
  }

  std::cout << "    added " << (n - initial_n) << " buffer particles" << std::endl;

  // finish timer and report
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    pse.add_new_boundary_parts:\t[%.4f] cpu seconds\n", (float)elapsed_seconds.count());

  return (n - initial_n);
}


//
// use a ring of sites to determine the location of a new particle
//
template <class ST, class CT>
std::array<ST,3> PSE<ST,CT>::fill_neighborhood_search(const int32_t idx,
                                                      const Vector<ST>& x,
                                                      const Vector<ST>& y,
                                                      const Vector<ST>& z,
                                                      const std::vector<int32_t>& inear,
                                                      const ST nom_sep) {

  // create array of potential sites
  std::vector<ST> tx,ty,tz,nearest;
  tx.resize(num_sites);
  ty.resize(num_sites);
  tz.resize(num_sites);
  nearest.resize(num_sites);
  for (size_t i=0; i<num_sites; ++i) {
    tx[i] = x[idx] + nom_sep * xsite[i];
    ty[i] = y[idx] + nom_sep * ysite[i];
    tz[i] = z[idx] + nom_sep * zsite[i];
  }

  // test all points vs. all sites
  size_t iopen = 0;
  ST maxmindist = 0.0;
  for (size_t i=0; i<num_sites; ++i) {
    // find the nearest particle to this site
    ST mindistsq = nom_sep * nom_sep;
    for (size_t j=0; j<inear.size(); ++j) {
      const ST distsq = std::pow(x[inear[j]]-tx[i], 2) + std::pow(y[inear[j]]-ty[i], 2) + std::pow(z[inear[j]]-tz[i], 2);
      if (distsq < mindistsq) mindistsq = distsq;
    }
    nearest[i] = mindistsq;
    if (mindistsq > maxmindist) {
      maxmindist = mindistsq;
      iopen = i;
    }
  }

  //std::cout << "    creating particle at " << tx[iopen] << " " << ty[iopen] << " " << tz[iopen]
  //          << ", now there are " << n << " particles" << std::endl;

  const std::array<ST,3> retval = {{tx[iopen], ty[iopen], tz[iopen]}};
  return retval;
}


//
// This is the Particle-Strength-Exchange scheme - simulate viscous
// diffusion by exchanging strengths between neighboring particles
// Find the change in strength and radius that would occur over one dt
//
template <class ST, class CT>
void PSE<ST,CT>::diffuse_all(std::array<Vector<ST>,3>& pos,
                             std::array<Vector<ST>,3>& str,
                             Vector<ST>& rad,
                             const ST h_nu,
                             const CoreType core_func,
                             const ST particle_overlap) {

  // make sure all vector sizes are identical
  assert(pos[0].size()==pos[1].size() && "Input arrays are not uniform size");
  assert(pos[0].size()==pos[2].size() && "Input arrays are not uniform size");
  assert(pos[0].size()==str[0].size() && "Input arrays are not uniform size");
  assert(pos[0].size()==str[1].size() && "Input arrays are not uniform size");
  assert(pos[0].size()==str[2].size() && "Input arrays are not uniform size");
  assert(pos[0].size()==rad.size() && "Input arrays are not uniform size");

  //std::cout << "  Running PSE with n " << n << std::endl;

  // start timer
  auto start = std::chrono::system_clock::now();

  // reference or generate the local set of vectors
  Vector<ST>& x = pos[0];
  Vector<ST>& y = pos[1];
  Vector<ST>& z = pos[2];
  Vector<ST>& r = rad;
  Vector<ST>& sx = str[0];
  Vector<ST>& sy = str[1];
  Vector<ST>& sz = str[2];


  //
  // first step is to add buffer particles where we need them
  //
  (void) add_new_boundary_parts(x,y,z,r,sx,sy,sz,particle_overlap);
  size_t n = x.size();

  // generate and zero out delta vector
  Vector<ST> dsx, dsy, dsz;
  dsx.resize(n);
  dsy.resize(n);
  dsz.resize(n);

  // zero out delta vector
  std::fill(dsx.begin(), dsx.end(), 0.0);
  std::fill(dsy.begin(), dsy.end(), 0.0);
  std::fill(dsz.begin(), dsz.end(), 0.0);

  // convert particle positions into something nanoflann can understand
  Eigen::Matrix<ST, Eigen::Dynamic, Dimensions> xp;
  xp.resize(n,Dimensions);
  xp.col(0) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(x.data(), n);
  xp.col(1) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(y.data(), n);
  xp.col(2) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(z.data(), n);
  
  // generate the searchable data structure
  typedef typename Eigen::Matrix<ST, Eigen::Dynamic, 3> EigenMatType;
  typedef typename EigenMatType::Index EigenIndexType;
  typedef nanoflann::KDTreeEigenMatrixAdaptor< EigenMatType >  my_kd_tree_t;
  my_kd_tree_t mat_index(Dimensions, std::cref(xp));
  if (use_tree) mat_index.index->buildIndex();

  std::vector<std::pair<EigenIndexType,ST> > ret_matches;
  ret_matches.reserve(max_near);
  nanoflann::SearchParams params;
  params.sorted = true;


  //
  // second step (optional) is to calculate particle volumes
  //

  Vector<float> vol;
  vol.resize(n);

  if (use_volumes) {
    std::cout << "  Find volumes with n " << n << std::endl;

    // loop over every particle
    for (size_t i=0; i<n; ++i) {

      // don't need to go out as far as for the derivative
      ST search_rad = r[i] * 2.4;
      if (core_func == CoreType::compactg) search_rad = r[i] * 2.0;

      // initialize vector of indexes of nearest particles
      std::vector<int32_t> inear;

      // switch on search method
      if (use_tree) {
        // tree-based search with nanoflann
        const ST distsq_thresh = std::pow(search_rad, 2);
        const ST query_pt[3] = { x[i], y[i], z[i] };
        (void) mat_index.index->radiusSearch(query_pt, distsq_thresh, ret_matches, params);

        // copy the indexes into my vector
        for (size_t j=0; j<ret_matches.size(); ++j) inear.push_back((int32_t)ret_matches[j].first);

      } else {
        // direct search: look for all neighboring particles
        const ST distsq_thresh = std::pow(search_rad, 2);
        for (size_t j=0; j<n; ++j) {
          const ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2) + std::pow(z[i]-z[j], 2);
          if (distsq < distsq_thresh) inear.push_back((int32_t)j);
        }
      }

      // zero - just to be sure
      vol[i] = 0.0;

      // the constant factor
      const CT factor = get_core_coefficient<CT>(core_func);

      // split on particle core function
      if (core_func == CoreType::gaussian) {

        // loop over all participating particles
        for (size_t j=0; j<inear.size(); ++j) {
          const int32_t jdx = inear[j];
          // all distances are unnormalized
          const CT dx = x[i] - x[jdx];
          const CT dy = y[i] - y[jdx];
          const CT dz = z[i] - z[jdx];
          // for symmetry, use product of radii
          const CT distsq = (dx*dx + dy*dy + dz*dz) / (r[i]*r[jdx]);
          // compute the volume component
          vol[i] += factor * std::exp(-distsq);
        }

        // scale the volume by the proper constant factor
        vol[i] *= 1.0 / std::pow(r[i], 3);

      } else /* compact gaussian */ {

        // loop over all participating particles
        for (size_t j=0; j<inear.size(); ++j) {
          const int32_t jdx = inear[j];
          // all distances are unnormalized
          const CT dx = x[i] - x[jdx];
          const CT dy = y[i] - y[jdx];
          const CT dz = z[i] - z[jdx];
          // for symmetry, use product of radii
          const CT distsq = (dx*dx + dy*dy + dz*dz) / (r[i]*r[jdx]);
          const CT dist = std::sqrt(distsq);
          const CT distcub = distsq * dist;
          // compute the volume component
          vol[i] += factor * std::exp(-distcub);
        }

        // scale the volume by the proper factor
        vol[i] *= 1.0 / std::pow(r[i], 3);
      }

      // invert weight to get particle volume
      vol[i] = 1.0 / vol[i];

      //std::cout << "particle " << i << " has volume " << vol[i] << std::endl;
    } // end loop over all particles
  }


  //
  // finally, perform the PSE operation
  //

  std::cout << "  Running PSE with n " << n << std::endl;

  // for each particle (can parallelize this part)
  size_t nneibs = 0;
  size_t minneibs = 999999;
  size_t maxneibs = 0;
  for (size_t i=0; i<n; ++i) {

    // find the nearest neighbor particles
    //std::cout << "\nDiffusing particle " << i << " with strength " << s[i] << std::endl;

    // do not check vs. threshold - just do all of them
    //   (this particle still core-spreads somewhat)
    //if ((thresholds_are_relative && (std::abs(s[i]) < maxAbsStr * ignore_thresh)) or
    //    (!thresholds_are_relative && (std::abs(s[i]) < ignore_thresh))) continue;

    // nominal separation for this particle (insertion distance)
    const ST nom_sep = r[i] / particle_overlap;

    // must go out to 2-3 core radii to get this right
    // these are set to introduce no more than 1% error vs. very large search radii
    ST search_rad = r[i] * 2.65;
    if (core_func == CoreType::compactg) search_rad = r[i] * 2.0;

    // initialize vector of indexes of nearest particles
    std::vector<int32_t> inear;

    // switch on search method
    if (use_tree) {
      // tree-based search with nanoflann
      const ST distsq_thresh = std::pow(search_rad, 2);
      const ST query_pt[3] = { x[i], y[i], z[i] };
      (void) mat_index.index->radiusSearch(query_pt, distsq_thresh, ret_matches, params);

      // copy the indexes into my vector
      for (size_t j=0; j<ret_matches.size(); ++j) inear.push_back((int32_t)ret_matches[j].first);

    } else {
      // direct search: look for all neighboring particles
      const ST distsq_thresh = std::pow(search_rad, 2);
      for (size_t j=0; j<n; ++j) {
        const ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2) + std::pow(z[i]-z[j], 2);
        if (distsq < distsq_thresh) inear.push_back((int32_t)j);
      }

    }

    // core of the PSE algorithm is here

    // zero - just to be sure
    dsx[i] = 0.0;
    dsy[i] = 0.0;
    dsz[i] = 0.0;

    // split on particle core function
    if (core_func == CoreType::gaussian) {

      // loop over all participating particles
      for (size_t j=0; j<inear.size(); ++j) {
        const int32_t jdx = inear[j];
        // interparticle distance is unnormalized
        const CT dx = x[i] - x[jdx];
        const CT dy = y[i] - y[jdx];
        const CT dz = z[i] - z[jdx];
        const CT distsq = (dx*dx + dy*dy + dz*dz) / std::pow(r[i], 2);
        const CT dist = std::sqrt(distsq);
        // eta_eps in PSE terminology is corefunc / -dist   but apply the constant later
        const CT eta = dist * std::exp(-distsq);
        //printf("  part %d w str %g has neib %d %g away, eta is %g and change is %g\n", i, s[i], jdx, std::sqrt(distsq), eta, 2.0*h_nu*str_diff*eta);
        //if (i==525) printf("  part %d is %g from part %ld eta is %g and str is %g\n", jdx, std::sqrt(distsq), i, eta, s[jdx]);
        // compute and apply the strength change
        if (use_volumes) {
          // scale strengths by volumes
          dsx[i] += eta * (vol[i]*sx[jdx] - vol[jdx]*sx[i]);
          dsy[i] += eta * (vol[i]*sy[jdx] - vol[jdx]*sy[i]);
          dsz[i] += eta * (vol[i]*sz[jdx] - vol[jdx]*sz[i]);
        } else {
          // no volumes - just compare strengths
          dsx[i] += eta * (sx[jdx] - sx[i]);
          dsy[i] += eta * (sy[jdx] - sy[i]);
          dsz[i] += eta * (sz[jdx] - sz[i]);
        }
      }

      // scale by the proper constant factor
      //ds[i] *= 4.0 * std::pow(r[i], -4) / M_PI; // this should be correct for 2D, according to C&K VM pg 145
      const ST factor = 2.0 * std::pow(r[i], -5) / M_PI;
      dsx[i] *= factor;
      dsy[i] *= factor;
      dsz[i] *= factor;

    // compute as if every particle were a compact Gaussian
    } else {

      // loop over all participating particles
      for (size_t j=0; j<inear.size(); ++j) {
        const int32_t jdx = inear[j];
        // all distances are unnormalized
        const CT dx = x[i] - x[jdx];
        const CT dy = y[i] - y[jdx];
        const CT dz = z[i] - z[jdx];
        const CT distsq = (dx*dx + dy*dy + dz*dz) / std::pow(r[i], 2);
        const CT dist = std::sqrt(distsq);
        const CT distcub = distsq * dist;
        // eta_eps in PSE terminology is corefunc / -dist
        const CT eta = distsq * std::exp(-distcub);
        //printf("  part %d w str %g has neib %d %g away, eta is %g and change is %g\n", i, s[i], jdx, dist*r[i], eta, 2.0*h_nu*str_diff*eta);

        // compute and apply the strength change
        if (use_volumes) {
          // scale strengths by volumes
          dsx[i] += eta * (vol[i]*sx[jdx] - vol[jdx]*sx[i]);
          dsy[i] += eta * (vol[i]*sy[jdx] - vol[jdx]*sy[i]);
          dsz[i] += eta * (vol[i]*sz[jdx] - vol[jdx]*sz[i]);
        } else {
          // no volumes - just compare strengths
          dsx[i] += eta * (sx[jdx] - sx[i]);
          dsy[i] += eta * (sy[jdx] - sy[i]);
          dsz[i] += eta * (sz[jdx] - sz[i]);
        }
      }

      // scale by the proper constant factor
      // 0.716197244 is from OmegaFlow-2, its -9/4pi
      // but 0.6 produces something closer to the VRM answer
      const ST factor = 2.0 * 0.716197244 * std::pow(r[i], -5);
      dsx[i] *= factor;
      dsy[i] *= factor;
      dsz[i] *= factor;
    }

    if (not use_volumes) {
      // and correct for overlap
      const ST factor = std::pow(nom_sep, 3);
      dsx[i] *= factor;
      dsy[i] *= factor;
      dsz[i] *= factor;
    }

    // always scale the ds by the proper constant factor
    const ST factor = std::pow(h_nu, 2);
    dsx[i] *= factor;
    dsy[i] *= factor;
    dsz[i] *= factor;

    // tally neighbor statistics
    nneibs += inear.size();
    if (inear.size() < minneibs) minneibs = inear.size();
    if (inear.size() > maxneibs) maxneibs = inear.size();

    // is ds significantly larger than s?
    //if (std::abs(ds[i]) > 3.0 * std::abs(s[i]) and s[i] != 0.0) {
    //  std::cout << "    part " << i << " at dist " << std::sqrt(x[i]*x[i]+y[i]*y[i]) << " had str before " << s[i] << " str after " << (s[i]+ds[i]) << std::endl;
    //}

  } // end loop over all current particles

  std::cout << "    neighbors: min/avg/max " << minneibs << "/" << ((ST)nneibs / (ST)n) << "/" << maxneibs << std::endl;
  std::cout << "    after PSE, n is " << x.size() << std::endl;

  // apply the changes to the master vectors
  assert(n==sx.size() and dsx.size()==sx.size() && "Array size mismatch in PSE");
  assert(n==r.size() && "Array size mismatch in PSE");
  for (size_t i=0; i<n; ++i) {
    sx[i] += dsx[i];
    sy[i] += dsy[i];
    sz[i] += dsz[i];
  }

  // finish timer and report
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    pse.diffuse_all:\t[%.4f] seconds\n", (float)elapsed_seconds.count());
}

//
// read/write parameters to json
//

// create and write a json object for all diffusion parameters
template <class ST, class CT>
void PSE<ST,CT>::from_json(const nlohmann::json simj) {

  if (simj.find("PSE") != simj.end()) {
    nlohmann::json j = simj["PSE"];

    if (j.find("ignoreBelow") != j.end()) {
      ignore_thresh = j["ignoreBelow"];
      std::cout << "  setting ignore_thresh= " << ignore_thresh << std::endl;
    }

    if (j.find("relativeThresholds") != j.end()) {
      thresholds_are_relative = j["relativeThresholds"];
      std::cout << "  setting thresholds_are_relative= " << thresholds_are_relative << std::endl;
    }

    if (j.find("useVolumes") != j.end()) {
      use_volumes = j["useVolumes"];
      std::cout << "  setting use_volumes= " << use_volumes << std::endl;
    }
  }
}

// create and write a json object for all diffusion parameters
template <class ST, class CT>
void PSE<ST,CT>::add_to_json(nlohmann::json& simj) const {

  // set vrm-specific parameters
  nlohmann::json j;
  j["ignoreBelow"] = ignore_thresh;
  j["relativeThresholds"] = thresholds_are_relative;
  j["useVolumes"] = use_volumes;
  simj["PSE"] = j;
}

