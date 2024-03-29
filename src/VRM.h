/*
 * VRM.h - the Vorticity Redistribution Method for 3D vortex particles
 *
 * (c)2017-20 Applied Scientific Research, Inc.
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
#include "Icosahedron.h"
#include "VectorHelper.h"
#ifdef PLUGIN_SIMPLEX
#include "simplex.h"
#endif

#include "eigen-nnls/nnls.h"
#include "nanoflann/nanoflann.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

enum SolverType { nnls, simplex };

//
// Class to hold VRM parameters and temporaries
//
// templatized on storage type ST, compute type CT, and max number of moments to solve for
//
template <class ST, class CT, uint8_t MAXMOM>
class VRM {
public:
  VRM();

  void set_adaptive_radii(const bool);
  void set_relative(const bool _in) { thresholds_are_relative = _in; }
  void set_ignore(const float _in) { ignore_thresh = _in; }
  void set_simplex(const bool _in) { use_solver = (_in ? simplex : nnls); }

  const bool get_adaptive_radii() const { return adapt_radii; }
  const bool get_relative() const { return thresholds_are_relative; }
  const float get_ignore() const { return ignore_thresh; }
  const bool get_simplex() const { return (use_solver==simplex); }

  // all-to-all diffuse; can change array sizes
  void diffuse_all(std::array<Vector<ST>,3>&,
                   std::array<Vector<ST>,3>&,
                   Vector<ST>&,
                   const ST,
                   const CoreType,
                   const ST);

  // other functions to eventually support:
  // two-to-one merge (when particles are close to each other)
  // one-to-many elongate (re-sphericalize a stretched particle)

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

  // set up and call the solver
  bool attempt_solution(const int32_t,
                        std::vector<int32_t>&,
                        Vector<ST>&,
                        Vector<ST>&,
                        Vector<ST>&,
                        Vector<ST>&,
                        Vector<ST>&,
                        const ST,
                        const CoreType,
                        Eigen::Matrix<CT, Eigen::Dynamic, 1>&);

private:
  // solve VRM to how many moments?
  static const int32_t num_moments = MAXMOM;
  static constexpr int32_t num_rows = (num_moments+1) * (num_moments+2) * (num_moments+3) / 6;
  // we needed 16 here for static solutions, 20 for dynamic, and 24 for dynamic with adaptivity
  static constexpr int32_t max_near = 48 * num_moments;

  // new point insertion sites (normalized to h_nu and centered around origin)
  size_t num_sites;
  std::vector<ST> xsite,ysite,zsite;
  void initialize_sites();

  // for standard VRM

  // do not perform VRM if source particle strength is less than
  //   this fraction of max particle strength
  ST ignore_thresh = 1.e-4;
  // are thresholds absolute or relative to strongest particle?
  bool thresholds_are_relative = true;

  // for any adaptive particle size VRM
  bool adapt_radii = false;
  // would have more here

  // use nanoflann for nearest-neighbor searching? false uses direct search
  const bool use_tree = true;

  SolverType use_solver = nnls;
  //SolverType use_solver = simplex;
};

// primary constructor
template <class ST, class CT, uint8_t MAXMOM>
VRM<ST,CT,MAXMOM>::VRM() {
  initialize_sites();
}

//
// initialize the particle placement sites array
//
template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::initialize_sites() {

  const std::vector<ST>& ico = ico2;

  if (MAXMOM <= 2) {
    //std::cout << "Creating one layer of insertion sites" << std::endl;

    // for first and second-moment VRM, one layer only
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

  } else {
    //std::cout << "Creating two layers of insertion sites" << std::endl;

    if (true) {
      // pull data from a multi-level refinement of an r=0.5 icosahedron
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

    // generate a second layer farther away
    // to solve for 3rd-4th moments stably, we need more points
    if (true) {
      // pull data from a multi-level refinement of an r=0.5 icosahedron
      const size_t add_num_sites = ico.size() / Dimensions;
      xsite.resize(num_sites+add_num_sites);
      ysite.resize(num_sites+add_num_sites);
      zsite.resize(num_sites+add_num_sites);
      for (size_t i=0; i<add_num_sites; ++i) {
        const size_t idx = i+num_sites;
        xsite[idx] = 3.7 * ico[3*i+0];
        ysite[idx] = 3.7 * ico[3*i+1];
        zsite[idx] = 3.7 * ico[3*i+2];
      }
      num_sites += add_num_sites;
    }
  }
}

template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::set_adaptive_radii(const bool _doamr) {
  //if (!adapt_radii and _doamr) std::cout << "Particle radii will adapt to solution" << std::endl;
  //if (adapt_radii and !_doamr) std::cout << "Particle radii will not adapt to solution" << std::endl;
  adapt_radii = _doamr;
}

//
// use a ring of sites to determine the location of a new particle
//
template <class ST, class CT, uint8_t MAXMOM>
std::array<ST,3> VRM<ST,CT,MAXMOM>::fill_neighborhood_search(const int32_t idx,
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
      ST distsq = std::pow(x[inear[j]]-tx[i], 2) + std::pow(y[inear[j]]-ty[i], 2) + std::pow(z[inear[j]]-tz[i], 2);
      if (distsq < mindistsq) mindistsq = distsq;
    }
    nearest[i] = mindistsq;
    if (mindistsq > maxmindist) {
      maxmindist = mindistsq;
      iopen = i;
    }
  }

  //std::cout << "    creating particle at " << tx[iopen] << " " << ty[iopen]
  //          << ", now there are " << n << " particles" << std::endl;

  const std::array<ST,3> retval = {{tx[iopen], ty[iopen], tz[iopen]}};
  return retval;
  //return std::array<ST,3>(tx[iopen], ty[iopen], tz[iopen]);
}


//
// Find the change in strength and radius that would occur over one dt and apply it
//
template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::diffuse_all(std::array<Vector<ST>,3>& pos,
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
  size_t n = rad.size();

  std::cout << "  Running VRM with n " << n << std::endl;

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

  Vector<ST> newr, dsx, dsy, dsz;
  newr.resize(n);
  dsx.resize(n);
  dsy.resize(n);
  dsz.resize(n);

  // zero out delta vector
  std::fill(dsx.begin(), dsx.end(), 0.0);
  std::fill(dsy.begin(), dsy.end(), 0.0);
  std::fill(dsz.begin(), dsz.end(), 0.0);

  // and copy the new radius
  std::copy(r.begin(), r.end(), newr.begin());

  // create the matrix elements, reusable, and dynamically allocated
  Eigen::Matrix<CT, Eigen::Dynamic, 1> fractions;

  const size_t minNearby = 11;
  const size_t maxNewParts = num_moments*12 - 4;

  // what is maximum strength of all particles?
  ST maxStr = 0.0;
  for (size_t i=0; i<n; ++i) {
    const ST thisstr = sx[i]*sx[i] + sy[i]*sy[i] + sz[i]*sz[i];
    if (thisstr > maxStr) maxStr = thisstr;
  }
  const ST maxStrSqrd = maxStr;
  //std::cout << "    maxStrSqrd " << maxStrSqrd << std::endl;

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

  // do not adapt particle radii -- copy current to new
  newr = r;

  // for each particle (can parallelize this part)
  const size_t initial_n = n;
  size_t nsolved = 0;
  size_t nneibs = 0;
  //size_t ntooclose = 0;
  size_t minneibs = 999999;
  size_t maxneibs = 0;
  // note that an OpenMP loop here will need to use int32_t as the counter variable type
  for (size_t i=0; i<initial_n; ++i) {

    // find the nearest neighbor particles
    //std::cout << "\nDiffusing particle " << i << " with strength " << sx[i] << " " << sy[i] << " " << sz[i] << std::endl;
    //std::cout << "\nDiffusing particle " << i << " at " << x[i] << " " << y[i] << " " << z[i] << std::endl;

    // if current particle strength is very small, skip out
    //   (this particle could still core-spread if adaptive particle size is on)
    const ST thisstr = sx[i]*sx[i] + sy[i]*sy[i] + sz[i]*sz[i];
    if ((thresholds_are_relative && (thisstr < maxStrSqrd * std::pow(ignore_thresh,2))) or
        (!thresholds_are_relative && (thisstr < std::pow(ignore_thresh,2)))) continue;

    nsolved++;

    // nominal separation for this particle (insertion distance)
    const ST nom_sep = r[i] / particle_overlap;

    // what is search radius?
    const ST search_rad = nom_sep * ((num_moments > 2) ? 2.5 : 1.6);

    // initialize vector of indexes of nearest particles
    std::vector<int32_t> inear;

    // switch on search method
    if (use_tree) {
      // tree-based search with nanoflann
      const ST distsq_thresh = std::pow(search_rad, 2);
      const ST query_pt[3] = { x[i], y[i], z[i] };
      (void) mat_index.index->radiusSearch(query_pt, distsq_thresh, ret_matches, params);
      //if (ret_matches.size() > 100) std::cout << "part " << i << " at " << x[i] << " " << y[i] << " " << z[i] << " has " << ret_matches.size() << " matches" << std::endl;

      // copy the indexes into my vector
      for (size_t j=0; j<ret_matches.size(); ++j) inear.push_back((int32_t)ret_matches[j].first);

      // now direct search over all newer particles
      for (size_t j=initial_n; j<n; ++j) {
        ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2) + std::pow(z[i]-z[j], 2);
        if (distsq < distsq_thresh) inear.push_back((int32_t)j);
      }
      //std::cout << " and " << (inear.size()-nMatches) << " matches";

      // is the closest of these "too close"? Start at 0.5 if you want to catch the most egregious.
      //if (std::sqrt(ret_matches[1].second) < 0.5*nom_sep) ntooclose++;

      //std::cout << "radiusSearch(): radius " << search_rad << " found " << inear.size();
      //std::cout << std::endl;
      //for (size_t j=0; j<ret_matches.size(); ++j) std::cout << "   " << ret_matches[j].first << "\t" << ret_matches[j].second << std::endl;
    } else {
      // direct search: look for all neighboring particles, include newly-created ones
      //   ideally this would be a tree search
      const ST distsq_thresh = std::pow(search_rad, 2);
      for (size_t j=0; j<n; ++j) {
        ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2) + std::pow(z[i]-z[j], 2);
        if (distsq < distsq_thresh) inear.push_back((int32_t)j);
      }
      //std::cout << "radiusSearch(): radius " << search_rad << " found " << inear.size() << " matches";

    }
    //std::cout << "  found " << inear.size() << " particles close to particle " << i << std::endl;
    //std::cout << " :";
    //for (size_t j=0; j<inear.size(); ++j) std::cout << " " << inear[j];
    //std::cout << std::endl;

    // if there are less than, say, 6, we should just add some now
    while (inear.size() < minNearby) {
      auto newpt = fill_neighborhood_search(i, x, y, z, inear, nom_sep);
      // add the index to the near list
      inear.push_back(n);
      // add it to the master list
      x.push_back(newpt[0]);
      y.push_back(newpt[1]);
      z.push_back(newpt[2]);
      r.push_back(newr[i]);
      newr.push_back(newr[i]);
      sx.push_back(0.0);
      sy.push_back(0.0);
      sz.push_back(0.0);
      dsx.push_back(0.0);
      dsy.push_back(0.0);
      dsz.push_back(0.0);
      n++;
      //std::cout << "  inear is";
      //for (int32_t j=0; j<inear.size(); ++j) std::cout << " " << inear[j];
      //std::cout << std::endl;
      //std::cout << "    fill neib with new part at " << newpt[0] << " " << newpt[1] << " " << newpt[2] << std::endl;
    }

    // now remove close parts if we have more than max_near
    while (inear.size() > max_near) {
      // look for the part closest to the diffusing particle
      int32_t jclose = 0;
      ST distnear = std::numeric_limits<ST>::max();
      for (size_t j=0; j<inear.size(); ++j) {
        const size_t inearj = static_cast<size_t>(inear[j]);
        if (inearj != i) {
          const ST distsq = std::pow(x[i]-x[inearj], 2) + std::pow(y[i]-y[inearj], 2) + std::pow(z[i]-z[inearj], 2);
          if (distsq < distnear) {
            distnear = distsq;
            jclose = j;
          }
        }
      }
      // and remove it
      //std::cout << "removing pt near " << x[i] << " " << y[i] << std::endl;
      inear.erase(inear.begin()+jclose);
    }

    bool haveSolution = false;
    size_t numNewParts = 0;

    // assemble the underdetermined system
    while (not haveSolution and ++numNewParts < maxNewParts) {
      //std::cout << "  attempt solution with " << inear.size() << " close particles" << std::endl;

      // this does the heavy lifting - assemble and solve the VRM equations for the 
      //   diffusion from particle i to particles in inear
      haveSolution = attempt_solution(i, inear, x, y, z, r, newr, h_nu, core_func, fractions);

      // if that didn't work, add a particle and try again
      if (not haveSolution) {
        auto newpt = fill_neighborhood_search(i, x, y, z, inear, nom_sep);

        if (inear.size() == max_near) {
          // replace an old particle with this new one
          size_t ireplace = 1;	// default is 1 because diffusing particle is probably position 0
          for (size_t j=0; j<inear.size(); ++j) {
            const size_t inj = static_cast<size_t>(inear[j]);
            if (inj != i and inj < initial_n) {
              ireplace = j;
              break;
            }
          }
          // we are moving an original particle from the near list, but not the global list
          // but we are adding a new particle to the global list, hence the n
          inear[ireplace] = n;
          //std::cout << "  replacing pt " << ireplace << " in the list of " << inear.size() << std::endl;
        } else {
          // add a new one to the inear list and the master particle lists
          inear.push_back(n);
        }

        // add it to the master list
        x.push_back(newpt[0]);
        y.push_back(newpt[1]);
        z.push_back(newpt[2]);
        r.push_back(newr[i]);
        newr.push_back(newr[i]);
        sx.push_back(0.0);
        sy.push_back(0.0);
        sz.push_back(0.0);
        dsx.push_back(0.0);
        dsy.push_back(0.0);
        dsz.push_back(0.0);
        n++;
        //std::cout << "    no solution, added part at " << newpt[0] << " " << newpt[1] << " " << newpt[2] << std::endl;
      }
    }

    // did we eventually reach a solution?
    if (numNewParts >= maxNewParts) {
      std::cout << "Something went wrong" << std::endl;
      std::cout << "  at " << x[i] << " " << y[i] << " " << z[i] << std::endl;
      std::cout << "  with " << inear.size() << " near neibs" << std::endl;
      std::cout << "  needed numNewParts= " << numNewParts << std::endl;
      // ideally, in this situation, we would create 6+ new particles around the original particle with optimal fractions,
      //   ignoring every other nearby particle - let merge take care of the higher density later
      exit(0);
    }

    nneibs += inear.size();
    if (inear.size() < minneibs) minneibs = inear.size();
    if (inear.size() > maxneibs) maxneibs = inear.size();

    // apply those fractions to the delta vector
    //std::cout << "Added strengths" << std::endl;
    for (size_t j=0; j<inear.size(); ++j) {
      size_t idx = inear[j];
      if (idx == i) {
        // self-influence
        dsx[idx] += sx[i] * (fractions(j) - 1.0);
        dsy[idx] += sy[i] * (fractions(j) - 1.0);
        dsz[idx] += sz[i] * (fractions(j) - 1.0);
      } else {
        dsx[idx] += sx[i] * fractions(j);
        dsy[idx] += sy[i] * fractions(j);
        dsz[idx] += sz[i] * fractions(j);
      }
      //std::cout << "  " << (s[i]*fractions(j)) << " to particle " << idx << std::endl;
    }

  } // end loop over all current particles

  std::cout << "    neighbors: min/avg/max " << minneibs << "/" << ((ST)nneibs / (ST)nsolved) << "/" << maxneibs << std::endl;
  //std::cout << "  number of close pairs " << (ntooclose/2) << std::endl;
  std::cout << "    after VRM, n is " << n << std::endl;

  // apply the changes to the master vectors
  for (size_t i=0; i<n; ++i) {
    r[i] = newr[i];
  }
  assert(n==sx.size() and dsx.size()==sx.size() && "Array size mismatch in VRM");
  for (size_t i=0; i<n; ++i) {
    sx[i] += dsx[i];
  }
  assert(n==sy.size() and dsy.size()==sy.size() && "Array size mismatch in VRM");
  for (size_t i=0; i<n; ++i) {
    sy[i] += dsy[i];
  }
  assert(n==sz.size() and dsz.size()==sz.size() && "Array size mismatch in VRM");
  for (size_t i=0; i<n; ++i) {
    sz[i] += dsz[i];
  }

  // finish timer and report
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    vrm.diffuse_all:\t[%.4f] seconds\n", (float)elapsed_seconds.count());
}

//
// Set up and solve the VRM equations
//
template <class ST, class CT, uint8_t MAXMOM>
bool VRM<ST,CT,MAXMOM>::attempt_solution(const int32_t idiff,
                                         std::vector<int32_t>& inear,
                                         Vector<ST>& x,
                                         Vector<ST>& y,
                                         Vector<ST>& z,
                                         Vector<ST>& r,
                                         Vector<ST>& newr,
                                         const ST h_nu,
                                         const CoreType core_func,
                                         Eigen::Matrix<CT, Eigen::Dynamic, 1>& fracout) {

  bool haveSolution = false;

  // the matricies that we will repeatedly work on
  static Eigen::Matrix<CT, num_rows, Eigen::Dynamic, 0, num_rows, max_near> A;
  static Eigen::Matrix<CT, num_rows, 1> b;
  static Eigen::Matrix<CT, Eigen::Dynamic, 1, 0, max_near, 1> fractions;

  // second moment in each direction
  // one dt should generate 6 hnu^2 of second moment, or when distances
  //   are normalized by hnu, 2.0 in each direction
  // core moments are the coefficients on the moments based on core type
  static const CT second_moment = 2.0;
  static const CT fourth_moment = 12.0;

  const CT oohnu = 1.0 / h_nu;

  // the Ixx and Iyy moments of these core functions is half of the 2nd radial moment
  //static const CT core_second_mom = get_core_second_mom<CT>(core_func);
  // the Ixxxx and Iyyyy moments of these core functions is 3/8th of the 4th radial moment
  //static const CT core_fourth_mom = get_core_fourth_mom<CT>(core_func);

  // for non-adaptive method and floats, 1e-6 fails immediately, 1e-5 fails quickly, 3e-5 seems to work
  // for doubles, can use 1e-6, will increase accuracy for slight performance hit (see vrm3d)
  static const CT nnls_eps = 1.e-5;
  // default to 1e-6, but drop to 1e-4 for adaptive with high overlap?
  static const CT nnls_thresh = 1.e-6;
#ifdef PLUGIN_SIMPLEX
  // default to 1e-6, but drop to 1e-4 for adaptive with high overlap?
  static const CT simplex_thresh = 1.e-6;
#endif

  // reset the arrays
  //std::cout << "\nSetting up Ax=b least-squares problem" << std::endl;
  assert(inear.size() <= static_cast<size_t>(max_near) && "Too many neighbors in VRM");
  b.setZero();
  A.resize(num_rows, inear.size());
  A.setZero();
  fractions.resize(inear.size());
  fractions.setZero();

  // fill it in
  for (size_t j=0; j<inear.size(); ++j) {
    const int32_t jdx = inear[j];
    // all distances are normalized to h_nu
    CT dx = (x[idiff]-x[jdx]) * oohnu;
    CT dy = (y[idiff]-y[jdx]) * oohnu;
    CT dz = (z[idiff]-z[jdx]) * oohnu;
    A(0,j) = 1.0;
    if (num_moments > 0) {
      A(1,j) = dx;
      A(2,j) = dy;
      A(3,j) = dz;
    }
    if (num_moments > 1) {
      A(4,j) = dx*dx;
      A(5,j) = dx*dy;
      A(6,j) = dx*dz;
      A(7,j) = dy*dy;
      A(8,j) = dy*dz;
      A(9,j) = dz*dz;
    }
    if (num_moments > 2) {
      A(10,j) = dx*dx*dx;
      A(11,j) = dx*dx*dy;
      A(12,j) = dx*dx*dz;
      A(13,j) = dx*dy*dy;
      A(14,j) = dx*dy*dz;
      A(15,j) = dx*dz*dz;
      A(16,j) = dy*dy*dy;
      A(17,j) = dy*dy*dz;
      A(18,j) = dy*dz*dz;
      A(19,j) = dz*dz*dz;
    }
    // fourth moments should be ?
    if (num_moments > 3) {
      A(20,j) = dx*dx*dx*dx;
      A(21,j) = dx*dx*dx*dy;
      A(22,j) = dx*dx*dx*dz;
      A(23,j) = dx*dx*dy*dy;
      A(24,j) = dx*dx*dy*dz;
      A(25,j) = dx*dx*dz*dz;
      A(26,j) = dx*dy*dy*dy;
      A(27,j) = dx*dy*dy*dz;
      A(28,j) = dx*dy*dz*dz;
      A(29,j) = dx*dz*dz*dz;
      A(30,j) = dy*dy*dy*dy;
      A(31,j) = dy*dy*dy*dz;
      A(32,j) = dy*dy*dz*dz;
      A(33,j) = dy*dz*dz*dz;
      A(34,j) = dz*dz*dz*dz;
      // now 20, 30, 35, cross at 23, 25, 32
    }
  }

  b(0) = 1.f;
  if (num_moments > 1) {
    b(4) = second_moment;
    b(7) = second_moment;
    b(9) = second_moment;
  }
  if (num_moments > 3) {
    b(20) = fourth_moment;
    b(30) = fourth_moment;
    b(34) = fourth_moment;
    b(23) = fourth_moment / 3.0;
    b(25) = fourth_moment / 3.0;
    b(32) = fourth_moment / 3.0;
  }
  
  if (VERBOSE and false) {
    std::cout << "  Here is the matrix A^T:\n" << A.transpose() << std::endl;
    std::cout << "  Here is the right hand side b:\n\t" << b.transpose() << std::endl;
    std::cout << "  Here is the solution vector:\n\t" << fractions.transpose() << std::endl;
  }

  if (use_solver == nnls) {
    //std::cout << "    using NNLS solver\n" << std::endl;

    // solve with non-negative least-squares
    Eigen::NNLS<Eigen::Matrix<CT,Eigen::Dynamic,Eigen::Dynamic> > nnls_solver(A, 100, nnls_eps);

    //std::cout << "A is" << std::endl << A << std::endl;
    //std::cout << "b is" << std::endl << b.transpose() << std::endl;
    if (nnls_solver.solve(b)) {
      fractions = nnls_solver.x();
      //std::cout << "  success! required " << nnls_solver.numLS() << " LS problems" << std::endl;
      //std::cout << "  check says " << nnls_solver.check(b) << std::endl;
    } else {
      for (size_t j=0; j<inear.size(); ++j) fractions(j) = 0.f;
      //std::cout << "  fail!" << std::endl;
    }

    //std::cout << "  fractions are:\n\t" << fractions.transpose() << std::endl;

    // measure the results
    Eigen::Matrix<CT,Eigen::Dynamic,1> err = A*fractions - b;
    //std::cout << "  error is:\n" << err.transpose() << std::endl;
    //std::cout << "  error magnitude is " << std::sqrt(err.dot(err)) << std::endl;

    // was this solution successful?
    if (err.dot(err) < nnls_thresh) {
      // this is good enough!
      haveSolution = true;
    }

  } else {
    //std::cout << "    using Simplex solver\n" << std::endl;

#ifdef PLUGIN_SIMPLEX
    // solve with simplex solver

    // first, adjust the RHS
    for (size_t j=0; j<num_rows; ++j) {
      b(j) -= 0.5 * A.row(j).sum();
    }
    //std::cout << "  New right hand side b:\n\t" << b.transpose() << std::endl;

    // finally call the solver
    CT LInfNorm = 1.0;
    int retval = undr_dtrmn_solvr<CT,num_rows,max_near>(A, b, fractions, LInfNorm);
    //std::cout << "  undr_dtrmn_solvr returned " << retval << " " << LInfNorm << std::endl;

    // if we used the simplex solver, adjust the fractions here
    fractions.array() += 0.5;
    //std::cout << "  final fractions are:\n\t" << fractions.transpose() << std::endl;

    // was this solution successful?
    if (retval == 0 and LInfNorm < 0.5 + simplex_thresh) {
      // this is good enough!
      haveSolution = true;
    }
#else
    // we should never get here
    throw "Simplex solver is not available.";
#endif
  }

  // set output fractions and result
  fracout = fractions;
  return haveSolution;
}

//
// read/write parameters to json
//

// read a json object and retrieve all diffusion parameters
template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::from_json(const nlohmann::json simj) {

  if (simj.find("VRM") != simj.end()) {
    nlohmann::json j = simj["VRM"];

    if (j.find("ignoreBelow") != j.end()) {
      ignore_thresh = j["ignoreBelow"];
      std::cout << "  setting ignore_thresh= " << ignore_thresh << std::endl;
    }

    if (j.find("relativeThresholds") != j.end()) {
      thresholds_are_relative = j["relativeThresholds"];
      std::cout << "  setting thresholds_are_relative= " << thresholds_are_relative << std::endl;
    }

    if (j.find("solver") != j.end()) {
      std::string solverstr = j["solver"];
      if (solverstr == "simplex") {
        use_solver = simplex;
      } else {
        // default is nnls
        use_solver = nnls;
      }
    } else {
      // default is nnls
      use_solver = nnls;
    }
#ifdef PLUGIN_SIMPLEX
    std::cout << "  setting VRM solver= " << (use_solver ? "simplex" : "nnls") << std::endl;
#else
    if (use_solver == simplex) {
      std::cout << "  setting VRM solver= nnls because PLUGIN_SIMPLEX is not set" << std::endl;
      use_solver = nnls;
    } else {
      std::cout << "  setting VRM solver= nnls" << std::endl;
    }
#endif
  }
}

// create and write a json object for all diffusion parameters
template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::add_to_json(nlohmann::json& simj) const {

  // set vrm-specific parameters
  nlohmann::json j;
  j["ignoreBelow"] = ignore_thresh;
  j["relativeThresholds"] = thresholds_are_relative;
  j["solver"] = use_solver ? "simplex" : "nnls";
  simj["VRM"] = j;
}

