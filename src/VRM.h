/*
 * VRM.h - the Vorticity Redistribution Method for 3D vortex particles
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Core.h"
#include "Icosahedron.h"
#include "VectorHelper.h"
#include "nanoflann.hpp"
#include "nnls.h"

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

//
// Class to hold VRM parameters and temporaries
//
// templatized on storage type ST, compute type CT, and max number of moments to solve for
//
template <class ST, class CT, uint8_t MAXMOM>
class VRM {
public:
  VRM();
  VRM(const CT);
  void set_hnu(const ST);
  ST get_hnu();
  void set_adaptive_radii(const bool);

  // all-to-all diffuse; can change array sizes
  void diffuse_all(Vector<ST>&, Vector<ST>&, Vector<ST>&,
                   Vector<ST>&, Vector<ST>&, Vector<ST>&,
                   Vector<ST>&, Vector<ST>&, Vector<ST>&,
                   Vector<ST>&, Vector<ST>&,
                   const CoreType, const ST);

protected:
  // search for new target location
  std::array<ST,3> fill_neighborhood_search(const int32_t,
                                            const Vector<ST>&,
                                            const Vector<ST>&,
                                            const Vector<ST>&,
                                            const std::vector<int32_t>&,
                                            const ST);

  bool attempt_solution(const int32_t,
                        std::vector<int32_t>&,
                        Vector<ST>&,
                        Vector<ST>&,
                        Vector<ST>&,
                        Vector<ST>&,
                        Vector<ST>&,
                        const CoreType,
                        Eigen::Matrix<CT, Eigen::Dynamic, 1>&);

private:
  // solve VRM to how many moments?
  static const int32_t num_moments = MAXMOM;
  static constexpr int32_t num_rows = (num_moments+1) * (num_moments+2) * (num_moments+3) / 6;
  // we needed 16 here for static solutions, 20 for dynamic, and 24 for dynamic with adaptivity
  static constexpr int32_t max_near = 48 * num_moments;

  // h_nu is sqrt(dt*nu) or sqrt(dt/Re)
  CT h_nu;

  // new point insertion sites (normalized to h_nu and centered around origin)
  size_t num_sites;
  std::vector<ST> xsite,ysite,zsite;
  void initialize_sites();

  // for adaptive particle size VRM
  bool adapt_radii = false;
  //const ST radius_lapse = 0.3;
  // only adapt particles if their strength is less than this
  //   fraction of max particle strength
  //const ST adapt_thresh = 1.e-2;

  // do not perform VRM if source particle strength is less than
  //   this fraction of max particle strength
  const ST ignore_thresh = 1.e-5;
  // are thresholds absolute or relative to strongest particle?
  const bool thresholds_are_relative = true;

  // use nanoflann for nearest-neighbor searching? false uses direct search
  const bool use_tree = true;
};

// delegating ctor
template <class ST, class CT, uint8_t MAXMOM>
VRM<ST,CT,MAXMOM>::VRM()
  : VRM(1.0)
  {}

// primary constructor
template <class ST, class CT, uint8_t MAXMOM>
VRM<ST,CT,MAXMOM>::VRM(const CT _hnu)
  : h_nu(_hnu) {

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
void VRM<ST,CT,MAXMOM>::set_hnu(const ST _newhnu) {
  h_nu = _newhnu;
}

template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::set_adaptive_radii(const bool _doamr) {
  //if (!adapt_radii and _doamr) std::cout << "Particle radii will adapt to solution" << std::endl;
  //if (adapt_radii and !_doamr) std::cout << "Particle radii will not adapt to solution" << std::endl;
  adapt_radii = _doamr;
}

template <class ST, class CT, uint8_t MAXMOM>
ST VRM<ST,CT,MAXMOM>::get_hnu() {
  return (ST)h_nu;
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
// Find the change in strength and radius that would occur over one dt
//
template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::diffuse_all(Vector<ST>& x, Vector<ST>& y, Vector<ST>& z,
                                    Vector<ST>& r, Vector<ST>& newr,
                                    Vector<ST>& sx, Vector<ST>& sy, Vector<ST>& sz,
                                    Vector<ST>& dsx, Vector<ST>& dsy, Vector<ST>& dsz,
                                    const CoreType core_func,
                                    const ST particle_overlap) {

  // make sure all vector sizes are identical
  assert(x.size()==y.size() && "Input array sizes do not match");
  assert(x.size()==z.size() && "Input array sizes do not match");
  assert(x.size()==r.size() && "Input array sizes do not match");
  assert(x.size()==newr.size() && "Input array sizes do not match");
  assert(x.size()==sx.size() && "Input array sizes do not match");
  assert(x.size()==sy.size() && "Input array sizes do not match");
  assert(x.size()==sz.size() && "Input array sizes do not match");
  assert(x.size()==dsx.size() && "Input array sizes do not match");
  assert(x.size()==dsy.size() && "Input array sizes do not match");
  assert(x.size()==dsz.size() && "Input array sizes do not match");
  size_t n = x.size();

  // start timers
  std::chrono::system_clock::time_point start, end;
  std::chrono::duration<double> elapsed_seconds;
  start = std::chrono::system_clock::now();

  std::cout << "  Running VRM with n " << n << std::endl;

  // zero out delta vector
  std::fill(dsx.begin(), dsx.end(), 0.0);
  std::fill(dsy.begin(), dsy.end(), 0.0);
  std::fill(dsz.begin(), dsz.end(), 0.0);

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
  std::cout << "    maxStrSqrd " << maxStrSqrd << std::endl;

  // convert particle positions into something nanoflann can understand
  Eigen::Matrix<ST, Eigen::Dynamic, 3> xp;
  xp.resize(n,3);
  xp.col(0) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(x.data(), n);
  xp.col(1) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(y.data(), n);
  xp.col(2) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(z.data(), n);
  
  // generate the searchable data structure
  typedef nanoflann::KDTreeEigenMatrixAdaptor< Eigen::Matrix<ST, Eigen::Dynamic, 3> >  my_kd_tree_t;
  my_kd_tree_t mat_index(xp, 20);
  if (use_tree) mat_index.index->buildIndex();
  std::vector<std::pair<long int,ST> > ret_matches;
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
  for (size_t i=0; i<initial_n; ++i) {

    // find the nearest neighbor particles
    //std::cout << "\nDiffusing particle " << i << " with strength " << sx[i] << " " << sy[i] << " " << sz[i] << std::endl;
    //std::cout << "\nDiffusing particle " << i << " at " << x[i] << " " << y[i] << " " << z[i] << std::endl;

    // if current particle strength is very small, skip out
    //   (this particle still core-spreads somewhat)
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
      //const size_t nMatches = mat_index.index->radiusSearch(query_pt, distsq_thresh, ret_matches, params);
      //std::cout << "radiusSearch(): radius " << search_rad << " found " << nMatches;

      // copy the indexes into my vector
      for (size_t j=0; j<ret_matches.size(); ++j) inear.push_back(ret_matches[j].first);

      // now direct search over all newer particles
      for (size_t j=initial_n; j<n; ++j) {
        ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2) + std::pow(z[i]-z[j], 2);
        if (distsq < distsq_thresh) inear.push_back(j);
      }
      //std::cout << " and " << (inear.size()-nMatches) << " matches";

      // is the closest of these "too close"? Start at 0.5 if you want to catch the most egregious.
      //if (std::sqrt(ret_matches[1].second) < 0.5*nom_sep) ntooclose++;

      //std::cout << std::endl;
      //for (size_t j=0; j<ret_matches.size(); ++j) std::cout << "   " << ret_matches[j].first << "\t" << ret_matches[j].second << std::endl;
    } else {
      // direct search: look for all neighboring particles, include newly-created ones
      //   ideally this would be a tree search
      const ST distsq_thresh = std::pow(search_rad, 2);
      for (size_t j=0; j<n; ++j) {
        ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2) + std::pow(z[i]-z[j], 2);
        if (distsq < distsq_thresh) inear.push_back(j);
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

    bool haveSolution = false;
    size_t numNewParts = 0;

    // assemble the underdetermined system
    while (not haveSolution and ++numNewParts < maxNewParts) {
      //std::cout << "  attempt solution with " << inear.size() << " close particles" << std::endl;

      // this does the heavy lifting - assemble and solve the VRM equations for the 
      //   diffusion from particle i to particles in inear
      haveSolution = attempt_solution(i, inear, x, y, z, r, newr, core_func, fractions);

      // if that didn't work, add a particle and try again
      if (not haveSolution) {
        // solution is bad, add a particle and try again
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
        //std::cout << "    no solution, added part at " << newpt[0] << " " << newpt[1] << " " << newpt[2] << std::endl;
      }
    }

    // did we eventually reach a solution?
    if (numNewParts >= maxNewParts) {
      std::cout << "Something went wrong" << std::endl;
      std::cout << "  needed numNewParts= " << numNewParts << std::endl;
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

  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  //const float gflops = (flops / 1.e+9) / (float)elapsed_seconds.count();
  //printf("    diffuse_all:\t[%.6f] seconds at [%.3f] GFlop/s\n", (float)elapsed_seconds.count(), gflops);
  printf("    diffuse_all:\t[%.4f] seconds\n", (float)elapsed_seconds.count());
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

  // for non-adaptive method and floats, 1e-6 fails immediately, 1e-5 fails quickly, 3e-5 seems to work
  // for doubles, can use 1e-6, will increase accuracy for slight performance hit (see vrm3d)
  static const CT nnls_eps = 3.e-5;
  // default to 1e-6, but drop to 1e-4 for adaptive with high overlap
  static const CT nnls_thresh = 1.e-6;

  // reset the arrays
  //std::cout << "\nSetting up Ax=b least-squares problem" << std::endl;
  b.setZero();
  A.resize(num_rows, inear.size());
  A.setZero();
  fractions.resize(inear.size());
  fractions.setZero();

  // fill it in
  for (size_t j=0; j<inear.size(); ++j) {
    const int32_t jdx = inear[j];
    // all distances are normalized to h_nu
    CT dx = (x[idiff]-x[jdx]) / h_nu;
    CT dy = (y[idiff]-y[jdx]) / h_nu;
    CT dz = (z[idiff]-z[jdx]) / h_nu;
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
  //std::cout << "  Here is the matrix A^T:\n" << A.transpose() << std::endl;
  //std::cout << "  Here is the right hand side b:\n\t" << b.transpose() << std::endl;
  //std::cout << "  Here is the solution vector:\n\t" << fractions.transpose() << std::endl;

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

  // set output fractions and result
  fracout = fractions;
  return haveSolution;
}

