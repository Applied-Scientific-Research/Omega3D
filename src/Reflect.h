/*
 * Reflect.h - Non-class particle-panel reflecting operation
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega3D.h"
#include "Points.h"
#include "Surfaces.h"
#include "MathHelper.h"

#include <cstdlib>
#include <limits>
#include <vector>
#include <cmath>

enum ClosestType { panel, edge, node };

//
// on output of this routine, jidx is:
//   -1..-3 if an edge is the closest
//    0 if the panel face is closest
//    1..3 if a node is closest
//
template <class S> 
struct ClosestReturn { 
  int jidx; 
  S distsq; 
  S cpx, cpy, cpz;
  ClosestType disttype; 
}; 

//
// find closest distance from a point to a triangle
// logic taken from pointElemDistance3d
//
// most likely 147 flops
//
template <class S>
ClosestReturn<S> panel_point_distance(const S sx0, const S sy0, const S sz0,
                                      const S sx1, const S sy1, const S sz1,
                                      const S sx2, const S sy2, const S sz2,
                                      const S sxn, const S syn, const S szn,
                                      const S tx,  const S ty,  const S tz) {

  ClosestReturn<S> retval;
  retval.jidx = 0;
  retval.distsq = 9.9e+9;

  // test against each of the three corners
  const std::array<S,3> dt0 = {tx-sx0, ty-sy0, tz-sz0};
  const S dt0ds = dot_product(dt0, dt0);
  if (dt0ds < retval.distsq) {
    retval.jidx = 1;
    retval.distsq = dt0ds;
    retval.disttype = node;
    retval.cpx = sx0;
    retval.cpy = sy0;
    retval.cpz = sz0;
  }

  const std::array<S,3> dt1 = {tx-sx1, ty-sy1, tz-sz1};
  const S dt1ds = dot_product(dt1, dt1);
  if (dt1ds < retval.distsq) {
    retval.jidx = 2;
    retval.distsq = dt1ds;
    retval.disttype = node;
    retval.cpx = sx1;
    retval.cpy = sy1;
    retval.cpz = sz1;
  }

  const std::array<S,3> dt2 = {tx-sx2, ty-sy2, tz-sz2};
  const S dt2ds = dot_product(dt2, dt2);
  if (dt2ds < retval.distsq) {
    retval.jidx = 3;
    retval.distsq = dt2ds;
    retval.disttype = node;
    retval.cpx = sx2;
    retval.cpy = sy2;
    retval.cpz = sz2;
  }
  // 27 flops to here

  //std::cout << "test point " << tx << " " << ty << " " << tz << std::endl;
  //std::cout << "  current distsq " << retval.distsq << std::endl;

  // compare to the edges
  const std::array<S,3> ds01 = {sx1-sx0, sy1-sy0, sz1-sz0};
  //std::cout << "  ds01 " << ds01[0] << " " << ds01[1] << " " << ds01[2] << std::endl;
  const S ds01ds = 1.0 / dot_product(ds01, ds01);
  //std::cout << "  ds01ds " << ds01ds << std::endl;
  std::array<S,3> rx;
  cross_product(ds01, dt0, rx);
  //std::cout << "  rx   " << rx[0] << " " << rx[1] << " " << rx[2] << std::endl;
  const S dt01ds = dot_product(rx, rx) * ds01ds;
  //std::cout << "  dt01ds " << dt01ds << std::endl;
  if (dt01ds < retval.distsq) {
    const S t = dot_product(ds01, dt0) * ds01ds;
    //std::cout << "  t " << t << std::endl;
    if (0.0 < t and t < 1.0) {
      retval.jidx = -1;
      retval.distsq = dt01ds;
      retval.disttype = edge;
      retval.cpx = sx0 + t*ds01[0];
      retval.cpy = sy0 + t*ds01[1];
      retval.cpz = sz0 + t*ds01[2];
    }
  }
  // 25 minimum, 33 maybe, 39 possibly

  const std::array<S,3> ds12 = {sx2-sx1, sy2-sy1, sz2-sz1};
  const S ds12ds = 1.0 / dot_product(ds12, ds12);
  cross_product(ds12, dt1, rx);
  const S dt12ds = dot_product(rx, rx) * ds12ds;
  //std::cout << "  dt12ds " << dt12ds << std::endl;
  if (dt12ds < retval.distsq) {
    const S t = dot_product(ds12, dt1) * ds12ds;
    //std::cout << "  t " << t << std::endl;
    if (0.0 < t and t < 1.0) {
      retval.jidx = -2;
      retval.distsq = dt12ds;
      retval.disttype = edge;
      retval.cpx = sx1 + t*ds12[0];
      retval.cpy = sy1 + t*ds12[1];
      retval.cpz = sz1 + t*ds12[2];
    }
  }

  const std::array<S,3> ds20 = {sx0-sx2, sy0-sy2, sz0-sz2};
  const S ds20ds = 1.0 / dot_product(ds20, ds20);
  cross_product(ds20, dt2, rx);
  const S dt20ds = dot_product(rx, rx) * ds20ds;
  //std::cout << "  dt20ds " << dt20ds << std::endl;
  if (dt20ds < retval.distsq) {
    const S t = dot_product(ds20, dt2) * ds20ds;
    //std::cout << "  t " << t << std::endl;
    if (0.0 < t and t < 1.0) {
      retval.jidx = -3;
      retval.distsq = dt20ds;
      retval.disttype = edge;
      retval.cpx = sx2 + t*ds20[0];
      retval.cpy = sy2 + t*ds20[1];
      retval.cpz = sz2 + t*ds20[2];
    }
  }
  // 102 minimum flops to here

  // finally, check vs. the panel prism
  // test to see if the test point is in the positive half-space of each edge
  const std::array<S,3> norm = {sxn, syn, szn};
  std::array<S,3> interior;
  cross_product(norm, ds01, interior);
  const S in01 = dot_product(dt0, interior);
  cross_product(norm, ds12, interior);
  const S in12 = dot_product(dt1, interior);
  cross_product(norm, ds20, interior);
  const S in20 = dot_product(dt2, interior);
  //std::cout << "  interior distances are " << in01 << " " << in12 << " " << in20 << std::endl;
  if (in01 > 0.0 and in12 > 0.0 and in20 > 0.0) {
    retval.jidx = 0;
    retval.disttype = panel;
    const S truedist = dot_product(dt0, norm);
    retval.distsq = truedist*truedist;
    retval.cpx = tx - sxn*truedist;
    retval.cpy = ty - syn*truedist;
    retval.cpz = tz - szn*truedist;
  }
  // minumum 45 more flops here, possibly 57 if match

  return retval;
}


//
// naive caller for the O(N^2) panel-particle reflection kernel
//
template <class S>
void reflect_panp2 (Surfaces<S> const& _src, Points<S>& _targ) {
  //std::cout << "  inside reflect(Surfaces, Points)" << std::endl;
  std::cout << "  Reflecting" << _targ.to_string() << " from near" << _src.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  // get handles for the vectors
  std::array<Vector<S>,Dimensions> const& sx = _src.get_pos();
  std::vector<Int> const&                 si = _src.get_idx();
  std::array<Vector<S>,3> const&          sn = _src.get_norm();
  std::array<Vector<S>,Dimensions>&       tx = _targ.get_pos();

  size_t num_reflected = 0;
  const S eps = 10.0*std::numeric_limits<S>::epsilon();

  #pragma omp parallel for reduction(+:num_reflected)
  for (size_t i=0; i<_targ.get_n(); ++i) {

    S mindist = std::numeric_limits<S>::max();
    std::vector<ClosestReturn<S>> hits;

    // iterate and search for closest panel
    for (size_t j=0; j<_src.get_npanels(); ++j) {
      const Int jp0 = si[3*j+0];
      const Int jp1 = si[3*j+1];
      const Int jp2 = si[3*j+2];
      ClosestReturn<S> result = panel_point_distance<S>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                                                        sx[0][jp1], sx[1][jp1], sx[2][jp1],
                                                        sx[0][jp2], sx[1][jp2], sx[2][jp2],
                                                        sn[0][j],   sn[1][j],   sn[2][j],
                                                        tx[0][i],   tx[1][i],   tx[2][i]);

      if (result.distsq < mindist - eps) {
        // we blew the old one away
        mindist = result.distsq;
        // rewrite jidx as the panel index
        result.jidx = j;
        hits.clear();
        hits.push_back(result);
        //std::cout << "  THIS BLOWS AWAY THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;

      } else if (result.distsq < mindist + eps) {
        // we are effectively the same as the old closest
        // rewrite jidx as the panel index
        result.jidx = j;
        hits.push_back(result);
        //std::cout << "  THIS TIES THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;
      }
    }

    // dump out the hits
    if (false) {
      std::cout << "point " << i << " is " << tx[0][i] << " " << tx[1][i] << std::endl;
      for (auto & ahit: hits) {
        if (ahit.disttype == node) {
          std::cout << "  node on panel " << ahit.jidx << " is " << std::sqrt(ahit.distsq) << std::endl;
        } else if (ahit.disttype == edge) {
          std::cout << "  edge on panel " << ahit.jidx << " is " << std::sqrt(ahit.distsq) << std::endl;
        } else {
          std::cout << "  panel " << ahit.jidx << " is " << std::sqrt(ahit.distsq) << std::endl;
        }
      }
    }

    // if no hits, then something is wrong
    assert(hits.size() > 0 && "No nearest neighbors");

    //std::cout << "  REFLECTING pt at " << tx[0][i] << " " << tx[1][i] << " " << tx[2][i] << std::endl;

    // now look at the vector of return values and decide if we're under or above the panel!
    // init the mean normal and the mean contact point
    std::array<S,3> mnorm = {0.0, 0.0, 0.0};
    std::array<S,3> mcp = {0.0, 0.0, 0.0};

    // accumulate mean normal and the mean contact point
    for (size_t k=0; k<hits.size(); ++k) {
      //std::cout << "    cp at " << hits[k].cpx << " " << hits[k].cpy << std::endl;

      const size_t j = hits[k].jidx;
      // doesn't matter if it hits a panel, edge, or node, just add the participating panel's norm
      for (size_t d=0; d<3; ++d) mnorm[d] += sn[d][j];

      mcp[0] += hits[k].cpx;
      mcp[1] += hits[k].cpy;
      mcp[2] += hits[k].cpz;
    }

    normalizeVec(mnorm);
    for (size_t d=0; d<3; ++d) mcp[d] /= (S)hits.size();

    // compare this mean norm to the vector from the contact point to the particle
    std::array<S,3> dx = {tx[0][i]-mcp[0], tx[1][i]-mcp[1], tx[2][i]-mcp[2]};
    const S dotp = dot_product(mnorm, dx);

    if (dotp < 0.0) {
      // this point is under the panel - reflect it off entry 0
      // this is reasonable for most cases, except very sharp angles between adjacent panels
      const S dist = std::sqrt(hits[0].distsq);
      //std::cout << "  REFLECTING pt at rad " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]+tx[2][i]*tx[2][i]);
      for (size_t d=0; d<3; ++d) tx[d][i] = mcp[d] + dist*mnorm[d];
      //std::cout << "  to " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]+tx[2][i]*tx[2][i]) << std::endl;
      num_reflected++;
    }
  }

  std::cout << "    reflected " << num_reflected << " particles" << std::endl;
  const S flops = 149.0 * (float)_targ.get_n() * (float)_src.get_npanels();

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    reflect_panp2:\t[%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


//
// reflect interior particles to exterior because VRM only works in free space
//
template <class S>
void reflect_interior(std::vector<Collection> const & _bdry,
                      std::vector<Collection>&        _vort) {

  // may need to do this multiple times to clear out concave zones!
  // this should only function when _vort is Points and _bdry is Surfaces
  for (auto &targ : _vort) {
    if (std::holds_alternative<Points<S>>(targ)) {
      Points<S>& pts = std::get<Points<S>>(targ);

      for (auto &src : _bdry) {
        if (std::holds_alternative<Surfaces<S>>(src)) {
          Surfaces<S> const & surf = std::get<Surfaces<S>>(src);

          // call the specific panels-affect-points routine
          (void) reflect_panp2<S>(surf, pts);
        }
      }
    }
  }
}


//
// generate the cut tables
//
template <class S>
std::vector<std::tuple<S,S,S>> init_cut_tables (const S _dx) {
  // get the vector ready
  std::vector<std::tuple<S,S,S>> ct;

  // all distances are normalized to default vdelta (not nominal separation)
  // going out to 3 vdeltas catches all but 1e-4 of the strength
  const S max_rad = 2.5;
  // how many entries?
  const int nx = (int)(0.5 + max_rad / _dx);
  const S dx = max_rad / (float)nx;
  //std::cout << "Making cut tables with nx " << nx << " and dx " << dx << std::endl;

  // add the first entry (remove all strength, set dshift later)
  ct.push_back(std::make_tuple((S)(-nx-0.5)*dx, 0.0, 0.0));

  // generate the entries
  S twgt = 0.0;
  S tmom = 0.0;
  for (int i=-nx; i<nx+1; ++i) {
    const S distxs = std::pow((S)i * dx, 2);

    // compute the weight of this plane
    S rwgt = 0.0;
    for (int j=-nx; j<nx+1; ++j) {
      const S distxys = distxs + std::pow((S)j * dx, 2);
      for (int k=-nx; k<nx+1; ++k) {
        const S distxyzs = distxys + std::pow((S)k * dx, 2);
        // this is a pure 2D Gaussian
        rwgt += std::exp(-distxyzs);
        // this is for a 2D compact Gaussian
        //rwgt += std::exp(-std::pow(distxys, 1.5));
      }
    }
    // scale by the cell size
    rwgt *= std::pow(dx, 3);
    // this is a pure 3D Gaussian ( pi^(-3/2) )
    rwgt *= 0.179587122;
    // this is for a 3D compact Gaussian
    //rwgt *= 0.23873241;

    //std::cout << "  row " << i << " at " << ((S)i*dx) << " has weight " << rwgt << std::endl;

    // total weight accumulated is the fraction of strength to keep
    twgt += rwgt;

    // total moment is the amount to shift
    tmom += ((S)i*dx) * rwgt;

    // add an entry
    ct.push_back(std::make_tuple(((S)i+0.5)*dx, twgt, -tmom/twgt));
  }
  //std::cout << "  total weight " << twgt << std::endl;

  // add the last entry (keep all strength, set dshift to zero)
  ct.push_back(std::make_tuple((S)(nx+1.5)*dx, 1.0, 0.0));

  //std::cout << "Cut table is" << std::endl;
  //for (auto &entry : ct) {
    //std::cout << "  " << std::get<0>(entry) << " " << std::get<1>(entry) << " " << std::get<2>(entry) << std::endl;
  //}

  return ct;
}

//
// use the cut tables - assume _pos is normalized by vdelta
//
template <class S>
std::pair<S,S> get_cut_entry (std::vector<std::tuple<S,S,S>>& ct, const S _pos) {
  // set defaults (change nothing)
  S smult = 1.0;
  S dshift = 0.0;

  const size_t inum = ct.size();
  assert(inum > 0 && "Cut tables not initialized");

  // check vs. low and high bounds
  if (_pos < std::get<0>(ct[0])) {
    // particle is too far inside to survive
    smult = 0.0;
    dshift = 1.0 - _pos;
  } else if (_pos > std::get<0>(ct[inum-1])) {
    // particle is too far outside to be affected
    smult = 1.0;
    dshift = 0.0;
  } else {
    // linearly interpolate between two values
    for (size_t i=1; i<inum; ++i) {
      if (_pos < std::get<0>(ct[i])) {
        // point lies between this and the previous entry
        const S frac = (_pos - std::get<0>(ct[i-1])) / (std::get<0>(ct[i]) - std::get<0>(ct[i-1]));
        smult = frac*std::get<1>(ct[i]) + (1.0-frac)*std::get<1>(ct[i-1]);
        dshift = frac*std::get<2>(ct[i]) + (1.0-frac)*std::get<2>(ct[i-1]);
        //std::cout << "new dist " << _pos << " found smult " << smult << " and dshift " << dshift << std::endl;
        break;
      }
    }
  }

  return std::pair<S,S>(smult,dshift);
}


//
// naive caller for the O(N^2) panel-particle clear-inner-layer kernel
//
template <class S>
void clear_inner_panp2 (Surfaces<S> const & _src,
                        Points<S>& _targ,
                        const S _cutoff_mult,
                        const S _ips) {

  //std::cout << "  inside clear_inner_layer(Surfaces, Points)" << std::endl;
  std::cout << "  Clearing" << _targ.to_string() << " from near" << _src.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  static bool made_cut_tables = false;
  static std::vector<std::tuple<S,S,S>> ct;
  if (not made_cut_tables) {
    ct = init_cut_tables<S>(0.1);
    made_cut_tables = true;
  }

  // get handles for the vectors
  std::array<Vector<S>,Dimensions> const& sx = _src.get_pos();
  std::vector<Int> const&                 si = _src.get_idx();
  std::array<Vector<S>,3> const&          sn = _src.get_norm();

  std::array<Vector<S>,Dimensions>&       tx = _targ.get_pos();
  std::array<Vector<S>,Dimensions>&       ts = _targ.get_str();
  // if called on field points, there is no tr
  Vector<S>&                              tr = _targ.get_rad();
  const bool are_fldpts = tr.empty();

  size_t num_cropped = 0;
  const S eps = 10.0*std::numeric_limits<S>::epsilon();

  #pragma omp parallel for reduction(+:num_cropped)
  for (size_t i=0; i<_targ.get_n(); ++i) {

    S mindist = std::numeric_limits<S>::max();
    std::vector<ClosestReturn<S>> hits;

    // iterate and search for closest panel
    for (size_t j=0; j<_src.get_npanels(); ++j) {
      const Int jp0 = si[3*j+0];
      const Int jp1 = si[3*j+1];
      const Int jp2 = si[3*j+2];
      ClosestReturn<S> result = panel_point_distance<S>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                                                        sx[0][jp1], sx[1][jp1], sx[2][jp1],
                                                        sx[0][jp2], sx[1][jp2], sx[2][jp2],
                                                        sn[0][j],   sn[1][j],   sn[2][j],
                                                        tx[0][i],   tx[1][i],   tx[2][i]);

      if (result.distsq < mindist - eps) {
        // we blew the old one away
        mindist = result.distsq;
        // rewrite jidx as the panel index
        result.jidx = j;
        hits.clear();
        hits.push_back(result);
        //std::cout << "  THIS BLOWS AWAY THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;

      } else if (result.distsq < mindist + eps) {
        // we are effectively the same as the old closest
        // rewrite jidx as the panel index
        result.jidx = j;
        hits.push_back(result);
        //std::cout << "  THIS TIES THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;
      }
    } // end of loop over panels

    // dump out the hits
    if (false) {
      std::cout << "point " << i << " is " << tx[0][i] << " " << tx[1][i] << std::endl;
      for (auto & ahit: hits) {
        if (ahit.disttype == node) {
          std::cout << "  node on panel " << ahit.jidx << " is " << std::sqrt(ahit.distsq) << std::endl;
        } else if (ahit.disttype == edge) {
          std::cout << "  edge on panel " << ahit.jidx << " is " << std::sqrt(ahit.distsq) << std::endl;
        } else {
          std::cout << "  panel " << ahit.jidx << " is " << std::sqrt(ahit.distsq) << std::endl;
        }
      }
    }

    // if no hits, then something is wrong
    assert(hits.size() > 0 && "No nearest neighbors");

    //std::cout << "  CLEARING pt at " << tx[0][i] << " " << tx[1][i] << std::endl;

    // init the mean normal and the mean contact point
    std::array<S,3> mnorm = {0.0, 0.0, 0.0};
    std::array<S,3> mcp = {0.0, 0.0, 0.0};

    // accumulate mean normal and the mean contact point
    for (size_t k=0; k<hits.size(); ++k) {
      //std::cout << "    cp at " << hits[k].cpx << " " << hits[k].cpy << std::endl;

      const size_t j = hits[k].jidx;
      // doesn't matter if it hits a panel, edge, or node, just add the participating panel's norm
      for (size_t d=0; d<3; ++d) mnorm[d] += sn[d][j];

      mcp[0] += hits[k].cpx;
      mcp[1] += hits[k].cpy;
      mcp[2] += hits[k].cpz;
    }

    normalizeVec(mnorm);
    for (size_t d=0; d<3; ++d) mcp[d] /= (S)hits.size();

    // compare this mean norm to the vector from the contact point to the particle
    std::array<S,3> dx = {tx[0][i]-mcp[0], tx[1][i]-mcp[1], tx[2][i]-mcp[2]};
    const S dotp = dot_product(mnorm, dx) - _cutoff_mult*_ips;
    // now dotp is how far this point is above the cutoff layer
    // if dotp == 0.0 then the point is exactly on the cutoff layer, and it loses half of its strength
    // if dotp < -vdelta then the point loses all of its strength
    // if dotp > vdelta then the point passes unmodified

    const S this_radius = are_fldpts ? _ips : tr[i];

    if (dotp < this_radius) {
      //std::cout << "  CLEARING pt at " << tx[0][i] << " " << tx[1][i] << " because dotp " << dotp << " and norm " << norm[0] << " " << norm[1] << std::endl;

      // use precomputed table lookups for new position and remaining strength
      std::pair<S,S> entry = get_cut_entry(ct, dotp/this_radius);
      if (not are_fldpts) for (size_t d=0; d<3; ++d) ts[d][i] *= std::get<0>(entry);
      for (size_t d=0; d<3; ++d) tx[d][i] += std::get<1>(entry) * this_radius * mnorm[d];

      //std::cout << "  PUSHING " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]);
      //std::cout << "    to " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]) << " and scale str by " << sfac << std::endl;
      // do not change radius yet
      //std::cout << "    TO " << tx[0][i] << " " << tx[1][i] << " and weaken by " << sfac << std::endl;

      num_cropped++;
    }

  } // end loop over particles

  // we did not resize the x array, so we don't need to touch the u array

  std::cout << "    cropped " << num_cropped << " particles" << std::endl;
  // flops count here is taken from reflect - might be different here
  const S flops = _targ.get_n() * (0.0 + 0.0*_src.get_npanels());

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    clear_inner_panp2:\t[%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


//
// clean up by removing the innermost layer - the one that will be represented by boundary strengths
//
template <class S>
void clear_inner_layer(std::vector<Collection> const & _bdry,
                       std::vector<Collection>&        _vort,
                       const S                         _cutoff_factor,
                       const S                         _ips) {

  // may need to do this multiple times to clear out concave zones!
  // this should only function when _vort is Points and _bdry is Surfaces
  for (auto &targ : _vort) {
    if (std::holds_alternative<Points<S>>(targ)) {
      Points<S>& pts = std::get<Points<S>>(targ);

      for (auto &src : _bdry) {
        if (std::holds_alternative<Surfaces<S>>(src)) {
          Surfaces<S> const & surf = std::get<Surfaces<S>>(src);

          // call the specific panels-affect-points routine
          (void) clear_inner_panp2<S>(surf, pts, _cutoff_factor, _ips);
        }
      }
    }
  }
}

