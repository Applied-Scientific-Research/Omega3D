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

template <class S> 
struct ClosestReturn { 
  size_t jpanel; 
  S distsq; 
  S cpx, cpy, cpz;
  ClosestType disttype; 
}; 

//
// find closest distance from a point to a triangle
// logic taken from pointElemDistance3d
//
template <class S>
ClosestReturn<S> panel_point_distance(const S sx0, const S sy0, const S sz0,
                                      const S sx1, const S sy1, const S sz1,
                                      const S sx2, const S sy2, const S sz2,
                                      const S sxn, const S syn, const S szn,
                                      const S tx,  const S ty,  const S tz) {

  ClosestReturn<S> retval;
  retval.jpanel = 0;
  retval.distsq = 9.9e+9;

  // test against each of the three corners
  const std::array<S,3> dt0 = {tx-sx0, ty-sy0, tz-sz0};
  const S dt0ds = dot_product(dt0, dt0);
  if (dt0ds < retval.distsq) {
    retval.distsq = dt0ds;
    retval.disttype = node;
    retval.cpx = sx0;
    retval.cpy = sy0;
    retval.cpz = sz0;
  }

  const std::array<S,3> dt1 = {tx-sx1, ty-sy1, tz-sz1};
  const S dt1ds = dot_product(dt1, dt1);
  if (dt1ds < retval.distsq) {
    retval.distsq = dt1ds;
    retval.disttype = node;
    retval.cpx = sx1;
    retval.cpy = sy1;
    retval.cpz = sz1;
  }

  const std::array<S,3> dt2 = {tx-sx2, ty-sy2, tz-sz2};
  const S dt2ds = dot_product(dt2, dt2);
  if (dt2ds < retval.distsq) {
    retval.distsq = dt2ds;
    retval.disttype = node;
    retval.cpx = sx2;
    retval.cpy = sy2;
    retval.cpz = sz2;
  }

  std::cout << "test point " << tx << " " << ty << " " << tz << std::endl;
  std::cout << "  current distsq " << retval.distsq << std::endl;

  // compare to the edges
  const std::array<S,3> ds01 = {sx1-sx0, sy1-sy0, sz1-sz0};
  std::cout << "  ds01 " << ds01[0] << " " << ds01[1] << " " << ds01[2] << std::endl;
  const S ds01ds = 1.0 / dot_product(ds01, ds01);
  std::cout << "  ds01ds " << ds01ds << std::endl;
  std::array<S,3> rx;
  cross_product(ds01, dt0, rx);
  std::cout << "  rx   " << rx[0] << " " << rx[1] << " " << rx[2] << std::endl;
  const S dt01ds = dot_product(rx, rx) * ds01ds;
  std::cout << "  dt01ds " << dt01ds << std::endl;
  if (dt01ds < retval.distsq) {
    const S t = dot_product(ds01, dt0) * ds01ds;
    std::cout << "  t " << t << std::endl;
    if (0.0 < t and t < 1.0) {
      retval.distsq = dt01ds;
      retval.disttype = edge;
      retval.cpx = sx0 + t*ds01[0];
      retval.cpy = sy0 + t*ds01[1];
      retval.cpz = sz0 + t*ds01[2];
    }
  }

  const std::array<S,3> ds12 = {sx2-sx1, sy2-sy1, sz2-sz1};
  const S ds12ds = 1.0 / dot_product(ds12, ds12);
  cross_product(ds12, dt1, rx);
  const S dt12ds = dot_product(rx, rx) * ds12ds;
  std::cout << "  dt12ds " << dt12ds << std::endl;
  if (dt12ds < retval.distsq) {
    const S t = dot_product(ds12, dt1) * ds12ds;
    std::cout << "  t " << t << std::endl;
    if (0.0 < t and t < 1.0) {
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
  std::cout << "  dt20ds " << dt20ds << std::endl;
  if (dt20ds < retval.distsq) {
    const S t = dot_product(ds20, dt2) * ds20ds;
    std::cout << "  t " << t << std::endl;
    if (0.0 < t and t < 1.0) {
      retval.distsq = dt20ds;
      retval.disttype = edge;
      retval.cpx = sx2 + t*ds20[0];
      retval.cpy = sy2 + t*ds20[1];
      retval.cpz = sz2 + t*ds20[2];
    }
  }

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
  std::cout << "  interior distances are " << in01 << " " << in12 << " " << in20 << std::endl;
  if (in01 > 0.0 and in12 > 0.0 and in20 > 0.0) {
    retval.disttype = panel;
    const S truedist = dot_product(dt0, norm);
    retval.distsq = truedist*truedist
    retval.cpx = tx - sxn*truedist;
    retval.cpy = ty - syn*truedist;
    retval.cpz = tz - szn*truedist;
  }

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

  // run some tests first
  ClosestReturn<S> r = panel_point_distance<S>(0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, 20.0,0.001,0.0);
  //std::cout << "  type " << r.disttype << " dist " << std::sqrt(r.distsq) << " cp " << r.cpx << " " << r.cpy << " " << r.cpz << std::endl;
  r = panel_point_distance<S>(0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, 20.0,-0.001,0.0);
  //std::cout << "  type " << r.disttype << " dist " << std::sqrt(r.distsq) << " cp " << r.cpx << " " << r.cpy << " " << r.cpz << std::endl;
  r = panel_point_distance<S>(0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, 0.001,-20.0,1.0);
  std::cout << "  type " << r.disttype << " dist " << std::sqrt(r.distsq) << " cp " << r.cpx << " " << r.cpy << " " << r.cpz << std::endl;
  r = panel_point_distance<S>(0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, -0.001,-20.0,1.0);
  std::cout << "  type " << r.disttype << " dist " << std::sqrt(r.distsq) << " cp " << r.cpx << " " << r.cpy << " " << r.cpz << std::endl;
  r = panel_point_distance<S>(0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, 0.1,-20.0,1.0);
  std::cout << "  type " << r.disttype << " dist " << std::sqrt(r.distsq) << " cp " << r.cpx << " " << r.cpy << " " << r.cpz << std::endl;
  r = panel_point_distance<S>(0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, 0.1,0.001,20.0);
  std::cout << "  type " << r.disttype << " dist " << std::sqrt(r.distsq) << " cp " << r.cpx << " " << r.cpy << " " << r.cpz << std::endl;
  r = panel_point_distance<S>(0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, 0.1,-0.001,20.0);
  std::cout << "  type " << r.disttype << " dist " << std::sqrt(r.distsq) << " cp " << r.cpx << " " << r.cpy << " " << r.cpz << std::endl;

  assert(false && "quit");

  #pragma omp parallel for
  for (size_t i=0; i<_targ.get_n(); ++i) {
    S mindist = std::numeric_limits<S>::max();
    std::vector<ClosestReturn<S>> hits;

    // iterate and search for closest panel
    for (size_t j=0; j<_src.get_npanels(); ++j) {
      const Int jp0 = 3*j+0;
      const Int jp1 = 3*j+1;
      const Int jp2 = 3*j+2;
      ClosestReturn<S> result = panel_point_distance<S>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                                                        sx[0][jp1], sx[1][jp1], sx[2][jp1],
                                                        sx[0][jp2], sx[1][jp2], sx[2][jp2],
                                                        sn[0][j],   sn[1][j],   sn[2][j],
                                                        tx[0][i],   tx[1][i],   tx[2][i]);

      if (result.distsq < mindist - eps) {
        // we blew the old one away
        mindist = result.distsq;
        result.jpanel = j;
        hits.clear();
        hits.push_back(result);
        //std::cout << "  THIS BLOWS AWAY THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;

      } else if (result.distsq < mindist + eps) {
        // we are effectively the same as the old closest
        result.jpanel = j;
        hits.push_back(result);
        //std::cout << "  THIS TIES THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;
      }
    }

    // dump out the hits
    if (false) {
    std::cout << "point " << i << " is " << tx[0][i] << " " << tx[1][i] << std::endl;
    for (auto & ahit: hits) {
      if (ahit.disttype == node) {
        std::cout << "  node " << ahit.jpanel << " is " << std::sqrt(ahit.distsq) << std::endl;
      } else {
        std::cout << "  panel " << ahit.jpanel << " is " << std::sqrt(ahit.distsq) << std::endl;
      }
    }
    }

    //std::cout << "  REFLECTING pt at " << tx[0][i] << " " << tx[1][i] << std::endl;

    // now look at the vector of return values and decide if we're under or above the panel!
    if (hits.size() == 1) {
      // this is easy if the closest is a panel!
      if (hits[0].disttype == panel) {
        // compare vector to normal vector
        const size_t j = hits[0].jpanel;
        const S bx = sx[0][si[2*j+1]] - sx[0][si[2*j]];
        const S by = sx[1][si[2*j+1]] - sx[1][si[2*j]];
        const S blen  = 1.0 / std::sqrt(bx*bx + by*by);
        const S normx = -by * blen;
        const S normy =  bx * blen;
        const S dist  = normx*(tx[0][i]-hits[0].cpx) + normy*(tx[1][i]-hits[0].cpy);
        //std::cout << "  dist is actually " << dist << std::endl;
        if (dist < 0.0) {
          // this point is under the panel - reflect it
          tx[0][i] -= 2.0*dist*normx;
          tx[1][i] -= 2.0*dist*normy;
          //std::cout << "    TYPE 1 TO " << tx[0][i] << " " << tx[1][i] << std::endl;
          num_reflected++;
        }
      } else {
        // this should never happen - any single hit must be a panel because every
        //   node effectively gets checked twice
        const S dist = std::sqrt(hits[0].distsq);
        // if the distance is large enough, we don't need to care
        if (dist < 0.1) {
          std::cout << "WARNING: point at " << tx[0][i] << " " << tx[1][i] << std::endl;
          std::cout << "  only hits one node on panel " << hits[0].jpanel << " with dist " << dist << std::endl;
        }
        // we cannot define a distance!
      }
    } else if (hits.size() == 2) {
      // if two hits, they should both be nodes, but in rare cases can be anything
      S normx = 0.0;
      S normy = 0.0;
      S cpx = 0.0;
      S cpy = 0.0;
      // find the mean normal and the mean contact point
      for (size_t k=0; k<hits.size(); ++k) {
        //std::cout << "    cp at " << hits[k].cpx << " " << hits[k].cpy << std::endl;
        const size_t j = hits[k].jpanel;
        const S bx = sx[0][si[2*j+1]] - sx[0][si[2*j]];
        const S by = sx[1][si[2*j+1]] - sx[1][si[2*j]];
        const S blen = 1.0 / std::sqrt(bx*bx + by*by);
        normx += -by * blen;
        normy +=  bx * blen;
        cpx += hits[k].cpx;
        cpy += hits[k].cpy;
      }
      const S normilen = 1.0 / std::sqrt(normx*normx + normy*normy);
      normx *= normilen;
      normy *= normilen;
      cpx /= (S)hits.size();
      cpy /= (S)hits.size();
      // compare this mean norm to the vector from the contact point to the particle
      const S dotp = normx*(tx[0][i]-cpx) + normy*(tx[1][i]-cpy);
      if (dotp < 0.0) {
        // this point is under the panel - reflect it off entry 0
        // this is reasonable for most cases, except very sharp angles between adjacent panels
        const S dist = std::sqrt(hits[0].distsq);
        //std::cout << "  REFLECTING " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]);
        tx[0][i] = hits[0].cpx + dist*normx;
        tx[1][i] = hits[0].cpy + dist*normy;
        //std::cout << " to " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]) << std::endl;
        //std::cout << "    TYPE 2 TO " << tx[0][i] << " " << tx[1][i] << std::endl;
        num_reflected++;
      }
    } else {
      // this should never happen!
      std::cout << "WARNING: point at " << tx[0][i] << " " << tx[1][i] << " hits " << hits.size() << " nodes/panels!" << std::endl;
    }
  }

  std::cout << "    reflected " << num_reflected << " particles" << std::endl;

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    reflect_panp2:\t[%.4f] seconds\n", (float)elapsed_seconds.count());
  //const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  //printf("    panels_affect_points: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


//
// naive caller for the O(N^2) panel-particle reflection kernel
//
template <class S>
void clear_inner_panp2 (Surfaces<S> const & _src, Points<S>& _targ, const S _cutoff) {
  //std::cout << "  inside clear_inner_layer(Surfaces, Points)" << std::endl;
  std::cout << "  Clearing" << _targ.to_string() << " from near" << _src.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  // get handles for the vectors
  std::array<Vector<S>,Dimensions> const& sx = _src.get_pos();
  std::vector<Int> const&                 si = _src.get_idx();
  std::array<Vector<S>,Dimensions>&       tx = _targ.get_pos();
  std::array<Vector<S>,Dimensions>&       ts = _targ.get_str();
  Vector<S>&                              tr = _targ.get_rad();

  size_t num_cropped = 0;
  const S eps = 5.0*std::numeric_limits<S>::epsilon();
  //const S opeps = 1.0 + eps;
  //const S omeps = 1.0 - eps;

  // accumulate results into targvel
  #pragma omp parallel for reduction(+:num_cropped)
  for (size_t i=0; i<_targ.get_n(); ++i) {

    // reserve space for oft-used temporaries
    S along[2], norm[2], thisnorm[2];
    S oopanlen, dotp, dist_sqrd, x0, y0, x1, y1;

    S near_dist_sqrd = std::numeric_limits<S>::max();
    // inear is the NODE id that the particle is closest to
    Int inear = std::numeric_limits<Int>::max();
    norm[0] = 0.0; norm[1] = 0.0;

    // first check - is the point far enough away from the body to ensure that it is outside/inside?
    // implement this later

    // iterate and search for closest panel/node
    for (size_t j=0; j<_src.get_n(); ++j) {

      x0 = sx[0][si[2*j]];
      y0 = sx[1][si[2*j]];
      x1 = sx[0][si[2*j+1]];
      y1 = sx[1][si[2*j+1]];
      // is the point closest to the panel line?

      // find vector along panel segment
      along[0] = x1 - x0;
      along[1] = y1 - y0;
      // one over the panel length is useful
      oopanlen = 1.0 / std::sqrt(along[0]*along[0] + along[1]*along[1]);
      // normalize along
      along[0] *= oopanlen;
      along[1] *= oopanlen;
      // and normal
      thisnorm[0] = -along[1];
      thisnorm[1] =  along[0];
      // find projection
      dotp = (tx[0][i]-x0)*along[0] + (tx[1][i]-y0)*along[1];
      dotp *= oopanlen;
      // dotp is now 0.0 to 1.0 for points ON the panel

      if (dotp < 0.0 + eps) {
        // particle is closer to first node
        dist_sqrd = std::pow(tx[0][i]-x0, 2) + pow(tx[1][i]-y0, 2);
        const S oodist = 1.0 / std::sqrt(dist_sqrd);
        if (dist_sqrd < near_dist_sqrd - eps) {
          // point is clearly the closest
          near_dist_sqrd = dist_sqrd;
          inear = si[2*j];
          // replace the running normal
          norm[0] = (tx[0][i]-x0) * oodist;
          norm[1] = (tx[1][i]-y0) * oodist;
        } else if (dist_sqrd < near_dist_sqrd + eps) {
          // point is just as close as another point
          // add the normal to the running sum
          norm[0] += (tx[0][i]-x0) * oodist;
          norm[1] += (tx[1][i]-y0) * oodist;
        }

      } else if (dotp > 1.0 - eps) {
        // particle is closer to second node
        dist_sqrd = std::pow(tx[0][i]-x1, 2) + std::pow(tx[1][i]-y1, 2);
        const S oodist = 1.0 / std::sqrt(dist_sqrd);
        if (dist_sqrd < near_dist_sqrd - eps) {
          // point is clearly the closest
          near_dist_sqrd = dist_sqrd;
          inear = si[2*j+1];
          // replace the running normal
          norm[0] = (tx[0][i]-x1) * oodist;
          norm[1] = (tx[1][i]-y1) * oodist;
        } else if (dist_sqrd < near_dist_sqrd + eps) {
          // point is just as close as another point
          // add the normal to the running sum
          norm[0] += (tx[0][i]-x1) * oodist;
          norm[1] += (tx[1][i]-y1) * oodist;
        }

      } else {
        // particle is closest to segment, not ends
        // find actual distance
        // first find closest point along segment
        const S cpx = x0 + dotp*along[0]/oopanlen;
        const S cpy = y0 + dotp*along[1]/oopanlen;
        // then compute distance
        dist_sqrd = std::pow(tx[0][i]-cpx, 2) + std::pow(tx[1][i]-cpy, 2);
        // and compare to saved
        if (dist_sqrd < near_dist_sqrd - eps) {
          // point is clearly the closest
          near_dist_sqrd = dist_sqrd;
          // and it doesn't matter which node we pick, the normal calculation is always the same
          inear = si[2*j];
          // replace the normal with this normal
          norm[0] = thisnorm[0];
          norm[1] = thisnorm[1];
        } else if (dist_sqrd < near_dist_sqrd + eps) {
          // point is just as close as another point
          // add the normal to the running sum
          norm[0] += thisnorm[0];
          norm[1] += thisnorm[1];
        }
      }

    } // end of loop over panels

    //std::cout << "  CLEARING pt at " << tx[0][i] << " " << tx[1][i] << std::endl;

    // compare the particle to the normal to see which side of the object it is on
    // then, find out which side the point is on (dot with normal)
    oopanlen = 1.0 / std::sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
    norm[0] *= oopanlen;
    norm[1] *= oopanlen;
    dotp = (tx[0][i]-sx[0][inear])*norm[0] + (tx[1][i]-sx[1][inear])*norm[1];
    if (dotp > 0.0) {
      // particle is above panel
    } else {
      // particle is beneath panel
    }
    //std::cout << "    part is " << dotp << " away from node " << inear << std::endl;

    // HACK, weaken and push it out
    // eventually precompute table lookups for new position and remaining strength
    dotp -= _cutoff;
    if (dotp < tr[i]) {
      //std::cout << "  CLEARING pt at " << tx[0][i] << " " << tx[1][i] << " because dotp " << dotp << " and norm " << norm[0] << " " << norm[1] << std::endl;

      // smoothly vary the strength scaling factor from 0..1 over range dotp/vdelta = -1..1
      const S sfac = 0.5 * (1.0 + sin(0.5*M_PI*std::max((S)-1.0, dotp/tr[i])));
      for (size_t d=0; d<3; ++d) ts[d][i] *= sfac;

      // move particle up to above the cutoff
      const S shiftd = tr[i] * 0.5 * (1.0 - dotp/tr[i]);
      //std::cout << "  PUSHING " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]);
      tx[0][i] += shiftd * norm[0];
      tx[1][i] += shiftd * norm[1];
      //std::cout << "    to " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]) << " and weaken by " << sfac << std::endl;
      // do not change radius yet
      //std::cout << "    TO " << tx[0][i] << " " << tx[1][i] << " and weaken by " << sfac << std::endl;

      num_cropped++;
    }

  } // end loop over particles

  // we did not resize the x array, so we don't need to touch the u array

  std::cout << "    cropped " << num_cropped << " particles" << std::endl;

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    clear_inner_panp2:\t[%.4f] seconds\n", (float)elapsed_seconds.count());
}

