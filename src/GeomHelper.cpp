/*
 * GeomHelper.cpp - Generators for basic solid shapes
 *
 * (c)2021 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#include "Omega3D.h"
#include "Icosahedron.h"
#include "Cube.h"
#include "ElementPacket.h"
#include "IglRefine.h"
#include "IglMergeDups.h"
#include "IglDecimate.h"

#include <cmath>
#include <iostream>
#include <algorithm>	// for std::fill

//
// the goal here is to create canonical solid geometries at the proper size
//   and resolution, but not rotated or translated to their starting positions
//

//
// Create a closed spherical/ovoid object from a icosahedron
//
ElementPacket<float> generate_ovoid(const float _rx, const float _ry, const float _rz, const float _ips) {

  std::cout << "Creating ovoid..." << std::endl;

  // first, make an icosahedron
  std::vector<float>   x = ico0;
  std::vector<Int>   idx = ico0idx;
  std::vector<float> val;
  ElementPacket<float> epack {x, idx, val, x.size()/Dimensions, 2};

  // estimate the triangle spacing for the scaled ovoid
  const float maxscale = std::max(_rx, std::max(_ry, _rz));
  float meansize = 0.35 * maxscale;

  // then, iteratively refine it
  std::cout << "  sphere is icosahedron with 20";
  while (meansize > 1.2*_ips) {

    // refine once
    refine_geometry(epack);

    // re-sphericalize (r=0.5)
    for (size_t i=0; i<epack.x.size()/Dimensions; ++i) {
      const float rad = std::sqrt( std::pow(epack.x[Dimensions*i], 2) +
                                   std::pow(epack.x[Dimensions*i+1], 2) +
                                   std::pow(epack.x[Dimensions*i+2], 2));
      const float scale = 0.5 / rad;
      epack.x[Dimensions*i]   *= scale;
      epack.x[Dimensions*i+1] *= scale;
      epack.x[Dimensions*i+2] *= scale;
    }

    meansize *= 0.5;
  }
  std::cout << " panels" << std::endl;

  // final scaling
  for (size_t i=0; i<epack.x.size()/Dimensions; ++i) {
    epack.x[Dimensions*i]   *= _rx;
    epack.x[Dimensions*i+1] *= _ry;
    epack.x[Dimensions*i+2] *= _rz;
  }

  // keep val empty
  epack.nelem = epack.x.size()/Dimensions;

  return epack;
}


//
// Create a closed rectangular object
//
ElementPacket<float> generate_cuboid(const float _dx, const float _dy, const float _dz, const float _ips) {

  std::cout << "Creating solid rectangle..." << std::endl;

  // first, make 12-triangle rectangle
  std::vector<float>   x = cube0;
  std::vector<Int>   idx = cube0idx;
  std::vector<float> val;
  ElementPacket<float> epack {x, idx, val, x.size()/Dimensions, 2};

  // estimate the triangle spacing for the scaled rectangle
  float maxscale = std::max(_dx, std::max(_dy, _dz));
  float meansize = 0.7 * maxscale;

  // then, iteratively refine it
  std::cout << "  rectangle is cube with 12";
  while (meansize > 1.2*_ips) {
    refine_geometry(epack);
    meansize *= 0.5;
  }
  std::cout << " panels" << std::endl;

  // final scaling here
  for (size_t i=0; i<epack.x.size()/Dimensions; ++i) {
    epack.x[Dimensions*i]   *= _dx;
    epack.x[Dimensions*i+1] *= _dy;
    epack.x[Dimensions*i+2] *= _dz;
  }

  // keep val empty
  epack.nelem = epack.x.size()/Dimensions;

  return epack;
}


//
// Create a closed torus
//
ElementPacket<float> generate_torus(const float _rmax, const float _rmin, const float _ips) {

  std::cout << "Creating solid torus..." << std::endl;

  // set up the element packet
  std::vector<float>   x;
  std::vector<Int>   idx;
  std::vector<float> val;

  // Number of steps around the torus
  const int ts = std::max(77,int(M_PI*_rmax/_ips));
  // Number of steps around the circle perpendicular to the table
  const int cs = std::max(7,int(M_PI*_rmin/_ips));

  // We first make a circle we will then drag in a circular motion to create the torus
  std::vector<float> circle;
  for (int i=0; i<cs; i++) {
    const float alpha = 2.0*M_PI*i/(float)cs;
    circle.emplace_back(_rmin*std::sin(alpha));
    circle.emplace_back(_rmin*std::cos(alpha));
  }

  // loop over integer indices
  for (int i=0; i<ts; ++i) {
    const float theta = 2.0 * M_PI * (float)i / (float)ts;
    const float st = std::sin(theta);
    const float ct = std::cos(theta);

    for (int j=0; j<cs; ++j) {
      // helper indices for the disk particles
      const size_t ix = 2*j;
      const size_t iy = 2*j+1;

      // create a point here
      x.emplace_back((_rmax + circle[ix]) * ct);
      x.emplace_back((_rmax + circle[ix]) * st);
      x.emplace_back(circle[iy]);
    }
  }

  // Create triangles from points
  for (int i=0; i<ts; i++) {
    const Int ir1 = cs*i;
    const Int ir2 = (i==ts-1 ? 0 : cs*(i+1));
    for (int j=0; j<cs-1; j++) {
      // set 1
      idx.emplace_back(ir1+j);
      idx.emplace_back(ir2+j);
      idx.emplace_back(ir1+j+1);
      // set 2
      idx.emplace_back(ir2+j);
      idx.emplace_back(ir2+j+1);
      idx.emplace_back(ir1+j+1);
    }
    // set 1
    idx.emplace_back(ir1+cs-1);
    idx.emplace_back(ir2+cs-1);
    idx.emplace_back(ir1);
    // set 2
    idx.emplace_back(ir2+cs-1);
    idx.emplace_back(ir2);
    idx.emplace_back(ir1);
  }

  // keep val empty
  ElementPacket<float> epack {x, idx, val, x.size()/Dimensions, 2};

  return epack;
}


//
// Create a closed discoid
//
ElementPacket<float> generate_discoid(const float _r, const float _h, const float _ips) {

  std::cout << "Creating solid discoid..." << std::endl;

  // set up the element packet
  std::vector<float>   x;
  std::vector<Int>   idx;
  std::vector<float> val;

  // the strategy here is to create a simple discoid, oriented in the xy plane,
  //   with many polygons sharing the center node; then use igl::decimate to 
  //   make better triangles in the middles of the two disc sides

  // num radial nodes and num radial elements
  const size_t nre = 1 + size_t(_r/_ips);
  const size_t nrn = 1 + nre;

  // num circumferential nodes and elements (min is 6)
  const size_t nce = size_t(std::max(6,int(0.5 + M_PI*_r/_ips)));
  const size_t ncn = nce;

  // same for the length of the cylinder (min is 3)
  const size_t nle = size_t(std::max(3,int(0.5 + _h/_ips)));
  const size_t nln = 1 + nle;

  // we do not care that the nodes on the lip (sharp edge) are repeated twice
  //   that is a requirement for watertightness in some software, but not here

  // create nodes for top surface
  x.emplace_back(0.0);
  x.emplace_back(0.0);
  x.emplace_back(_h);
  for (size_t i=1; i<nrn; ++i) {
    const float thisr = _r * (float)i / (float)nre;
    for (size_t j=0; j<ncn; ++j) {
      const float thisa = 2.0*M_PI * (float)j / (float)ncn;
      x.emplace_back(thisr * std::cos(thisa));
      x.emplace_back(thisr * std::sin(thisa));
      x.emplace_back(_h);
    }
  }

  // create elements for top surface
  // inner loop
    for (size_t j=0; j<nce-1; ++j) {
      idx.emplace_back(0);
      idx.emplace_back(j+1);
      idx.emplace_back(j+2);
    }
      // last element
      idx.emplace_back(0);
      idx.emplace_back(nce);
      idx.emplace_back(1);
  // the second to final loops
  for (size_t i=1; i<nre; ++i) {
    const Int iin = 1 + (i-1)*nce;
    const Int iout = 1 + i*nce;
    for (size_t j=0; j<nce-1; ++j) {
      idx.emplace_back(iin+j);
      idx.emplace_back(iout+j);
      idx.emplace_back(iout+1+j);
      idx.emplace_back(iin+j);
      idx.emplace_back(iout+1+j);
      idx.emplace_back(iin+1+j);
    }
      // last elements
      idx.emplace_back(iout-1);
      idx.emplace_back(iout-1+nce);
      idx.emplace_back(iout);
      idx.emplace_back(iout-1);
      idx.emplace_back(iout);
      idx.emplace_back(iin);
  }

  // create nodes for side surface
  Int ifirst = x.size()/Dimensions;
  for (size_t i=0; i<nln; ++i) {
    const float thish = _h * (float)i / (float)nle;
    for (size_t j=0; j<ncn; ++j) {
      const float thisa = 2.0*M_PI * (float)j / (float)ncn;
      x.emplace_back(_r * std::cos(thisa));
      x.emplace_back(_r * std::sin(thisa));
      x.emplace_back(thish);
    }
  }
  // create elements for side surface
  for (size_t i=0; i<nle; ++i) {
    const Int ilo = ifirst + i*nce;
    const Int ihi = ilo + nce;
    for (size_t j=0; j<nce-1; ++j) {
      idx.emplace_back(ilo+j);
      idx.emplace_back(ilo+1+j);
      idx.emplace_back(ihi+1+j);
      idx.emplace_back(ilo+j);
      idx.emplace_back(ihi+1+j);
      idx.emplace_back(ihi+j);
    }
      // last elements
      idx.emplace_back(ihi-1);
      idx.emplace_back(ilo);
      idx.emplace_back(ihi);
      idx.emplace_back(ihi-1);
      idx.emplace_back(ihi);
      idx.emplace_back(ihi+nce-1);
  }

  // create nodes for bottom surface
  ifirst = x.size()/Dimensions;
  x.emplace_back(0.0);
  x.emplace_back(0.0);
  x.emplace_back(0.0);
  for (size_t i=1; i<nrn; ++i) {
    const float thisr = _r * (float)i / (float)nre;
    for (size_t j=0; j<ncn; ++j) {
      // NOTE this is negative, so that they go the other way
      const float thisa = -2.0*M_PI * (float)j / (float)ncn;
      x.emplace_back(thisr * std::cos(thisa));
      x.emplace_back(thisr * std::sin(thisa));
      x.emplace_back(0.0);
    }
  }

  // create elements for bottom surface
  // inner loop
    for (size_t j=0; j<nce-1; ++j) {
      idx.emplace_back(ifirst+0);
      idx.emplace_back(ifirst+j+1);
      idx.emplace_back(ifirst+j+2);
    }
      // last element
      idx.emplace_back(ifirst+0);
      idx.emplace_back(ifirst+nce);
      idx.emplace_back(ifirst+1);
  // the second to final loops
  for (size_t i=1; i<nre; ++i) {
    const Int iin = ifirst + 1 + (i-1)*nce;
    const Int iout = ifirst + 1 + i*nce;
    for (size_t j=0; j<nce-1; ++j) {
      idx.emplace_back(iin+j);
      idx.emplace_back(iout+j);
      idx.emplace_back(iout+1+j);
      idx.emplace_back(iin+j);
      idx.emplace_back(iout+1+j);
      idx.emplace_back(iin+1+j);
    }
      // last elements
      idx.emplace_back(iout-1);
      idx.emplace_back(iout-1+nce);
      idx.emplace_back(iout);
      idx.emplace_back(iout-1);
      idx.emplace_back(iout);
      idx.emplace_back(iin);
  }


  // keep val empty
  ElementPacket<float> epack {x, idx, val, x.size()/Dimensions, 2};

  // now - very important - use decimation to limit the element count!
  mergeduplicates_geometry(epack,1.e-6);
  const float area = 2.0*M_PI*_r * (_r+_h);
  const size_t maxelems = std::max(size_t(1500), size_t(1.4*area/(_ips*_ips)));
  decimate_geometry(epack, maxelems);

  return epack;
}

