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



  // keep val empty
  ElementPacket<float> epack {x, idx, val, x.size()/Dimensions, 2};

  return epack;
}

