/*
 * IglReadGeom.cpp - Call into igl to read a triangle mesh
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
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

#include "IglReadGeom.h"

#include "igl/read_triangle_mesh.h"

#include <vector>
#include <iostream>

// templatized on float-like storage class
ElementPacket<float> read_geometry_file(const std::string _infile) {

  // generate temporary vectors to accept the mesh data
  std::vector<std::vector<float>> v;
  std::vector<std::vector<Int>> f;

  // call the reader
  std::cout << std::endl << "Reading " << _infile << std::endl;
  if (!igl::read_triangle_mesh<float,Int>(_infile, v, f)) {
    // failed
    std::cout << "  Geometry file is unreadable, abandoning!" << std::endl;
    throw std::runtime_error("  Geometry file is unreadable, abandoning");
    return ElementPacket<float>({std::vector<float>(), std::vector<Int>(), std::vector<float>()});
  }

  // are these all triangles?

  // should we refine these?

  //std::cout << "v is " << v.size() << std::endl;
  //std::cout << "v[0] is " << v[0].size() << std::endl;
  //std::cout << "f is " << f.size() << std::endl;
  //std::cout << "f[0] is " << f[0].size() << std::endl;

  // create the flattened arrays
  const size_t num_nodes = v.size();
  const size_t num_panels = f.size();
  std::cout << "  Igl read " << num_nodes << " nodes and " << num_panels << " triangles" << std::endl;

  std::vector<float> x(num_nodes*3);
  std::vector<Int>   idx(num_panels*3);
  std::vector<float> val;

  // and fill them up
  for (size_t i=0; i<v.size(); ++i) {
    x[3*i+0] = v[i][0];
    x[3*i+1] = v[i][1];
    x[3*i+2] = v[i][2];
    //std::cout << "node " << i << " at " << x[3*i+0] << " " << x[3*i+1] << " " << x[3*i+2] << std::endl;
  }
  for (size_t i=0; i<f.size(); ++i) {
    idx[3*i+0] = f[i][0];
    idx[3*i+1] = f[i][1];
    idx[3*i+2] = f[i][2];
    //std::cout << "tri " << i << " at " << idx[3*i+0] << " " << idx[3*i+1] << " " << idx[3*i+2] << std::endl;
  }

  // do not set val! let caller specialize
  //std::fill(val.begin(), val.end(), 0.0);

  return ElementPacket<float>({x, idx, val, (size_t)(num_panels), (uint8_t)2});
}
