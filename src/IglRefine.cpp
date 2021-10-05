/*
 * IglRefine.cpp - Call into igl to refine a triangle mesh
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
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

#include "IglRefine.h"

#define IGL_PARALLEL_FOR_FORCE_SERIAL
#include "igl/upsample.h"

#include <Eigen/Core>

#include <vector>
#include <iostream>

// templatized on float-like storage class
void refine_geometry(ElementPacket<float>& _mesh) {

  // convert to igl representation
  Eigen::MatrixXd V(_mesh.x.size()/3, 3);
  for (size_t i=0; i<_mesh.x.size()/3; ++i) {
    V(i,0) = _mesh.x[3*i+0];
    V(i,1) = _mesh.x[3*i+1];
    V(i,2) = _mesh.x[3*i+2];
  }
  //std::cout << "  input V is " << V.rows() << " by " << V.cols() << std::endl;

  Eigen::MatrixXi F(_mesh.idx.size()/3, 3);
  for (size_t i=0; i<_mesh.idx.size()/3; ++i) {
    F(i,0) = _mesh.idx[3*i+0];
    F(i,1) = _mesh.idx[3*i+1];
    F(i,2) = _mesh.idx[3*i+2];
  }
  //std::cout << "  input F is " << F.rows() << " by " << F.cols() << std::endl;

  // refine (in-plane subdivision)
  igl::upsample( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);

  // check results
  //std::cout << "  output V is " << V.rows() << " by " << V.cols() << std::endl;
  //std::cout << "  output F is " << F.rows() << " by " << F.cols() << std::endl;
  std::cout << " ...refined to " << F.rows();

  // create the flattened arrays
  const size_t num_nodes = V.rows();
  const size_t num_panels = F.rows();
  _mesh.x.resize(num_nodes*3);
  _mesh.idx.resize(num_panels*3);

  // and fill them up
  for (size_t i=0; i<num_nodes; ++i) {
    _mesh.x[3*i+0] = V(i,0);
    _mesh.x[3*i+1] = V(i,1);
    _mesh.x[3*i+2] = V(i,2);
  }
  for (size_t i=0; i<num_panels; ++i) {
    _mesh.idx[3*i+0] = F(i,0);
    _mesh.idx[3*i+1] = F(i,1);
    _mesh.idx[3*i+2] = F(i,2);
  }

  return;
}
