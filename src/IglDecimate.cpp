/*
 * IglDecimate.cpp - Call into igl to decimate a triangle mesh
 *
 * (c)2019,21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#include "IglDecimate.h"

#define IGL_PARALLEL_FOR_FORCE_SERIAL
#include "igl/decimate.h"

#include <Eigen/Core>

#include <vector>
#include <iostream>

// templatized on float-like storage class
void decimate_geometry(ElementPacket<float>& _mesh, const size_t _nend) {

  const bool debug = true;

  // convert to igl representation
  Eigen::MatrixXd V(_mesh.x.size()/3, 3);
  for (size_t i=0; i<_mesh.x.size()/3; ++i) {
    V(i,0) = _mesh.x[3*i+0];
    V(i,1) = _mesh.x[3*i+1];
    V(i,2) = _mesh.x[3*i+2];
  }
  if (debug) std::cout << "  input V is " << V.rows() << " by " << V.cols() << std::endl;

  Eigen::MatrixXi F(_mesh.idx.size()/3, 3);
  for (size_t i=0; i<_mesh.idx.size()/3; ++i) {
    F(i,0) = _mesh.idx[3*i+0];
    F(i,1) = _mesh.idx[3*i+1];
    F(i,2) = _mesh.idx[3*i+2];
  }
  if (debug) std::cout << "  input F is " << F.rows() << " by " << F.cols() << std::endl;

  Eigen::VectorXi J;
  Eigen::VectorXi I;

  // decimate (simple)
  if (debug) std::cout << "  DECIMATING\n";
  igl::decimate( Eigen::MatrixXd(V), Eigen::MatrixXi(F), _nend, V,F, J,I);

  // check results
  if (debug) std::cout << "  output V is " << V.rows() << " by " << V.cols() << std::endl;
  if (debug) std::cout << "  output F is " << F.rows() << " by " << F.cols() << std::endl;
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
