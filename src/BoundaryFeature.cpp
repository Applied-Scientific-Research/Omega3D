/*
 * BoundaryFeature.cpp - GUI-side descriptions of boundary features
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "BoundaryFeature.h"
#include "IglReadGeom.h"
#include "IglRefine.h"
#include "Icosahedron.h"
#include "Cube.h"

#include <cmath>
#include <iostream>
#include <sstream>

// write out any object of parent type BoundaryFeature by dispatching to appropriate "debug" method
std::ostream& operator<<(std::ostream& os, BoundaryFeature const& ff) { 
  ff.debug(os);
  return os;
}


//
// parse the json and dispatch the constructors
//
void parse_boundary_json(std::vector<std::unique_ptr<BoundaryFeature>>& _flist,
                         std::shared_ptr<Body> _bp,
                         const nlohmann::json _jin) {

  // must have one and only one type
  if (_jin.count("geometry") != 1) return;

  const std::string ftype = _jin["geometry"];

  // if it's a recognized keyword, generate the geometry
  if      (ftype == "sphere") { _flist.emplace_back(std::make_unique<Ovoid>(_bp)); }
  else if (ftype == "rect")   { _flist.emplace_back(std::make_unique<SolidRect>(_bp)); }
  // otherwise it's a file to read
  else                        { _flist.emplace_back(std::make_unique<ExteriorFromFile>(_bp)); }

  // and pass the json object to the specific parser
  _flist.back()->from_json(_jin);

  std::cout << "  found " << _flist.back()->to_string() << std::endl;
}


//
// Create a closed spherical/ovoid object from a icosahedron
//
ElementPacket<float>
Ovoid::init_elements(const float _ips) const {

  if (not this->is_enabled()) return ElementPacket<float>();

  // first, make an icosahedron
  std::vector<float>   x = ico0;
  std::vector<Int>   idx = ico0idx;
  std::vector<float> val;
  ElementPacket<float> epack {x, idx, val};

  // estimate the triangle spacing for the scaled ovoid
  float maxscale = std::max(m_sx, std::max(m_sy, m_sz));
  float meansize = 0.35 * maxscale;

  // then, iteratively refine it
  std::cout << "  sphere is icosahedron with 20";
  while (meansize > 1.2*_ips) {

    // refine once
    refine_geometry(epack);

    // re-sphericalize (r=0.5)
    for (size_t i=0; i<epack.x.size()/3; ++i) {
      const float rad = std::sqrt( std::pow(epack.x[3*i], 2) +
                                   std::pow(epack.x[3*i+1], 2) +
                                   std::pow(epack.x[3*i+2], 2));
      const float scale = 0.5 / rad;
      epack.x[3*i]   *= scale;
      epack.x[3*i+1] *= scale;
      epack.x[3*i+2] *= scale;
    }

    meansize *= 0.5;
  }
  std::cout << " panels" << std::endl;

  // scale and translate here
  for (size_t i=0; i<epack.x.size()/3; ++i) {
    const float in_x = epack.x[3*i];
    const float in_y = epack.x[3*i+1];
    const float in_z = epack.x[3*i+2];

    // first scale, then translate
    epack.x[3*i]   = in_x * m_sx + m_x;
    epack.x[3*i+1] = in_y * m_sy + m_y;
    epack.x[3*i+2] = in_z * m_sz + m_z;
  }

  // finally, assume standard behavior: reactive, zero-flow panels
  const size_t nsurfs = epack.idx.size() / 3;
  epack.val.resize(3*nsurfs);
  std::fill(epack.val.begin(), epack.val.end(), 0.0);

  return epack;
}

void
Ovoid::debug(std::ostream& os) const {
  os << to_string();
}

std::string
Ovoid::to_string() const {
  std::stringstream ss;

  if (m_sx == m_sy and m_sy == m_sz) {
    ss << "sphere at " << m_x << " " << m_y << " " << m_z << " scaled by " << m_sx;
  } else {
    ss << "ovoid at " << m_x << " " << m_y << " " << m_z << " scaled by " << m_sx << " " << m_sy << " " << m_sz;
  }

  return ss.str();
}

void
Ovoid::from_json(const nlohmann::json j) {
  const std::vector<float> tr = j["translation"];
  m_x = tr[0];
  m_y = tr[1];
  m_z = tr[2];
  if (j["scale"].is_array()) {
    const std::vector<float> sc = j["scale"].get<std::vector<float>>();
    m_sx = sc[0];
    m_sy = sc[1];
    m_sz = sc[2];
  } else if (j["scale"].is_number()) {
    const float sc = j["scale"].get<float>();
    m_sx = sc;
    m_sy = sc;
    m_sz = sc;
  }
  m_external = j.value("external", true);
}

nlohmann::json
Ovoid::to_json() const {
  // make an object for the mesh
  nlohmann::json mesh = nlohmann::json::object();
  mesh["geometry"] = "sphere";
  mesh["translation"] = {m_x, m_y, m_z};
  mesh["scale"] = {m_sx, m_sy, m_sz};
  mesh["external"] = m_external;
  return mesh;
}


//
// Create a closed rectangular object
//
ElementPacket<float>
SolidRect::init_elements(const float _ips) const {

  if (not this->is_enabled()) return ElementPacket<float>();

  // first, make 12-triangle rectangle
  std::vector<float>   x = cube0;
  std::vector<Int>   idx = cube0idx;
  std::vector<float> val;
  ElementPacket<float> epack {x, idx, val};

  // estimate the triangle spacing for the scaled rectangle
  float maxscale = std::max(m_sx, std::max(m_sy, m_sz));
  float meansize = 0.7 * maxscale;

  // then, iteratively refine it
  std::cout << "  rectangle is cube with 12";
  while (meansize > 1.2*_ips) {
    refine_geometry(epack);
    meansize *= 0.5;
  }
  std::cout << " panels" << std::endl;

  // scale and translate here
  for (size_t i=0; i<epack.x.size()/3; ++i) {
    const float in_x = epack.x[3*i];
    const float in_y = epack.x[3*i+1];
    const float in_z = epack.x[3*i+2];

    // first scale, then translate
    epack.x[3*i]   = in_x * m_sx + m_x;
    epack.x[3*i+1] = in_y * m_sy + m_y;
    epack.x[3*i+2] = in_z * m_sz + m_z;
  }

  // finally, assume standard behavior: reactive, zero-flow panels
  const size_t nsurfs = epack.idx.size() / 3;
  epack.val.resize(3*nsurfs);
  std::fill(epack.val.begin(), epack.val.end(), 0.0);

  return epack;
}

void
SolidRect::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SolidRect::to_string() const {
  std::stringstream ss;

  if (m_sx == m_sy and m_sy == m_sz) {
    ss << "cube at " << m_x << " " << m_y << " " << m_z << " scaled by " << m_sx;
  } else {
    ss << "rectangle at " << m_x << " " << m_y << " " << m_z << " scaled by " << m_sx << " " << m_sy << " " << m_sz;
  }

  return ss.str();
}

void
SolidRect::from_json(const nlohmann::json j) {
  const std::vector<float> tr = j["translation"];
  m_x = tr[0];
  m_y = tr[1];
  m_z = tr[2];
  if (j["scale"].is_array()) {
    const std::vector<float> sc = j["scale"].get<std::vector<float>>();
    m_sx = sc[0];
    m_sy = sc[1];
    m_sz = sc[2];
  } else if (j["scale"].is_number()) {
    const float sc = j["scale"].get<float>();
    m_sx = sc;
    m_sy = sc;
    m_sz = sc;
  }
  m_external = j.value("external", true);
}

nlohmann::json
SolidRect::to_json() const {
  // make an object for the mesh
  nlohmann::json mesh = nlohmann::json::object();
  mesh["geometry"] = "rect";
  mesh["translation"] = {m_x, m_y, m_z};
  mesh["scale"] = {m_sx, m_sy, m_sz};
  mesh["external"] = m_external;
  return mesh;
}


//
// Create a triangulated quad of a solid boundary
//
ElementPacket<float>
BoundaryQuad::init_elements(const float _ips) const {

  if (not this->is_enabled()) return ElementPacket<float>();

  // how many panels?
  const float side1 = 0.5 * (std::sqrt(std::pow(m_x1-m_x,2)  + std::pow(m_y1-m_y,2)  + std::pow(m_z1-m_z,2))
                            +std::sqrt(std::pow(m_x3-m_x2,2) + std::pow(m_y3-m_y2,2) + std::pow(m_z3-m_z2,2)));
  const float side2 = 0.5 * (std::sqrt(std::pow(m_x2-m_x1,2) + std::pow(m_y2-m_y1,2) + std::pow(m_z2-m_z1,2))
                            +std::sqrt(std::pow(m_x-m_x3,2)  + std::pow(m_y-m_y3,2)  + std::pow(m_z-m_z3,2)));
  const size_t num1 = std::max(1, (int)(side1 / _ips));
  const size_t num2 = std::max(1, (int)(side2 / _ips));
  const size_t num_panels = 2 * num1 * num2;

  std::cout << "Creating quad with " << num_panels << " panels" << std::endl;
  std::cout << "  " << to_string() << std::endl;

  // created once
  std::vector<float>   x(3*(num1+1)*(num2+1));
  std::vector<Int>   idx(num_panels*3);
  std::vector<float> val(num_panels*3);

  // when quad nodes appear in CCW order, the viewer is "outside" the object (inside the fluid)
  // so go CW around the body
  size_t icnt = 0;
  for (size_t i=0; i<num1+1; ++i) {
    const float s = (float)i / (float)num1;
    for (size_t j=0; j<num2+1; ++j) {
      const float t = (float)j / (float)num2;
      const float w0 = (1.0-s)*(1.0-t);
      const float w1 =      s *(1.0-t);
      const float w2 =      s *     t ;
      const float w3 = (1.0-s)*     t ;
      x[icnt++] = w0*m_x + w1*m_x1 + w2*m_x2 + w3*m_x3;
      x[icnt++] = w0*m_y + w1*m_y1 + w2*m_y2 + w3*m_y3;
      x[icnt++] = w0*m_z + w1*m_z1 + w2*m_z2 + w3*m_z3;
    }
  }

  // outside is to the left walking from one point to the next
  icnt = 0;
  for (size_t i=0; i<num1; ++i) {
    const size_t idx1 = (num2+1) * i;
    const size_t idx2 = (num2+1) * (i+1);
    for (size_t j=0; j<num2; ++j) {
      // make the first of 2 tris
      idx[3*icnt]   = idx1 + j + 0;
      idx[3*icnt+1] = idx2 + j + 1;
      idx[3*icnt+2] = idx2 + j + 0;
      icnt++;
      // do the other triangle
      idx[3*icnt]   = idx1 + j + 0;
      idx[3*icnt+1] = idx1 + j + 1;
      idx[3*icnt+2] = idx2 + j + 1;
      icnt++;
    }
  }

  // all triangles share the same boundary condition velocity (in global coords)
  for (size_t i=0; i<num_panels; ++i) {
    val[3*i+0] = m_bcx;
    val[3*i+1] = m_bcy;
    val[3*i+2] = m_bcz;
  }

  return ElementPacket<float>({x, idx, val});
}

void
BoundaryQuad::debug(std::ostream& os) const {
  os << to_string();
}

std::string
BoundaryQuad::to_string() const {
  std::stringstream ss;
  ss << "quad at " << m_x << " " << m_y << " " << m_x;
  if (std::abs(m_bcx)+std::abs(m_bcy)+std::abs(m_bcz) > std::numeric_limits<float>::epsilon()) {
    ss << " with vel " << m_bcx << " " << m_bcy << " " << m_bcz;
  }
  return ss.str();
}

void
BoundaryQuad::from_json(const nlohmann::json j) {
  const std::vector<float> tr = j["startpt"];
  m_x = tr[0];
  m_y = tr[1];
  m_z = tr[2];
  const std::vector<float> p1 = j["pt1"];
  m_x1 = p1[0];
  m_y1 = p1[1];
  m_z1 = p1[2];
  const std::vector<float> p2 = j["pt2"];
  m_x2 = p2[0];
  m_y2 = p2[1];
  m_z2 = p2[2];
  const std::vector<float> p3 = j["pt3"];
  m_x3 = p3[0];
  m_y3 = p3[1];
  m_z3 = p3[2];
  const std::vector<float> bv = j["vel"];
  m_bcx = bv[0];
  m_bcy = bv[1];
  m_bcz = bv[2];
  m_external = true;//j.value("external", true);
}

nlohmann::json
BoundaryQuad::to_json() const {
  // make an object for the mesh
  nlohmann::json mesh = nlohmann::json::object();
  mesh["geometry"] = "segment";
  mesh["startpt"] = {m_x, m_y, m_z};
  mesh["p1"] = {m_x1, m_y1, m_z1};
  mesh["p2"] = {m_x2, m_y2, m_z2};
  mesh["p3"] = {m_x3, m_y3, m_z3};
  mesh["vel"] = {m_bcx, m_bcy, m_bcz};
  //mesh["external"] = m_external;
  return mesh;
}


//
// Create a closed object from a geometry file (fluid is outside)
//
ElementPacket<float>
ExteriorFromFile::init_elements(const float _ips) const {

  if (not this->is_enabled()) return ElementPacket<float>();

  ElementPacket<float> epack = read_geometry_file(m_infile);

  // assume standard behavior: reactive, zero-flow panels
  const size_t nsurfs = epack.idx.size() / 3;
  epack.val.resize(3*nsurfs);
  std::fill(epack.val.begin(), epack.val.end(), 0.0);

  // scale and translate here
  for (size_t i=0; i<epack.x.size()/3; ++i) {
    const float in_x = epack.x[3*i];
    const float in_y = epack.x[3*i+1];
    const float in_z = epack.x[3*i+2];

    // first scale, then translate
    epack.x[3*i]   = in_x * m_sx + m_x;
    epack.x[3*i+1] = in_y * m_sy + m_y;
    epack.x[3*i+2] = in_z * m_sz + m_z;
  }

  return epack;
}

void
ExteriorFromFile::debug(std::ostream& os) const {
  os << to_string();
}

std::string
ExteriorFromFile::to_string() const {
  std::stringstream ss;

  // shorten the filename
  const size_t lastchar = m_infile.find_last_of("/\\");

  ss << m_infile.substr(lastchar+1) << " at " << m_x << " " << m_y << " " << m_z << " scaled by " << m_sx << " " << m_sy << " " << m_sz;
  return ss.str();
}

void
ExteriorFromFile::from_json(const nlohmann::json j) {

  // be smarter about the string - this allows the code to compile with GCC 7.2.0
  nlohmann::json infile_object = j["geometry"];
  if (infile_object.is_string()) {
    m_infile = infile_object.get<std::string>();
  }
  // this used to be:
  //m_infile = j["geometry"];

  const std::vector<float> tr = j["translation"];
  m_x = tr[0];
  m_y = tr[1];
  m_z = tr[2];
  if (j["scale"].is_array()) {
    const std::vector<float> sc = j["scale"].get<std::vector<float>>();
    m_sx = sc[0];
    m_sy = sc[1];
    m_sz = sc[2];
  } else if (j["scale"].is_number()) {
    const float sc = j["scale"].get<float>();
    m_sx = sc;
    m_sy = sc;
    m_sz = sc;
  }
  m_external = j.value("external", true);
}

nlohmann::json
ExteriorFromFile::to_json() const {
  // make an object for the mesh
  nlohmann::json mesh = nlohmann::json::object();
  mesh["geometry"] = m_infile;
  mesh["translation"] = {m_x, m_y, m_z};
  mesh["scale"] = {m_sx, m_sy, m_sz};
  mesh["external"] = m_external;
  return mesh;
}

