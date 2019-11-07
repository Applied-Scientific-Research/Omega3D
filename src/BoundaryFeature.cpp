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
  std::cout << "  found " << ftype << std::endl;

  // if it's a recognized keyword, generate the geometry
  if      (ftype == "sphere") { _flist.emplace_back(std::make_unique<Ovoid>(_bp)); }
  else if (ftype == "rect")   { _flist.emplace_back(std::make_unique<SolidRect>(_bp)); }
  // otherwise it's a file to read
  else                        { _flist.emplace_back(std::make_unique<ExteriorFromFile>(_bp)); }

  // and pass the json object to the specific parser
  _flist.back()->from_json(_jin);
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
  //mesh["external"] = m_external;
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
  m_infile = j["geometry"];
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
  //mesh["external"] = m_external;
  return mesh;
}

