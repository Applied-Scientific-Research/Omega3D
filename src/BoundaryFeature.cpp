/*
 * BoundaryFeature.cpp - GUI-side descriptions of boundary features
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "BoundaryFeature.h"

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
  //if      (ftype == "sphere") { _flist.emplace_back(std::make_unique<SolidSphere>(_bp)); }
  //else if (ftype == "cube")   { _flist.emplace_back(std::make_unique<SolidCube>(_bp)); }
  // otherwise it's a file to read
  //else                        { _flist.emplace_back(std::make_unique<ExteriorFromFile>(_bp)); }

  // but until those are supported, always call the reader
  _flist.emplace_back(std::make_unique<ExteriorFromFile>(_bp));

  // and pass the json object to the specific parser
  _flist.back()->from_json(_jin);
}


//
// Create a closed object from a geometry file (fluid is outside)
//
ElementPacket<float>
ExteriorFromFile::init_elements(const float _ips) const {

  // how many panels?
  const size_t num_panels = 1000;

  std::cout << "Reading file with " << num_panels << " panels" << std::endl;

  // created once
  std::vector<float>   x(num_panels*3);
  std::vector<Int>   idx(num_panels*3);
  std::vector<float> val(num_panels);

  // outside is to the left walking from one point to the next
  // so go CW around the circle starting at theta=0 (+x axis)
  for (size_t i=0; i<num_panels; i++) {
    x[3*i]     = m_x + std::cos(2.0 * M_PI * (float)i / (float)num_panels);
    x[3*i+1]   = m_y - std::sin(2.0 * M_PI * (float)i / (float)num_panels);
    x[3*i+2]   = m_z - std::sin(2.0 * M_PI * (float)i / (float)num_panels);
    idx[3*i]   = i;
    idx[3*i+1] = i+1;
    idx[3*i+2] = i+2;
    val[i]     = 0.0;
  }

  return ElementPacket<float>({x, idx, val});
}

void
ExteriorFromFile::debug(std::ostream& os) const {
  os << to_string();
}

std::string
ExteriorFromFile::to_string() const {
  std::stringstream ss;
  ss << m_infile << " at " << m_x << " " << m_y << " " << m_z << " scaled by " << m_sx << " " << m_sy << " " << m_sz;
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
}

nlohmann::json
ExteriorFromFile::to_json() const {
  // make an object for the mesh
  nlohmann::json mesh = nlohmann::json::object();
  mesh["geometry"] = m_infile;
  mesh["translation"] = {m_x, m_y, m_z};
  mesh["scale"] = {m_sx, m_sy, m_sz};
  return mesh;
}

