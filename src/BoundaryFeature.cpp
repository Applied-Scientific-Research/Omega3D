/*
 * BoundaryFeature.cpp - GUI-side descriptions of boundary features
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#include "BoundaryFeature.h"
#include "GuiHelper.h"
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
  else if (ftype == "quad")   { _flist.emplace_back(std::make_unique<BoundaryQuad>(_bp)); }
  // otherwise it's a file to read
  else                        { _flist.emplace_back(std::make_unique<ExteriorFromFile>(_bp)); }

  // and pass the json object to the specific parser
  _flist.back()->from_json(_jin);

  std::cout << "  finished " << _flist.back()->to_string() << std::endl;
}

#ifdef USE_IMGUI
int BoundaryFeature::obj_movement_gui(int &mitem, char* strx, char* stry, char* strz,
                                       char* strrx, char* strry, char* strrz) {
  const char* mitems[] = { "fixed to ground", "attached to previous", "according to formula" };
  int changed = 0;
  static int tmp = -1;
  //const char* mitems [] = { "fixed", "attached to previous", "according to formula", "dynamic" };
  ImGui::Combo("movement", &mitem, mitems, 3);
  if (tmp != mitem) {
    tmp = mitem;
    changed += 1;
  }

  if (mitem == 2) {
    ImGui::InputText("x position", strx, 512);
    ImGui::SameLine();
    ShowHelpMarker("Position as constant or function of time t. Use C-style expressions: + - / * % ^ ( ) pi e abs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
    ImGui::InputText("y position", stry, 512);
    ImGui::InputText("z position", strz, 512);
    ImGui::InputText("x rotation", strrx, 512);
    ImGui::SameLine();
    ShowHelpMarker("Rotation as axis times angle (in radians) as constant or function of time t. Use C-style expressions: + - / * % ^ ( ) pi e abs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
    ImGui::InputText("y rotation", strry, 512);
    ImGui::InputText("z rotation", strrz, 512);
  }

  return changed;
}

// 0 means keep open, 1 means create, 2 means cancel
int BoundaryFeature::draw_creation_gui(std::vector<std::unique_ptr<BoundaryFeature>>& _bfs, Simulation& _sim) {
  // define movement first
  static int mitem = 0;
  static char strx[512] = "0.0*t";
  static char stry[512] = "0.0*t";
  static char strz[512] = "0.0*t";
  static char strrx[512] = "0.0*t";
  static char strry[512] = "0.0*t";
  static char strrz[512] = "0.0*t";
  int changed = BoundaryFeature::obj_movement_gui(mitem, strx, stry, strz, strrx, strry, strrz);

  // static bp prevents a bunch of pointers from being created during the same boundary creation
  // The switch prevents constant assignment (mainly to prevent the terminal from being flooded from messages)
  static std::shared_ptr<Body> bp = nullptr;
  if (changed) {
    switch(mitem) {
      case 0:
        // this geometry is fixed (attached to inertial)
        bp = _sim.get_pointer_to_body("ground");
        break;
      case 1:
        // this geometry is attached to the previous geometry (or ground)
        bp = _sim.get_last_body();
        break;
      case 2:
        // this geometry is attached to a new moving body
        bp = std::make_shared<Body>();
        bp->set_pos(0, std::string(strx));
        bp->set_pos(1, std::string(stry));
        bp->set_pos(2, std::string(strz));
        bp->set_rot(0, std::string(strrx));
        bp->set_rot(1, std::string(strry));
        bp->set_rot(2, std::string(strrz));
        break;
    }
  }

  // define geometry second
  static int item = 0;
  static int oldItem = -1;
  const char* items[] = { "sphere", "rectangle", "boundary quad", "from file" };
  ImGui::Spacing();
  ImGui::Combo("geometry type", &item, items, 4);


  // show different inputs based on what is selected
  static std::unique_ptr<BoundaryFeature> bf = nullptr;
  if (oldItem != item) {
    switch(item) {
      case 0: {
        bf = std::make_unique<Ovoid>(bp);
      } break;
      case 1: {
        bf = std::make_unique<SolidRect>(bp);
      } break;
      case 2: {
        bf = std::make_unique<BoundaryQuad>(bp);
      } break;
      case 3: {
        bf = std::make_unique<ExteriorFromFile>(bp);
      } break;
    } // end switch for geometry
    oldItem = item;
  }

  int created = 0;
  if (bf->draw_info_gui("Add")) {
    if (!bp) { abort(); }
    if (mitem == 2) {
      bp->set_name(bf->to_short_string());
      _sim.add_body(bp);
    }
    bf->set_body(bp);
    bf->generate_draw_geom();
    _bfs.emplace_back(std::move(bf));
    bf = nullptr;
    oldItem = -1;
    created = 1;
  }

  ImGui::SameLine();
  if (ImGui::Button("Cancel", ImVec2(120,0))) {
    oldItem = -1;
    created = 2;
    bf = nullptr;
  }

  return created;
}
#endif

//
// Create a closed spherical/ovoid object from a icosahedron
//
ElementPacket<float>
Ovoid::init_elements(const float _ips) const {

  std::cout << "Creating ovoid..." << std::endl;

  // first, make an icosahedron
  std::vector<float>   x = ico0;
  std::vector<Int>   idx = ico0idx;
  std::vector<float> val;
  ElementPacket<float> epack {x, idx, val, val.size()/Dimensions, 2};

  // estimate the triangle spacing for the scaled ovoid
  float maxscale = std::max(m_sx, std::max(m_sy, m_sz));
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

  // scale and translate here
  for (size_t i=0; i<epack.x.size()/Dimensions; ++i) {
    const float in_x = epack.x[Dimensions*i];
    const float in_y = epack.x[Dimensions*i+1];
    const float in_z = epack.x[Dimensions*i+2];

    // first scale, then translate
    epack.x[Dimensions*i]   = in_x * m_sx + m_x;
    epack.x[Dimensions*i+1] = in_y * m_sy + m_y;
    epack.x[Dimensions*i+2] = in_z * m_sz + m_z;
  }

  // finally, assume standard behavior: reactive, zero-flow panels
  const size_t nsurfs = epack.idx.size() / Dimensions;
  epack.val.resize(Dimensions*nsurfs);
  std::fill(epack.val.begin(), epack.val.end(), 0.0);
  epack.nelem = epack.val.size()/Dimensions;
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

void Ovoid::generate_draw_geom() {
  m_draw = init_elements(0.125);
  m_draw.val.resize(m_draw.val.size()/Dimensions);
}

#ifdef USE_IMGUI
bool Ovoid::draw_info_gui(const std::string _action) {
  float xc[3] = {m_x, m_y, m_z};
  float scale = m_sx;
  std::string buttonText = _action+" spherical body";
  bool add = false;

  // create a solid sphere
  ImGui::InputFloat3("center", xc);
  ImGui::SliderFloat("diameter", &scale, 0.01f, 10.0f, "%.4f", 2.0);
  ImGui::Spacing();
  ImGui::TextWrapped("This feature will add a solid spherical body centered at the given coordinates");
  ImGui::Spacing();
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  m_z = xc[2];
  m_sx = m_sy = m_sz = scale;
  return add;
}
#endif

//
// Create a closed rectangular object
//
ElementPacket<float>
SolidRect::init_elements(const float _ips) const {

  std::cout << "Creating solid rectangle..." << std::endl;

  // first, make 12-triangle rectangle
  std::vector<float>   x = cube0;
  std::vector<Int>   idx = cube0idx;
  std::vector<float> val;
  ElementPacket<float> epack {x, idx, val, val.size()/Dimensions, 2};

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
  for (size_t i=0; i<epack.x.size()/Dimensions; ++i) {
    const float in_x = epack.x[Dimensions*i];
    const float in_y = epack.x[Dimensions*i+1];
    const float in_z = epack.x[Dimensions*i+2];

    // first scale, then translate
    epack.x[Dimensions*i]   = in_x * m_sx + m_x;
    epack.x[Dimensions*i+1] = in_y * m_sy + m_y;
    epack.x[Dimensions*i+2] = in_z * m_sz + m_z;
  }

  // finally, assume standard behavior: reactive, zero-flow panels
  const size_t nsurfs = epack.idx.size() / Dimensions;
  epack.val.resize(Dimensions*nsurfs);
  std::fill(epack.val.begin(), epack.val.end(), 0.0);
  epack.nelem = epack.val.size()/Dimensions;

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

void SolidRect::generate_draw_geom() {
  m_draw = init_elements(1);
  m_draw.val.resize(m_draw.val.size()/Dimensions);
}

#ifdef USE_IMGUI
bool SolidRect::draw_info_gui(const std::string _action) {
  //static bool external_flow = true;
  static float xc[3] = {m_x, m_y, m_z};
  static float xs[3] = {m_sx, m_sy, m_sz};
  //static float rotdeg = 0.0f;
  std::string buttonText = _action+" rectangular body";
  bool add = false;
  
  // create a solid rectangle
  ImGui::InputFloat3("center", xc);
  ImGui::InputFloat3("side lengths", xs);
  //ImGui::SliderFloat("orientation", &rotdeg, 0.0f, 89.0f, "%.0f");
  //ImGui::SliderAngle("orientation", &rotdeg);
  ImGui::Spacing();
  ImGui::TextWrapped("This feature will add a solid rectangular body centered at the given coordinates");
  ImGui::Spacing();
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  m_z = xc[2];
  m_sx = xs[0];
  m_sy = xs[1];
  m_sz = xs[2];
  return add;
}
#endif

//
// Create a triangulated quad of a solid boundary
//
ElementPacket<float>
BoundaryQuad::init_elements(const float _ips) const {

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
  std::vector<float>   x(Dimensions*(num1+1)*(num2+1));
  std::vector<Int>   idx(num_panels*Dimensions);
  std::vector<float> val(num_panels*Dimensions);

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
      idx[Dimensions*icnt]   = idx1 + j + 0;
      idx[Dimensions*icnt+1] = idx2 + j + 0;
      idx[Dimensions*icnt+2] = idx2 + j + 1;
      icnt++;
      // do the other triangle
      idx[Dimensions*icnt]   = idx1 + j + 0;
      idx[Dimensions*icnt+1] = idx2 + j + 1;
      idx[Dimensions*icnt+2] = idx1 + j + 1;
      icnt++;
    }
  }

  // all triangles share the same boundary condition velocity (in global coords)
  for (size_t i=0; i<num_panels; ++i) {
    val[Dimensions*i+0] = m_bcx;
    val[Dimensions*i+1] = m_bcy;
    val[Dimensions*i+2] = m_bcz;
  }

  return ElementPacket<float>({x, idx, val, val.size()/Dimensions, 2});
}

void
BoundaryQuad::debug(std::ostream& os) const {
  os << to_string();
}

std::string
BoundaryQuad::to_string() const {
  std::stringstream ss;
  ss << "quad at " << m_x << " " << m_y << " " << m_z;
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
  //if (j["scale"].is_array()) {
    //const std::vector<float> sc = j["scale"].get<std::vector<float>>();
  const std::vector<float> p1 = j["p1"];
  m_x1 = p1[0];
  m_y1 = p1[1];
  m_z1 = p1[2];
  const std::vector<float> p2 = j["p2"];
  m_x2 = p2[0];
  m_y2 = p2[1];
  m_z2 = p2[2];
  const std::vector<float> p3 = j["p3"];
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
  mesh["geometry"] = "quad";
  mesh["startpt"] = {m_x, m_y, m_z};
  mesh["p1"] = {m_x1, m_y1, m_z1};
  mesh["p2"] = {m_x2, m_y2, m_z2};
  mesh["p3"] = {m_x3, m_y3, m_z3};
  mesh["vel"] = {m_bcx, m_bcy, m_bcz};
  //mesh["external"] = m_external;
  return mesh;
}

void BoundaryQuad::generate_draw_geom() {
  m_draw = init_elements(1);
  m_draw.val.resize(m_draw.val.size()/Dimensions);
}

#ifdef USE_IMGUI
bool BoundaryQuad::draw_info_gui(const std::string _action) {
  float xc[3] = {m_x, m_y, m_z};
  float x1[3] = {m_x1, m_y1, m_z1};
  float x2[3] = {m_x2, m_y2, m_z2};
  float x3[3] = {m_x3, m_y3, m_z3};
  float v[3] = {m_bcx, m_bcy, m_bcz};
  std::string buttonText = _action+" boundary quad";
  bool add = false;
  
  // create a solid rectangle
  ImGui::InputFloat3("Starting point", xc);
  ImGui::InputFloat3("Point 1", x1);
  ImGui::InputFloat3("Point 2", x2);
  ImGui::InputFloat3("Point 3", x3);
  ImGui::InputFloat3("Velocities", v);
  ImGui::Spacing();
  ImGui::TextWrapped("This feature will add a boundary quad");
  ImGui::Spacing();
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  m_z = xc[2];
  m_x1 = x1[0];
  m_y1 = x1[1];
  m_z1 = x1[2];
  m_x2 = x2[0];
  m_y2 = x2[1];
  m_z2 = x2[2];
  m_x3 = x3[0];
  m_y3 = x3[1];
  m_z3 = x3[2];
  m_bcx = v[0];
  m_bcy = v[1];
  m_bcz = v[2];

  return add;
}
#endif

//
// Create a closed object from a geometry file (fluid is outside)
//
ElementPacket<float>
ExteriorFromFile::init_elements(const float _ips) const {

  std::cout << "Creating object from file..." << std::endl;

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

void ExteriorFromFile::generate_draw_geom() {
  m_draw = init_elements(1);
  m_draw.val.resize(m_draw.val.size()/Dimensions);
}

#ifdef USE_IMGUI
bool ExteriorFromFile::draw_info_gui(const std::string _action) {
  // load a geometry file
  static std::string infile = m_infile;
  static std::string shortname = infile;
  static std::vector<std::string> recent_geom_files;
  static bool show_geom_input_window = false;
  if (ImGui::Button("Geometry file", ImVec2(200,0))) show_geom_input_window = true;

  if (show_geom_input_window) {
    bool try_it = false;

    if (fileIOWindow( try_it, infile, recent_geom_files, "Open", {"*.obj", "*.stl", "*.ply", "*.*"}, true, ImVec2(500,250))) {
      show_geom_input_window = false;

      if (try_it and !infile.empty()) {
        // remember
        recent_geom_files.push_back( infile );

        // now remove the leading directories from the string
        const size_t lastchar = infile.find_last_of("/\\");
        shortname = infile.substr(lastchar+1);
      }
    }
  }
  
  //static bool external_flow = true;
  static float xc[3] = {m_x, m_y, m_z};
  //static float rotdeg = 0.0f;
  static float scale = m_sx;
  std::string buttonText = _action+" geometry from file";
  bool add = false;
  
  ImGui::SameLine();
  ImGui::Text(shortname.c_str());
  ImGui::InputFloat3("center", xc);
  ImGui::SliderFloat("scale", &scale, 0.01f, 10.0f, "%.4f", 2.0);
  ImGui::Spacing();
  ImGui::TextWrapped("This feature will add a solid body centered at the given coordinates");
  ImGui::Spacing();
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  m_z = xc[2];
  m_sx = m_sy = m_sz = scale;
  return add;
}
#endif
