/*
 * FlowFeature.cpp - GUI-side descriptions of flow features
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#include "BoundaryFeature.h"
#include "FlowFeature.h"
#include "MathHelper.h"
#include "imgui/imgui.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <random>

// write out any object of parent type FlowFeature by dispatching to appropriate "debug" method
std::ostream& operator<<(std::ostream& os, FlowFeature const& ff) {
  ff.debug(os);
  return os;
}


//
// parse the json and dispatch the constructors
//
void parse_flow_json(std::vector<std::unique_ptr<FlowFeature>>& _flist,
                     const nlohmann::json _jin) {

  // must have one and only one type
  if (_jin.count("type") != 1) return;

  const std::string ftype = _jin["type"];

  if      (ftype == "single particle") {  _flist.emplace_back(std::make_unique<SingleParticle>()); }
  else if (ftype == "vortex blob") {      _flist.emplace_back(std::make_unique<VortexBlob>()); }
  else if (ftype == "block of random") {  _flist.emplace_back(std::make_unique<BlockOfRandom>()); }
  else if (ftype == "particle emitter") { _flist.emplace_back(std::make_unique<ParticleEmitter>()); }
  else if (ftype == "singular ring") {    _flist.emplace_back(std::make_unique<SingularRing>()); }
  else if (ftype == "thick ring") {       _flist.emplace_back(std::make_unique<ThickRing>()); }
  else {
    std::cout << "  type " << ftype << " does not name an available flow feature, ignoring" << std::endl;
    return;
  }

  // and pass the json object to the specific parser
  _flist.back()->from_json(_jin);

  std::cout << "  found " << ftype << std::endl;
}

#ifdef USE_IMGUI
bool FlowFeature::draw_creation_gui(std::vector<std::unique_ptr<FlowFeature>> &ffs, const float ips) {
  static int item = 1;
  static int oldItem = -1;
  const char* items[] = { "vortex blob", "random particles", "singular vortex ring", "thick vortex ring" };
  ImGui::Combo("type", &item, items, 4);

  // show different inputs based on what is selected
  static std::unique_ptr<FlowFeature> ff = nullptr;
  if (oldItem != item) {
    switch(item) {
      case 0: {
        ff = std::make_unique<VortexBlob>();
      } break;
      case 1: {
        ff = std::make_unique<BlockOfRandom>();
      } break;
      case 2: {
        ff = std::make_unique<SingularRing>();
      } break;
      case 3: {
        ff = std::make_unique<ThickRing>();
      } break;
    }
    oldItem = item;
  }

  bool created = false;
  if (ff->draw_info_gui("Add", ips)) {
    ff->generate_draw_geom();
    ffs.emplace_back(std::move(ff));
    ff = nullptr;
    created = true;
    oldItem = -1;
    ImGui::CloseCurrentPopup();
  }

  ImGui::SameLine();
  if (ImGui::Button("Cancel", ImVec2(120,0))) {
    oldItem = -1;
    ImGui::CloseCurrentPopup();
  }

  ImGui::EndPopup();
  return created;
}
#endif

//
// important feature: convert flow feature definition into 1-D vector of values
//
// each 7 floats is one particle's: [xyz] location, [xyz] strength, vdelta (radius)
//

//
// drop a single particle
//
ElementPacket<float>
SingleParticle::init_elements(float _ips) const {
  //if (this->is_enabled()) return std::vector<float>({m_x, m_y, m_z, m_sx, m_sy, m_sz, 0.0});
  //else return std::vector<float>();
  std::vector<float> x = {m_x, m_y, m_z};
  std::vector<Int> idx = {};
  std::vector<float> vals = {m_sx, m_sy, m_sz};
  ElementPacket<float> packet({x, idx, vals, (size_t)1, 0});
  if (packet.verify(packet.x.size()+packet.val.size(), 6)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

ElementPacket<float>
SingleParticle::step_elements(float _ips) const {
  return ElementPacket<float>();
}

void
SingleParticle::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SingleParticle::to_string() const {
  std::stringstream ss;
  ss << "single particle at " << m_x << " " << m_y << " " << m_z;
  ss << " with strength " << m_sx << " " << m_sy << " " << m_sz;
  return ss.str();
}

void
SingleParticle::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
  const std::vector<float> s = j["strength"];
  m_sx = s[0];
  m_sy = s[1];
  m_sz = s[2];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
SingleParticle::to_json() const {
  nlohmann::json j;
  j["type"] = "single particle";
  j["center"] = {m_x, m_y, m_z};
  j["strength"] = {m_sx, m_sy, m_sz};
  j["enabled"] = m_enabled;
  return j;
}

//Single Particles cant be made by user
void SingleParticle::generate_draw_geom() {
  //const float diam = 0.01;
  //std::unique_ptr<Ovoid> tmp = std::make_unique<SolidCircle>(nullptr, true, m_x, m_y, m_z,
  //                                                           diam, diam, diam);
  //m_draw = tmp->init_elements(diam/25.0);
  //std::fill(m_draw.val.begin(), m_draw.val.end(), m_str);
}

#ifdef USE_IMGUI
// User can't actually create this
bool SingleParticle::draw_info_gui(const std::string action, const float ips) {
  return false;
}
#endif

//
// make a circular vortex blob with soft transition
//
ElementPacket<float>
VortexBlob::init_elements(float _ips) const {
  // create a new vector to pass on
  std::vector<float> x;
  std::vector<Int> idx;
  std::vector<float> vals;

  // what size 2D integer array will we loop over
  int irad = 1 + (m_rad + 0.5*m_softness) / _ips;
  std::cout << "blob needs " << (-irad) << " to " << irad << " spaces" << std::endl;

  // and a counter for the total circulation
  double tot_wgt = 0.0;

  // loop over integer indices
  for (int i=-irad; i<=irad; ++i) {
  for (int j=-irad; j<=irad; ++j) {
  for (int k=-irad; k<=irad; ++k) {

    // how far from the center are we?
    float dr = std::sqrt((float)(i*i+j*j+k*k)) * _ips;
    if (dr < m_rad + 0.5*m_softness) {

      // create a particle here
      x.emplace_back(m_x + _ips*(float)i);
      x.emplace_back(m_y + _ips*(float)j);
      x.emplace_back(m_z + _ips*(float)k);

      // figure out the strength from another check
      double this_wgt = 1.0;
      if (dr > m_rad - 0.5*m_softness) {
        // create a weaker particle
        this_wgt = 0.5 - 0.5*std::sin(M_PI * (dr - m_rad) / m_softness);
      }
      vals.emplace_back(m_sx * (float)this_wgt);
      vals.emplace_back(m_sy * (float)this_wgt);
      vals.emplace_back(m_sz * (float)this_wgt);
      tot_wgt += this_wgt;

      // this is the radius - still zero for now
      vals.emplace_back(0.0f);
    }
  }
  }
  }

  // finally, normalize all particle strengths so that the whole blob
  //   has exactly the right strength
  std::cout << "blob had " << tot_wgt << " initial circulation" << std::endl;
  double str_scale = 1.0 / tot_wgt;
  for (size_t i=0; i<x.size(); i+=4) {
    vals[i+0] = (float)((double)vals[i+0] * str_scale);
    vals[i+1] = (float)((double)vals[i+1] * str_scale);
    vals[i+2] = (float)((double)vals[i+2] * str_scale);
  }

  ElementPacket<float> packet({x, idx, vals, x.size()/3, 0});
  if (packet.verify(packet.x.size()+packet.val.size(), 7)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

ElementPacket<float>
VortexBlob::step_elements(float _ips) const {
  return std::vector<float>();
}

void
VortexBlob::debug(std::ostream& os) const {
  os << to_string();
}

std::string
VortexBlob::to_string() const {
  std::stringstream ss;
  ss << "vortex blob at " << m_x << " " << m_y << " " << m_z << ", radius " << m_rad << ", softness " << m_softness;
  ss << ", and strength " << m_sx << " " << m_sy << " " << m_sz;
  return ss.str();
}

void
VortexBlob::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
  const std::vector<float> s = j["strength"];
  m_sx = s[0];
  m_sy = s[1];
  m_sz = s[2];
  m_rad = j["rad"];
  m_softness = j["softness"];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
VortexBlob::to_json() const {
  nlohmann::json j;
  j["type"] = "vortex blob";
  j["center"] = {m_x, m_y, m_z};
  j["radius"] = m_rad;
  j["softness"] = m_softness;
  j["strength"] = {m_sx, m_sy, m_sz};
  j["enabled"] = m_enabled;
  return j;
}

void VortexBlob::generate_draw_geom() {
  // Based on irad in init_elems
  const float rad = 1+2*m_rad+m_softness;
  std::unique_ptr<Ovoid> tmp = std::make_unique<Ovoid>(nullptr, true, m_x, m_y, m_z,
                                                       2*rad, 2*rad, 2*rad);
  m_draw = tmp->init_elements(0.125);
  std::fill(m_draw.val.begin(), m_draw.val.end(), m_sz);
}

#ifdef USE_IMGUI
bool VortexBlob::draw_info_gui(const std::string _action, const float _ips) {
  static float xc[3] = {m_x, m_y, m_z};
  static float vstr[3] = {m_sx, m_sy, m_sz};
  static float guess_n = 4.1888f * std::pow((2.0f*m_rad+m_softness)/_ips, 3);
  std::string buttonText = _action+" vortex blob";
  bool add = false;

  ImGui::InputFloat3("center", xc);
  ImGui::InputFloat3("strength", vstr);
  ImGui::SliderFloat("radius", &m_rad, _ips, 10.0f*_ips, "%.4f");
  ImGui::SliderFloat("softness", &m_softness, _ips, m_rad, "%.4f");
  ImGui::Spacing();
  ImGui::TextWrapped("This feature will add about %f particles", guess_n);
  ImGui::Spacing();
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  m_z = xc[2];
  m_sx = vstr[0];
  m_sy = vstr[1];
  m_sz = vstr[2];
  return add;
}
#endif

//
// make the block of randomly-placed and random-strength particles
//
ElementPacket<float>
BlockOfRandom::init_elements(float _ips) const {

  // set up the random number generator
  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> zmean_dist(-0.5, 0.5);
  static std::uniform_real_distribution<> zo_dist(0.0, 1.0);

  // create a new vector to pass on
  std::vector<float> x;
  std::vector<Int> idx;
  std::vector<float> vals;
  x.resize(Dimensions*m_num);
  idx.resize(4*m_num);

  // initialize the particles' locations and strengths, leave radius zero for now
  for (size_t i=0; i<(size_t)m_num; ++i) {
    size_t idx = 3*i;
    // positions
    x[idx+0] = m_x + m_xsize*zmean_dist(gen);
    x[idx+1] = m_y + m_ysize*zmean_dist(gen);
    x[idx+2] = m_z + m_zsize*zmean_dist(gen);
  }

  for (size_t i=0; i<(size_t)m_num; ++i) {
    size_t idx = 4*i;
    // strengths
    vals[idx+0] = m_maxstr * zmean_dist(gen) / (float)m_num;
    vals[idx+1] = m_maxstr * zmean_dist(gen) / (float)m_num;
    vals[idx+2] = m_maxstr * zmean_dist(gen) / (float)m_num;
    // radius will get set later
    vals[idx+3] = 0.0f;
  }

  ElementPacket<float> packet({x, idx, vals, x.size()/3, 0});
  if (packet.verify(packet.x.size()+packet.val.size(), 1)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

ElementPacket<float>
BlockOfRandom::step_elements(float _ips) const {
  return ElementPacket<float>();
}

void
BlockOfRandom::debug(std::ostream& os) const {
  os << to_string();
}

std::string
BlockOfRandom::to_string() const {
  std::stringstream ss;
  ss << "block of " << m_num << " particles in [" << (m_x-0.5*m_xsize) << " " << (m_x+0.5*m_xsize) << "] ["
                                                  << (m_y-0.5*m_ysize) << " " << (m_y+0.5*m_ysize) << "] ["
                                                  << (m_z-0.5*m_zsize) << " " << (m_z+0.5*m_zsize) <<
                                               "] with max str mag " << m_maxstr;
  return ss.str();
}

void
BlockOfRandom::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
  const std::vector<float> s = j["size"];
  m_xsize = s[0];
  m_ysize = s[1];
  m_zsize = s[2];
  m_maxstr = j["max strength"];
  m_num = j["num"];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
BlockOfRandom::to_json() const {
  nlohmann::json j;
  j["type"] = "block of random";
  j["center"] = {m_x, m_y, m_z};
  j["size"] = {m_xsize, m_ysize, m_zsize};
  j["max strength"] = m_maxstr;
  j["num"] = m_num;
  j["enabled"] = m_enabled;
  return j;
}

void BlockOfRandom::generate_draw_geom() {
  std::unique_ptr<SolidRect> tmp = std::make_unique<SolidRect>(nullptr, true, m_x, m_y, m_z, m_xsize, m_ysize, m_zsize);
  m_draw = tmp->init_elements(1);
  std::fill(m_draw.val.begin(), m_draw.val.end(), m_maxstr);
}

#ifdef USE_IMGUI
bool BlockOfRandom::draw_info_gui(const std::string action, const float ips) {
  static float xs[3] = {m_xsize, m_ysize, m_zsize};
  static float xc[3] = {m_x, m_y, m_z};
  std::string buttonText = action+" random vorticies";
  bool add = false;

  ImGui::SliderInt("number", &m_num, 10, 100000);
  ImGui::SliderFloat3("box size", xs, 0.01f, 10.0f, "%.4f", 2.0f);
  ImGui::InputFloat3("center", xc);
  ImGui::SliderFloat("strength magnitude", &m_maxstr, 0.01f, 10.0f, "%.3f", 2.0f);
  ImGui::Spacing();
  ImGui::TextWrapped("This feature will add %d particles", m_num);
  ImGui::Spacing();
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_xsize = xs[0];
  m_ysize = xs[1];
  m_zsize = xs[2];
  m_x = xc[0];
  m_y = xc[1];
  m_z = xc[2];
  return add;
}
#endif

//
// drop a single particle from the emitter
//
ElementPacket<float>
ParticleEmitter::init_elements(float _ips) const {
  return ElementPacket<float>();
}

ElementPacket<float>
ParticleEmitter::step_elements(float _ips) const {
  std::vector<float> x = {m_x, m_y, m_z};
  std::vector<Int> idx;
  std::vector<float> vals = {m_sx, m_sy, m_sz, 0.0};
  ElementPacket<float> packet(x, idx, vals, x.size()/3, 0);

  if (packet.verify(packet.x.size()+packet.val.size(), 7)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

void
ParticleEmitter::debug(std::ostream& os) const {
  os << to_string();
}

std::string
ParticleEmitter::to_string() const {
  std::stringstream ss;
  ss << "particle emitter at " << m_x << " " << m_y << " " << m_z << " spawning particles";
  ss << " with strength " << m_sx << " " << m_sy << " " << m_sz;
  return ss.str();
}

void
ParticleEmitter::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
  const std::vector<float> s = j["strength"];
  m_sx = s[0];
  m_sy = s[1];
  m_sz = s[2];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
ParticleEmitter::to_json() const {
  nlohmann::json j;
  j["type"] = "particle emitter";
  j["center"] = {m_x, m_y, m_z};
  j["strength"] = {m_sx, m_sy, m_sz};
  j["enabled"] = m_enabled;
  return j;
}

void ParticleEmitter::generate_draw_geom() {
  const float diam = 0.01;
  std::unique_ptr<Ovoid> tmp = std::make_unique<Ovoid>(nullptr, true, m_x, m_y, m_z, m_sx, m_sy, m_sz);
  m_draw = tmp->init_elements(diam/25.0);
  //std::fill(m_draw.val.begin(), m_draw.val.end(), m_str);
}

#ifdef USE_IMGUI
bool ParticleEmitter::draw_info_gui(const std::string action, const float ips) {
  return false;
}
#endif

//
// make a singular (one row) vortex ring
//
ElementPacket<float>
SingularRing::init_elements(float _ips) const {

// create a new vector to pass on
  std::vector<float> x;
  std::vector<Int> idx;
  std::vector<float> vals;

  // what size 2D integer array will we loop over
  const int ndiam = 1 + (2.0 * M_PI * m_majrad) / _ips;
  std::cout << "  ring needs " << ndiam << " particles" << std::endl;
  const float this_ips = (2.0 * M_PI * m_majrad) / (float)ndiam;

  // generate a set of orthogonal basis vectors for the given normal
  std::array<float,3> norm = {m_nx, m_ny, m_nz};
  normalizeVec(norm);
  std::array<float,3> b1, b2;
  branchlessONB<float>(norm, b1, b2);

  // loop over integer indices
  for (int i=0; i<ndiam; ++i) {
    const float theta = 2.0 * M_PI * (float)i / (float)ndiam;

    // create a particle here
    x.emplace_back(m_x + m_majrad * (b1[0]*std::cos(theta) + b2[0]*std::sin(theta)));
    x.emplace_back(m_y + m_majrad * (b1[1]*std::cos(theta) + b2[1]*std::sin(theta)));
    x.emplace_back(m_z + m_majrad * (b1[2]*std::cos(theta) + b2[2]*std::sin(theta)));

    // set the strength
    vals.emplace_back(this_ips * m_circ * (b2[0]*std::cos(theta) - b1[0]*std::sin(theta)));
    vals.emplace_back(this_ips * m_circ * (b2[1]*std::cos(theta) - b1[1]*std::sin(theta)));
    vals.emplace_back(this_ips * m_circ * (b2[2]*std::cos(theta) - b1[2]*std::sin(theta)));
    // this is the radius - still zero for now
    vals.emplace_back(0.0f);
  }

  ElementPacket<float> packet({x, idx, vals, (size_t)ndiam, 0});
  if (packet.verify(packet.x.size()+packet.val.size(), 7)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

ElementPacket<float>
SingularRing::step_elements(float _ips) const {
  return ElementPacket<float>();
}

void
SingularRing::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SingularRing::to_string() const {
  std::stringstream ss;
  ss << "singular vortex ring at " << m_x << " " << m_y << " " << m_z << ", radius " << m_majrad << ", circulation " << m_circ;
  ss << ", aimed along " << m_nx << " " << m_ny << " " << m_nz;
  return ss.str();
}

void
SingularRing::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
  const std::vector<float> n = j["normal"];
  m_nx = n[0];
  m_ny = n[1];
  m_nz = n[2];
  m_majrad = j["major radius"];
  m_circ = j["circulation"];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
SingularRing::to_json() const {
  nlohmann::json j;
  j["type"] = "singular ring";
  j["center"] = {m_x, m_y, m_z};
  j["normal"] = {m_nx, m_ny, m_nz};
  j["major radius"] = m_majrad;
  j["circulation"] = m_circ;
  j["enabled"] = m_enabled;
  return j;
}

//void create
void SingularRing::generate_draw_geom() {
  // For sake of visualization, imagine the torus sitting on your table
  // like a doughnut
  // Number of steps around the torus
  const int ts = 12;
  // Number of steps around the circle perpendicular to the table
  const int cs = 12;
  const float m_minrad = 0.01;
  /*std::vector<float> x;
  std::vector<Int> idx;
  for (int i=0; i<ts; i++) {
    const float alpha = (2*M_PI)*i/ts;
    for (int j=0; j<cs; j++) {
      const float beta = (2*M_PI)*j/cs;
      x.emplace_back((m_majrad+m_minrad*std::cos(beta))*std::cos(alpha));
      x.emplace_back((m_majrad+m_minrad*std::cos(beta))*std::sin(alpha));
      x.emplace_back(m_minrad*sin(beta));
    }
  }*/
  // We first make a circle we will then drag in a circular motion to create the torus
  std::vector<float> circle;
  for (int i=0; i<cs; i++) {
    const float alpha = 2*M_PI*i/cs;
    circle.emplace_back(m_minrad*std::sin(alpha));
    circle.emplace_back(m_minrad*std::cos(alpha));
  }
  // generate a set of orthogonal basis vectors for the given normal
  std::array<float,3> norm = {m_nx, m_ny, m_nz};
  normalizeVec(norm);
  std::array<float,3> b1, b2;
  branchlessONB<float>(norm, b1, b2);
  
  std::vector<float> x;
  // loop over integer indices
  for (int i=0; i<ts; ++i) {
    const float theta = 2.0 * M_PI * (float)i / (float)ts;
    const float st = std::sin(theta);
    const float ct = std::cos(theta);

    for (int j=0; j<cs; ++j) {
      // helper indices for the disk particles
      const size_t ix = 3*j;
      const size_t iy = 3*j+1;

      // create a particle here
      x.emplace_back(m_x + (m_majrad + circle[ix]) * (b1[0]*ct + b2[0]*st) + circle[iy]*norm[0]);
      x.emplace_back(m_y + (m_majrad + circle[ix]) * (b1[1]*ct + b2[1]*st) + circle[iy]*norm[1]);
      x.emplace_back(m_z + (m_majrad + circle[ix]) * (b1[2]*ct + b2[2]*st) + circle[iy]*norm[2]);
    }
  }
      
  // Create triangles from points
  std::vector<Int> idx;
  for (int i=0; i<ts; i++) {
    const Int ir1 = cs*i;
    const Int ir2 = cs*(i+1);
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
    idx.emplace_back(ir2);
    // set 2
    idx.emplace_back(ir2-1);
    idx.emplace_back(ir2);
    idx.emplace_back(ir1);
  }
  const int numPoints = cs*ts;
  for (size_t i=idx.size()-cs; i<idx.size(); i++) {
    idx[i] = idx[i]%numPoints;
  }

  std::vector<float> val;
  val.resize(idx.size());
  std::fill(val.begin(), val.end(), 1.0);
  ElementPacket epack {x, idx, val, val.size()/Dimensions, 2};
  m_draw = epack;
}

#ifdef USE_IMGUI
bool SingularRing::draw_info_gui(const std::string _action, const float _ips) {
  float xc[3] = {m_x, m_y, m_z};
  float vstr[3] = {m_nx, m_ny, m_nz};
  float guess_n = 1 + (2.0f * 3.1416f * m_majrad / _ips);
  std::string buttonText = _action+" singular ring";
  bool add = false;

  ImGui::InputFloat3("center", xc);
  ImGui::InputFloat3("direction", vstr);
  ImGui::SliderFloat("circulation", &m_circ, 0.001f, 10.0f, "%.3f");
  ImGui::SliderFloat("radius", &m_majrad, 3.0f*_ips, 10.0f, "%.3f");
  ImGui::Spacing();
  ImGui::TextWrapped("This feature will add about %f particles", guess_n);
  ImGui::Spacing();
  if (ImGui::Button("Add singular vortex ring")) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  m_z = xc[2];
  m_nx = vstr[0];
  m_ny = vstr[1];
  m_nz = vstr[2];
  return add;
}
#endif

//
// make a thick-cored vortex ring
//
ElementPacket<float>
ThickRing::init_elements(float _ips) const {

  if (not this->is_enabled()) return std::vector<float>();

  // make a temporary array of the particles at one station around the ring (a disk)
  //   for each particle in the disk, this is the local x,y, and a length scale
  //   where +x is pointing directly away from the center of the vortex ring,
  //   and +y is along the ring's normal vector
  std::vector<float> disk;
  const int nlayers = 1 + m_minrad / _ips;
  // the central particle
  disk.emplace_back(0.0);
  disk.emplace_back(0.0);
  disk.emplace_back(1.0);
  int nthisdisk = 1;
  // and each ring of particles
  for (int l=1; l<nlayers; ++l) {
    const float thisrad = (float)l * _ips;
    const int nthislayer = 1 + (2.0 * M_PI * thisrad) / _ips;
    for (int i=0; i<nthislayer; ++i) {
      const float phi = 2.0 * M_PI * (float)i / (float)nthislayer;
      disk.emplace_back(thisrad * std::cos(phi));
      disk.emplace_back(thisrad * std::sin(phi));
      disk.emplace_back((m_majrad + thisrad*std::cos(phi)) / m_majrad);
    }
    nthisdisk += nthislayer;
  }
  std::cout << "  ring needs " << nlayers << " layers and " << nthisdisk << " particles per azimuthal station" << std::endl;

  // And how many stations around the ring?
  const int ndiam = 1 + (2.0 * M_PI * m_majrad) / _ips;
  //std::cout << "  ring needs " << ndiam << " azimuthal stations" << std::endl;
  const float this_ips = (2.0 * M_PI * m_majrad) / (float)ndiam;

  // generate a set of orthogonal basis vectors for the given normal
  std::array<float,3> norm = {m_nx, m_ny, m_nz};
  normalizeVec(norm);
  std::array<float,3> b1, b2;
  branchlessONB<float>(norm, b1, b2);

  std::vector<float> x;
  std::vector<Int> idx;
  std::vector<float> vals;

  // loop over integer indices
  for (int i=0; i<ndiam; ++i) {
    const float theta = 2.0 * M_PI * (float)i / (float)ndiam;
    const float st = std::sin(theta);
    const float ct = std::cos(theta);

    for (int j=0; j<nthisdisk; ++j) {
      // helper indices for the disk particles
      const size_t ix = 3*j;
      const size_t iy = 3*j+1;
      const size_t il = 3*j+2;

      // create a particle here
      x.emplace_back(m_x + (m_majrad + disk[ix]) * (b1[0]*ct + b2[0]*st) + disk[iy]*norm[0]);
      x.emplace_back(m_y + (m_majrad + disk[ix]) * (b1[1]*ct + b2[1]*st) + disk[iy]*norm[1]);
      x.emplace_back(m_z + (m_majrad + disk[ix]) * (b1[2]*ct + b2[2]*st) + disk[iy]*norm[2]);

      // set the strength
      const float sscale = disk[il] * this_ips * m_circ / (float)nthisdisk;
      vals.emplace_back(sscale * (b2[0]*ct - b1[0]*st));
      vals.emplace_back(sscale * (b2[1]*ct - b1[1]*st));
      vals.emplace_back(sscale * (b2[2]*ct - b1[2]*st));

      // this is the radius - still zero for now
      vals.emplace_back(0.0f);
    }
  }

  ElementPacket<float> packet({x, idx, vals, (size_t)x.size()/3, 0});
  if (packet.verify(packet.x.size()+packet.val.size(), 1)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

ElementPacket<float>
ThickRing::step_elements(float _ips) const {
  return std::vector<float>();
}

void
ThickRing::debug(std::ostream& os) const {
  os << to_string();
}

std::string
ThickRing::to_string() const {
  std::stringstream ss;
  ss << "thick vortex ring at " << m_x << " " << m_y << " " << m_z << ", radii " << m_majrad << " " << m_minrad << ", circulation " << m_circ;
  ss << ", aimed along " << m_nx << " " << m_ny << " " << m_nz;
  return ss.str();
}

void
ThickRing::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
  const std::vector<float> n = j["normal"];
  m_nx = n[0];
  m_ny = n[1];
  m_nz = n[2];
  m_majrad = j["major radius"];
  m_minrad = j["minor radius"];
  m_circ = j["circulation"];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
ThickRing::to_json() const {
  nlohmann::json j;
  j["type"] = "thick ring";
  j["center"] = {m_x, m_y, m_z};
  j["normal"] = {m_nx, m_ny, m_nz};
  j["major radius"] = m_majrad;
  j["minor radius"] = m_minrad;
  j["circulation"] = m_circ;
  j["enabled"] = m_enabled;
  return j;
}

void ThickRing::generate_draw_geom() {

}

#ifdef USE_IMGUI
bool ThickRing::draw_info_gui(const std::string action, const float ips) {
  float xc[3] = {m_x, m_y, m_z};
  float vstr[3] = {m_nx, m_ny, m_nz};
  float guess_n = (1 + (2.0f * 3.1416f * m_majrad / ips) * std::pow(m_minrad/ips, 2));
  std::string buttonText = action+" thick vortex ring";
  bool add = false;

  ImGui::InputFloat3("center", xc);
  ImGui::InputFloat3("direction", vstr);
  ImGui::SliderFloat("circulation", &m_circ, 0.001f, 10.0f, "%.4f");
  ImGui::SliderFloat("radius", &m_majrad, 3.0f*ips, 10.0f, "%.3f");
  ImGui::SliderFloat("thickness", &m_minrad, ips, 10.0f*ips, "%.4f");
  ImGui::Spacing();
  ImGui::TextWrapped("This feature will add about %f particles", guess_n);
  ImGui::Spacing();
  if (ImGui::Button("Add thick vortex ring")) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  m_z = xc[2];
  m_nx = vstr[0];
  m_ny = vstr[1];
  m_nz = vstr[2];
  return add;
}
#endif
