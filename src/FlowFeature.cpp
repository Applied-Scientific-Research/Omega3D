/*
 * FlowFeature.cpp - GUI-side descriptions of flow features
 *
 * (c)2017-21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
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

#include "FlowFeature.h"
#include "MathHelper.h"
#include "GeomHelper.h"

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
  else if (ftype == "uniform block") {    _flist.emplace_back(std::make_unique<UniformBlock>()); }
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

  std::cout << "  finished " << _flist.back()->to_string() << std::endl;
}

#ifdef USE_IMGUI
// 0 means keep open, 1 means create, 2 means cancel
int FlowFeature::draw_creation_gui(std::vector<std::unique_ptr<FlowFeature>> &ffs, const float ips) {
  static int item = 1;
  static int oldItem = -1;
  const char* items[] = { "vortex particle", "vortex blob", "uniform block", "random particles",  
                          "singular vortex ring", "thick vortex ring", "particle emitter" };
  ImGui::Combo("type", &item, items, 7);

  // show different inputs based on what is selected
  static std::unique_ptr<FlowFeature> ff = nullptr;
  if (oldItem != item) {
    switch(item) {
      case 0: {
        ff = std::make_unique<SingleParticle>();
      } break;
      case 1: {
        ff = std::make_unique<VortexBlob>();
      } break;
      case 2: {
        ff = std::make_unique<UniformBlock>();
      } break;
      case 3: {
        ff = std::make_unique<BlockOfRandom>();
      } break;
      case 4: {
        ff = std::make_unique<SingularRing>();
      } break;
      case 5: {
        ff = std::make_unique<ThickRing>();
      } break;
      case 6: {
        ff = std::make_unique<ParticleEmitter>();
      } break;
    }
    oldItem = item;
  }

  int created = 0;
  if (ff->draw_info_gui("Add", ips)) {
    ff->generate_draw_geom();
    ffs.emplace_back(std::move(ff));
    ff = nullptr;
    created = 1;
    oldItem = -1;
  }

  ImGui::SameLine();
  if (ImGui::Button("Cancel", ImVec2(120,0))) {
    oldItem = -1;
    created = 2;
    ff = nullptr;
  }

  return created;
}

void FlowFeature::draw_feature_list(std::vector<std::unique_ptr<FlowFeature>> &feat,
                                    std::unique_ptr<FlowFeature> &editingFeat, int &edit_feat_index,
                                    int &del_feat_index, bool &redraw, int &buttonIDs) {
  for (int i=0; i<(int)feat.size(); ++i) {
    ImGui::PushID(++buttonIDs);
    if (ImGui::Checkbox("", feat[i]->addr_enabled())) { redraw = true; }
    ImGui::PopID();
    
    // add an "edit" button after the checkbox (so it's not easy to accidentally hit remove)
    ImGui::SameLine();
    ImGui::PushID(++buttonIDs);
    if (ImGui::SmallButton("edit")) {
      editingFeat = std::unique_ptr<FlowFeature>(feat[i]->copy());
      edit_feat_index = i;
    }
    ImGui::PopID();
    
    if (feat[i]->is_enabled()) {
      ImGui::SameLine();
      ImGui::Text("%s", feat[i]->to_string().c_str());
    } else {
      ImGui::SameLine();
      ImGui::TextColored(ImVec4(0.5f,0.5f,0.5f,1.0f), "%s", feat[i]->to_string().c_str());
    }

    // add a "remove" button at the end of the line (so it's not easy to accidentally hit)
    ImGui::SameLine();
    ImGui::PushID(++buttonIDs);
    if (ImGui::SmallButton("remove")) { del_feat_index = i; }
    ImGui::PopID();
  }
}
#endif

//
// important feature: convert flow feature definition into an ElementPacket
//

//
// drop a single particle
//
ElementPacket<float>
SingleParticle::init_elements(float _ips) const {
  std::cout << "Creating single particle" << std::endl;

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

// Single Particles cant be made by user
void SingleParticle::generate_draw_geom() {
  const float diam = 0.01;
  m_draw = generate_ovoid(diam, diam, diam, 0.05*diam);
  m_draw.translate(m_x, m_y, m_z);

  // OpenGL expects a val for every point (3x's)
  const int numPts = m_draw.x.size()/Dimensions;
  m_draw.val.resize(numPts);
  std::fill(m_draw.val.begin(), m_draw.val.end(), length(std::array<float,3>{m_sx, m_sy, m_sz}));
}

#ifdef USE_IMGUI
// User can't actually create this
bool SingleParticle::draw_info_gui(const std::string _action, const float _ips) {
  float xc[3] = {m_x, m_y, m_z};
  float xs[3] = {m_sx, m_sy, m_sz};
  std::string buttonText = _action+" single particle";
  bool add = false;

  ImGui::InputFloat3("center", xc);
  ImGui::InputFloat3("strengths", xs);
  ImGui::Spacing();
  ImGui::TextWrapped("This feature will add a single particle");
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

  std::cout << "Creating vortex blob with up to " << std::pow(2*irad+1,3) << " particles" << std::endl;

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
      //vals.emplace_back(0.0f);
    }
  }
  }
  }

  // finally, normalize all particle strengths so that the whole blob
  //   has exactly the right strength
  std::cout << "blob had " << tot_wgt << " initial circulation" << std::endl;
  double str_scale = 1.0 / tot_wgt;
  for (size_t i=0; i<vals.size(); i+=3) {
    vals[i+0] = (float)((double)vals[i+0] * str_scale);
    vals[i+1] = (float)((double)vals[i+1] * str_scale);
    vals[i+2] = (float)((double)vals[i+2] * str_scale);
  }

  ElementPacket<float> packet({x, idx, vals, x.size()/Dimensions, 0});
  if (packet.verify(packet.x.size()+packet.val.size(), 6)) {
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
  m_rad = j["radius"];
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
  const float diam = 2.0*m_rad + m_softness;
  m_draw = generate_ovoid(diam, diam, diam, 0.05*diam);
  m_draw.translate(m_x, m_y, m_z);

  // OpenGL expects a val for every point (3x's)
  const int numPts = m_draw.x.size()/Dimensions;
  m_draw.val.resize(numPts);
  std::fill(m_draw.val.begin(), m_draw.val.end(), length(std::array<float,3>{m_sx, m_sy, m_sz}));

  //for (size_t i=0; i<m_draw.val.size()/Dimensions; i+=Dimensions) {
  //  m_draw.val[i] = m_sx;
  //  m_draw.val[i+1] = m_sy;
  //  m_draw.val[i+2] = m_sz;
  //}
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
  ImGui::TextWrapped("This feature will add about %d particles", (int)guess_n);
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
// make the block of regular, and uniform-strength particles
//
ElementPacket<float>
UniformBlock::init_elements(float _ips) const {

  // what size 2D integer array will we loop over
  const int isize = 1 + m_xsize / _ips;
  const int jsize = 1 + m_ysize / _ips;
  const int ksize = 1 + m_zsize / _ips;
  const size_t totn = isize*jsize*ksize;
  std::cout << "Creating block with " << totn << " particles" << std::endl;
  //std::cout << "block needs " << isize << " by " << jsize << " particles" << std::endl;

  // create a new vector to pass on
  std::vector<float> x(Dimensions*totn);
  std::vector<Int> idx;
  std::vector<float> vals(numStrPerNode*totn);

  const float each_str = 1.0 / (float)(totn);

  // initialize the particles' locations and strengths, leave radius zero for now
  size_t ix = 0;
  size_t iv = 0;
  for (int i=0; i<isize; ++i) {
  for (int j=0; j<jsize; ++j) {
  for (int k=0; k<ksize; ++k) {
    x[ix++] = m_x + m_xsize * (((float)i + 0.5)/(float)isize - 0.5);
    x[ix++] = m_y + m_ysize * (((float)j + 0.5)/(float)jsize - 0.5);
    x[ix++] = m_z + m_zsize * (((float)k + 0.5)/(float)ksize - 0.5);
    vals[iv++] = m_sx * each_str;
    vals[iv++] = m_sy * each_str;
    vals[iv++] = m_sz * each_str;
  }
  }
  }

  ElementPacket<float> packet({x, idx, vals, totn, 0});
  if (packet.verify(packet.x.size()+packet.val.size(), 6)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

ElementPacket<float>
UniformBlock::step_elements(float _ips) const {
  return ElementPacket<float>();
}

void
UniformBlock::debug(std::ostream& os) const {
  os << to_string();
}

std::string
UniformBlock::to_string() const {
  std::stringstream ss;
  ss << "block of uniform particles in [" << (m_x-0.5*m_xsize) << " " << (m_x+0.5*m_xsize) << "] ["
                                          << (m_y-0.5*m_ysize) << " " << (m_y+0.5*m_ysize) << "] ["
                                          << (m_z-0.5*m_zsize) << " " << (m_z+0.5*m_zsize) <<
                                        "] with str " << m_sx << " " << m_sy << " " << m_sz;
  return ss.str();
}

void
UniformBlock::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
  const std::vector<float> z = j["size"];
  m_xsize = z[0];
  m_ysize = z[1];
  m_zsize = z[2];
  const std::vector<float> s = j["strength"];
  m_sx = s[0];
  m_sy = s[1];
  m_sz = s[2];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
UniformBlock::to_json() const {
  nlohmann::json j;
  j["type"] = "uniform block";
  j["center"] = {m_x, m_y, m_z};
  j["size"] = {m_xsize, m_ysize, m_zsize};
  j["strength"] = {m_sx, m_sy, m_sz};
  j["enabled"] = m_enabled;
  return j;
}

void UniformBlock::generate_draw_geom() {

  // generate draw geometry, OK to use more than 12 triangles
  m_draw = generate_cuboid(m_xsize, m_ysize, m_zsize, 0.1*(m_xsize+m_ysize+m_zsize));
  m_draw.translate(m_x-0.5*m_xsize, m_y-0.5*m_ysize, m_z-0.5*m_zsize);

  // OpenGL expects a val for every point (3x's)
  const size_t numPts = m_draw.x.size()/Dimensions;
  m_draw.val.resize(numPts);
  const float maxval = std::sqrt(m_sx*m_sx + m_sy*m_sy + m_sz*m_sz);
  std::fill(m_draw.val.begin(), m_draw.val.end(), maxval);
  //for (size_t i = 0; i < numPts; i++) {
  //  m_draw.val[i] = maxval;
  //}
}

#ifdef USE_IMGUI
bool UniformBlock::draw_info_gui(const std::string action, const float ips) {
  bool add = false;
  static float xc[3] = {m_x, m_y, m_z};
  static float xs[3] = {m_xsize, m_ysize, m_zsize};
  static float vstr[3] = {m_sx, m_sy, m_sz};
  const std::string buttonText = action+" block of vorticies";

  ImGui::InputFloat3("center", xc);
  ImGui::InputFloat3("strength", vstr);
  ImGui::SliderFloat3("box size", xs, 0.01f, 10.0f, "%.4f", 2.0f);
  ImGui::TextWrapped("This feature will add about %d particles", (int)(xs[0]*xs[1]*xs[2]/std::pow(ips,3)));
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  m_z = xc[2];
  m_xsize = xs[0];
  m_ysize = xs[1];
  m_zsize = xs[2];
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

  std::cout << "Creating random block with " << m_num << " particles" << std::endl;

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
  vals.resize(numStrPerNode*m_num);

  // initialize the particles' locations and strengths, leave radius zero for now
  for (size_t i=0; i<(size_t)m_num; ++i) {
    size_t idx = 3*i;
    // positions
    x[idx+0] = m_x + m_xsize*zmean_dist(gen);
    x[idx+1] = m_y + m_ysize*zmean_dist(gen);
    x[idx+2] = m_z + m_zsize*zmean_dist(gen);
  }

  for (size_t i=0; i<(size_t)m_num; ++i) {
    size_t idx = 3*i;
    // strengths
    vals[idx+0] = m_maxstr * zmean_dist(gen) / (float)m_num;
    vals[idx+1] = m_maxstr * zmean_dist(gen) / (float)m_num;
    vals[idx+2] = m_maxstr * zmean_dist(gen) / (float)m_num;
    // radius will get set later
    //vals[idx+3] = 0.0f;
  }

  ElementPacket<float> packet({x, idx, vals, x.size()/Dimensions, 0});
  if (packet.verify(packet.x.size()+packet.val.size(), 6)) {
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

  // generate draw geometry, OK to use more than 12 triangles
  m_draw = generate_cuboid(m_xsize, m_ysize, m_zsize, 0.04*(m_xsize+m_ysize+m_zsize));
  m_draw.translate(m_x-0.5*m_xsize, m_y-0.5*m_ysize, m_z-0.5*m_zsize);

  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> zmean_dist(-0.5, 0.5);

  // OpenGL expects a val for every point (3x's)
  const size_t numPts = m_draw.x.size()/Dimensions;
  m_draw.val.resize(numPts);
  for (size_t i = 0; i < numPts; i++) {
    m_draw.val[i] = m_maxstr*zmean_dist(gen);
  }
}

#ifdef USE_IMGUI
bool BlockOfRandom::draw_info_gui(const std::string _action, const float _ips) {
  static float xs[3] = {m_xsize, m_ysize, m_zsize};
  static float xc[3] = {m_x, m_y, m_z};
  std::string buttonText = _action+" random vorticies";
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
  std::cout << "Creating particle emitter" << std::endl;
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
  m_draw = generate_ovoid(diam, diam, diam, 0.05*diam);
  m_draw.translate(m_x, m_y, m_z);

  // OpenGL expects a val for every point (3x's)
  const int numPts = m_draw.val.size()/Dimensions;
  m_draw.val.resize(numPts);
  const float sign = std::copysign(1.0, m_sx+m_sy+m_sz);
  std::fill(m_draw.val.begin(), m_draw.val.end(), sign*length(std::array<float,3>{m_sx, m_sy, m_sz}));
}

#ifdef USE_IMGUI
bool ParticleEmitter::draw_info_gui(const std::string _action, const float _ips) {
  float xc[3] = {m_x, m_y, m_z};
  float xs[3] = {m_sx, m_sy, m_sz};
  std::string buttonText = _action+" particle emitter";
  bool add = false;

  ImGui::InputFloat3("center", xc);
  ImGui::InputFloat3("strengths", xs);
  ImGui::Spacing();
  ImGui::TextWrapped("This feature will add a particle emitter");
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
  const float this_ips = (2.0 * M_PI * m_majrad) / (float)ndiam;

  std::cout << "Creating singular ring with " << ndiam << " particles" << std::endl;

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
    //vals.emplace_back(0.0f);
  }

  ElementPacket<float> packet({x, idx, vals, (size_t)ndiam, 0});
  if (packet.verify(packet.x.size()+packet.val.size(), 6)) {
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

void SingularRing::generate_draw_geom() {
  // For sake of visualization, imagine the torus sitting on your table like a doughnut

  // generate the ElementPacket here - scaled and along +z
  m_draw = generate_torus(m_majrad, 0.02*m_majrad, 0.013*m_majrad);

  // generate a set of orthogonal basis vectors for the given normal
  std::array<float,3> norm = {m_nx, m_ny, m_nz};
  normalizeVec(norm);
  std::array<float,3> b1, b2;
  branchlessONB<float>(norm, b1, b2);
  
  // rotate the points to the new basis
  for (size_t i=0; i<m_draw.x.size()/Dimensions; ++i) {
    const float px = m_draw.x[i*Dimensions+0];
    const float py = m_draw.x[i*Dimensions+1];
    const float pz = m_draw.x[i*Dimensions+2];

    m_draw.x[i*Dimensions+0] = b1[0]*px + b2[0]*py + norm[0]*pz;
    m_draw.x[i*Dimensions+1] = b1[1]*px + b2[1]*py + norm[1]*pz;
    m_draw.x[i*Dimensions+2] = b1[2]*px + b2[2]*py + norm[2]*pz;
  }

  // finally, translate to the desired center
  m_draw.translate(m_x, m_y, m_z);
      
  // OpenGL expects a val for every point (3x's)
  const int numPts = m_draw.x.size()/Dimensions;
  m_draw.val.resize(numPts);
  std::fill(m_draw.val.begin(), m_draw.val.end(), length(std::array<float,3>{m_nx, m_ny, m_nz}));
}

#ifdef USE_IMGUI
bool SingularRing::draw_info_gui(const std::string _action, const float _ips) {

  float xc[3] = {m_x, m_y, m_z};
  float vstr[3] = {m_nx, m_ny, m_nz};
  float guess_n = 1 + (2.0f * M_PI * m_majrad / _ips);
  std::string buttonText = _action+" singular ring";
  bool add = false;

  ImGui::InputFloat3("center", xc);
  ImGui::InputFloat3("direction", vstr);
  ImGui::SliderFloat("circulation", &m_circ, 0.001f, 10.0f, "%.3f");
  ImGui::SliderFloat("radius", &m_majrad, 3.0f*_ips, 10.0f, "%.3f");
  ImGui::Spacing();
  ImGui::TextWrapped("This feature will add about %d particles", (int)guess_n);
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

  std::cout << "Creating thick vortex ring" << std::endl;

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
      //vals.emplace_back(0.0f);
    }
  }

  ElementPacket<float> packet({x, idx, vals, (size_t)x.size()/Dimensions, 0});
  if (packet.verify(packet.x.size()+packet.val.size(), 6)) {
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
  // For sake of visualization, imagine the torus sitting on your table like a doughnut

  // generate the ElementPacket here - scaled and along +z
  m_draw = generate_torus(m_majrad, m_minrad, 0.013*m_majrad);

  // generate a set of orthogonal basis vectors for the given normal
  std::array<float,3> norm = {m_nx, m_ny, m_nz};
  normalizeVec(norm);
  std::array<float,3> b1, b2;
  branchlessONB<float>(norm, b1, b2);

  // rotate the points to the new basis
  for (size_t i=0; i<m_draw.x.size()/Dimensions; ++i) {
    const float px = m_draw.x[i*Dimensions+0];
    const float py = m_draw.x[i*Dimensions+1];
    const float pz = m_draw.x[i*Dimensions+2];

    m_draw.x[i*Dimensions+0] = b1[0]*px + b2[0]*py + norm[0]*pz;
    m_draw.x[i*Dimensions+1] = b1[1]*px + b2[1]*py + norm[1]*pz;
    m_draw.x[i*Dimensions+2] = b1[2]*px + b2[2]*py + norm[2]*pz;
  }

  // finally, translate to the desired center
  m_draw.translate(m_x, m_y, m_z);
      
  // OpenGL expects a val for every point (3x's)
  const int numPts = m_draw.x.size()/Dimensions;
  m_draw.val.resize(numPts);
  std::fill(m_draw.val.begin(), m_draw.val.end(), length(std::array<float,3>{m_nx, m_ny, m_nz}));
}

#ifdef USE_IMGUI
bool ThickRing::draw_info_gui(const std::string _action, const float _ips) {
  float xc[3] = {m_x, m_y, m_z};
  float vstr[3] = {m_nx, m_ny, m_nz};
  float guess_n = (1 + (2.0f * M_PI * m_majrad / _ips) * std::pow(m_minrad/_ips, 2));
  std::string buttonText = _action+" thick vortex ring";
  bool add = false;

  ImGui::InputFloat3("center", xc);
  ImGui::InputFloat3("direction", vstr);
  ImGui::SliderFloat("circulation", &m_circ, 0.001f, 10.0f, "%.4f");
  ImGui::SliderFloat("radius", &m_majrad, 3.0f*_ips, 10.0f, "%.3f");
  ImGui::SliderFloat("thickness", &m_minrad, _ips, 10.0f*_ips, "%.4f");
  ImGui::Spacing();
  ImGui::TextWrapped("This feature will add about %d particles", (int)guess_n);
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
