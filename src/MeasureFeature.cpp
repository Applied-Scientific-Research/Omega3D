/*
 * MeasureFeature.cpp - GUI-side descriptions of flow features
 *
 * (c)2018-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#include "BoundaryFeature.h"
#include "MeasureFeature.h"
#include "imgui/imgui.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <random>

// write out any object of parent type MeasureFeature by dispatching to appropriate "debug" method
std::ostream& operator<<(std::ostream& os, MeasureFeature const& ff) {
  ff.debug(os);
  return os;
}


//
// parse the json and dispatch the constructors
//
void parse_measure_json(std::vector<std::unique_ptr<MeasureFeature>>& _flist,
                        const nlohmann::json _jin) {

  // must have one and only one type
  if (_jin.count("type") != 1) return;

  const std::string ftype = _jin["type"];

  if      ((ftype == "tracer") || (ftype == "point")) {
    _flist.emplace_back(std::make_unique<SinglePoint>());
  } else if (ftype == "tracer emitter") {   _flist.emplace_back(std::make_unique<SinglePoint>(0.0, 0.0, false, true)); }
  else if (ftype == "tracer blob") {      _flist.emplace_back(std::make_unique<MeasurementBlob>()); }
  else if (ftype == "tracer line") {      _flist.emplace_back(std::make_unique<MeasurementLine>(0.0, 0.0, false, true)); }
  else if (ftype == "measurement line") { _flist.emplace_back(std::make_unique<MeasurementLine>()); }
  else if (ftype == "measurement grid") { _flist.emplace_back(std::make_unique<Grid2dPoints>()); }
  else {
    std::cout << "  type " << ftype << " does not name an available measurement feature, ignoring" << std::endl;
    return;
  }

  // and pass the json object to the specific parser
  _flist.back()->from_json(_jin);

  std::cout << "  found " << ftype << std::endl;
}

#ifdef USE_IMGUI
bool MeasureFeature::draw_creation_gui(std::vector<std::unique_ptr<MeasureFeature>> &mfs, const float _ips, const float &_tracerScale) {
  static int item = 0;
  static int oldItem = -1;
  const char* items[] = { "single point", "measurement circle", "measurement line", "2d grid" };
  ImGui::Combo("type", &item, items, 4);

  // show different inputs based on what is selected
  static std::unique_ptr<MeasureFeature> mf = nullptr;
  if (oldItem != item) {
    switch(item) {
      case 0: {
        mf = std::make_unique<SinglePoint>();
      } break;
      case 1: {
        mf = std::make_unique<MeasurementBlob>();
      } break;
      case 2: {
        mf = std::make_unique<MeasurementLine>();
      } break;
      case 3: {
        mf = std::make_unique<Grid2dPoints>();
      } break;
    }
    oldItem = item;
    }

  bool created = false;  
  if (mf->draw_info_gui("Add", _tracerScale, _ips)) {
    mf->generate_draw_geom();
    mfs.emplace_back(std::move(mf));
    mf = nullptr;
    oldItem = -1;
    created = true;
    ImGui::CloseCurrentPopup();
  }

  ImGui::SameLine();
  if (ImGui::Button("Cancel", ImVec2(120,0))) {
    oldItem = -1;
    mf = nullptr;
    ImGui::CloseCurrentPopup();
  }

  ImGui::EndPopup();
  return created;
}
#endif

float MeasureFeature::jitter(const float _z, const float _ips) const {
  // set up the random number generator
  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<float> dist(-0.5, 0.5);
  // emits one per step, jittered slightly
  return _z+_ips*dist(gen);
}

//
// Create a single measurement point
//
ElementPacket<float>
SinglePoint::init_elements(float _ips) const {
  // created once
  std::vector<float> x = {m_x, m_y, m_z};
  std::vector<Int> idx;
  std::vector<float> vals;
  ElementPacket<float> packet({x, idx, vals, (size_t)1, (uint8_t)0});
  if (packet.verify(packet.x.size(), Dimensions)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

ElementPacket<float>
SinglePoint::step_elements(float _ips) const {
  if ((m_enabled) && (m_emits)) {
    std::vector<float> x = {jitter(m_x, _ips), jitter(m_y, _ips), jitter(m_z, _ips)};
    std::vector<Int> idx;
    std::vector<float> vals;
    ElementPacket<float> packet({x, idx, vals, (size_t)1, (uint8_t)0});
    if (packet.verify(packet.x.size(), Dimensions)) {
      return packet;
    } else {
      return ElementPacket<float>();
    }
  } else {
    return ElementPacket<float>();
  }
}

void
SinglePoint::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SinglePoint::to_string() const {
  std::stringstream ss;
  if (m_emits) {
    ss << "emiter";
  } else if (m_is_lagrangian) {
    ss << "tracer";
  } else {
    ss << "stationary";
  }
  ss << " point at " << m_x << " " << m_y << " " << m_z;
  return ss.str();
}

void
SinglePoint::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
  m_enabled = j.value("enabled", true);
  m_is_lagrangian = j.value("lagrangian", m_is_lagrangian);
  m_emits = j.value("emits", m_emits);
}

nlohmann::json
SinglePoint::to_json() const {
  nlohmann::json j;
  j["type"] = "point";
  j["center"] = {m_x, m_y, m_z};
  j["enabled"] = m_enabled;
  j["lagrangian"] = m_is_lagrangian;
  j["emits"] = m_emits;
  return j;
}

void SinglePoint::generate_draw_geom() {
  const float diam = 0.005;
  std::unique_ptr<Ovoid> tmp = std::make_unique<Ovoid>(nullptr, true, m_x, m_y, m_z, diam, diam, diam);
  m_draw = tmp->init_elements(diam/25.0);

  const int numPts = m_draw.val.size()/Dimensions;
  m_draw.val.resize(numPts);
}

#ifdef USE_IMGUI
bool SinglePoint::draw_info_gui(const std::string _action, const float &_tracerScale,
                                const float _ips) {
  float xc[3] = {m_x, m_y, m_z};
  const std::string buttonText = _action+" single point";
  bool add = false;

  ImGui::InputFloat3("position", xc);
  if (!m_emits) {
    ImGui::Checkbox("Point follows flow", &m_is_lagrangian);
  }
  if (!m_is_lagrangian) {
    ImGui::Checkbox("Point emits particles", &m_emits);
  }
  ImGui::TextWrapped("\nThis feature will add 1 point");
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  
  return add;
}
#endif

//
// Create a circle of tracer points
//
ElementPacket<float>
MeasurementBlob::init_elements(float _ips) const {

  // set up the random number generator
  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<float> zmean_dist(-0.5, 0.5);

  // create a new vector to pass on
  std::vector<float> x;
  std::vector<Int> idx;
  std::vector<float> vals;

  // what size 2D integer array will we loop over
  int irad = 1 + m_rad / _ips;
  //std::cout << "blob needs " << (-irad) << " to " << irad << " spaces" << std::endl;

  // loop over integer indices
  for (int i=-irad; i<=irad; ++i) {
  for (int j=-irad; j<=irad; ++j) {
  for (int k=-irad; k<=irad; ++k) {

    // how far from the center are we?
    float dr = sqrt((float)(i*i+j*j+k*k)) * _ips;
    if (dr < m_rad) {
      // create a particle here
      x.emplace_back(m_x + _ips*((float)i+zmean_dist(gen)));
      x.emplace_back(m_y + _ips*((float)j+zmean_dist(gen)));
      x.emplace_back(m_z + _ips*((float)k+zmean_dist(gen)));
    }
  }
  }
  }
  
  ElementPacket<float> packet({x, idx, vals, (size_t)(x.size()/2), (uint8_t)0});
  if (packet.verify(packet.x.size(), Dimensions)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

ElementPacket<float>
MeasurementBlob::step_elements(float _ips) const {
  if ((m_enabled) && (m_emits)) {
    ElementPacket<float> packet = init_elements(_ips);
    for (size_t i=0; i<packet.x.size(); i++) {
      packet.x[i] = jitter(packet.x[i], _ips);
    }
    if (packet.verify(packet.x.size(), Dimensions)) {
      return packet;
    } else {
      return ElementPacket<float>();
    }
  } else {
    return ElementPacket<float>();
  }
}

void
MeasurementBlob::debug(std::ostream& os) const {
  os << to_string();
}

std::string
MeasurementBlob::to_string() const {
  std::stringstream ss;
  if (m_emits) {
    ss << "emiter";
  } else if (m_is_lagrangian) {
    ss << "tracer";
  } else {
    ss << "stationary";
  }
  ss << " blob at " << m_x << " " << m_y << " " << m_z << " with radius " << m_rad;
  return ss.str();
}

void
MeasurementBlob::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
  m_rad = j["rad"];
  m_enabled = j.value("enabled", true);
  m_is_lagrangian = j.value("lagrangian", m_is_lagrangian);
  m_emits = j.value("emits", m_emits);
}

nlohmann::json
MeasurementBlob::to_json() const {
  nlohmann::json j;
  j["type"] = "tracer blob";
  j["center"] = {m_x, m_y, m_z};
  j["rad"] = m_rad;
  j["enabled"] = m_enabled;
  j["lagrangian"] = m_is_lagrangian;
  j["emits"] = m_emits;
  return j;
}

void MeasurementBlob::generate_draw_geom() {
  std::unique_ptr<Ovoid> tmp = std::make_unique<Ovoid>(nullptr, true, m_x, m_y, m_z, m_rad*2.0, m_rad*2.0, m_rad*2.0);
  m_draw = tmp->init_elements(m_rad/12.5);

  const int numPts = m_draw.val.size()/Dimensions;
  m_draw.val.resize(numPts);
}

#ifdef USE_IMGUI
bool MeasurementBlob::draw_info_gui(const std::string _action, const float &_tracerScale, float _ips) {
  float xc[3] = {m_x, m_y, m_z};
  const std::string buttonText = _action+" circle of tracers";
  bool add = false;

  ImGui::InputFloat3("center", xc);
  ImGui::SliderFloat("radius", &m_rad, 0.5f*_ips, 0.5f, "%.4f");
  if (!m_emits) {
    ImGui::Checkbox("Point follows flow", &m_is_lagrangian);
  }
  if (!m_is_lagrangian) {
    ImGui::Checkbox("Point emits particles", &m_emits);
  }
  ImGui::TextWrapped("This feature will add about %d field points",
                     (int)(0.6*std::pow(2*m_rad/(_tracerScale*_ips), 3)));
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  m_z = xc[2];

  return add;
}
#endif

//
// Create a line of static measurement points
//
ElementPacket<float>
MeasurementLine::init_elements(float _ips) const {

  // create a new vector to pass on
  std::vector<float> x;
  std::vector<Int> idx;
  std::vector<float> vals;

  // how many points do we need?
  float llen = std::sqrt( std::pow(m_xf-m_x, 2) + std::pow(m_yf-m_y, 2) + std::pow(m_zf-m_z, 2));
  int ilen = 1 + llen / _ips;

  // loop over integer indices
  for (int i=0; i<ilen; ++i) {
    // how far along the line?
    float frac = (float)i / (float)(ilen-1);

    // create a particle here
    x.emplace_back((1.0-frac)*m_x + frac*m_xf);
    x.emplace_back((1.0-frac)*m_y + frac*m_yf);
    x.emplace_back((1.0-frac)*m_z + frac*m_zf);
  }

  ElementPacket<float> packet({x, idx, vals, (size_t)(Dimensions*ilen), (uint8_t)0});
  if (packet.verify(packet.x.size(), Dimensions)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

ElementPacket<float>
MeasurementLine::step_elements(float _ips) const {
  if ((m_enabled) && (m_emits)) {
    ElementPacket<float> packet = init_elements(_ips);
    for (size_t i=0; i<packet.x.size(); i++) {
      packet.x[i] = jitter(packet.x[i], _ips);
    }
    if (packet.verify(packet.x.size(), Dimensions)) {
      return packet;
    } else {
      return ElementPacket<float>();
    }
  } else {
    return ElementPacket<float>();
  }
}

void
MeasurementLine::debug(std::ostream& os) const {
  os << to_string();
}

std::string
MeasurementLine::to_string() const {
  std::stringstream ss;
  if (m_emits) {
    ss << "emiter";
  } else if (m_is_lagrangian) {
    ss << "tracer";
  } else {
    ss << "stationary";
  }
  ss << " line from " << m_x << " " << m_y << " " << m_z << " to " << m_xf << " " << m_yf << " " << m_zf << " with dx " << m_dx;
  return ss.str();
}

void
MeasurementLine::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
  const std::vector<float> e = j["end"];
  m_xf = e[0];
  m_yf = e[1];
  m_zf = e[2];
  m_dx = j.value("dx", 0.1);
  m_enabled = j.value("enabled", true);
  m_is_lagrangian = j.value("lagrangian", m_is_lagrangian);
  m_emits = j.value("emits", m_emits);
}

nlohmann::json
MeasurementLine::to_json() const {
  nlohmann::json j;
  j["type"] = "measurement line";
  j["center"] = {m_x, m_y, m_z};
  j["end"] = {m_xf, m_yf, m_zf};
  j["dx"] = m_dx;
  j["enabled"] = m_enabled;
  j["lagrangian"] = m_is_lagrangian;
  j["emits"] = m_emits;
  return j;
}

void MeasurementLine::generate_draw_geom() {
  const float minS = 0.01;
  std::unique_ptr<SolidRect> tmp = std::make_unique<SolidRect>(nullptr, true, m_x, m_y, m_z,
                                                               std::max(minS, m_xf-m_x), std::max(minS, m_yf-m_y),
                                                               std::max(minS, m_zf-m_z));
  m_draw = tmp->init_elements(1.0);

  const int numPts = m_draw.val.size()/Dimensions;
  m_draw.val.resize(numPts);
}

#ifdef USE_IMGUI
bool MeasurementLine::draw_info_gui(const std::string _action, const float &_tracerScale, float _ips) {
  float xc[3] = {m_x, m_y, m_z};
  float xf[3] = {m_xf, m_yf, m_zf};
  const std::string buttonText = _action+" line of measurement points";
  bool add = false;
  
  ImGui::InputFloat3("start", xc);
  ImGui::InputFloat3("finish", xf);
  ImGui::TextWrapped("This feature will add about %d field points",
                     1+(int)(std::sqrt(std::pow(xf[0]-xc[0],2)+std::pow(xf[1]-xc[1],2)+std::pow(xf[2]-xc[2],2))/(_tracerScale*_ips)));
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  m_z = xc[2];
  m_xf = xf[0];
  m_yf = xf[1];
  m_zf = xf[2];

  return add;
}
#endif

//
// Create a 2D grid of static measurement points
//
ElementPacket<float>
Grid2dPoints::init_elements(float _ips) const {

  // create a new vector to pass on
  std::vector<float> x;
  std::vector<Int> idx;
  std::vector<float> vals;

  // ignore _ips and use m_dx to define grid density

  // calculate length of two axes
  const float dist_s = std::sqrt( std::pow(m_xs, 2) + std::pow(m_ys, 2) + std::pow(m_zs, 2));
  const float dist_f = std::sqrt( std::pow(m_xf, 2) + std::pow(m_yf, 2) + std::pow(m_zf, 2));

  // loop over integer indices
  for (float sp=0.5*m_ds/dist_s; sp<1.0+0.01*m_ds/dist_s; sp+=m_ds/dist_s) {
    for (float fp=0.5*m_df/dist_f; fp<1.0+0.01*m_df/dist_f; fp+=m_df/dist_f) {
      // create a field point here
      x.emplace_back(m_x + m_xs*sp + m_xf*fp);
      x.emplace_back(m_y + m_ys*sp + m_yf*fp);
      x.emplace_back(m_z + m_zs*sp + m_zf*fp);
    }
  }

  ElementPacket<float> packet({x, idx, vals, (size_t)(x.size()/Dimensions), (uint8_t)0});
  if (packet.verify(packet.x.size(), Dimensions)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

ElementPacket<float>
Grid2dPoints::step_elements(float _ips) const {
  // does not emit
  return ElementPacket<float>();
}

void
Grid2dPoints::debug(std::ostream& os) const {
  os << to_string();
}

std::string
Grid2dPoints::to_string() const {
  std::stringstream ss;
  ss << "measurement plane at " << m_x << " " << m_y << " " << m_z << " with ds,df " << m_ds << " " << m_df;
  return ss.str();
}

void
Grid2dPoints::from_json(const nlohmann::json j) {
  const std::vector<float> s = j["start"];
  m_x = s[0];
  m_y = s[1];
  m_z = s[2];
  const std::vector<float> e = j["axis1"];
  m_xs = e[0];
  m_ys = e[1];
  m_zs = e[2];
  const std::vector<float> f = j["axis2"];
  m_xf = f[0];
  m_yf = f[1];
  m_zf = f[2];
  const std::vector<float> d = j["dx"];
  m_ds = d[0];
  m_df = d[1];
  m_enabled = j.value("enabled", true);
  m_is_lagrangian = j.value("lagrangian", m_is_lagrangian);
  m_emits= j.value("emits", m_emits);
}

nlohmann::json
Grid2dPoints::to_json() const {
  nlohmann::json j;
  j["type"] = "measurement plane";
  j["start"] = {m_x, m_y, m_z};
  j["axis1"] = {m_xs, m_ys, m_zs};
  j["axis2"] = {m_xf, m_yf, m_zf};
  j["dx"] = {m_ds, m_df};
  j["enabled"] = m_enabled;
  j["lagrangian"] = m_is_lagrangian;
  j["emits"] = m_emits;
  return j;
}

void Grid2dPoints::generate_draw_geom() {
  const float normS = std::sqrt( std::pow(m_xs, 2) + std::pow(m_ys, 2) + std::pow(m_zs, 2));
  const float normF = std::sqrt( std::pow(m_xf, 2) + std::pow(m_yf, 2) + std::pow(m_zf, 2));
  const float p0[3] = { m_x, m_y, m_z };
  const float p1[3] = { m_x+(m_xs*0.5f*m_ds/normS)+(m_xf),
                        m_y+(m_ys*0.5f*m_ds/normS)+(m_yf),
                        m_z+(m_zs*0.5f*m_ds/normS)+(m_zf) };
  const float p2[3] = { m_x+(m_xf*0.5f*m_df/normF)+(m_xs),
                        m_y+(m_yf*0.5f*m_df/normF)+(m_ys),
                        m_z+(m_zf*0.5f*m_df/normF)+(m_zs) };
  const float p3[3] = { m_x+(m_xs)+(m_xf),
                        m_y+(m_ys)+(m_yf),
                        m_z+(m_zs)+(m_zf) };

  std::unique_ptr<BoundaryQuad> tmp = std::make_unique<BoundaryQuad>(nullptr, p0[0], p0[1], p0[2],
                                                                              p2[0], p2[1], p2[2], 
                                                                              p3[0], p3[1], p3[2], 
                                                                              p1[0], p1[1], p1[2],
                                                                              0.0, 0.0, 0.0);
  m_draw = tmp->init_elements(1.0);
  
  const int numPts = m_draw.val.size()/Dimensions;
  m_draw.val.resize(numPts);
}

#ifdef USE_IMGUI
bool Grid2dPoints::draw_info_gui(const std::string _action, const float &tracer_scale, const float _ips) {
  float xc[3] = {m_x, m_y, m_z};
  float xs[3] = {m_xs, m_ys, m_zs};
  float xf[3] = {m_xf, m_yf, m_zf};
  float dx[2] = {m_ds, m_df};
  const std::string buttonText = _action+" 2D grid of measurement points";
  bool add = false;
  
  ImGui::InputFloat3("corner", xc);
  ImGui::InputFloat3("axis 1", xs);
  ImGui::InputFloat3("axis 2", xf);
  ImGui::InputFloat2("dx", dx);
  ImGui::TextWrapped("This feature will add about %d field points",
                     1+(int)(std::sqrt(xs[0]*xs[0]+xs[1]*xs[1]+xs[2]*xs[2])*
                             std::sqrt(xf[0]*xf[0]+xf[1]*xf[1]+xf[2]*xf[2])/
                             (dx[0]*dx[1])));
  if (ImGui::Button("Add 2D grid of measurement points")) { add = true;}
  m_x = xc[0];
  m_y = xc[1];
  m_z = xc[2];
  m_xs = xs[0];
  m_ys = xs[1];
  m_zs = xs[2];
  m_xf = xf[0];
  m_yf = xf[1];
  m_zf = xf[2];
  m_ds = dx[0];
  m_df = dx[1];

  return add;
}
#endif
