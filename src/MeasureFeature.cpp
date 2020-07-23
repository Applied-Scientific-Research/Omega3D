/*
 * MeasureFeature.cpp - GUI-side descriptions of flow features
 *
 * (c)2018-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

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

  if      (ftype == "tracer") {           _flist.emplace_back(std::make_unique<SinglePoint>()); }
  else if (ftype == "tracer emitter") {   _flist.emplace_back(std::make_unique<TracerEmitter>()); }
  else if (ftype == "tracer blob") {      _flist.emplace_back(std::make_unique<TracerBlob>()); }
  else if (ftype == "tracer line") {      _flist.emplace_back(std::make_unique<TracerLine>()); }
  else if (ftype == "measurement line") { _flist.emplace_back(std::make_unique<MeasurementLine>()); }
  else if (ftype == "measurement plane") { _flist.emplace_back(std::make_unique<Grid2dPoints>()); }
  else {
    std::cout << "  type " << ftype << " does not name an available measurement feature, ignoring" << std::endl;
    return;
  }

  // and pass the json object to the specific parser
  _flist.back()->from_json(_jin);

  std::cout << "  found " << ftype << std::endl;
}

#ifdef USE_IMGUI
void MeasureFeature::draw_creation_gui(std::vector<std::unique_ptr<MeasureFeature>> &mfs, const float ips, const float &tracerScale) {
  static int item = 0;
  const char* items[] = { "single point/tracer", "streakline", "circle of tracers", "line of tracers", "measurement line", "2D grid" };
  ImGui::Combo("type", &item, items, 6);

  // show different inputs based on what is selected
  std::unique_ptr<MeasureFeature> mf = nullptr;
  switch(item) {
    case 0: {
      mf = std::make_unique<SinglePoint>();
    } break;
    case 1: {
      mf = std::make_unique<TracerEmitter>();
    } break;
    case 2: {
      mf = std::make_unique<TracerBlob>();
    } break;
    case 3: {
      mf = std::make_unique<TracerLine>();
    } break;
    case 4: {
      mf = std::make_unique<MeasurementLine>();
    } break;
    case 5: {
      mf = std::make_unique<Grid2dPoints>();
    } break;
  }

  if (mf->draw_info_gui("Add", tracerScale, ips)) {
    //std::cout << "Added " << mf << std::endl;
    mfs.emplace_back(std::move(mf));
    mf = nullptr;
    ImGui::CloseCurrentPopup();
  }

  ImGui::SameLine();
  if (ImGui::Button("Cancel", ImVec2(120,0))) {
    mf = nullptr;
    ImGui::CloseCurrentPopup();
  }

  ImGui::EndPopup();
}
#endif

//
// Create a single measurement point
//
std::vector<float>
SinglePoint::init_particles(float _ips) const {
  // created once
  if (this->is_enabled()) return std::vector<float>({m_x, m_y, m_z});
  else return std::vector<float>();
}

std::vector<float>
SinglePoint::step_particles(float _ips) const {
  // does not emit
  return std::vector<float>();
}

void
SinglePoint::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SinglePoint::to_string() const {
  std::stringstream ss;
  ss << "single field point at " << m_x << " " << m_y << " " << m_z;
  return ss.str();
}

void
SinglePoint::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
}

nlohmann::json
SinglePoint::to_json() const {
  nlohmann::json j;
  j["type"] = "tracer";
  j["center"] = {m_x, m_y, m_z};
  return j;
}

#ifdef USE_IMGUI
bool SinglePoint::draw_info_gui(const std::string action, const float ips, const float &tracerScale) {
  static float xc[3] = {m_x, m_y, m_z};
  static bool lagrangian = m_is_lagrangian;
  const std::string buttonText = action+" single point";
  bool add = false;

  ImGui::InputFloat3("position", xc);
  ImGui::Checkbox("Point follows flow", &lagrangian);
  ImGui::TextWrapped("This feature will add 1 point");
  if (ImGui::Button(buttonText.c_str())) {
    m_x = xc[0];
    m_y = xc[1];
    m_z = xc[2];
    add = true;
  }

  return add;
}
#endif

// Create a single, stable point which emits Lagrangian points
//
std::vector<float>
TracerEmitter::init_particles(float _ips) const {
  // is not a measurement point in itself
  // but if it was, we could use the local velocity to help generate points at any given time
  return std::vector<float>();
}

std::vector<float>
TracerEmitter::step_particles(float _ips) const {

  if (not this->is_enabled()) return std::vector<float>();

  // set up the random number generator
  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<float> zmean_dist(-0.5, 0.5);

  // emits one per step, jittered slightly
  return std::vector<float>({m_x + _ips*zmean_dist(gen),
                             m_y + _ips*zmean_dist(gen),
                             m_z + _ips*zmean_dist(gen)});
}

void
TracerEmitter::debug(std::ostream& os) const {
  os << to_string();
}

std::string
TracerEmitter::to_string() const {
  std::stringstream ss;
  ss << "tracer emitter at " << m_x << " " << m_y << " " << m_z << " spawning tracers every step";
  return ss.str();
}

void
TracerEmitter::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
}

nlohmann::json
TracerEmitter::to_json() const {
  nlohmann::json j;
  j["type"] = "tracer emitter";
  j["center"] = {m_x, m_y, m_z};
  return j;
}

#ifdef USE_IMGUI
bool TracerEmitter::draw_info_gui(const std::string action, const float ips,
                                      const float &tracer_scale) {
  static float xc[3] = {m_x, m_y, m_z};
  const std::string buttonText = action+" streakline";
  bool add = false;
  
  ImGui::InputFloat3("position", xc);
  ImGui::TextWrapped("This feature will add 1 tracer emitter");
  if (ImGui::Button("Add streakline")) {
    m_x = xc[0];
    m_y = xc[1];
    m_z = xc[2];
    add = true;
  }

  return add;
}
#endif

//
// Create a blob of tracer points
//
std::vector<float>
TracerBlob::init_particles(float _ips) const {

  if (not this->is_enabled()) return std::vector<float>();

  // set up the random number generator
  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<float> zmean_dist(-0.5, 0.5);

  // create a new vector to pass on
  std::vector<float> x;

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

  return x;
}

std::vector<float>
TracerBlob::step_particles(float _ips) const {
  // does not emit
  return std::vector<float>();
}

void
TracerBlob::debug(std::ostream& os) const {
  os << to_string();
}

std::string
TracerBlob::to_string() const {
  std::stringstream ss;
  ss << "tracer blob at " << m_x << " " << m_y << " " << m_z << " with radius " << m_rad;
  return ss.str();
}

void
TracerBlob::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
  m_rad = j["rad"];
}

nlohmann::json
TracerBlob::to_json() const {
  nlohmann::json j;
  j["type"] = "tracer blob";
  j["center"] = {m_x, m_y, m_z};
  j["rad"] = m_rad;
  return j;
}

#ifdef USE_IMGUI
bool TracerBlob::draw_info_gui(const std::string action, const float ips,
                               const float &tracerScale) {
  static float xc[3] = {m_x, m_y, m_z};
  static float rad = m_rad;
  const std::string buttonText = action+" circle of tracers";
  bool add = false;

  ImGui::InputFloat3("center", xc);
  ImGui::SliderFloat("radius", &rad, 0.5f*ips, 0.5f, "%.4f");
  ImGui::TextWrapped("This feature will add about %d field points",
                     (int)(0.6*std::pow(2*rad/(tracerScale*ips), 3)));
  if (ImGui::Button("Add circle of tracers")) {
    m_x = xc[0];
    m_y = xc[1];
    m_z = xc[2];
    m_rad = rad;
    add = true;
  }

  return add;
}
#endif

//
// Create a line of tracer points
//
std::vector<float>
TracerLine::init_particles(float _ips) const {

  // create a new vector to pass on
  std::vector<float> x;

  if (not this->is_enabled()) return x;

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

  return x;
}

std::vector<float>
TracerLine::step_particles(float _ips) const {
  // does not emit
  return std::vector<float>();
}

void
TracerLine::debug(std::ostream& os) const {
  os << to_string();
}

std::string
TracerLine::to_string() const {
  std::stringstream ss;
  ss << "tracer line from " << m_x << " " << m_y << " " << m_z << " to " << m_xf << " " << m_yf << " " << m_zf;
  return ss.str();
}

void
TracerLine::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_z = c[2];
  const std::vector<float> e = j["end"];
  m_xf = e[0];
  m_yf = e[1];
  m_zf = e[2];
}

nlohmann::json
TracerLine::to_json() const {
  nlohmann::json j;
  j["type"] = "tracer line";
  j["center"] = {m_x, m_y, m_z};
  j["end"] = {m_xf, m_yf, m_zf};
  return j;
}

#ifdef USE_IMGUI
bool TracerLine::draw_info_gui(const std::string action, const float ips, const float &tracerScale) {
  static float xc[3] = {m_x, m_y, m_z};
  static float xf[3] = {m_xf, m_yf, m_zf};
  const std::string buttonText = action+" line of tracers";
  bool add = false;

  ImGui::InputFloat3("start", xc);
  ImGui::InputFloat3("finish", xf);
  ImGui::TextWrapped("This feature will add about %d field points",
                     1+(int)(std::sqrt(std::pow(xf[0]-xc[0],2)+std::pow(xf[1]-xc[1],2)+std::pow(xf[2]-xc[2],2))/(tracerScale*ips)));
  if (ImGui::Button(buttonText.c_str())) {
    m_x = xc[0];
    m_y = xc[1];
    m_z = xc[2];
    m_xf = xf[0];
    m_yf = xf[1];
    m_zf = xf[2];
    add = true;
  }

  return add;
}
#endif

//
// Create a line of static measurement points
//
std::vector<float>
MeasurementLine::init_particles(float _ips) const {

  // create a new vector to pass on
  std::vector<float> x;

  if (not this->is_enabled()) return x;

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

  return x;
}

std::vector<float>
MeasurementLine::step_particles(float _ips) const {
  // does not emit
  return std::vector<float>();
}

void
MeasurementLine::debug(std::ostream& os) const {
  os << to_string();
}

std::string
MeasurementLine::to_string() const {
  std::stringstream ss;
  ss << "measurement line from " << m_x << " " << m_y << " " << m_z << " to " << m_xf << " " << m_yf << " " << m_zf;
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
}

nlohmann::json
MeasurementLine::to_json() const {
  nlohmann::json j;
  j["type"] = "measurement line";
  j["center"] = {m_x, m_y, m_z};
  j["end"] = {m_xf, m_yf, m_zf};
  return j;
}

#ifdef USE_IMGUI
bool MeasurementLine::draw_info_gui(const std::string action, const float ips, const float &tracerScale) {
  static float xc[3] = {m_x, m_y, m_z};
  static float xf[3] = {m_xf, m_yf, m_zf};
  const std::string buttonText = action+" line of measurement points";
  bool add = false;
  
  ImGui::InputFloat3("start", xc);
  ImGui::InputFloat3("finish", xf);
  ImGui::TextWrapped("This feature will add about %d field points",
                     1+(int)(std::sqrt(std::pow(xf[0]-xc[0],2)+std::pow(xf[1]-xc[1],2)+std::pow(xf[2]-xc[2],2))/(tracerScale*ips)));
  if (ImGui::Button(buttonText.c_str())) {
    m_x = xc[0];
    m_y = xc[1];
    m_z = xc[2];
    m_xf = xf[0];
    m_yf = xf[1];
    m_zf = xf[2];
    add = true;
  }

  return add;
}
#endif

//
// Create a 2D grid of static measurement points
//
std::vector<float>
Grid2dPoints::init_particles(float _ips) const {

  // create a new vector to pass on
  std::vector<float> x;

  if (not this->is_enabled()) return x;

  // ignore _ips and use m_dx to define grid density

  // calculate length of two axes
  const float dist_s = std::sqrt( std::pow(m_xs, 2) + std::pow(m_ys, 2) + std::pow(m_zs, 2));
  const float dist_t = std::sqrt( std::pow(m_xt, 2) + std::pow(m_yt, 2) + std::pow(m_zt, 2));

  // loop over integer indices
  for (float sp=0.5*m_ds/dist_s; sp<1.0+0.01*m_ds/dist_s; sp+=m_ds/dist_s) {
    for (float tp=0.5*m_dt/dist_t; tp<1.0+0.01*m_dt/dist_t; tp+=m_dt/dist_t) {
      // create a field point here
      x.emplace_back(m_x + m_xs*sp + m_xt*tp);
      x.emplace_back(m_y + m_ys*sp + m_yt*tp);
      x.emplace_back(m_z + m_zs*sp + m_zt*tp);
    }
  }

  return x;
}

std::vector<float>
Grid2dPoints::step_particles(float _ips) const {
  // does not emit
  return std::vector<float>();
}

void
Grid2dPoints::debug(std::ostream& os) const {
  os << to_string();
}

std::string
Grid2dPoints::to_string() const {
  std::stringstream ss;
  ss << "measurement plane at " << m_x << " " << m_y << " " << m_z << " with ds,dt " << m_ds << " " << m_dt;
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
  m_xt = f[0];
  m_yt = f[1];
  m_zt = f[2];
  const std::vector<float> d = j["dx"];
  m_ds = d[0];
  m_dt = d[1];
}

nlohmann::json
Grid2dPoints::to_json() const {
  nlohmann::json j;
  j["type"] = "measurement plane";
  j["start"] = {m_x, m_y, m_z};
  j["axis1"] = {m_xs, m_ys, m_zs};
  j["axis2"] = {m_xt, m_yt, m_zt};
  j["dx"] = {m_ds, m_dt};
  return j;
}

#ifdef USE_IMGUI
bool Grid2dPoints::draw_info_gui(const std::string action, const float ips, const float &tracerScale) {
  static float xc[3] = {m_x, m_y, m_z};
  static float xs[3] = {m_xs, m_ys, m_zs};
  static float xt[3] = {m_xt, m_yt, m_zt};
  static float dx[2] = {m_ds, m_dt};
  const std::string buttonText = action+" 2D grid of measurement points";
  bool add = false;
  
  ImGui::InputFloat3("corner", xc);
  ImGui::InputFloat3("axis 1", xs);
  ImGui::InputFloat3("axis 2", xt);
  ImGui::InputFloat2("dx", dx);
  ImGui::TextWrapped("This feature will add about %d field points",
                     1+(int)(std::sqrt(xs[0]*xs[0]+xs[1]*xs[1]+xs[2]*xs[2])*
                             std::sqrt(xt[0]*xt[0]+xt[1]*xt[1]+xt[2]*xt[2])/
                             (dx[0]*dx[1])));
  if (ImGui::Button("Add 2D grid of measurement points")) {
    m_x = xc[0];
    m_y = xc[1];
    m_z = xc[2];
    m_xs = xs[0];
    m_ys = xs[1];
    m_zs = xs[2];
    m_ds = dx[0];
    m_dt = dx[1];
    add = true;
  }

  return add;
}
#endif
