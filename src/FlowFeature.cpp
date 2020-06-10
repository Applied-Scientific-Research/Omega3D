/*
 * FlowFeature.cpp - GUI-side descriptions of flow features
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#include "FlowFeature.h"
#include "MathHelper.h"

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


//
// important feature: convert flow feature definition into 1-D vector of values
//
// each 7 floats is one particle's: [xyz] location, [xyz] strength, vdelta (radius)
//

//
// drop a single particle
//
std::vector<float>
SingleParticle::init_particles(float _ips) const {
  if (this->is_enabled()) return std::vector<float>({m_x, m_y, m_z, m_sx, m_sy, m_sz, 0.0});
  else return std::vector<float>();
}

std::vector<float>
SingleParticle::step_particles(float _ips) const {
  return std::vector<float>();
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
}

nlohmann::json
SingleParticle::to_json() const {
  nlohmann::json j;
  j["type"] = "single particle";
  j["center"] = {m_x, m_y, m_z};
  j["strength"] = {m_sx, m_sy, m_sz};
  return j;
}


//
// make a circular vortex blob with soft transition
//
std::vector<float>
VortexBlob::init_particles(float _ips) const {
  // create a new vector to pass on
  std::vector<float> x;

  if (not this->is_enabled()) return x;

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
      tot_wgt += this_wgt;
      x.emplace_back(m_sx * (float)this_wgt);
      x.emplace_back(m_sy * (float)this_wgt);
      x.emplace_back(m_sz * (float)this_wgt);

      // this is the radius - still zero for now
      x.emplace_back(0.0f);
    }
  }
  }
  }

  // finally, normalize all particle strengths so that the whole blob
  //   has exactly the right strength
  std::cout << "blob had " << tot_wgt << " initial circulation" << std::endl;
  double str_scale = 1.0 / tot_wgt;
  for (size_t i=3; i<x.size(); i+=7) {
    x[i+0] = (float)((double)x[i+0] * str_scale);
    x[i+1] = (float)((double)x[i+1] * str_scale);
    x[i+2] = (float)((double)x[i+2] * str_scale);
  }

  return x;
}

std::vector<float>
VortexBlob::step_particles(float _ips) const {
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
}

nlohmann::json
VortexBlob::to_json() const {
  nlohmann::json j;
  j["type"] = "vortex blob";
  j["center"] = {m_x, m_y, m_z};
  j["radius"] = m_rad;
  j["softness"] = m_softness;
  j["strength"] = {m_sx, m_sy, m_sz};
  return j;
}


//
// make the block of randomly-placed and random-strength particles
//
std::vector<float>
BlockOfRandom::init_particles(float _ips) const {

  if (not this->is_enabled()) return std::vector<float>();

  // set up the random number generator
  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> zmean_dist(-0.5, 0.5);
  static std::uniform_real_distribution<> zo_dist(0.0, 1.0);

  std::vector<float> x(7*m_num);
  // initialize the particles' locations and strengths, leave radius zero for now
  for (size_t i=0; i<(size_t)m_num; ++i) {
    size_t idx = 7*i;
    // positions
    x[idx+0] = m_x + m_xsize*zmean_dist(gen);
    x[idx+1] = m_y + m_ysize*zmean_dist(gen);
    x[idx+2] = m_z + m_zsize*zmean_dist(gen);
    // strengths
    x[idx+3] = m_maxstr * zmean_dist(gen) / (float)m_num;
    x[idx+4] = m_maxstr * zmean_dist(gen) / (float)m_num;
    x[idx+5] = m_maxstr * zmean_dist(gen) / (float)m_num;
    // radius will get set later
    x[idx+6] = 0.0f;
  }
  return x;
}

std::vector<float>
BlockOfRandom::step_particles(float _ips) const {
  return std::vector<float>();
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
}

nlohmann::json
BlockOfRandom::to_json() const {
  nlohmann::json j;
  j["type"] = "block of random";
  j["center"] = {m_x, m_y, m_z};
  j["size"] = {m_xsize, m_ysize, m_zsize};
  j["max strength"] = m_maxstr;
  j["num"] = m_num;
  return j;
}


//
// drop a single particle from the emitter
//
std::vector<float>
ParticleEmitter::init_particles(float _ips) const {
  return std::vector<float>();
}

std::vector<float>
ParticleEmitter::step_particles(float _ips) const {
  if (this->is_enabled()) return std::vector<float>({m_x, m_y, m_z, m_sx, m_sy, m_sz, 0.0});
  else return std::vector<float>();
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
}

nlohmann::json
ParticleEmitter::to_json() const {
  nlohmann::json j;
  j["type"] = "particle emitter";
  j["center"] = {m_x, m_y, m_z};
  j["strength"] = {m_sx, m_sy, m_sz};
  return j;
}


//
// make a singular (one row) vortex ring
//
std::vector<float>
SingularRing::init_particles(float _ips) const {
  // create a new vector to pass on
  std::vector<float> x;

  if (not this->is_enabled()) return x;

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
    x.emplace_back(this_ips * m_circ * (b2[0]*std::cos(theta) - b1[0]*std::sin(theta)));
    x.emplace_back(this_ips * m_circ * (b2[1]*std::cos(theta) - b1[1]*std::sin(theta)));
    x.emplace_back(this_ips * m_circ * (b2[2]*std::cos(theta) - b1[2]*std::sin(theta)));

    // this is the radius - still zero for now
    x.emplace_back(0.0f);
  }

  return x;
}

std::vector<float>
SingularRing::step_particles(float _ips) const {
  return std::vector<float>();
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
}

nlohmann::json
SingularRing::to_json() const {
  nlohmann::json j;
  j["type"] = "singular ring";
  j["center"] = {m_x, m_y, m_z};
  j["normal"] = {m_nx, m_ny, m_nz};
  j["major radius"] = m_majrad;
  j["circulation"] = m_circ;
  return j;
}


//
// make a thick-cored vortex ring
//
std::vector<float>
ThickRing::init_particles(float _ips) const {

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

  // create a new vector with the particle values
  std::vector<float> x;

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
      x.emplace_back(sscale * (b2[0]*ct - b1[0]*st));
      x.emplace_back(sscale * (b2[1]*ct - b1[1]*st));
      x.emplace_back(sscale * (b2[2]*ct - b1[2]*st));

      // this is the radius - still zero for now
      x.emplace_back(0.0f);
    }
  }

  return x;
}

std::vector<float>
ThickRing::step_particles(float _ips) const {
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
  return j;
}


//
// various GUI draw methods for the subclasses
//
/*
void
SingleParticle::draw_creation_gui(std::vector< std::unique_ptr<FlowFeature> >& features) {
  static float xc[2] = {0.0f, 0.0f};
  static float str = 1.0f;
  ImGui::InputFloat2("center", xc);
  ImGui::SliderFloat("strength", &str, -1.0f, 1.0f, "%.4f");
  ImGui::TextWrapped("This feature will add 1 particle");
  if (ImGui::Button("Add single particle")) {
    // this is C++14
    //features.emplace_back(std::make_unique<SingleParticle>(xc[0], xc[1], str));
    // this is C++11
    features.emplace_back(std::unique_ptr<SingleParticle>(new SingleParticle(xc[0], xc[1], str)));
    std::cout << "Added " << (*features.back()) << std::endl;
    ImGui::CloseCurrentPopup();
  }
  ImGui::SameLine();
}
*/

