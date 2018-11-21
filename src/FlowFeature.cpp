#include "FlowFeature.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <random>

// write out any object of parent type FlowFeature by dispatching to appropriate "debug" method
std::ostream& operator<<(std::ostream& os, FlowFeature const& ff) {
  ff.debug(os);
  return os;
}

// various debug print methods for the subclasses
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
  return std::vector<float>({m_x, m_y, m_z, m_sx, m_sy, m_sz, 0.0});
}

std::vector<float>
SingleParticle::step_particles(float _ips) const {
  return std::vector<float>();
}


//
// make a circular vortex blob with soft transition
//
std::vector<float>
VortexBlob::init_particles(float _ips) const {
  // create a new vector to pass on
  std::vector<float> x;

  // what size 2D integer array will we loop over
  int irad = 1 + (m_rad + 0.5*m_softness) / _ips;
  std::cout << "blob needs " << (-irad) << " to " << irad << " spaces" << std::endl;

  // and a counter for the total circulation
  double tot_wgt = 0.0;

  // will need this
  const double pi = std::acos(-1);

  // loop over integer indices
  for (int i=-irad; i<=irad; ++i) {
  for (int j=-irad; j<=irad; ++j) {
  for (int k=-irad; k<=irad; ++k) {

    // how far from the center are we?
    float dr = sqrt((float)(i*i+j*j+k*k)) * _ips;
    if (dr < m_rad + 0.5*m_softness) {

      // create a particle here
      x.emplace_back(m_x + _ips*(float)i);
      x.emplace_back(m_y + _ips*(float)j);
      x.emplace_back(m_z + _ips*(float)k);

      // figure out the strength from another check
      double this_wgt = 1.0;
      if (dr > m_rad - 0.5*m_softness) {
        // create a weaker particle
        this_wgt = 0.5 - 0.5*std::sin(pi * (dr - m_rad) / m_softness);
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


//
// make the block of randomly-placed and random-strength particles
//
std::vector<float>
BlockOfRandom::init_particles(float _ips) const {
  // set up the random number generator
  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> zmean_dist(-1.0, 1.0);
  static std::uniform_real_distribution<> zo_dist(0.0, 1.0);

  std::vector<float> x(4*m_num);
  // initialize the particles' locations and strengths, leave radius zero for now
  for (size_t i=0; i<(size_t)m_num; ++i) {
    size_t idx = 7*i;
    // positions
    x[idx+0] = m_x + m_xsize*zmean_dist(gen);
    x[idx+1] = m_y + m_ysize*zmean_dist(gen);
    x[idx+2] = m_z + m_zsize*zmean_dist(gen);
    // strengths
    x[idx+3] = m_maxstr * zmean_dist(gen);
    x[idx+4] = m_maxstr * zmean_dist(gen);
    x[idx+5] = m_maxstr * zmean_dist(gen);
    // radius is zero still
    x[idx+6] = 0.0f;
  }
  return x;
}

std::vector<float>
BlockOfRandom::step_particles(float _ips) const {
  return std::vector<float>();
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
  return std::vector<float>({m_x, m_y, m_z, m_sx, m_sy, m_sz, 0.0});
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

