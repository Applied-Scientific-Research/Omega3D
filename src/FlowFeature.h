#pragma once

#include <iostream>
#include <vector>

//
// Abstract class for any flow feature present initially
//
class FlowFeature {
public:
  explicit
  FlowFeature(float _x, float _y, float _z)
    : m_x(_x),
      m_y(_y),
      m_z(_z)
    {}

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual std::vector<float> init_particles(float) const = 0;
  virtual std::vector<float> step_particles(float) const = 0;

  // emit particles as vector of float4

protected:
  float m_x;
  float m_y;
  float m_z;
};

std::ostream& operator<<(std::ostream& os, FlowFeature const& ff);


//
// make intermediate abstract classes for flow and solid elements?
//


//
// Concrete class for a single particle
//
class SingleParticle : public FlowFeature {
public:
  SingleParticle(float _x, float _y, float _z, float _sx, float _sy, float _sz)
    : FlowFeature(_x, _y, _z),
      m_sx(_sx), m_sy(_sy), m_sz(_sz)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

protected:
  float m_sx, m_sy, m_sz;
};


//
// Concrete class for a vortex blob
//
class VortexBlob : public SingleParticle {
public:
  VortexBlob(float _x, float _y, float _z, float _sx, float _sy, float _sz, float _rad, float _soft)
    : SingleParticle(_x, _y, _z, _sx, _sy, _sz),
      m_rad(_rad),
      m_softness(_soft)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

private:
  float m_rad;
  float m_softness;
};


//
// Concrete class for a rectangle of randomly-placed particles
//
class BlockOfRandom : public FlowFeature {
public:
  BlockOfRandom(float _x,
                float _y,
                float _z,
                float _xsize,
                float _ysize,
                float _zsize,
                float _maxstr,
                int   _num)
    : FlowFeature(_x, _y, _z),
      m_xsize(_xsize),
      m_ysize(_ysize),
      m_zsize(_zsize),
      m_maxstr(_maxstr),
      m_num(_num)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

private:
  float m_xsize;
  float m_ysize;
  float m_zsize;
  float m_maxstr;
  int m_num;
};


//
// Concrete class for a particle emitter (one per frame)
//
class ParticleEmitter : public SingleParticle {
public:
  ParticleEmitter(float _x, float _y, float _z, float _sx, float _sy, float _sz)
    : SingleParticle(_x, _y, _z, _sx, _sy, _sz)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

private:
};


// vortex ring (singular)
// vortex ring (thick)
// vortex ring emitter (singular)

// uniformly-spaced particles

// particles from file

// panels: circle, rectangle, from file


//
// Concrete class for a solid circle
//
/*
class SolidCircle : public FlowFeature {
public:
  SolidCircle(float _x, float _y, float _diam)
    : FlowFeature(_x, _y),
      m_diam(_diam)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

protected:
  float m_diam;
};
*/

