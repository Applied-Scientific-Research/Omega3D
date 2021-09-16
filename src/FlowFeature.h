/*
 * FlowFeature.h - GUI-side descriptions of flow features
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#pragma once

#include "Body.h"
#include "Feature.h"
#include "ElementPacket.h"

#include "json/json.hpp"

#include <memory>
#include <iostream>
#include <vector>
#include <string>

//
// Abstract class for any flow feature present initially
//
class FlowFeature : public Feature {
public:
  explicit
  FlowFeature(float _x,
              float _y,
              float _z,
              std::shared_ptr<Body> _bp)
    : Feature(true),
      m_x(_x),
      m_y(_y),
      m_z(_z),
      m_bp(_bp)
    {}
  virtual ~FlowFeature() = default; 
  virtual FlowFeature* copy() const = 0;

  std::shared_ptr<Body> get_body() { return m_bp; }

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual void from_json(const nlohmann::json) = 0;
  virtual nlohmann::json to_json() const = 0;
  virtual ElementPacket<float> init_elements(float) const = 0;
  virtual ElementPacket<float> step_elements(float) const = 0;
  virtual void generate_draw_geom() = 0;
  virtual ElementPacket<float> get_draw_packet() { return m_draw; }

#ifdef USE_IMGUI
  virtual bool draw_info_gui(const std::string, const float) = 0;
  static int draw_creation_gui(std::vector<std::unique_ptr<FlowFeature>> &, const float);
#endif

protected:
  float m_x;
  float m_y;
  float m_z;
  ElementPacket<float> m_draw;
  std::shared_ptr<Body> m_bp;
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
  SingleParticle(float _x = 0.0,
                 float _y = 0.0,
                 float _z = 0.0,
                 float _sx = 0.0,
                 float _sy = 0.0,
                 float _sz = 1.0,
                 std::shared_ptr<Body> _bp = nullptr)
    : FlowFeature(_x, _y, _z, _bp),
      m_sx(_sx),
      m_sy(_sy),
      m_sz(_sz)
    {}
  SingleParticle* copy() const override 
                  { return new SingleParticle(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
#endif

protected:
  float m_sx, m_sy, m_sz;
};


//
// Concrete class for a vortex blob
//
class VortexBlob : public SingleParticle {
public:
  VortexBlob(float _x = 0.0,
             float _y = 0.0,
             float _z = 0.0,
             float _sx = 0.0,
             float _sy = 0.0,
             float _sz = 1.0,
             float _rad = 0.1,
             float _soft = 0.1,
             std::shared_ptr<Body> _bp = nullptr)
    : SingleParticle(_x, _y, _z, _sx, _sy, _sz, _bp),
      m_rad(_rad),
      m_softness(_soft)
    {}
  VortexBlob* copy() const override 
              { return new VortexBlob(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
#endif

protected:
  float m_rad;
  float m_softness;
};


//
// Concrete class for a rectangle of randomly-placed particles
//
class BlockOfRandom : public FlowFeature {
public:
  BlockOfRandom(float _x = 0.0,
                float _y = 0.0,
                float _z = 0.0,
                float _xsize = 1.0,
                float _ysize = 1.0,
                float _zsize = 1.0,
                float _maxstr = 0.01,
                int   _num = 100,
                std::shared_ptr<Body> _bp = nullptr)
    : FlowFeature(_x, _y, _z, _bp),
      m_xsize(_xsize),
      m_ysize(_ysize),
      m_zsize(_zsize),
      m_maxstr(_maxstr),
      m_num(_num)
    {}
  BlockOfRandom* copy() const override 
                 { return new BlockOfRandom(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
#endif

protected:
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
  ParticleEmitter(float _x = 0.0,
                  float _y = 0.0,
                  float _z = 0.0,
                  float _sx = 0.0,
                  float _sy = 0.0,
                  float _sz = 0.01,
                  std::shared_ptr<Body> _bp = nullptr)
    : SingleParticle(_x, _y, _z, _sx, _sy, _sz, _bp)
    {}
  ParticleEmitter* copy() const override 
                   { return new ParticleEmitter(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
#endif

protected:
};


//
// Concrete class for a singular vortex ring
//
class SingularRing : public FlowFeature {
public:
  SingularRing(float _x = 0.0,
               float _y = 0.0,
               float _z = 0.0,
               float _nx = 1.0,
               float _ny = 0.0,
               float _nz = 0.0,
               float _majrad = 0.5,
               float _circ = 1.0,
               std::shared_ptr<Body> _bp = nullptr)
    : FlowFeature(_x, _y, _z, _bp),
      m_nx(_nx),
      m_ny(_ny),
      m_nz(_nz),
      m_majrad(_majrad),
      m_circ(_circ)
    {}
  SingularRing* copy() const override 
                { return new SingularRing(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
#endif

protected:
  float m_nx;
  float m_ny;
  float m_nz;
  float m_majrad;
  float m_circ;
};


//
// Concrete class for a thick vortex ring
//
class ThickRing : public SingularRing {
public:
  ThickRing(float _x = 0.0,
            float _y = 0.0,
            float _z = 0.0,
            float _nx = 1.0,
            float _ny = 0.0,
            float _nz = 0.0,
            float _majrad = 0.5,
            float _minrad = 0.05,
            float _circ = 1.0,
            std::shared_ptr<Body> _bp = nullptr)
    : SingularRing(_x, _y, _z, _nx, _ny, _nz, _majrad, _circ),
      m_minrad(_minrad)
    {}
  ThickRing* copy() const override 
             { return new ThickRing(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
#endif

protected:
  float m_minrad;
};

// uniformly-spaced particles

// how about an oval ring? requires no radii, but two basis vectors: long axis and short axis
// vortex ring emitter (singular)

// read particles from file

//
// Parser for converting json object to new feature
//
void parse_flow_json(std::vector<std::unique_ptr<FlowFeature>>&, const nlohmann::json);

