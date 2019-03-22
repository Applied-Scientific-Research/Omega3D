/*
 * BoundaryFeature.h - GUI-side descriptions of boundary features
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega3D.h"
#include "Body.h"

#include "json/json.hpp"

#include <iostream>
#include <vector>
#include <memory>

//
// Abstract class for any boundary feature present initially
//
class BoundaryFeature {
public:
  explicit
  BoundaryFeature(std::shared_ptr<Body> _bp,
                  float _x,
                  float _y,
                  float _z)
    : m_bp(_bp),
      m_x(_x),
      m_y(_y),
      m_z(_z)
    {}
  virtual ~BoundaryFeature() {}

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual void from_json(const nlohmann::json) = 0;
  virtual nlohmann::json to_json() const = 0;
  virtual ElementPacket<float> init_elements(const float) const = 0;
  //virtual std::vector<float> step_elements(const float) const = 0;
  std::shared_ptr<Body> get_body() { return m_bp; }

protected:
  std::shared_ptr<Body> m_bp;
  float m_x;
  float m_y;
  float m_z;
};

std::ostream& operator<<(std::ostream& os, BoundaryFeature const& ff);


//
// Concrete class for geometry from a file (fluid is outside)
//
class ExteriorFromFile : public BoundaryFeature {
public:
  ExteriorFromFile(std::shared_ptr<Body> _bp = nullptr,
                   float _x = 0.0,
                   float _y = 0.0,
                   float _z = 0.0,
                   float _sx = 1.0,
                   float _sy = 1.0,
                   float _sz = 1.0,
                   std::string _infile = "")
    : BoundaryFeature(_bp, _x, _y, _z),
      m_sx(_sx),
      m_sy(_sy),
      m_sz(_sz),
      m_infile(_infile)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(const float) const override;

protected:
  float m_sx;
  float m_sy;
  float m_sz;
  std::string m_infile;
};


//
// Parser for converting json object to new feature
//
void parse_boundary_json(std::vector<std::unique_ptr<BoundaryFeature>>&,
                         std::shared_ptr<Body>,
                         const nlohmann::json);

