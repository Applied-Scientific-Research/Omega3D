/*
 * MeasureFeature.h - GUI-side descriptions of flow measurement features
 *
 * (c)2018-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Feature.h"

#include "json/json.hpp"

#include <iostream>
#include <vector>

//
// Abstract class for any measurement feature (streamlines, rakes, tracers, etc.) present initially
//
class MeasureFeature : public Feature {
public:
  explicit
  MeasureFeature(float _x,
                 float _y,
                 float _z,
                 bool _moves)
    : Feature(true),
      m_x(_x),
      m_y(_y),
      m_z(_z),
      m_is_lagrangian(_moves)
    {}

  bool moves() const { return m_is_lagrangian; }
  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual void from_json(const nlohmann::json) = 0;
  virtual nlohmann::json to_json() const = 0;
  virtual std::vector<float> init_particles(float) const = 0;
  virtual std::vector<float> step_particles(float) const = 0;
#ifdef USE_IMGUI
  static void draw_creation_gui(std::vector<std::unique_ptr<MeasureFeature>> &, const float, const float &);
  virtual bool draw_info_gui(const std::string, const float, const float &) = 0;
#endif

protected:
  float m_x;
  float m_y;
  float m_z;
  bool  m_is_lagrangian;
};

std::ostream& operator<<(std::ostream& os, MeasureFeature const& ff);


//
// types of measurement features:
//
// single origin point, continuous tracer emitter
// single set of tracer particles
// fixed set of field points
// periodic rake tracer emitter
// grid of fixed field points
// solid block (square, circle) of tracers
// single streamline (save all positions of a single point, draw as a line)
//



//
// Concrete class for a single measurement point
//
class SinglePoint : public MeasureFeature {
public:
  SinglePoint(float _x = 0.0,
              float _y = 0.0,
              float _z = 0.0,
              bool _moves = true)
    : MeasureFeature(_x, _y, _z, _moves)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float, const float &) override;
#endif

protected:
  //float m_str;
};


//
// Concrete class for an immobile particle emitter (one per frame)
//
class TracerEmitter : public SinglePoint {
public:
  TracerEmitter(float _x = 0.0,
                float _y = 0.0,
                float _z = 0.0)
    : SinglePoint(_x, _y, _z, false)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float, const float &) override;
#endif

protected:
  // eventually implement frequency, but for now, once per step
  //float m_frequency;
};


//
// Concrete class for a circle of tracer points
//
class TracerBlob : public SinglePoint {
public:
  TracerBlob(float _x = 0.0,
             float _y = 0.0,
             float _z = 0.0,
             float _rad = 0.1)
    : SinglePoint(_x, _y, _z, true),
      m_rad(_rad)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float, const float &) override;
#endif

protected:
  float m_rad;
};


//
// Concrete class for a tracer line
//
class TracerLine : public SinglePoint {
public:
  TracerLine(float _x = 0.0,
             float _y = 0.0,
             float _z = 0.0,
             float _xf = 1.0,
             float _yf = 0.0,
             float _zf = 0.0)
    : SinglePoint(_x, _y, _z, true),
      m_xf(_xf),
      m_yf(_yf),
      m_zf(_zf)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float, const float &) override;
#endif

protected:
  float m_xf, m_yf, m_zf;
};


//
// Concrete class for a line of static measurement points
//
class MeasurementLine : public SinglePoint {
public:
  MeasurementLine(float _x = 0.0,
                  float _y = 0.0,
                  float _z = 0.0,
                  float _xf = 1.0,
                  float _yf = 0.0,
                  float _zf = 0.0)
    : SinglePoint(_x, _y, _z, false),
      m_xf(_xf),
      m_yf(_yf),
      m_zf(_zf)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float, const float &) override;
#endif

protected:
  float m_xf, m_yf, m_zf;
};


//
// Concrete class for a 2D grid of measurement points
//
class Grid2dPoints : public MeasureFeature {
public:
  Grid2dPoints(float _x = -1.0,
               float _y = -1.0,
               float _z =  0.0,
               float _xs = 2.0,
               float _ys = 0.0,
               float _zs = 0.0,
               float _xt = 0.0,
               float _yt = 2.0,
               float _zt = 0.0,
               float _ds = 0.1,
               float _dt = 0.1)
    : MeasureFeature(_x, _y, _z, false),
      m_xs(_xs),
      m_ys(_ys),
      m_zs(_zs),
      m_xt(_xt),
      m_yt(_yt),
      m_zt(_zt),
      m_ds(_ds),
      m_dt(_dt)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float, const float &) override;
#endif

protected:
  float m_xs, m_ys, m_zs;
  float m_xt, m_yt, m_zt;
  float m_ds, m_dt;
};

//
// Parser for converting json object to new feature
//
void parse_measure_json(std::vector<std::unique_ptr<MeasureFeature>>&, const nlohmann::json);

