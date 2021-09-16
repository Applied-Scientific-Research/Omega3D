/*
 * MeasureFeature.h - GUI-side descriptions of flow measurement features
 *
 * (c)2018-9 Applied Scientific Research, Inc.
 *           Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Body.h"
#include "Omega3D.h"
#include "Feature.h"
#include "ElementPacket.h"

#include "json/json.hpp"

#include <iostream>
#include <vector>
#include <string>

//
// Abstract class for any measurement feature (streamlines, rakes, tracers, etc.) present initially
//
class MeasureFeature : public Feature {
public:
  explicit
  MeasureFeature(float _x,
                 float _y,
                 float _z,
                 bool _moves,
                 bool _emits,
                 std::shared_ptr<Body> _bp)
    : Feature(true),
      m_x(_x),
      m_y(_y),
      m_z(_z),
      m_is_lagrangian(_moves),
      m_emits(_emits),
      m_bp(_bp)
    {}
  virtual ~MeasureFeature() {}
  virtual MeasureFeature* copy() const = 0;

  bool moves() const { return m_is_lagrangian; }
  bool emits() const { return m_emits; }
  float jitter(const float, const float) const;
  ElementPacket<float> get_draw_packet() const { return m_draw; }
  bool get_is_lagrangian() { return m_is_lagrangian; }
  std::shared_ptr<Body> get_body() { return m_bp; }

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual void from_json(const nlohmann::json) = 0;
  virtual nlohmann::json to_json() const = 0;
  virtual ElementPacket<float> init_elements(float) const = 0;
  virtual ElementPacket<float> step_elements(float) const = 0;
  virtual void generate_draw_geom() = 0;
#ifdef USE_IMGUI
  static int draw_creation_gui(std::vector<std::unique_ptr<MeasureFeature>> &, const float, const float &);
  virtual bool draw_info_gui(const std::string, const float &, const float) = 0;
#endif

protected:
  float m_x;
  float m_y;
  float m_z;
  bool  m_is_lagrangian;
  bool m_emits;
  std::shared_ptr<Body> m_bp;
  ElementPacket<float> m_draw;
};

std::ostream& operator<<(std::ostream& os, MeasureFeature const& ff);

//
// types of measurement features:
//
// single origin point, continuous tracer emitter
// fixed set of field points
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
              bool _moves = false,
              bool _emits = false,
              std::shared_ptr<Body> _bp = nullptr)
    : MeasureFeature(_x, _y, _z, _moves, _emits, _bp)
    {}
  SinglePoint* copy() const override { return new SinglePoint(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float&, const float) override;
#endif

protected:
  //float m_str;
};

//
// Concrete class for a circle of tracer points
//
class MeasurementBlob : public SinglePoint {
public:
  MeasurementBlob(float _x = 0.0,
                  float _y = 0.0,
                  float _z = 0.0,
                  bool _moves = false,
                  bool _emits = false,
                  float _rad = 0.1,
                  std::shared_ptr<Body> _bp = nullptr)
    : SinglePoint(_x, _y, _z, _moves, _emits, _bp),
      m_rad(_rad)
    {}
  MeasurementBlob* copy() const override { return new MeasurementBlob(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float&, const float) override;
#endif

protected:
  float m_rad;
};

//
// Concrete class for a line of measurement points
//
class MeasurementLine : public SinglePoint {
public:
  MeasurementLine(float _x = 0.0,
                  float _y = 0.0,
                  float _z = 0.0,
                  bool _moves = false,
                  bool _emits = false,
                  float _xf = 1.0,
                  float _yf = 0.0,
                  float _zf = 0.0,
                  float _dx = 0.1,
                  std::shared_ptr<Body> _bp = nullptr)
    : SinglePoint(_x, _y, _z, _moves, _emits, _bp),
      m_xf(_xf),
      m_yf(_yf),
      m_zf(_zf),
      m_dx(_dx)
    {}
  MeasurementLine* copy() const override { return new MeasurementLine(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float&, const float) override;
#endif

protected:
  float m_xf, m_yf, m_zf;
  float m_dx;
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
               float _xf = 0.0,
               float _yf = 2.0,
               float _zf = 0.0,
               float _ds = 0.1,
               float _df = 0.1)
    : MeasureFeature(_x, _y, _z, false, false, nullptr),
      m_xs(_xs),
      m_ys(_ys),
      m_zs(_zs),
      m_xf(_xf),
      m_yf(_yf),
      m_zf(_zf),
      m_ds(_ds),
      m_df(_df)
    {}
  Grid2dPoints* copy() const override { return new Grid2dPoints(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float&, const float) override;
#endif

protected:
  float m_xs, m_ys, m_zs;
  float m_xf, m_yf, m_zf;
  float m_ds, m_df;
};

//
// Parser for converting json object to new feature
//
void parse_measure_json(std::vector<std::unique_ptr<MeasureFeature>>&, const nlohmann::json);

