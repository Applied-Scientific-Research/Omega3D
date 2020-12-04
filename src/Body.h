/*
 * Body.h - class for an independent solid boundary
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega3D.h"

#define TE_NAT_LOG
#include <tinyexpr/tinyexpr.h>
#include <json/json.hpp>

#include <Eigen/Geometry>

#include <string>
#include <memory>
#include <cmath>
#include <array>

using Vec = std::array<double,Dimensions>;
using Trans = Eigen::Transform<double,Dimensions,Eigen::Affine>;

//------------------------------------------------------------------------
//
// A single rigid body
//
// Holds and reports only motions
//
// how would we support inside-out bodies? say we're simulating flow inside a circle or box?
//   those collections would be attached to the ground body - the default 0th body
//
// for use in 2D and 3D codes, should really use a custom Rotation object instead of a double
//

class Body {
public:
  Body();
  Body(const double, const double, const double);
  // destructor needs to call te_free on all te_expr pointers
  ~Body();// = default;

  // dump a json object for writing
  nlohmann::json to_json() const;

  // setters, as we may not construct the class at once
  void set_name(const std::string);
  void set_parent_name(const std::string);
  void set_pos(const size_t, const double);
  void set_pos(const size_t, const std::string);
  void set_rot(const size_t, const double);
  void set_rot(const size_t, const std::string);

  // getters
  std::string get_name();

  Vec get_pos();
  Vec get_pos(const double);
  Vec get_vel();
  Vec get_vel(const double);

  // six ways to get the orientation (3 formats, 2 times)
  Vec get_orient_vec();
  Eigen::AngleAxis<double> get_orient_aa();
  Eigen::Quaternion<double> get_orient_quat();
  Vec get_orient_vec(const double);
  Eigen::AngleAxis<double> get_orient_aa(const double);
  Eigen::Quaternion<double> get_orient_quat(const double);

  Vec get_rotvel_vec();
  Eigen::AngleAxis<double> get_rotvel_aa();
  Vec get_rotvel_vec(const double);
  Eigen::AngleAxis<double> get_rotvel_aa(const double);

  // set and get the transform for a given time
  void transform(const double);
  Trans get_transform_mat();
  Trans get_transform_mat(const double);

  // compare motion vs another Body
  bool relative_motion_vs(std::shared_ptr<Body>, const double, const double);

private:
  // a name to refer to this body and echo when asked
  std::string name;

  // string containing expression to be parsed when needed
  std::array<std::string,Dimensions> pos_expr;
  std::array<std::string,Dimensions> apos_expr;
  // why not std::variant<double, std::string> for these?

  // needed by tinyexpr
  double this_time;
  std::vector<te_variable> func_vars;
  std::vector<te_expr*> pos_func;
  std::vector<te_expr*> apos_func;

  // 3D position and velocity (initial, or constant)
  Vec pos;
  Vec vel;

  // angular position and velocity (initial, or constant)
  // first as scaled axis-angle, then as quaternion
  Vec apos;
  Eigen::Quaternion<double> qpos;
  Vec avel;

  // enclosed volume (needed for total circulation of rotating body)
  double vol;

  // name of parent
  std::string parent;
  // pointer to parent
  std::shared_ptr<Body> pp;
};

