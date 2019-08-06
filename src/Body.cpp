/*
 * Body.cpp - class for an independent solid boundary
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "Body.h"

#include <cassert>
#include <iostream>

// use this differential to evaluate velocities (in the absence of a differentiator function)
const double DT = 1.e-6;

//
// A single rigid body
//
// Note that Vec = std::array<double,Dimensions>
//

// delegating ctor
Body::Body() :
  Body(0.0, 0.0, 0.0)
  {}

// primary constructor
Body::Body(const double _x, const double _y, const double _z) :
  this_time(0.0),
  pos(Vec({{_x, _y, _z}})),
  vel(Vec({{0.0, 0.0, 0.0}})),
  apos(Vec({{0.0, 0.0, 0.0}})),
  avel(Vec({{0.0, 0.0, 0.0}})),
  vol(0.0)
{
  // time (t) is the only variable allowed in the equations
  func_vars.push_back({"t", &this_time});
  // quaternion needs special initialization
  qpos.setIdentity();
  // ane make space for the compiled functions
  pos_func.resize(Dimensions);
  apos_func.resize(Dimensions);
}

Body::~Body() {
  // free the internal memory used by tinyexpr
  for (size_t i=0; i<pos_func.size(); ++i) te_free(pos_func[i]);
  for (size_t i=0; i<apos_func.size(); ++i) te_free(apos_func[i]);
}


// create and write a json object to which to add geometries
nlohmann::json
Body::to_json() const {
  nlohmann::json j;

  if (not name.empty()) j["name"] = name;
  if (not parent.empty()) j["parent"] = parent;

  // translation has to be an array
  nlohmann::json jpos = nlohmann::json::array();
  for (size_t i=0; i<Dimensions; ++i) {
    if (pos_func[i]) {
      jpos.push_back(pos_expr[i]);
    } else {
      jpos.push_back(pos[i]);
    }
  }
  j["translation"] = jpos;

  // HACK - remember to tell if the rotation is "euler", "vector" (scaled axis-angle), or "quat" !
  // here we assume "vector"
  nlohmann::json jrot = nlohmann::json::array();
  for (size_t i=0; i<Dimensions; ++i) {
    if (apos_func[i]) {
      jrot.push_back(apos_expr[i]);
    } else {
      jrot.push_back(apos[i]);
    }
  }
  j["rotation"] = jrot;

  return j;
}


// getters/setters

void Body::set_name(const std::string _name) { name = _name; }
void Body::set_parent_name(const std::string _name) { parent = _name; }
std::string Body::get_name() { return name; }

void Body::set_pos(const size_t _i, const double _val) {
  assert(_i>=0 and _i<Dimensions && "Invalid index into array");
  pos[_i] = _val;
}

void Body::set_pos(const size_t _i, const std::string _val) {
  assert(_i>=0 and _i<Dimensions && "Invalid index into array");
  // store the expression locally
  pos_expr[_i] = _val;
  // compile it
  int ierr = 0;
  pos_func[_i] = te_compile(_val.c_str(), func_vars.data(), 1, &ierr);
  if (pos_func[_i]) {
    std::cout << "  read expression (" << pos_expr[_i] << ")" << std::endl;
    this_time = 0.0;
    std::cout << "  testing parsed expression, with t=0, value is " << te_eval(pos_func[_i]) << std::endl;
    this_time = 1.0;
    std::cout << "                                  t=1, value is " << te_eval(pos_func[_i]) << std::endl;
    //this_time = 2.0;
    //std::cout << "                                  t=2, value is " << te_eval(pos_func[_i]) << std::endl;
  } else {
    std::cout << "  Error parsing expression (" << _val << "), near character " << ierr << std::endl;
  }
}

void Body::set_rot(const size_t _i, const double _val) {
  assert(_i>=0 and _i<Dimensions && "Invalid index into array");
  apos[_i] = _val;
}

void Body::set_rot(const size_t _i, const std::string _val) {
  assert(_i>=0 and _i<Dimensions && "Invalid index into array");
  // store the expression locally
  apos_expr[_i] = _val;
  // compile it
  int ierr = 0;
  apos_func[_i] = te_compile(_val.c_str(), func_vars.data(), 1, &ierr);
  if (apos_func[_i]) {
    std::cout << "  read expression (" << apos_expr[_i] << ")" << std::endl;
    this_time = 0.0;
    std::cout << "  testing parsed expression, with t=0, value is " << te_eval(apos_func[_i]) << std::endl;
    this_time = 1.0;
    std::cout << "                                  t=1, value is " << te_eval(apos_func[_i]) << std::endl;
    //this_time = 2.0;
    //std::cout << "                                  t=2, value is " << te_eval(apos_func) << std::endl;
  } else {
    std::cout << "  Error parsing expression (" << _val << "), near character " << ierr << std::endl;
  }
}

Vec Body::get_pos() {
  return pos;
}
Vec Body::get_pos(const double _time) {
  Vec newpos;

  // get the value or evaluate the expression
  this_time = _time;
  //std::cout << "  MOVING BODY (" << get_name() << ") at time " << _time << std::endl;
  for (size_t i=0; i<Dimensions; ++i) {
    if (pos_func[i]) {
      newpos[i] = te_eval(pos_func[i]);
      //std::cout << "IDEAL POS " << (0.5 * (1.0-cos(2.0*_time))) << "  AND ACTUAL " << pos[i] << std::endl;
      //std::cout << "  MOVED BODY pos[" << i << "] to " << pos[i] << std::endl;
    } else {
      newpos[i] = pos[i];
    }
  }

  return newpos;
}

Vec Body::get_vel() {
  return vel;
}
Vec Body::get_vel(const double _time) {
  Vec newvel;

  // use 2-point first derivative estimate
  for (size_t i=0; i<Dimensions; ++i) {
    if (pos_func[i]) {
      this_time = _time + DT;
      const double pplus = te_eval(pos_func[i]);
      this_time = _time - DT;
      const double pminus = te_eval(pos_func[i]);
      newvel[i] = (pplus - pminus) / (2.0*DT);
      //std::cout << "      VEL " << (0.5 * 2.0*sin(2.0*_time)) << "  AND ACTUAL " << vel[i] << std::endl;
      //std::cout << "        used " << pplus << " and " << pminus << std::endl;
    } else {
      newvel[i] = vel[i];
    }
  }

  return newvel;
}
  

// three ways to get the current orientation - assume it's already set
Vec Body::get_orient_vec() {
  return apos;
}
Eigen::AngleAxis<double> Body::get_orient_aa() {
  return Eigen::AngleAxis<double>(qpos);
}
Eigen::Quaternion<double> Body::get_orient_quat() {
  return qpos;
}

// get the orientation at a given time
Vec Body::get_orient_vec(const double _time) {
  Vec orient;
  for (size_t i=0; i<Dimensions; ++i) {
    if (apos_func[i]) {
      this_time = _time;
      orient[i] = te_eval(apos_func[i]);
    } else {
      orient[i] = apos[i];
    }
  }
  return orient;
}
Eigen::AngleAxis<double> Body::get_orient_aa(const double _time) {
  const Vec orient = get_orient_vec(_time);
  Eigen::Vector3d axis(orient[0], orient[1], orient[2]);
  const double angle_in_radians = axis.norm();
  if (std::abs(angle_in_radians) < std::numeric_limits<double>::epsilon()) {
    // no rotation, but force a dummy axis
    axis(0) = 1.0;
  } else {
    // scale by non-zero angle
    axis *= 1.0/angle_in_radians;
  }
  return Eigen::AngleAxis<double>(angle_in_radians, axis);
}
Eigen::Quaternion<double> Body::get_orient_quat(const double _time) {
  const Eigen::AngleAxis<double> aa = get_orient_aa(_time);
  //std::cout << "AngleAxis is now" << std::endl;
  //std::cout << aa.toRotationMatrix() << std::endl;
  return Eigen::Quaternion<double>(aa);
}


// two ways to get current rotational velocity
Vec Body::get_rotvel_vec() {
  return avel;
}
Eigen::AngleAxis<double> Body::get_rotvel_aa() {
  // must convert from the stored Vec representation
  Eigen::Vector3d axis(avel[0], avel[1], avel[2]);
  const double angle_in_radians = axis.norm();
  axis *= 1.0/angle_in_radians;
  return Eigen::AngleAxis<double>(angle_in_radians, axis);
}
Vec Body::get_rotvel_vec(const double _time) {
  Vec newavel;

  // use 2-point first derivative estimate
  for (size_t i=0; i<Dimensions; ++i) {
    if (apos_func[i]) {
      this_time = _time + DT;
      const double pplus = te_eval(apos_func[i]);
      this_time = _time - DT;
      const double pminus = te_eval(apos_func[i]);
      newavel[i] = (pplus - pminus) / (2.0*DT);
      //std::cout << "      VEL " << (0.5 * 2.0*sin(2.0*_time)) << "  AND ACTUAL " << vel[i] << std::endl;
      //std::cout << "        used " << pplus << " and " << pminus << std::endl;
    } else {
      newavel[i] = avel[i];
    }
  }

  return newavel;
}
Eigen::AngleAxis<double> Body::get_rotvel_aa(const double _time) {
  const Vec rotvel = get_rotvel_vec(_time);
  Eigen::Vector3d axis(rotvel[0], rotvel[1], rotvel[2]);
  const double angle_in_radians = axis.norm();
  axis *= 1.0/angle_in_radians;
  return Eigen::AngleAxis<double>(angle_in_radians, axis);
}


// using the above getters, set this object's transform to the given time
void Body::transform(const double _time) {
  // first do positions and velocities
  pos = get_pos(_time);
  vel = get_vel(_time);

  // rotation has to work differently, as we need to construct quaternions
  apos = get_orient_vec(_time);
  qpos = get_orient_quat(_time);
  avel = get_rotvel_vec(_time);

  //std::cout << "Quaternion is now" << std::endl;
  //std::cout << qpos.toRotationMatrix() << std::endl;
}

// get an Eigen-like Transform
Trans Body::get_transform_mat() {
  // this is a rotation followed by a movement (rotation is farthest to the right in the equation)
  Trans t;
  t.setIdentity();
  t *= Eigen::Translation<double,3>(pos[0], pos[1], pos[2]);
  t *= Eigen::AngleAxisd(qpos);
  return t;
}

Trans Body::get_transform_mat(const double _time) {
  Trans t;
  t.setIdentity();
  const Vec thispos = get_pos(_time);
  t *= Eigen::Translation<double,3>(thispos[0], thispos[1], thispos[2]);
  const Eigen::AngleAxisd thisorient = get_orient_aa(_time);
  t *= thisorient;
  return t;
}


// compare motion vs another Body
bool Body::relative_motion_vs(std::shared_ptr<Body> _other, const double _last, const double _current) {

  // transformation from one body to the other at time _last
  const Trans oldtrans = get_transform_mat(_last) * _other->get_transform_mat(_last).inverse();

  // transformation from one body to the other at time _current
  const Trans newtrans = get_transform_mat(_current) * _other->get_transform_mat(_current).inverse();

  // are these similar?
  return not oldtrans.isApprox(newtrans);
}

