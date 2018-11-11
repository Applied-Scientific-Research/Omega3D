#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#include <variant>
#define _USE_MATH_DEFINES
#include <cmath>


// the superclass

template <class S>
class ElementBase {
public:
  ElementBase<S>(const size_t _n, const elem_t _e, const move_t _m) :
      E(_e), M(_m), n(_n) {
  }

  size_t getn() const { return n; }
  const std::array<std::vector<S>,Dimensions>& get_pos() const { return x; }
  const std::vector<S>&                        get_rad() const { return r; }
  const std::vector<S>&                        get_str() const { return *s; }
  std::array<std::vector<S>,Dimensions>&       get_vel()       { return u; }

  void zero_vels() {
    for (size_t d=0; d<Dimensions; ++d) {
      for (size_t i=0; i<getn(); ++i) {
        u[d][i] = 0.0;
      }
    }
  }
  void finalize_vels(const std::array<double,Dimensions>& _fs) {
    for (size_t d=0; d<Dimensions; ++d) {
      for (size_t i=0; i<getn(); ++i) {
        u[d][i] = _fs[d] + u[d][i] * 0.25/M_PI;
      }
    }
  }
  void move(const double _dt) {
    if (M == lagrangian) {
      std::cout << "  Moving" << to_string() << std::endl;

      // update positions
      for (size_t d=0; d<Dimensions; ++d) {
        for (size_t i=0; i<n; ++i) {
          x[d][i] += (S)_dt * u[d][i];
        }
      }

      // update elongation
      // update strengths (in derived class)
    }
  }

  std::string to_string() const {
    std::string mystr = " " + std::to_string(n);
    if (E == active) {
      mystr += " Active";
    } else if (E == reactive) {
      mystr += " Reactive";
    } else {
      mystr += " Inert";
    }
    if (M == lagrangian) {
      mystr += " Lagrangian";
    } else if (M == bodybound) {
      mystr += " Body-fixed";
    } else {
      mystr += " Fixed";
    }
    return mystr;
  }

protected:
  // active, reactive, or inert?
  elem_t E;
  // how does it move? use move_t or Body*
  move_t M;
  //Body* b = nullptr;
  //Move_t get_move_type() {
  //}

  // common arrays for all derived types
  size_t n;

  // state vector
  std::array<std::vector<S>,Dimensions> x;   // position
  std::vector<S> r;                          // thickness/radius
  std::optional<std::vector<S>> s;           // strength

  // time derivative of state vector
  std::array<std::vector<S>,Dimensions> u;   // velocity
  std::optional<std::vector<S>> dsdt;        // strength change
};

