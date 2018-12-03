/*
 * ElementBase.h - abstract class for arrays of any computational elements
 *
 * (c)2018 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "VectorHelper.h"
#include "Omega3D.h"

#include <iostream>
#include <vector>
#include <memory>
#include <cassert>
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
  const std::array<Vector<S>,Dimensions>& get_pos() const { return x; }
  const Vector<S>&                        get_rad() const { return r; }
  const std::array<Vector<S>,Dimensions>& get_str() const { return *s; }
  std::array<Vector<S>,Dimensions>&       get_vel()       { return u; }

  void add_new(std::vector<float>& _in) {

    // check inputs
    if (_in.size() == 0) return;
    assert(_in.size() % 7 == 0);
    const size_t nnew = _in.size()/7;

    // this initialization is specific to Points - so should we do it there?
    for (size_t d=0; d<Dimensions; ++d) {
      // copy existing vector
      //Vector<S> new_x = x[d];
      // extend with more space for new values
      //new_x.resize(n+nnew);
      x[d].resize(n+nnew);
      // copy new values to end of vector
      for (size_t i=0; i<nnew; ++i) {
        //new_x[n+i] = _in[7*i+d];
        x[d][n+i] = _in[7*i+d];
      }
      // finally, replace the master with this new vector
      //x[d] = std::move(new_x);
    }

    // do radius now
    {
      //Vector<S> new_r = r;
      //new_r.resize(n+nnew);
      r.resize(n+nnew);
      for (size_t i=0; i<nnew; ++i) {
        //new_r[n+i] = _in[7*i+6];
        r[n+i] = _in[7*i+6];
      }
      //r = std::move(new_r);
    }

    // strength
    if (s) {
      for (size_t d=0; d<Dimensions; ++d) {
        (*s)[d].resize(n+nnew);
        for (size_t i=0; i<nnew; ++i) {
          (*s)[d][n+i] = _in[7*i+d+3];
        }
      }
    }

    // extend the other vectors as well
    for (size_t d=0; d<Dimensions; ++d) {
      u[d].resize(n+nnew);
    }
    if (dsdt) {
      for (size_t d=0; d<Dimensions; ++d) {
        (*dsdt)[d].resize(n+nnew);
      }
    }

    // finally, update n
    n += nnew;
  }

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

      // update strengths (in derived class)
    }
  }
  void move(const double _dt,
            const double _wt1, ElementBase<S> const & _u1,
            const double _wt2, ElementBase<S> const & _u2) {
    // must confirm that incoming time derivates include velocity
    // if this has vels, then lets advect it
    if (M == lagrangian) {
      std::cout << "  Moving" << to_string() << std::endl;

      // update positions
      for (size_t d=0; d<Dimensions; ++d) {
        for (size_t i=0; i<n; ++i) {
          x[d][i] += (S)_dt * (_wt1*_u1.u[d][i] + _wt2*_u2.u[d][i]);
        }
      }

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
  std::array<Vector<S>,Dimensions> x;                   // position
  Vector<S> r;                                          // thickness/radius
  std::optional<std::array<Vector<S>,Dimensions>> s;    // strength

  // time derivative of state vector
  std::array<Vector<S>,Dimensions> u;                   // velocity
  //Vector<S> dr;                                       // thickness/radius
  std::optional<std::array<Vector<S>,Dimensions>> dsdt; // strength change
};

