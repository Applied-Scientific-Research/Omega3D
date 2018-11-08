#pragma once

#include "ElementBase.h"

#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#include <variant>
#define _USE_MATH_DEFINES
#include <cmath>


// 1-D elements
template <class S>
class Panels: public ElementBase<S> {
public:
  Panels(const size_t _n, const elem_t _e, const move_t _m) :
      ElementBase<S>(_n, _e, _m) {
    // this initialization specific to Panels - a circle
    this->x.resize(2*_n);
    this->r.resize(_n);
    for (size_t i=0; i<_n; ++i) {
      this->x[2*i+0] = 0.5 * cos(2.0*i*M_PI/_n);
      this->x[2*i+1] = 0.5 * sin(2.0*i*M_PI/_n);
      this->r[i] = 0.01;
    }
    // initialize indices to nodes
    idx.resize(2*_n);
    for (size_t i=0; i<_n; ++i) {
      idx[2*i] = i;
      idx[2*i+1] = i+1;
    }
    idx[2*_n-1] = 0;
    // just size vels
    this->u.resize(2*_n);
    // and generate panel strengths
    if (_e != inert) {
      // need to assign it a vector first!
      std::vector<S> new_s;
      new_s.resize(_n);
      for (size_t i=0; i<_n; ++i) {
        new_s[i] = 2.0 * sin(2.0*(i+0.5)*M_PI/_n);
      }
      this->s = std::move(new_s);
    }
  }

  const std::vector<uint16_t>& get_idx() const { return idx; }

  std::string to_string() const {
    return ElementBase<S>::to_string() + " Panels";
  }

protected:
  // geometry
  //std::vector<S> x;
  //std::vector<S> r;
  //std::variant<std::vector<uint16_t>, std::vector<uint32_t>> idx;
  std::vector<uint16_t> idx;
  // velocity
  //std::vector<S> u;
  // movement
  //std::optional<Body&> b;
  // curved panels need: normals
  std::optional<std::vector<S>> norm;
};

