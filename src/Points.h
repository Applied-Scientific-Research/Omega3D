#pragma once

#include "ElementBase.h"

#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#define _USE_MATH_DEFINES
#include <cmath>


// 0-D elements
template <class S>
class Points: public ElementBase<S> {
public:
  Points(const size_t _n, const elem_t _e, const move_t _m)
    : ElementBase<S>(_n, _e, _m) {

    // this initialization specific to Points
    for (size_t d=0; d<Dimensions; ++d) {
      this->x[d].resize(_n);
      for (size_t i=0; i<_n; ++i) {
        this->x[d][i] = -1.0 + 2.0*(S)rand()/(S)RAND_MAX;
      }
    }
    this->r.resize(_n);
    this->elong.resize(_n);
    for (size_t i=0; i<_n; ++i) {
      this->r[i] = 0.01;
      this->elong[i] = 1.0;
    }

    // optional strength in base class
    if (_e != inert) {
      // need to assign it a vector first!
      std::vector<S> new_s;
      new_s.resize(3*_n);
      for (size_t i=0; i<3*_n; ++i) {
        new_s[i] = (-1.0 + 2.0*(S)rand()/(S)RAND_MAX) / (S)_n;
      }
      this->s = std::move(new_s);
    }

    // velocity in base class
    this->u.resize(3*_n);
    // optional velgrads here
    if (_m == lagrangian) {
      std::vector<S> new_ug;
      new_ug.resize(9*_n);
      ug = std::move(new_ug);
    }
  }

  std::optional<std::vector<S>>& get_velgrad() { return ug; }

  void zero_vels() {
    // must explicitly call the method in the base class
    ElementBase<S>::zero_vels();
    // and specialize
    if (ug) {
      for( S& val : *ug ) val = 0.0;
    }
  }
  void move(const double _dt) {
    // must explicitly call the method in the base class
    ElementBase<S>::move(_dt);
    // and specialize
    if (this->M == lagrangian and ug) {
      std::cout << "  Stretching" << to_string() << std::endl;

      // get pointers to the right part of the vectors
      std::vector<S> all_ug = *ug;
      std::vector<S> all_s = *this->s;

      for (size_t i=0; i<this->n; ++i) {
        std::array<S,3> wdu = {0.0};
        S* this_ug = &all_ug[9*i];
        S* this_s = &all_s[3*i];

        // compute stretch term
        // note that multiplying by the transpose may maintain linear impulse better, but
        //   severely underestimates stretch!
        wdu[0] = this_s[0]*this_ug[0] + this_s[1]*this_ug[3] + this_s[2]*this_ug[6];
        wdu[1] = this_s[0]*this_ug[1] + this_s[1]*this_ug[4] + this_s[2]*this_ug[7];
        wdu[2] = this_s[0]*this_ug[2] + this_s[1]*this_ug[5] + this_s[2]*this_ug[8];

        // update elongation
        const S circmag = std::sqrt(this_s[0]*this_s[0] + this_s[1]*this_s[1] + this_s[2]*this_s[2]);
        const S elongfactor = (S)_dt * (this_s[0]*wdu[0] + this_s[1]*wdu[1] + this_s[2]*wdu[2]) / circmag;
        elong[i] *= 1.0 + elongfactor;

        // add Cottet SFS into stretch term (after elongation)

        // update strengths
        this_s[0] += _dt * wdu[0];
        this_s[1] += _dt * wdu[1];
        this_s[2] += _dt * wdu[2];
      }
    }
  }

  std::string to_string() const {
    std::string retstr = ElementBase<S>::to_string() + " Points";
    if (ug) retstr += " with grads";
    return retstr;
  }

protected:
  //Geometry<0,0,S> x;    // num dimensions, num order, storage type
  // movement
  //std::optional<Body&> b;
  // state vector
  std::vector<S> elong;   // scalar elongation
  // time derivative of state vector
  std::optional<std::vector<S>> ug;   // velocity gradients
};

