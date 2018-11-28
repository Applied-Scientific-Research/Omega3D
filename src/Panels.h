#pragma once

#include "VectorHelper.h"
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
  Panels(const size_t _n, const elem_t _e, const move_t _m)
    : ElementBase<S>(_n, _e, _m) {

    // this initialization specific to Panels - a circle
    this->x[0].resize(_n);
    this->x[1].resize(_n);
    this->x[2].resize(_n);
    for (size_t i=0; i<_n; ++i) {
      this->x[0][i] = 0.5 * cos(2.0*i*M_PI/_n);
      this->x[1][i] = 0.5 * sin(2.0*i*M_PI/_n);
      this->x[2][i] = 0.5 * cos(2.0*i*M_PI/_n) * sin(2.0*i*M_PI/_n);
    }
    this->r.resize(_n);
    for (size_t i=0; i<_n; ++i) {
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
    //this->u.resize(2*_n);
    // velocity in base class
    //for (size_t d=0; d<Dimensions; ++d) {
    //  this->u[d].resize(_n);
    //}

    // and generate panel strengths
    if (_e != inert) {
      // need to assign it a vector first!
      std::array<Vector<S>,3> new_s;
      for (size_t d=0; d<3; ++d) {
        new_s[d].resize(_n);
        for (size_t i=0; i<_n; ++i) {
          new_s[d][i] = 2.0 * sin(2.0*(i+0.5)*M_PI/_n);
        }
      }
      this->s = std::move(new_s);
    }
  }

  const std::vector<uint16_t>& get_idx() const { return idx; }

  void add_new(std::vector<float>& _in) {
    // must explicitly call the method in the base class first?
    ElementBase<S>::add_new(_in);
  }

  //
  // OpenGL functions
  //

  // this gets done once - load the shaders, set up the vao
  void initGL(std::vector<float>& _projmat,
              float*              _poscolor,
              float*              _negcolor) {

    std::cout << "inside Panels.initGL" << std::endl;
  }

  // this gets done every time we change the size of the positions array
  void updateGL() {
    std::cout << "inside Panels.updateGL" << std::endl;
  }

  // stuff to display points, called once per frame
  void drawGL(std::vector<float>& _projmat,
              float*              _poscolor,
              float*              _negcolor) {

    std::cout << "inside Panels.drawGL" << std::endl;
  }

  std::string to_string() const {
    return ElementBase<S>::to_string() + " Panels";
  }

protected:
  // geometry
  // x, r, s are in ElementBase
  //std::variant<std::vector<uint16_t>, std::vector<uint32_t>> idx;
  std::vector<uint16_t> idx;
  // velocity
  //std::vector<S> u;
  // movement
  //std::optional<Body&> b;
  // curved panels need: normals
  //std::optional<std::vector<S>> norm;
};

