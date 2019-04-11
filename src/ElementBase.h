/*
 * ElementBase.h - abstract class for arrays of any computational elements
 *
 * (c)2018-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "VectorHelper.h"
#include "Omega3D.h"
#include "Body.h"

#include <iostream>
#include <vector>
#include <memory>
#include <cassert>
#include <optional>
#include <variant>
#include <cmath>


// the superclass

template <class S>
class ElementBase {
public:
  ElementBase<S>(const size_t _n,
                 const elem_t _e,
                 const move_t _m,
                 std::shared_ptr<Body> _bp) :
      E(_e), M(_m), B(_bp), n(_n) {
  }

  size_t get_n() const { return n; }
  bool is_inert() const { return E==inert; }
  elem_t                                  get_elemt() const { return E; }
  move_t                                  get_movet() const { return M; }
  const std::shared_ptr<Body>             get_body_ptr() const { return B; }
  std::shared_ptr<Body>                   get_body_ptr()  { return B; }
  const std::array<Vector<S>,Dimensions>& get_pos() const { return x; }
  std::array<Vector<S>,Dimensions>&       get_pos()       { return x; }
  const std::array<Vector<S>,Dimensions>& get_str() const { return *s; }
  std::array<Vector<S>,Dimensions>&       get_str()       { return *s; }
  const std::array<Vector<S>,Dimensions>& get_vel() const { return u; }
  std::array<Vector<S>,Dimensions>&       get_vel()       { return u; }

  void set_str(const size_t ioffset, const size_t icnt, Vector<S> _in) {
    // only Points use s, and only Surfaces are supported in BEM
    assert(false && "Should not be in here");

    //assert(s && "Strength array does not exist");
    //assert(_in.size() == (*s).size() && "Set strength array size does not match");
    //assert(ioffset == 0 && "Offset is not zero");

    // copy over the strengths
    //*s = _in;
  }

  void add_new(std::vector<float>& _in) {

    // check inputs
    if (_in.size() == 0) return;
    const size_t nper = (this->E == inert) ? 3 : 7;
    assert(_in.size() % nper == 0 && "Input vector not a multiple of 3 or 7");
    const size_t nnew = _in.size()/nper;

    // this initialization is specific to Points - so should we do it there?
    for (size_t d=0; d<Dimensions; ++d) {
      // extend with more space for new values
      x[d].resize(n+nnew);
      // copy new values to end of vector
      for (size_t i=0; i<nnew; ++i) {
        x[d][n+i] = _in[nper*i+d];
      }
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
    //if (dsdt) {
    //  for (size_t d=0; d<Dimensions; ++d) {
    //    (*dsdt)[d].resize(n+nnew);
    //  }
    //}

    // finally, update n
    n += nnew;
  }

  // up-size all arrays to the new size, filling with sane values
  // this only happens right after diffusion
  void resize(const size_t _nnew) {
    const size_t currn = n;
    //std::cout << "  inside ElementBase::resize with " << currn << " " << _nnew << std::endl;
    if (_nnew == currn) return;

    // positions first
    for (size_t d=0; d<Dimensions; ++d) {
      const size_t thisn = x[d].size();
      x[d].resize(_nnew);
      for (size_t i=thisn; i<_nnew; ++i) {
        x[d][i] = 0.0;
      }
    }

    // strength
    if (s) {
      for (size_t d=0; d<Dimensions; ++d) {
        const size_t thisn = (*s)[d].size();
        (*s)[d].resize(_nnew);
        for (size_t i=thisn; i<_nnew; ++i) {
          (*s)[d][i] = 0.0;
        }
      }
    }

    // and finally velocity (no need to set it)
    for (size_t d=0; d<Dimensions; ++d) {
      u[d].resize(_nnew);
    }

    // lastly, update n
    n = _nnew;
  }

  void zero_vels() {
    for (size_t d=0; d<Dimensions; ++d) {
      std::fill(u[d].begin(), u[d].end(), 0.0);
    }
  }

  void finalize_vels(const std::array<double,Dimensions>& _fs) {
    const double factor = 0.25/M_PI;
    for (size_t d=0; d<Dimensions; ++d) {
      for (size_t i=0; i<get_n(); ++i) {
        u[d][i] = _fs[d] + u[d][i] * factor;
      }
    }
  }

  void zero_strengths() {
    if (s) {
      for (size_t d=0; d<Dimensions; ++d) {
        std::fill((*s)[d].begin(), (*s)[d].end(), 0.0);
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

    } else if (B and M == bodybound) {
      //transform(_time+_dt);
    }
  }

  // find the new peak strength magnitude
  S get_max_str() {
    if (s) {
      // we have strengths, go through and check them
      S thismax = 0.0;
      for (size_t i=0; i<(*s)[0].size(); ++i) {
        S thisstr = (*s)[0][i]*(*s)[0][i] + (*s)[1][i]*(*s)[1][i] + (*s)[2][i]*(*s)[2][i];
        if (thisstr > thismax) thismax = thisstr;
      }
      return std::sqrt(thismax);
    } else {
      return 1.0;
    }
  }

  // add and return the total circulation of all elements
  std::array<S,3> get_total_circ(const double _time) {
    std::array<S,3> circ;
    circ.fill(0.0);

    if (s) {
      // we have strengths, add them up
      // this is the c++17 way
      //return std::reduce(std::execution::par, s->begin(), s->end());
      // this is the c++11 way
      for (size_t d=0; d<3; ++d) {
        circ[d] = std::accumulate((*s)[d].begin(), (*s)[d].end(), 0.0);
      }
    }

    return circ;
  }

  std::string to_string() const {
    std::string mystr;
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
  // if you add anything here, you need to wipe out all build files and run cmake again!

  // active, reactive, or inert?
  elem_t E;
  // how does it move? use move_t and Body*
  move_t M;
  // if attached to a body, which one?
  std::shared_ptr<Body> B;

  // common arrays for all derived types
  size_t n;

  // state vector
  std::array<Vector<S>,Dimensions> x;                   // position
  std::optional<std::array<Vector<S>,Dimensions>> s;    // strength

  // time derivative of state vector
  std::array<Vector<S>,Dimensions> u;                   // velocity
  //std::optional<std::array<Vector<S>,Dimensions>> dsdt; // strength change

  // for objects moving with a body
  std::optional<std::array<Vector<S>,Dimensions>> ux;   // untransformed position
};

