/*
 * Points.h - Specialized class for 3D stretchable points with optional vel gradients
 *
 * (c)2018-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega3D.h"
#include "VectorHelper.h"
#include "MathHelper.h"
#include "ElementBase.h"

#ifdef USE_GL
#include "GlState.h"
#include "RenderParams.h"
#include "OglHelper.h"
#include "ShaderHelper.h"
#include "glad.h"
#endif

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <optional>
#include <random>
#include <cassert>
#include <cmath>
#include <algorithm>


// 0-D elements
template <class S>
class Points: public ElementBase<S> {
public:
  // flexible constructor - use input 7*n vector (x, y, z, sx, sy, sz, r)
  //                         or input 3*n vector (x, y, z)
  Points(const std::vector<S>& _in,
         const elem_t _e,
         const move_t _m,
         std::shared_ptr<Body> _bp)
    : ElementBase<S>(0, _e, _m, _bp),
      max_strength(-1.0) {

    const size_t nper = (_e == inert) ? 3 : 7;
    std::cout << "  new collection with " << (_in.size()/nper);
    std::cout << ((_e == inert) ? " tracers" : " vortons") << std::endl;

    // need to reset the base class n
    this->n = _in.size()/nper;

    // make sure we have a complete input vector
    assert(_in.size() % nper == 0 && "Input array size is not a multiple of 3 or 7");

    // this initialization specific to Points
    for (size_t d=0; d<Dimensions; ++d) {
      this->x[d].resize(this->n);
      for (size_t i=0; i<this->n; ++i) {
        this->x[d][i] = _in[nper*i+d];
      }
    }

    // save untransformed positions if we are given a Body pointer
    if (_bp) {
      this->ux = this->x;
    }

    if (_e == inert) {
      // field points need neither radius, elong, nor strength

    } else {
      // active vortons need a radius
      r.resize(this->n);
      for (size_t i=0; i<this->n; ++i) {
        r[i] = _in[7*i+6];
      }

      // and elongation
      this->elong.resize(this->n);
      for (size_t i=0; i<this->n; ++i) {
        this->elong[i] = 1.0;
      }

      // optional strength in base class
      // need to assign it a vector first!
      std::array<Vector<S>,3> new_s;
      for (size_t d=0; d<3; ++d) {
        new_s[d].resize(this->n);
        for (size_t i=0; i<this->n; ++i) {
          new_s[d][i] = _in[7*i+d+3];
        }
      }
      this->s = std::move(new_s);
    }

    // velocity in base class
    for (size_t d=0; d<Dimensions; ++d) {
      this->u[d].resize(this->n);
    }

    // optional velgrads here
    if (_m == lagrangian) {
      std::array<Vector<S>,Dimensions*Dimensions> new_ug;
      for (size_t d=0; d<Dimensions*Dimensions; ++d) {
        new_ug[d].resize(this->n);
      }
      ug = std::move(new_ug);
    }
  }

  Vector<S>& get_elong() { return elong; }
  std::optional<std::array<Vector<S>,Dimensions*Dimensions>>& get_velgrad() { return ug; }
  const Vector<S>& get_rad() const { return r; }
  Vector<S>&       get_rad()       { return r; }

  // a little logic to see if we should augment the BEM equations for this object (see Surfaces.h)
  const bool is_augmented() const { return false; }

  // find out the next row index in the BEM after this collection
  // once we start supporting BEM unknowns on points, we'll have to change these
  void set_first_row(const Int _i) { return; }
  const Int get_first_row() const { return 0; }
  const Int get_num_rows()  const { return 0; }
  const Int get_next_row()  const { return 0; }

  const float get_max_bc_value() const { return 0.0; }

  void add_new(std::vector<float>& _in) {
    // remember old size and incoming size
    const size_t nold = this->n;

    const size_t nper = (this->E == inert) ? 3 : 7;
    assert(_in.size() % nper == 0 && "Input array size is not a multiple of 3 or 7");

    const size_t nnew = _in.size()/nper;
    std::cout << "  adding " << nnew << " particles to collection..." << std::endl;

    // must explicitly call the method in the base class first
    ElementBase<S>::add_new(_in);

    // then do local stuff
    if (this->E == inert) {
      // no radius needed

    } else {
      // active points need radius and elongation
      r.resize(nold+nnew);
      for (size_t i=0; i<nnew; ++i) {
        r[nold+i] = _in[7*i+6];
      }

      elong.resize(nold+nnew);
      for (size_t i=nold; i<nold+nnew; ++i) {
        elong[i] = 1.0;
      }
    }

    // save the new untransformed positions if we have a Body pointer
    if (this->B) {
      for (size_t d=0; d<Dimensions; ++d) {
        (*this->ux)[d].resize(nold+nnew);
        for (size_t i=nold; i<nold+nnew; ++i) {
          (*this->ux)[d][i] = this->x[d][i];
        }
      }
    }

    // optional vel grads don't need to be initialized
    if (ug) {
      for (size_t d=0; d<Dimensions*Dimensions; ++d) {
        (*ug)[d].resize(nold+nnew);
      }
    }
  }

  // up-size all arrays to the new size, filling with sane values
  void resize(const size_t _nnew) {
    const size_t currn = this->n;
    //std::cout << "  inside Points::resize with " << currn << " " << _nnew << std::endl;

    // must explicitly call the method in the base class - this sets n
    ElementBase<S>::resize(_nnew);

    if (_nnew == currn) return;

    // radii here
    if (this->E == inert) {
      // no radii or elongation
    } else {
      const size_t thisrn = r.size();
      r.resize(_nnew);
      for (size_t i=thisrn; i<_nnew; ++i) {
        r[i] = 1.0;
      }

      // elongations here
      const size_t thisen = elong.size();
      elong.resize(_nnew);
      for (size_t i=thisen; i<_nnew; ++i) {
        elong[i] = 1.0;
      }
    }

    // vel grads ((no need to set it)
    if (ug) {
      for (size_t d=0; d<Dimensions*Dimensions; ++d) {
        (*ug)[d].resize(_nnew);
      }
    }
  }

  void zero_vels() {
    // must explicitly call the method in the base class to zero the vels
    ElementBase<S>::zero_vels();
    // and specialize here for the vel grads
    if (ug) {
      for (size_t d=0; d<Dimensions*Dimensions; ++d) {
        for (size_t i=0; i<this->n; ++i) {
          (*ug)[d][i] = 0.0;
        }
      }
    }
  }

  void finalize_vels(const std::array<double,Dimensions>& _fs) {
    // must explicitly call the method in the base class, too
    ElementBase<S>::finalize_vels(_fs);
    // and specialize here for the vel grads
    if (ug) {
      const S factor = 0.25/M_PI;
      for (size_t d=0; d<Dimensions*Dimensions; ++d) {
        for (size_t i=0; i<this->n; ++i) {
          (*ug)[d][i] = (*ug)[d][i] * factor;
        }
      }
    }
  }

  void transform(const double _time) {
    // must explicitly call the method in the base class
    ElementBase<S>::transform(_time);
    // no other specialization required here
  }

  //
  // 1st order Euler advection and stretch
  //
  void move(const double _time, const double _dt) {
    // must explicitly call the method in the base class
    ElementBase<S>::move(_time, _dt);

    // and specialize
    if (this->M == lagrangian and ug and this->E != inert) {
      std::cout << "  Stretching" << to_string() << " using 1st order" << std::endl;
      S thismax = 0.0;

      for (size_t i=0; i<this->n; ++i) {
        std::array<S,Dimensions*Dimensions> this_ug = {0.0};
        for (size_t d=0; d<Dimensions*Dimensions; ++d) {
          this_ug[d] = (*ug)[d][i];
        }
        std::array<S,Dimensions> this_s = {0.0};
        for (size_t d=0; d<Dimensions; ++d) {
          this_s[d] = (*this->s)[d][i];
        }

        // compute stretch term
        // note that multiplying by the transpose may maintain linear impulse better, but
        //   severely underestimates stretch!
        std::array<S,3> wdu = {0.0};
        wdu[0] = this_s[0]*this_ug[0] + this_s[1]*this_ug[3] + this_s[2]*this_ug[6];
        wdu[1] = this_s[0]*this_ug[1] + this_s[1]*this_ug[4] + this_s[2]*this_ug[7];
        wdu[2] = this_s[0]*this_ug[2] + this_s[1]*this_ug[5] + this_s[2]*this_ug[8];

        // update elongation
        const S circmag = std::sqrt(this_s[0]*this_s[0] + this_s[1]*this_s[1] + this_s[2]*this_s[2]);
        const S elongfactor = (S)_dt * (this_s[0]*wdu[0] + this_s[1]*wdu[1] + this_s[2]*wdu[2]) / circmag;
        elong[i] *= 1.0 + elongfactor;

        // add Cottet SFS into stretch term (after elongation)

        // update strengths
        (*this->s)[0][i] = this_s[0] + _dt * wdu[0];
        (*this->s)[1][i] = this_s[1] + _dt * wdu[1];
        (*this->s)[2][i] = this_s[2] + _dt * wdu[2];

        // check for max strength
        S thisstr = std::pow((*this->s)[0][i], 2) + std::pow((*this->s)[1][i], 2) + std::pow((*this->s)[2][i], 2);
        if (thisstr > thismax) thismax = thisstr;

        if (false) {
        //if (i == 0) {
        //if (i < 10) {
          std::cout << "  x " << this->x[0][i] << " " << this->x[1][i] << " " << this->x[2][i];// << std::endl;
          std::cout << "  v " << this->u[0][i] << " " << this->u[1][i] << " " << this->u[2][i];// << std::endl;
          std::cout << "  ug " << this_ug[0] << " " << this_ug[1] << " " << this_ug[2];// << std::endl;
          //std::cout << "  wdu " << wdu[0] << " " << wdu[1] << " " << wdu[2];// << std::endl;
          std::cout << "  s " << (*this->s)[0][i] << " " << (*this->s)[1][i] << " " << (*this->s)[2][i];// << std::endl;
          //std::cout << "  elong " << elong[i];
          std::cout << std::endl;
        }
      }
      if (max_strength < 0.0) {
        max_strength = std::sqrt(thismax);
      } else {
        max_strength = 0.1*std::sqrt(thismax) + 0.9*max_strength;
      }
      //std::cout << "  New max_strength is " << max_strength << std::endl;
    } else {
      //std::cout << "  Not stretching" << to_string() << std::endl;
      max_strength = 1.0;
    }
  }

  //
  // 2nd order RK advection and stretch
  //
  void move(const double _time, const double _dt,
            const double _wt1, Points<S> const & _u1,
            const double _wt2, Points<S> const & _u2) {
    // must explicitly call the method in the base class
    ElementBase<S>::move(_time, _dt, _wt1, _u1, _wt2, _u2);

    // must confirm that incoming time derivates include velocity

    // and specialize
    if (this->M == lagrangian and this->E != inert and _u1.ug and _u2.ug) {
      std::cout << "  Stretching" << to_string() << " using 2nd order" << std::endl;
      S thismax = 0.0;

      for (size_t i=0; i<this->n; ++i) {

        // set up some convenient temporaries
        std::array<S,Dimensions> this_s = {0.0};
        std::array<S,Dimensions*Dimensions> this_ug = {0.0};
        for (size_t d=0; d<Dimensions; ++d) {
          this_s[d] = (*this->s)[d][i];
        }
        auto& optug1 = _u1.ug;
        for (size_t d=0; d<Dimensions*Dimensions; ++d) {
          this_ug[d] = (*optug1)[d][i];
        }

        // compute stretch term
        // note that multiplying by the transpose may maintain linear impulse better, but
        //   severely underestimates stretch!
        std::array<S,3> wdu1 = {0.0};
        wdu1[0] = this_s[0]*this_ug[0] + this_s[1]*this_ug[3] + this_s[2]*this_ug[6];
        wdu1[1] = this_s[0]*this_ug[1] + this_s[1]*this_ug[4] + this_s[2]*this_ug[7];
        wdu1[2] = this_s[0]*this_ug[2] + this_s[1]*this_ug[5] + this_s[2]*this_ug[8];

        auto& optug2 = _u2.ug;
        for (size_t d=0; d<Dimensions*Dimensions; ++d) {
          this_ug[d] = (*optug2)[d][i];
        }
        std::array<S,3> wdu2 = {0.0};
        wdu2[0] = this_s[0]*this_ug[0] + this_s[1]*this_ug[3] + this_s[2]*this_ug[6];
        wdu2[1] = this_s[0]*this_ug[1] + this_s[1]*this_ug[4] + this_s[2]*this_ug[7];
        wdu2[2] = this_s[0]*this_ug[2] + this_s[1]*this_ug[5] + this_s[2]*this_ug[8];

        std::array<S,3> wdu = {0.0};
        wdu[0] = _wt1*wdu1[0] + _wt2*wdu2[0];
        wdu[1] = _wt1*wdu1[1] + _wt2*wdu2[1];
        wdu[2] = _wt1*wdu1[2] + _wt2*wdu2[2];

        // update elongation
        const S circmagsqrd = this_s[0]*this_s[0] + this_s[1]*this_s[1] + this_s[2]*this_s[2];
        if (circmagsqrd > 0.0) {
          const S elongfactor = (S)_dt * (this_s[0]*wdu[0] + this_s[1]*wdu[1] + this_s[2]*wdu[2]) / circmagsqrd;
          elong[i] *= 1.0 + elongfactor;
        }

        // add Cottet SFS into stretch term (after elongation)

        // update strengths
        (*this->s)[0][i] = this_s[0] + _dt * wdu[0];
        (*this->s)[1][i] = this_s[1] + _dt * wdu[1];
        (*this->s)[2][i] = this_s[2] + _dt * wdu[2];

        // check for max strength
        S thisstr = std::pow((*this->s)[0][i], 2) + std::pow((*this->s)[1][i], 2) + std::pow((*this->s)[2][i], 2);
        if (thisstr > thismax) thismax = thisstr;

        if (false) {
          std::array<S,3> thisx = {this->x[0][i], this->x[1][i], this->x[2][i]};
          std::cout << "  p " << i << " at rad " << length(thisx) << " has circmag " << std::sqrt(thisstr) << " and new elong " << elong[i] << std::endl;
        }
      }

      if (max_strength < 0.0) {
        max_strength = std::sqrt(thismax);
      } else {
        max_strength = 0.05*std::sqrt(thismax) + 0.95*max_strength;
      }

    } else {
      //std::cout << "  Not stretching" << to_string() << std::endl;
      max_strength = 1.0;
    }
  }

  // find largest elongation in this collection
  S get_max_elong() {
    // max_element returns an iterator
    auto imax = std::max_element(this->elong.begin(), this->elong.end());
    //std::cout << "  max elong " << *imax << std::endl;
    return *imax;
  }

  // find the new peak strength magnitude
  void update_max_str() {
    S thismax = ElementBase<S>::get_max_str();

    // and slowly update the current saved value
    if (max_strength < 0.0) {
      max_strength = thismax;
    } else {
      max_strength = 0.1*thismax + 0.9*max_strength;
    }
  }

  // add and return the total impulse of all elements
  std::array<S,Dimensions> get_total_impulse() {

    // here is the return vector
    std::array<S,Dimensions> imp;
    imp.fill(0.0);

    if (this->s) {
      // accumulate impulse from each particle
      for (size_t i=0; i<this->n; ++i) {
        imp[0] += (*this->s)[1][i] * this->x[2][i] - (*this->s)[2][i] * this->x[1][i];
        imp[1] += (*this->s)[2][i] * this->x[0][i] - (*this->s)[0][i] * this->x[2][i];
        imp[2] += (*this->s)[0][i] * this->x[1][i] - (*this->s)[1][i] * this->x[0][i];
      }
    }

    return imp;
  }

  void add_body_motion(const S _factor, const double _time) {
    // no need to call base class now
    //ElementBase<S>::add_body_motion(_factor);
    // apply a factor times the body motion
    //for (size_t i=0; i<this->get_n(); ++i) {
    //  std::array<S,Dimensions> thisvel = B->get_vel(_time);
    //  // apply the velocity
    //  for (size_t d=0; d<Dimensions; ++d) {
    //    this->u[d][i] += _factor * thisvel[d];
    //  }
    //}
  }

  // augment the strengths with a value equal to that which accounts for
  //   the solid-body rotation of the object
  void add_unit_rot_strengths(const int _d) {}
  void add_solved_rot_strengths(const S _factor) {}

#ifdef USE_GL
  //
  // OpenGL functions
  //

  // helper function to clean up initGL
  void prepare_opengl_buffer(GLuint _prog, GLuint _idx, const GLchar* _name) {
    glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[_idx]);
    const GLint position_attribute = glGetAttribLocation(_prog, _name);
    // Specify how the data for position can be accessed
    glVertexAttribPointer(position_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);
    // Enable the attribute
    glEnableVertexAttribArray(position_attribute);
    // and tell it to advance two primitives per point (2 tris per quad)
    glVertexAttribDivisor(position_attribute, 1);
  }

  // this gets done once - load the shaders, set up the vao
  void initGL(std::vector<float>& _projmat,
              float*              _poscolor,
              float*              _negcolor,
              float*              _defcolor) {

    //std::cout << "inside Points.initGL" << std::endl;
    std::cout << "inside Points.initGL with E=" << this->E << " and M=" << this->M << std::endl;

    // generate the opengl state object with space for 7 vbos and 2 shader programs
    mgl = std::make_shared<GlState>(7,2);

    // Allocate space, but don't upload the data from CPU to GPU yet
    for (size_t i=0; i<Dimensions; ++i) {
      glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i]);
      glBufferData(GL_ARRAY_BUFFER, 0, this->x[i].data(), GL_STATIC_DRAW);
    }

    // here is where we split on element type: active/reactive vs. inert
    if (this->E == inert) {

      // Load and create the blob-drawing shader program
      mgl->spo[0] = create_draw_point_program();

      // Only send position arrays - no radius or strength
      prepare_opengl_buffer(mgl->spo[0], 0, "px");
      prepare_opengl_buffer(mgl->spo[0], 1, "py");
      prepare_opengl_buffer(mgl->spo[0], 2, "posz");

      // and for the compute shaders!

      // Get the location of the attributes that enters in the vertex shader
      mgl->projmat_attribute_pt = glGetUniformLocation(mgl->spo[0], "Projection");

      // upload the projection matrix
      glUniformMatrix4fv(mgl->projmat_attribute_pt, 1, GL_FALSE, _projmat.data());

      // locate where the colors and color scales go
      mgl->def_color_attribute = glGetUniformLocation(mgl->spo[0], "def_color");

      // send the current values
      glUniform4fv(mgl->def_color_attribute, 1, _defcolor);

      // locate where the point radius goes
      mgl->unif_rad_attribute = glGetUniformLocation(mgl->spo[0], "rad");

      // send the current values
      glUniform1f(mgl->unif_rad_attribute, (const GLfloat)0.01);

      // and indicate the fragment color output
      glBindFragDataLocation(mgl->spo[0], 0, "frag_color");

      // Initialize the quad attributes
      std::vector<float> quadverts = {-1,-1, 1,-1, 1,1, -1,1};
      glGenBuffers(1, &(mgl->qvbo));
      glBindBuffer(GL_ARRAY_BUFFER, mgl->qvbo);
      glBufferData(GL_ARRAY_BUFFER, quadverts.size()*sizeof(float), quadverts.data(), GL_STATIC_DRAW);

      mgl->quad_attribute_pt = glGetAttribLocation(mgl->spo[0], "quad_attr");
      glVertexAttribPointer(mgl->quad_attribute_pt, 2, GL_FLOAT, GL_FALSE, 0, 0);
      glEnableVertexAttribArray(mgl->quad_attribute_pt);

    } else { // this->E is active or reactive

      glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[6]);
      glBufferData(GL_ARRAY_BUFFER, 0, r.data(), GL_STATIC_DRAW);
      for (size_t i=0; i<Dimensions; ++i) {
        if (this->s) {
          glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i+3]);
          glBufferData(GL_ARRAY_BUFFER, 0, (*this->s)[i].data(), GL_STATIC_DRAW);
        }
      }

      // Load and create the blob-drawing shader program
      mgl->spo[1] = create_draw_blob_program();

      // Now do the seven arrays
      prepare_opengl_buffer(mgl->spo[1], 0, "px");
      prepare_opengl_buffer(mgl->spo[1], 1, "py");
      prepare_opengl_buffer(mgl->spo[1], 2, "posz");
      prepare_opengl_buffer(mgl->spo[1], 3, "sx");
      prepare_opengl_buffer(mgl->spo[1], 4, "sy");
      prepare_opengl_buffer(mgl->spo[1], 5, "sz");
      prepare_opengl_buffer(mgl->spo[1], 6, "r");

      // and for the compute shaders!

      // Get the location of the attributes that enters in the vertex shader
      mgl->projmat_attribute_bl = glGetUniformLocation(mgl->spo[1], "Projection");

      // upload the projection matrix
      glUniformMatrix4fv(mgl->projmat_attribute_bl, 1, GL_FALSE, _projmat.data());

      // locate where the colors and color scales go
      mgl->pos_color_attribute = glGetUniformLocation(mgl->spo[1], "pos_color");
      mgl->neg_color_attribute = glGetUniformLocation(mgl->spo[1], "neg_color");
      mgl->str_scale_attribute = glGetUniformLocation(mgl->spo[1], "str_scale");
      mgl->rad_scale_attribute = glGetUniformLocation(mgl->spo[1], "rad_scale");

      // send the current values
      glUniform4fv(mgl->pos_color_attribute, 1, (const GLfloat *)_poscolor);
      glUniform4fv(mgl->neg_color_attribute, 1, (const GLfloat *)_negcolor);
      glUniform1f (mgl->str_scale_attribute, (const GLfloat)1.0);
      glUniform1f (mgl->rad_scale_attribute, (const GLfloat)1.0);

      // and indicate the fragment color output
      glBindFragDataLocation(mgl->spo[1], 0, "frag_color");

      // Initialize the quad attributes
      std::vector<float> quadverts = {-1,-1, 1,-1, 1,1, -1,1};
      glGenBuffers(1, &(mgl->qvbo));
      glBindBuffer(GL_ARRAY_BUFFER, mgl->qvbo);
      glBufferData(GL_ARRAY_BUFFER, quadverts.size()*sizeof(float), quadverts.data(), GL_STATIC_DRAW);

      mgl->quad_attribute_bl = glGetAttribLocation(mgl->spo[1], "quad_attr");
      glVertexAttribPointer(mgl->quad_attribute_bl, 2, GL_FLOAT, GL_FALSE, 0, 0);
      glEnableVertexAttribArray(mgl->quad_attribute_bl);

    } // end this->E is active or reactive

    glBindVertexArray(0);
  }

  // this gets done every time we change the size of the positions array
  void updateGL() {
    //std::cout << "inside Points.updateGL" << std::endl;

    // has this been init'd yet?
    if (not mgl) return;
    if (glIsVertexArray(mgl->vao) == GL_FALSE) return;

    const size_t vlen = this->x[0].size()*sizeof(S);
    if (vlen > 0) {
      glBindVertexArray(mgl->vao);

      // Indicate and upload the data from CPU to GPU
      for (size_t i=0; i<Dimensions; ++i) {
        // the positions
        glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i]);
        glBufferData(GL_ARRAY_BUFFER, vlen, this->x[i].data(), GL_DYNAMIC_DRAW);
      }

      // here is where we split on element type: active/reactive vs. inert
      if (this->E == inert) {

        // no strengths or radii needed or present

      } else { // this->E is active or reactive

        // the strengths
        if (this->s) {
          for (size_t i=0; i<Dimensions; ++i) {
            glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i+3]);
            glBufferData(GL_ARRAY_BUFFER, vlen, (*this->s)[i].data(), GL_DYNAMIC_DRAW);
          }
        }
        // the radii
        glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[6]);
        glBufferData(GL_ARRAY_BUFFER, vlen, r.data(), GL_DYNAMIC_DRAW);
      }

      glBindVertexArray(0);

      // must tell draw call how many elements are there
      mgl->num_uploaded = this->x[0].size();
    }
  }

  // OpenGL3 stuff to display points, called once per frame
  void drawGL(std::vector<float>& _projmat,
              RenderParams&       _rparams,
              const float         _vdelta) {

    //std::cout << "inside Points.drawGL" << std::endl;

    // has this been init'd yet?
    if (not mgl) {
      initGL(_projmat, _rparams.pos_circ_color,
                       _rparams.neg_circ_color,
                       _rparams.default_color);
      updateGL();
    }

    if (mgl->num_uploaded > 0) {
      glBindVertexArray(mgl->vao);

      // get blending ready
      glDisable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);
      glBlendFunc(GL_ONE, GL_ONE);

      // here is where we split on element type: active/reactive vs. inert
      if (this->E == inert) {

        // draw as small dots
        glUseProgram(mgl->spo[0]);

        glEnableVertexAttribArray(mgl->quad_attribute_pt);

        // upload the current uniforms
        glUniformMatrix4fv(mgl->projmat_attribute_pt, 1, GL_FALSE, _projmat.data());
        glUniform4fv(mgl->def_color_attribute, 1, (const GLfloat *)_rparams.default_color);
        glUniform1f (mgl->unif_rad_attribute, (const GLfloat)(1.0f*_rparams.tracer_size));

        // the one draw call here
        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 4, mgl->num_uploaded);

      } else { // this->E is active or reactive

        // draw as colored clouds
        glUseProgram(mgl->spo[1]);

        glEnableVertexAttribArray(mgl->quad_attribute_bl);

        // upload the current projection matrix
        glUniformMatrix4fv(mgl->projmat_attribute_bl, 1, GL_FALSE, _projmat.data());

        // upload the current color values
        const float color_scaling = _rparams.circ_density * std::pow(_vdelta/_rparams.vorton_scale,2) / max_strength;
        glUniform4fv(mgl->pos_color_attribute, 1, (const GLfloat *)_rparams.pos_circ_color);
        glUniform4fv(mgl->neg_color_attribute, 1, (const GLfloat *)_rparams.neg_circ_color);
        glUniform1f (mgl->str_scale_attribute, (const GLfloat)color_scaling);
        glUniform1f (mgl->rad_scale_attribute, (const GLfloat)_rparams.vorton_scale);

        // the one draw call here
        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 4, mgl->num_uploaded);
      }

      // return state
      glEnable(GL_DEPTH_TEST);
      glDisable(GL_BLEND);
      glBindVertexArray(0);
    }
  }
#endif

  std::string to_string() const {
    std::string retstr = " " + std::to_string(this->n) + ElementBase<S>::to_string() + " Points";
    if (ug) retstr += " with grads";
    return retstr;
  }

protected:
  // additional state vectors
  Vector<S> r;		// thickness/radius
  Vector<S> elong;	// scalar elongation, does not require register alignment

  // derivatives of state vector
  std::optional<std::array<Vector<S>,Dimensions*Dimensions>> ug;   // velocity gradients

private:
#ifdef USE_GL
  std::shared_ptr<GlState> mgl;
#endif
  float max_strength;
};

