/*
 * Points.h - Specialized class for 3D stretchable points with optional vel gradients
 *
 * (c)2018 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "VectorHelper.h"
#include "ElementBase.h"

#ifdef USE_GL
#include "ShaderHelper.h"
#include "glad.h"
#endif

#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#include <random>
#include <cassert>
#define _USE_MATH_DEFINES
#include <cmath>


// 0-D elements
template <class S>
class Points: public ElementBase<S> {
public:
  // flexible constructor - use input 7*n vector
  Points(const std::vector<S>& _in, const elem_t _e, const move_t _m)
    : ElementBase<S>(_in.size()/7, _e, _m),
      max_strength(-1.0) {

    std::cout << "  new collection with " << (_in.size()/7) << " particles..." << std::endl;

    // make sure we have a complete input vector
    assert(_in.size() % 7 == 0);

    // this initialization specific to Points
    for (size_t d=0; d<Dimensions; ++d) {
      this->x[d].resize(this->n);
      for (size_t i=0; i<this->n; ++i) {
        this->x[d][i] = _in[7*i+d];
      }
    }
    this->r.resize(this->n);
    this->elong.resize(this->n);
    for (size_t i=0; i<this->n; ++i) {
      this->r[i] = _in[7*i+6];
      this->elong[i] = 1.0;
    }

    // optional strength in base class
    if (_e != inert) {
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

  std::optional<std::array<Vector<S>,Dimensions*Dimensions>>& get_velgrad() { return ug; }

  void add_new(std::vector<float>& _in) {
    // remember old size and incoming size
    const size_t nold = this->n;
    const size_t nnew = _in.size()/7;
    std::cout << "  adding " << nnew << " particles to collection..." << std::endl;

    // must explicitly call the method in the base class first
    ElementBase<S>::add_new(_in);

    // then do local stuff
    elong.resize(nold+nnew);
    for (size_t i=nold; i<nold+nnew; ++i) {
      elong[i] = 1.0;
    }

    // optional vel grads don't need to be initialized
    if (ug) {
      for (size_t d=0; d<Dimensions*Dimensions; ++d) {
        (*ug)[d].resize(nold+nnew);
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
      for (size_t d=0; d<Dimensions*Dimensions; ++d) {
        for (size_t i=0; i<this->n; ++i) {
          (*ug)[d][i] = (*ug)[d][i] * 0.25/M_PI;
        }
      }
    }
  }

  void move(const double _dt) {
    // must explicitly call the method in the base class
    ElementBase<S>::move(_dt);
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

  void move(const double _dt,
            const double _wt1, Points<S> const & _u1,
            const double _wt2, Points<S> const & _u2) {
    // must explicitly call the method in the base class
    ElementBase<S>::move(_dt, _wt1, _u1, _wt2, _u2);

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
      }

      if (max_strength < 0.0) {
        max_strength = std::sqrt(thismax);
      } else {
        max_strength = 0.1*std::sqrt(thismax) + 0.9*max_strength;
      }

    } else {
      //std::cout << "  Not stretching" << to_string() << std::endl;
      max_strength = 1.0;
    }
  }

#ifdef USE_GL
  //
  // OpenGL functions
  //

  // this gets done once - load the shaders, set up the vao
  void initGL(std::vector<float>& _projmat,
              float*              _poscolor,
              float*              _negcolor) {

    std::cout << "inside Points.initGL" << std::endl;

    // Use a Vertex Array Object
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Load and create the blob-drawing shader program
    blob_program = create_particle_program();

    // Create seven Vector Buffer Objects that will store the vertices on video memory
    glGenBuffers(7, vbo);

    // Allocate space, but don't upload the data from CPU to GPU yet
    for (size_t i=0; i<Dimensions; ++i) {
      glBindBuffer(GL_ARRAY_BUFFER, vbo[i]);
      glBufferData(GL_ARRAY_BUFFER, 0, this->x[i].data(), GL_STATIC_DRAW);
      if (this->s) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo[i+3]);
        glBufferData(GL_ARRAY_BUFFER, 0, (*this->s)[i].data(), GL_STATIC_DRAW);
      }
    }
    glBindBuffer(GL_ARRAY_BUFFER, vbo[6]);
    glBufferData(GL_ARRAY_BUFFER, 0, this->r.data(), GL_STATIC_DRAW);

    // Get the location of the attributes that enters in the vertex shader
    projmat_attribute = glGetUniformLocation(blob_program, "Projection");

    // Now do the seven arrays
    GLint position_attribute;
    glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
    position_attribute = glGetAttribLocation(blob_program, "px");

    // Specify how the data for position can be accessed
    glVertexAttribPointer(position_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);

    // Enable the attribute
    glEnableVertexAttribArray(position_attribute);

    // and tell it to advance two primitives per point
    glVertexAttribDivisor(position_attribute, 1);

    // do it for the rest
    glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
    position_attribute = glGetAttribLocation(blob_program, "py");
    glVertexAttribPointer(position_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(position_attribute);
    glVertexAttribDivisor(position_attribute, 1);

    glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
    position_attribute = glGetAttribLocation(blob_program, "posz");
    glVertexAttribPointer(position_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(position_attribute);
    glVertexAttribDivisor(position_attribute, 1);

    glBindBuffer(GL_ARRAY_BUFFER, vbo[3]);
    position_attribute = glGetAttribLocation(blob_program, "sx");
    glVertexAttribPointer(position_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(position_attribute);
    glVertexAttribDivisor(position_attribute, 1);

    glBindBuffer(GL_ARRAY_BUFFER, vbo[4]);
    position_attribute = glGetAttribLocation(blob_program, "sy");
    glVertexAttribPointer(position_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(position_attribute);
    glVertexAttribDivisor(position_attribute, 1);

    glBindBuffer(GL_ARRAY_BUFFER, vbo[5]);
    position_attribute = glGetAttribLocation(blob_program, "sz");
    glVertexAttribPointer(position_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(position_attribute);
    glVertexAttribDivisor(position_attribute, 1);

    glBindBuffer(GL_ARRAY_BUFFER, vbo[6]);
    position_attribute = glGetAttribLocation(blob_program, "r");
    glVertexAttribPointer(position_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(position_attribute);
    glVertexAttribDivisor(position_attribute, 1);

    // upload the projection matrix
    glUniformMatrix4fv(projmat_attribute, 1, GL_FALSE, _projmat.data());

    // locate where the colors and color scales go
    pos_color_attribute = glGetUniformLocation(blob_program, "pos_color");
    neg_color_attribute = glGetUniformLocation(blob_program, "neg_color");
    str_scale_attribute = glGetUniformLocation(blob_program, "str_scale");

    // send the current values
    glUniform4fv(pos_color_attribute, 1, (const GLfloat *)_poscolor);
    glUniform4fv(neg_color_attribute, 1, (const GLfloat *)_negcolor);
    glUniform1f (str_scale_attribute, (const GLfloat)1.0);

    // and indicate the fragment color output
    glBindFragDataLocation(blob_program, 0, "frag_color");

    // Initialize the quad attributes
    std::vector<float> quadverts = {-1,-1, 1,-1, 1,1, -1,1};
    GLuint qvbo;
    glGenBuffers(1, &qvbo);
    glBindBuffer(GL_ARRAY_BUFFER, qvbo);
    glBufferData(GL_ARRAY_BUFFER, quadverts.size()*sizeof(float), quadverts.data(), GL_STATIC_DRAW);

    quad_attribute = glGetAttribLocation(blob_program, "quad_attr");
    glVertexAttribPointer(quad_attribute, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(quad_attribute);
  }

  // this gets done every time we change the size of the positions array
  void updateGL() {
    //std::cout << "inside Points.updateGL" << std::endl;

    // has this been init'd yet?
    if (glIsVertexArray(vao) == GL_FALSE) return;

    const size_t vlen = this->x[0].size()*sizeof(S);
    if (vlen > 0) {
      // Indicate and upload the data from CPU to GPU
      for (size_t i=0; i<Dimensions; ++i) {
        // the positions
        glBindBuffer(GL_ARRAY_BUFFER, vbo[i]);
        glBufferData(GL_ARRAY_BUFFER, vlen, this->x[i].data(), GL_DYNAMIC_DRAW);
        // the strengths
        if (this->s) {
          glBindBuffer(GL_ARRAY_BUFFER, vbo[i+3]);
          glBufferData(GL_ARRAY_BUFFER, vlen, (*this->s)[i].data(), GL_DYNAMIC_DRAW);
        }
      }
      // the radii
      glBindBuffer(GL_ARRAY_BUFFER, vbo[6]);
      glBufferData(GL_ARRAY_BUFFER, vlen, this->r.data(), GL_DYNAMIC_DRAW);
    }
  }

  // OpenGL3 stuff to display points, called once per frame
  void drawGL(std::vector<float>& _projmat,
              float*              _poscolor,
              float*              _negcolor) {

    //std::cout << "inside Points.drawGL" << std::endl;

    // has this been init'd yet?
    if (glIsVertexArray(vao) == GL_FALSE) {
      initGL(_projmat, _poscolor, _negcolor);
      updateGL();
    }

    if (this->n > 0) {
      glBindVertexArray(vao);

      // get blending ready
      glDisable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);
      glBlendFunc(GL_ONE, GL_ONE);

      // draw as colored clouds
      if (true) {
        glUseProgram(blob_program);

        glEnableVertexAttribArray(quad_attribute);

        // upload the current projection matrix
        glUniformMatrix4fv(projmat_attribute, 1, GL_FALSE, _projmat.data());

        // upload the current color values
        glUniform4fv(pos_color_attribute, 1, (const GLfloat *)_poscolor);
        glUniform4fv(neg_color_attribute, 1, (const GLfloat *)_negcolor);
        glUniform1f (str_scale_attribute, (const GLfloat)(0.15f/max_strength));

        // the one draw call here
        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 4, this->n);
      }

      // return state
      glEnable(GL_DEPTH_TEST);
      glDisable(GL_BLEND);
    }
  }
#endif

  std::string to_string() const {
    std::string retstr = ElementBase<S>::to_string() + " Points";
    if (ug) retstr += " with grads";
    return retstr;
  }

protected:
  // movement
  //std::optional<Body&> b;
  // state vector
  std::vector<S> elong;   // scalar elongation, does not require register alignment
  // time derivative of state vector
  std::optional<std::array<Vector<S>,Dimensions*Dimensions>> ug;   // velocity gradients

private:
#ifdef USE_GL
  // OpenGL stuff
  GLuint vao, vbo[7];
  GLuint blob_program;
  GLint projmat_attribute, quad_attribute;
  GLint pos_color_attribute, neg_color_attribute, str_scale_attribute;
#endif
  float max_strength;
};

