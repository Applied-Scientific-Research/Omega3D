/*
 * GlState.h - Store general OpenGL state for drawing a Collection
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include "ComputeState.h"

#include "glad.h"

#include <iostream>
#include <vector>
#include <cassert>
#include <atomic>


// Hold state for OpenGL
class GlState {
public:
  GlState(const int _nvbo, const int _nspo) {

    assert(_nvbo>=0 && "Invalid number of VBOs requested");
    assert(_nspo>=0 && "Invalid number of shader program objects requested");
    std::cout << "  new GlState with " << _nvbo << " buffers and " << _nspo << " shader programs" << std::endl;

    // generate the vao
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // and the vbo array
    vbo.resize(_nvbo);
    glGenBuffers(_nvbo, vbo.data());

    // finally the shader program array
    spo.resize(_nspo);

    cstate.store(no_compute);
    num_uploaded = 0;
  }

  // must specifically destroy buffers
  ~GlState() {
    //std::cout << "In ~GlState glIsVertexArray(vao) == " << (glIsVertexArray(vao)==GL_TRUE) << std::endl;
    // closing out of the program dies here, but we put GlState in a unique_ptr, so wtf?
    if (glIsVertexArray(vao) == GL_FALSE) return;
    //std::cout << " finishing ~GlState" << std::endl;

    glBindVertexArray(vao);
    glDeleteBuffers(vbo.size(), vbo.data());
    for (size_t i=0; i<spo.size(); ++i) {
      glDeleteProgram(spo[i]);
    }
    glDeleteVertexArrays(1,&vao);
    glBindVertexArray(0);
  }

  //
  // member variables
  //

  // the one vertex array object per Collection
  GLuint vao;

  // some number of vertex buffer objects (usually 4)
  std::vector<GLuint> vbo;
  GLuint qvbo;

  // shader program objects (1 or 2)
  std::vector<GLuint> spo;

  // we probably don't need this
  GLsizei num_uploaded;

  // testing the interthread coordination idea
  std::atomic<compute_state_t> cstate;

  // drawing attributes
  GLint projmat_attribute, projmat_attribute_bl, projmat_attribute_pt;
  GLint mvmat_attribute_bl, mvmat_attribute_pt;
  GLint quad_attribute_bl, quad_attribute_pt;
  GLint def_color_attribute, pos_color_attribute, neg_color_attribute; //, back_color_attribute; 
  GLint use_def_attribute; //, use_back_attribute;
  GLint str_scale_attribute, unif_rad_attribute, rad_scale_attribute;

  // compute attributes
  GLint source_offset_attr, source_count_attr, target_offset_attr, target_count_attr;
};

