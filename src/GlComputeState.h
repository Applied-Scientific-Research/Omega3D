/*
 * GlComputeState.h - Store compute-specific OpenGL state
 *
 * (c)2019-20 Applied Scientific Research, Inc.
 *            Written by Mark J Stock <markjstock@gmail.com>
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
#include "Collection.h"

#include "glad.h"

#include <iostream>
#include <vector>
#include <cassert>
#include <atomic>


// 1-D elements
template <class S>
class GlComputeState {
public:
  GlComputeState<S>(const int _nvbo, const int _nspo) {

    assert(_nvbo>=0 && "Invalid number of VBOs requested");
    assert(_nspo>=0 && "Invalid number of shader program objects requested");
    std::cout << "new GlComputeState with " << _nvbo << " buffers and " << _nspo << " shader programs" << std::endl;

    // generate the vao
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // and the vbo array
    vbo.resize(_nvbo);
    glGenBuffers(_nvbo, vbo.data());

    // finally the shader program array
    spo.resize(_nspo);

    cstate.store(no_compute);
  }

  // must specifically destroy buffers
  ~GlComputeState() {
    //std::cout << "In ~GlComputeState glIsVertexArray(vao) == " << (glIsVertexArray(vao)==GL_TRUE) << std::endl;
    if (glIsVertexArray(vao) == GL_FALSE) return;
    //std::cout << " finishing ~GlComputeState" << std::endl;

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

  // must remember how many are uploaded
  GLuint nsrc, ntarg;

  // interthread coordination state
  std::atomic<compute_state_t> cstate;

  // attributes
  GLint source_offset_attr, source_count_attr, target_offset_attr, target_count_attr;
  GLint target_offset_attr_z, target_count_attr_z;

  // the current source and target Collections
  //std::shared_ptr<Collection> src, targ;
  // can't store collections because of circular declarations
  // store pointers directly to the arrays
  //std::shared_ptr<std::array<Vector<S>,3>> sx, ss, tx, tu;
  //std::shared_ptr<std::array<Vector<S>,9>> tug;
  //std::shared_ptr<Vector<S>> sr, tr;
  // or, just rely on whomever sets the state to set the vectors below...

  // the host-side compacted arrays
  // two for sources (position and strength), one for targets (position),
  //   and three for results (3*4 = 12 = 3 vels + 9 vel grads)
  std::vector<S> hsx, hss, htx, hr1, hr2, hr3;
};

