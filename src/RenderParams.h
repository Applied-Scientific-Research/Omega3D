/*
 * RenderParams.h - Structure to contain rendering/drawing parameters
 *
 * (c)2019-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
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

#include "json/json.hpp"

//
// Class-like struct for all Imgui and OpenGL render parameters
//
class RenderParams {
public:
  RenderParams() = default;
  ~RenderParams() = default;

  RenderParams(RenderParams const&) = default; //allow copy
  RenderParams(RenderParams&&) = default; //allow move
  RenderParams& operator=(RenderParams const&) = default; //allow copy
  RenderParams& operator=(RenderParams&&) = default; //allow move

  // read to and write from a json object
  void from_json(const nlohmann::json);
  nlohmann::json to_json() const;

  // public-equivalent data

  // window space (from glfwGetWindowSize, not Framebuffer)
  int width = 1280;
  int height = 720;

  // view position
  float vcx = 0.0f;	// "truck"
  float vcy = 0.0f;	// "pedestal"
  float vcz = -3.0f;	// "dolly"
  float vsize = 2.0f;
  float vfov = 35.0f;

  // view rotation (radians)
  float rx = 0.0f;
  float ry = 0.0f;

  // colors
  float pos_circ_color[4] = {239./255.,  54./255.,  54./255., 1.0};
  float neg_circ_color[4] = { 64./255.,  64./255., 255./255., 1.0};
  float default_color[4]  = { 68./255.,  68./255.,  68./255., 1.0};
  float clear_color[4]    = {  0./255.,   0./255.,   0./255., 1.0};

  // other
  float circ_density = 0.15;	// non-dimensional
  float vorton_scale = 0.1;	// non-dimensional
  float tracer_scale = 0.15;	// non-dimensional
  float tracer_size;		// dimensional
};

