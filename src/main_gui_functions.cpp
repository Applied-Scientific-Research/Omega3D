/*
 * main_gui_functions.cpp - contains helper, non-class related functions used in main_gui.cpp
 *                          to in increase readability
 *
 * (c) 2020 Applied Scientific Reasearch, Inc.
 *          Mark J Stock <markjstock@gmail.com>
 *          Blake B Hillier <blakehillier@mac.com>
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

#ifdef USE_IMGUI
#include "imgui/imgui.h"

#include <Eigen/Geometry>
#include <GLFW/glfw3.h>

#include <algorithm>

static void error_callback(int error, const char* description) {
  fprintf(stderr, "Error %d: %s\n", error, description);
}

//static void keyboard_callback(int key, int action) {
//  printf("%d %d\n", key, action);
//}

//
// this is NOT a GLFW callback, but my own, and it needs the
// ImGuiIO data structure for information on the mouse state
// There are no OpenGL calls here
//
void mouse_callback(GLFWwindow* /*_thiswin*/,
                    ImGuiIO& io,
                    float*   _cx,
                    float*   _cy,
                    float*   _cz,
                    float*   _rx,
                    float*   _ry,
                    float*   _size) {

  // use left-click and drag to rotate
  static bool lbutton_down = false;
  if (io.MouseClicked[0]) lbutton_down = true;
  if (io.MouseReleased[0]) lbutton_down = false;

  if (lbutton_down) {
    //std::cout << "free mouse moved " << io.MouseDelta.x << " " << io.MouseDelta.y << std::endl;
    // this works on a Retina display:
    (*_rx) += 3.0f * (float)io.MouseDelta.x / io.DisplaySize.x;
    (*_ry) += 3.0f * (float)io.MouseDelta.y / io.DisplaySize.x;
  }

  // use left button drag to pan
  static bool rbutton_down = false;
  if (io.MouseClicked[1]) rbutton_down = true;
  if (io.MouseReleased[1]) rbutton_down = false;
  if (rbutton_down) {
    // this worked on Linux:
    //int display_w, display_h;
    //glfwGetFramebufferSize(_thiswin, &display_w, &display_h);
    //(*_cx) -= 2.0f * (*_size) * (float)io.MouseDelta.x / (float)display_w;
    //(*_cy) += 2.0f * (*_size) * (float)io.MouseDelta.y / (float)display_w;

    // this works on a Retina display:
    (*_cx) -= 2.0f * (*_size) * (float)io.MouseDelta.x / io.DisplaySize.x;
    (*_cy) += 2.0f * (*_size) * (float)io.MouseDelta.y / io.DisplaySize.x;
  }

  // then, use scroll wheel to zoom!
  //std::cout << "free mouse wheel " << io.MouseWheel << std::endl;
  if (io.MouseWheel != 0 and false) {
    // do your drag here
    //int display_w, display_h;
    //glfwGetFramebufferSize(_thiswin, &display_w, &display_h);

    // change the size
    const float oldsize = (*_size);
    (*_size) *= std::pow(1.1f, io.MouseWheel);

    // and adjust the center such that the zoom occurs about the pointer!
    //const float ar = (float)display_h / (float)display_w;
    const float ar = io.DisplaySize.y / io.DisplaySize.x;

    // this only scales around world origin
    //(*_cx) += 2.0f * ((float)io.MousePos.x / (float)display_w - 0.5f) * (oldsize - (*_size));
    //(*_cy) += 2.0f * (0.5f - (float)io.MousePos.y / (float)display_h) * (oldsize - (*_size)) * ar;
    (*_cx) += 2.0f * ((float)io.MousePos.x / io.DisplaySize.x - 0.5f) * (oldsize - (*_size));
    (*_cy) += 2.0f * (0.5f - (float)io.MousePos.y / io.DisplaySize.y) * (oldsize - (*_size)) * ar;
  }

  // mouse wheel dollies
  if (io.MouseWheel != 0) {
    //(*_cz) *= std::pow(1.05f, io.MouseWheel);
    (*_cz) += 0.1f * io.MouseWheel;
  }
}

//
// Generate and overwrite the modelview matrix (model space to world space to eye space)
//
void compute_modelview_mat(const float         _cx,
                           const float         _cy,
                           const float         _cz,
                           const float         _rx,
                           const float         _ry,
                           std::vector<float>& _mvmat) {

  // make an affine transformation
  Eigen::Transform<float,3,Eigen::Affine> trans(Eigen::Transform<float,3,Eigen::Affine>::Identity());

  // back the camera up so we can see the simulation
  trans.translate(Eigen::Vector3f(-_cx, -_cy, _cz));

  // auto-rotate around the origin
  //static float theta = 0.0;
  //trans.rotate(Eigen::AngleAxisf(theta, Eigen::Vector3f::UnitY()));
  //theta += 0.01f;
  trans.rotate(Eigen::AngleAxisf(_rx, Eigen::Vector3f::UnitY()));
  trans.rotate(Eigen::AngleAxisf(_ry, Eigen::Vector3f::UnitX()));

  // or arcball rotate

  // and write it into the matrix
  Eigen::Map<Eigen::Matrix<float,4,4>> pmat(_mvmat.data());
  pmat = trans.matrix();
}

//
// Helper routine to determine orthographic projection matrix
// given coords at screen center and a measure of size
// Also changes overall pixels-to-length scale
//
void compute_ortho_proj_mat(GLFWwindow*         _thiswin,
                            float&              _size,
                            std::vector<float>& _projmat) {

  // track changes in window!
  static int last_w, last_h = -1;

  // get current window size
  int display_w, display_h;
  glfwGetFramebufferSize(_thiswin, &display_w, &display_h);

  // compare window size to previous call
  if (last_h != -1) {
    if (last_h != display_h or last_w != display_w) {
      // window aspect ratio changed, adjust _size
      _size *= sqrt(  ((float)last_h   /(float)last_w   )
                    / ((float)display_h/(float)display_w));
    }
  }

  // near and far clipping planes (positive distance)
  const float near = -50.0;
  const float far = 50.0;

  const float vsx = _size;
  const float vsy = _size * (float)display_h / (float)display_w;
  const float fac = 1.0f / (near-far);
  _projmat =
    { 1.0f/vsx, 0.0f,     0.0f,           0.0f,
      0.0f,     1.0f/vsy, 0.0f,           0.0f,
      0.0f,     0.0f,     2.0f*fac,       0.0f,
      0.0f,     0.0f,     (far+near)*fac, 1.0f };

  // save window size for next call
  last_w = display_w;
  last_h = display_h;
}

//
// Helper routine to determine perspective projection matrix
// given coords at screen center and a measure of size
// Also changes overall pixels-to-length scale
//
void compute_persp_proj_mat(GLFWwindow*         _thiswin,
                            float&              _fov,
                            std::vector<float>& _projmat) {

  // get current window size
  int display_w, display_h;
  glfwGetFramebufferSize(_thiswin, &display_w, &display_h);

  // validate fov
  const float thisfov = std::min(std::max(_fov,10.f),160.f);
  _fov = thisfov;

  // near and far clipping planes (positive distance)
  const float near = 0.1;
  const float far = 10.0;

  // off-axis projection
  if (false) {
    // compute precursors
    const float viewscale = std::tan(thisfov * 0.5 * M_PI / 180.0) * near;
    const float right = viewscale * (float)display_w / (float)display_h;
    const float left = -right;
    const float top = viewscale;
    const float bottom = -top;
    // generate matrix
    _projmat =
      { 2.0f*near/(right-left),    0.0f,                      0.0f,                     0.0f,
        0.0f,                      2.0f*near/(top-bottom),    0.0f,                     0.0f,
        (right+left)/(right-left), (top+bottom)/(top-bottom), (near+far)/(near-far),    -1.0f,
        0.0f,                      0.0f,                      2.0f*far*near/(near-far), 0.0f };

  } else {
    // this is the same, but for on-axis only
    const float viewscale = std::tan(thisfov * 0.5 * M_PI / 180.0);
    const float right = viewscale * (float)display_w / (float)display_h;
    const float top = viewscale;
    _projmat =
      { 1.0f/right, 0.0f,     0.0f,                     0.0f,
        0.0f,       1.0f/top, 0.0f,                     0.0f,
        0.0f,       0.0f,     (near+far)/(near-far),    -1.0f,
        0.0f,       0.0f,     2.0f*far*near/(near-far), 0.0f };
  }

  // alternatively
  //Eigen::Map<Eigen::Matrix<float,4,4>> pmat(_projmat.data());

  // translate it to position
  //Eigen::Transform<float,3,Eigen::Affine> trans(Eigen::Transform<float,3,Eigen::Affine>::Identity());
  //static float theta = 0.0;
  //trans.translate(Eigen::Vector3f(0.0f, 0.0f, -2.0f));
  //trans.rotate(Eigen::AngleAxisf(theta, Eigen::Vector3f::UnitY()));
  //theta += 0.01f;

  // apply it to the matrix
  //pmat *= trans.matrix();
}

//
// resize a window and framebuffer programmatically
//
void resize_to_resolution(GLFWwindow* window, const int new_w, const int new_h) {

  // get framebuffer size
  int fb_w, fb_h;
  glfwGetFramebufferSize(window, &fb_w, &fb_h);
  //std::cout << "Framebuffer size is " << fb_w << " x " << fb_h << std::endl;

  // get window size
  int ws_w, ws_h;
  glfwGetWindowSize(window, &ws_w, &ws_h);
  //std::cout << "Window size is " << ws_w << " x " << ws_h << std::endl;

  // on normal monitors, these numbers should be the same; on retina displays, they may not

  // check and resize
  if (fb_w != new_w or fb_h != new_h) {
    // on a retina display...do anything different?

    glfwSetWindowSize(window, new_w, new_h);
    std::cout << "Resizing window/framebuffer to " << new_w << " x " << new_h << std::endl;
  }
}

void LoadJsonSims(std::vector<nlohmann::json> &sims, std::vector<std::string> &descriptions, const std::string dirPath) {
  std::list<std::string> fileNames = MiniPath::listFiles(dirPath, "*.json");
  std::string sysDelim = MiniPath::getSystemDelim();
  std::cout << "Reading in from " << dirPath << std::endl;
  for(const std::string& s : fileNames) {
    sims.push_back(read_json(dirPath+sysDelim+s));
    descriptions.push_back(sims.back()["description"]);
  }
}

// execution starts here
void draw_render_gui(RenderParams &rp) {
  ImGui::ColorEdit3("positive circulation", rp.pos_circ_color);
  ImGui::ColorEdit3("negative circulation", rp.neg_circ_color);
  ImGui::ColorEdit3("feature color",        rp.default_color);
  ImGui::ColorEdit3("background color",     rp.clear_color);
  //ImGui::Checkbox("show origin", &show_origin);
  ImGui::SliderFloat("particle brightness", &(rp.circ_density), 0.001f, 1.0f, "%.3f", 2.0f);
  ImGui::SliderFloat("particle scale", &(rp.vorton_scale), 0.01f, 2.0f, "%.2f", 2.0f);

  if (ImGui::Button("Recenter")) {
    // put everything back to center
    rp.vcx = 0.0f;
    rp.vcy = 0.0f;
    rp.vcz = -3.0f;
    rp.vsize = 2.0f;
    rp.vfov = 35.0f;
  }
  // add button to recenter on all vorticity?
}

void draw_stats_window(const long int numPanels, const long int numFieldPts, const long int step, const float time,
                       const long int numParticles, bool* showStatsWindow, const int fontSize, const float displayH) {
  //std::cout << "Creating stats window: " << std::endl;
  // there's no way to have this appear in the output png without the rest of the GUI
  const int numrows = 4 + (numPanels > 0 ? 1 : 0) + (numFieldPts > 0 ? 1 : 0);
  // std::cout << "   fontSize: " << fontSize << "\n   numrows: " << numrows << "\n   display_h: " << display_h << std::endl;
#ifdef __APPLE__
  ImGui::SetNextWindowSize(ImVec2(10+fontSize*11,10+1.1*fontSize*numrows));
  //ImGui::SetNextWindowPos(ImVec2(10.0f, (displayH-34*numrows-1.1*fontSize*4)/2));
  ImGui::SetNextWindowPos(ImVec2(20, (displayH-fontSize*(1.1*numrows+1)-30*numrows)/2));
#else
  ImGui::SetNextWindowSize(ImVec2(10+fontSize*11, 10+1.1*fontSize*numrows));
  ImGui::SetNextWindowPos(ImVec2(20, displayH-fontSize*(1.1*numrows+1)));
#endif
  ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
  bool h = ImGui::Begin("Statistics", showStatsWindow, window_flags);
  if (!h) {
    std::cout << "Uh oh" << std::endl;
  }
  ImGui::Text("Step %13ld", step);
  ImGui::Text("Time %13.4f", time);
  if (numPanels > 0) { ImGui::Text("Panels %11ld", numPanels); }
  ImGui::Text("Particles %8ld", numParticles);
  if (numFieldPts > 0) { ImGui::Text("Field Pts %8ld", numFieldPts); }
  ImGui::End();
}

bool draw_welcome_window(const float displayW, const float displayH) {
  bool show = true;
  ImGui::OpenPopup("Welcome!");
  ImGui::SetNextWindowSize(ImVec2(500,300));
  ImGui::SetNextWindowPos(ImVec2(displayW * 0.5f, displayH * 0.5f), ImGuiCond_Always, ImVec2(0.5f,0.5f));
  ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
  if (ImGui::BeginPopupModal("Welcome!", NULL, window_flags)) {
    //ImGui::Begin("Welcome", &show_welcome_window);
    ImGui::TextWrapped("Welcome to Omega2D! Select a simulation from the drop-down, load from a file, or manually set your simulation global properites and add one or more flow, boundary, or measurement structures. Space bar starts and stops the run, Reset clears and loads new simulation properties. Left mouse button drags the frame around, mouse scroll wheel zooms. Save your flow set-up to json or your flow image to png or vtu. Have fun!");
    ImGui::Spacing();
    //if (ImGui::Button("Got it.", ImVec2(120,0))) { show_welcome_window = false; }
    //ImGui::End();
    // const float xwid = ImGui::GetWindowContentRegionWidth();
    if (ImGui::Button("Got it", ImVec2(120,0))) {
      ImGui::CloseCurrentPopup(); 
      show = false;
    }
    ImGui::EndPopup();
  }
  return show;
}
#endif
